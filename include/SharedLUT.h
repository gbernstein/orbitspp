#include <vector>
#include <memory>
#include <algorithm>
#include <mutex>

template <class T, class Allocator=std::allocator<T>>
class SharedLUT {
  /* A lookup-table indexed by integers 0,1,...,N which is thread-safe
   * and retains lock-free, O(1) access on reads.
   * A C++ mutex lock is used for the appending of new entries to the table.
   * It is assumed that appending is rare, and deleting entries is not allowed.
   * No iterators will be issued.
   */
 public:
  typedef std::shared_ptr<std::vector<T,Allocator>> Sptr;
  SharedLUT() {
    sptr.reset(new std::vector<T,Allocator>);
  }
  template <typename Iter>
  SharedLUT(Iter it, Iter end) {
    sptr.reset(new std::vector<T,Allocator>() );
    sptr->reserve(end-it);
    for ( ; it!=end; ++it)
      sptr->push_back(*it);
  }

  // Read element
  T operator[](size_t i) const {
    // Get a shared pointer to the current vector of data,
    // so it will persist until the operation is done, even
    // if another thread changes the current vector's location.
    Sptr tmp_ptr(sptr);
    return (*tmp_ptr)[i];
  }

  // Current element count
  const size_t size() const {return sptr->size();}

  // Write one or a sequence of elements to the table.
  // argument i is the index for first write.
  // It is an error to attempt to write that starts more than one
  // past the end of the table, and any writes
  // to table elements that already exist are IGNORED.
  void append(size_t i, const T& t) {
    if (i<sptr->size()) {
      // We already go this far
      return;
    } else if (i>sptr->size()+1) {
      // 
      throw std::runtime_error("ERROR: appending past the end of a SharedLUT");
    }
    // Acquire the lock
    std::lock_guard<std::mutex> g(writeLock);
    if (sptr->capacity() - sptr->size() > 0) {
      sptr->push_back(t);
    } else {
      // Make a fresh vector by hand
      auto vnew = new std::vector<T,Allocator>();
      vnew->resize(sptr->size() + 1);
      std::copy(sptr->begin(), sptr->end(), vnew->begin());
      // Add new element
      (*vnew)[sptr->size()] = t;
      // atomically point sptr to the new vector.
      // The smart pointer will destroy the old one if
      // no one else is using it.
      sptr.reset(vnew);
    }
  }
  template <typename Iter>
  void append(size_t i, Iter it, Iter end) {
    size_t nAdd = end-it;
    size_t s = sptr->size();
    if (i+nAdd <= s) {
      // Table already goes this far, nothing to do.
      return;
    } else if (i>s+1) {
      // 
      throw std::runtime_error("ERROR: appending past the end of a SharedLUT");
    }
    // Acquire the lock
    std::lock_guard<std::mutex> g(writeLock);
    s = sptr->size(); // Re-read size while locked
    // Advance input iterator for already-known elements
    for (int j=i; j<s && it!=end; ++j) {
      ++it;
      --nAdd;
    }
    if (sptr->capacity() - sptr->size() >= nAdd) {
      for (; it!=end; ++it)
	sptr->push_back(*it);
    } else {
      // Make a fresh vector by hand
      auto vnew = new std::vector<T,Allocator>();
      vnew->resize(s + nAdd);
      std::copy(sptr->begin(), sptr->end(), vnew->begin());
      // Add new elements (potentially rewriting some old ones)
      std::copy(it, end, vnew->begin()+s);
      // atomically point sptr to the new vector.
      // The smart pointer will destroy the old one if
      // no one else is using it.
      sptr.reset(vnew);
    }
  }

private:
  Sptr sptr;
  // The write lock for this realization
  std::mutex writeLock;
};
