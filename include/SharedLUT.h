#include <vector>
#include <memory>
#include <algorithm>
#include <mutex>

#ifndef SHARED_LUT_H
#define SHARED_LUT_H

template <class T, class Allocator=std::allocator<T>>
class SharedLUT {
  /* A lookup-table indexed by integers iBegin,...,iEnd which is thread-safe
   * and retains lock-free, O(1) access on reads.
   * A C++ mutex lock is used for the appending or prepending of new entries to the table.
   * It is assumed that pre/appending is rare, and deleting entries is not allowed.
   * No iterators will be issued.
   */
 public:
  SharedLUT() {}

  SharedLUT(const SharedLUT& rhs) {
    // temporarily share rhs LUT to copy it
    LUTptr tmp_ptr;
    std::atomic_store(&tmp_ptr,rhs.lutptr);
    lutptr.reset(new LUT(*tmp_ptr));
  }

  SharedLUT& operator=(const SharedLUT& rhs) {
    if (this == &rhs) return *this; // Self-assignment
    // temporarily share rhs LUT to copy it
    LUTptr tmp_ptr;
    std::atomic_store(&tmp_ptr,rhs.lutptr);
    //??LUTptr tmp_ptr(rhs.lutptr);
    lutptr.reset(new LUT(*tmp_ptr)); // ???
    return *this;
  }
    
  template <typename Iter>
  SharedLUT(int iStart, Iter it, Iter end) {
    extend(iStart,it,end);
  }

  // Read element
  T operator[](size_t i) const {
    // Get a shared pointer to the current vector of data,
    // so it will persist until the operation is done, even
    // if another thread changes the current vector's location.
    LUTptr tmp_ptr;
    std::atomic_store(&tmp_ptr,lutptr);
    return (*tmp_ptr)[i];
  }
  // Range-checked read
  T at(size_t i) const {
    LUTptr tmp_ptr;
    std::atomic_store(&tmp_ptr,lutptr);
    if (!tmp_ptr || i<tmp_ptr->iStart() || i>=tmp_ptr->iEnd())
      throw std::out_of_range("SharedLUT index out of bounds");
    return (*tmp_ptr)[i];
  }

  // Current element counts
  bool empty() const {return !lutptr;}
  const size_t size() const {return lutptr ? lutptr->size() : 0;}
  const int iStart() const {return lutptr->iStart();}
  const int iEnd() const {return lutptr->iEnd();}

  // Extend LUT by one element or a sequence.
  // argument i is the index for first write.
  // It is an error to attempt to write that would leave
  // a gap between extension and existing data.
  // Any writes
  // to table elements that already exist are IGNORED.
  void extend(int i, const T& t) {
    // Acquire the lock
    std::lock_guard<std::mutex> g(writeLock);
    if (empty()) {
      auto newlut = std::make_shared<LUT>(); //new LUT();
      newlut->i0 = i;
      newlut->data.push_back(t);
      // Install the LUT atomically
      std::atomic_store(&lutptr,newlut);
    } else if (i>iEnd() || i<iStart()-1) {
      throw std::runtime_error("ERROR: non-contiguous SharedLUT::extend()");
    } else if (i==iStart()-1) {
      // prepend - make new copy as all iterators become invalid
      auto newlut = std::make_shared<LUT>(); //new LUT();
      newlut->data.resize(size()+1);
      newlut->i0 = lutptr->i0-1;
      // Install new point first
      (newlut->data)[0] = t;
      // Then copy remainder over
      std::copy(lutptr->data.begin(), lutptr->data.end(), newlut->data.begin()+1);
      // And redirect lutptr to new one.
      std::atomic_store(&lutptr,newlut);
    } else if (i==iEnd()) {
      if (lutptr->data.capacity() > lutptr->size()) {
	// Can safely just append
	lutptr->data.push_back(t);
      } else {
	// Need to expand capacity, so get new LUT object
	auto newlut = std::make_shared<LUT>(); //new LUT();
	newlut->data.resize(size()+1);
	newlut->i0 = lutptr->i0;
	// Copy data
	std::copy(lutptr->data.begin(), lutptr->data.end(), newlut->data.begin());
	// Install new point 
	(newlut->data)[iEnd()] = t;
	// And redirect lutptr to new one.
	std::atomic_store(&lutptr,newlut);
      }
    }
    // Lock released on exit
  }
  
  template <typename Iter>
  void extend(int i, Iter it, Iter end) {
    int nAdd = end-it;
    // Freeze out other changes to lut
    std::lock_guard<std::mutex> g(writeLock);
    /***
#pragma omp critical (io)
    {
      cerr << "->Extending from " << i
             << " to " << i+nAdd
	       << " currently ";
      if (empty())
	cerr << "empty" << endl;
      else
	cerr << iStart() << "--" << iEnd()  << endl;
    }
    **/
    if (empty()) {
      // Create and install new LUT
      auto newlut = std::make_shared<LUT>(); //new LUT();
      newlut->i0 = i;
      newlut->data.resize(nAdd);
      std::copy(it,end,newlut->data.begin());
      std::atomic_store(&lutptr,newlut);
      /**
#pragma omp critical (io)
      cerr << "<-Init to  " << i
	       << " to " << i+nAdd
	       << " currently " << iStart() << "--" << iEnd()
	       << endl;
      **/
      return;
    } else if (i>iEnd() || i+nAdd<iStart()) {
      // Leaving a gap - not allowed
      throw std::runtime_error("ERROR: non-contiguous SharedLUT::extend()");
    } else if (i>=iStart() && i+nAdd<=iEnd()) {
      // "extension" already contained in LUT.  Nothing to do.
      /***
#pragma omp critical (io)
      cerr << "<-Nop for  " << i
	       << " to " << i+nAdd
	       << endl;
      ***/
      return;
    }

    // Need to append and/or prepend.
    int newStart = std::min(i,iStart());
    int newEnd = std::max(i+nAdd, iEnd());
    size_t newsize = newEnd - newStart;
    if (newStart<iStart()) {
      // Prepending to do.  Make new LUT
      auto newlut = std::make_shared<LUT>(); //new LUT();
      newlut->i0 = newStart;
      newlut->data.resize(newsize);
      // Copy in new data until we reach the old stuff
      auto newit = newlut->data.begin();
      for (int j=i; j<iStart(); ++it, ++j, ++newit) {
	*newit = *it;
      }
      // Then copy the old data
      newit = std::copy(lutptr->data.begin(), lutptr->data.end(), newit);
      it += size();
      
      // Also copy appended elements if any
      if (it < end)
	std::copy(it,end,newit);

      // Install new table
      std::atomic_store(&lutptr,newlut);
      /***
#pragma omp critical (io)
      cerr << "<-Prepended " << i
	       << " to " << i+nAdd
	       << " currently " << iStart() << "--" << iEnd()
	       << endl;
      ***/
      return;
    }

    // Need to append?
    if (newEnd > iEnd()) {
      if (false) { //??lutptr->data.capacity() >= newsize) {
	// Just append to existing array, will not invalidate iterators
	size_t oldsize = lutptr->data.size();
	lutptr->data.resize(newsize);
	std::copy( it + (lutptr->i0 + oldsize -i), end, lutptr->data.begin()+oldsize);
      } else {
	// Need to expand array, make & install new LUT
	auto newlut = std::make_shared<LUT>(); //new LUT();
	newlut->i0 = newStart;
	newlut->data.resize(newsize);
	auto newit = newlut->data.begin();
	// Copy the old data
	newit = std::copy(lutptr->data.begin(), lutptr->data.end(), newit);
	std::copy( it + (iEnd()-i), end, newit);

	// Install new table
	std::atomic_store(&lutptr,newlut);
      }
      
      /***
#pragma omp critical (io)
      cerr << "<-Appended " << i
	       << " to " << i+nAdd
	       << " currently " << iStart() << "--" << iEnd()
	       << endl;
      **/
    }
    return;
  }

private:
 private:
  // Define private structure holding the data and the starting index
  struct LUT {
  public:
    int i0;  // Index that accesses first element of vector
    std::vector<T,Allocator> data;
    size_t size() const {return data.size();}
    int iStart() const {return i0;}
    int iEnd() const {return i0 + size();}
    T operator[](int i) const {return data[i-i0];} // No range check
  };
  // The table structure itself: hold via shared pointer
  typedef std::shared_ptr<LUT> LUTptr;
  LUTptr lutptr;  
  // The write lock for this realization
  std::mutex writeLock;
};

#endif
