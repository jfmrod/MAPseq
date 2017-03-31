#ifndef EKMERHASHMAP_DEC_H
#define EKMERHASHMAP_DEC_H

#include <eutils/ehashfunc.h>

inline size_t hash_lookup3_cstr(const char* cstr,int len)
{
#ifdef __i386__
  return(hash_lookup3(cstr,len,0x75483842));
#else  // __X86_64__
  return(hash_lookup3_64(cstr,len,0x75483842,0xB5749368));
#endif
}

template <int KSIZE,class T>
class ekmerhashitem
{
 public:
  size_t hash;
  const char *key;
  T *value;
  ekmerhashitem<KSIZE,T>* next;
  ekmerhashitem<KSIZE,T>* prev;

  ekmerhashitem(const char *key,T* value,ekmerhashitem<KSIZE,T>* next,size_t hash);
  ~ekmerhashitem();
};

template <int KSIZE,class T>
ekmerhashitem<KSIZE,T>::ekmerhashitem(const char *_key,T* _value,ekmerhashitem<KSIZE,T>* _next,size_t _hash): key(_key),value(_value),next(_next),prev(0x00),hash(_hash)
{
//  memcpy(key,_key,KSIZE);
  if (next)
    next->prev=this;
}

template <int KSIZE,class T>
ekmerhashitem<KSIZE,T>::~ekmerhashitem()
{
}


template <int KSIZE,class T,size_t (*hashfunc)(const char*,int)=hash_lookup3_cstr>
class ekmerhashmap
{
  mutable ekmerhashitem<KSIZE,T> **_hashitems;
  size_t count;
 public:
  mutable size_t _hashmask;

  ekmerhashmap();
  ekmerhashmap(const ekmerhashmap<KSIZE,T,hashfunc>& oldhm);
  ~ekmerhashmap();

  inline size_t size() const { return(count); }

  T& add(const char *key,const T& value);
  T& addref(const char *key,T* value);

  inline const T& operator[](const char *key) const;
  T& operator[](const char *key);

  ekmerhashmap<KSIZE,T,hashfunc>& operator+=(const ekmerhashmap<KSIZE,T,hashfunc>& hm);

  T& values(const char *key);
  const T& values(const char *key) const;

  bool exists(const char *key) const ;
  long findkey(const char *key,size_t pos=0) const;

  void clear();
  void resizehash(size_t i=0) const;

  void reserve(size_t i);
  void erase(const char *key);

  ekmerhashitem<KSIZE,T>* gethashitem(size_t khash,const char *key) const;

  inline size_t hash(const char *key) const { return(hashfunc(key,KSIZE)); }

  class iter
  {
   public:
    const ekmerhashmap<KSIZE,T,hashfunc> *hashmap;
    ekmerhashitem<KSIZE,T> *hitem;
    size_t  bucket;

    iter();

    iter& operator++();
 
    size_t hash() const;
    const char *key() const;
    T& value() const;

    T& operator*() const;
    T* operator->() const;

    bool operator==(const iter& i) const;
    bool operator!=(const iter& i) const;

    iter& operator=(const iter& i);
  };

  void erase(const typename ekmerhashmap<KSIZE,T,hashfunc>::iter&);
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter get(const char *key) const;
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter gethash(size_t hash,const char *key) const;
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter addhash(size_t hash,const char *key,const T&);
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter begin() const;
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter end() const;
};

#endif

