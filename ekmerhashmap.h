#ifndef EKMERHASHMAP_H
#define EKMERHASHMAP_H

#include "ekmerhashmap_dec.h"

#include <eutils/logger.h>

const size_t CSTRHASH_INIT_MASK=0xFF;

template <int KSIZE,class T,size_t (*hashfunc)(const char*,int)>
ekmerhashmap<KSIZE,T,hashfunc>::iter::iter(): hashmap(0x00), hitem(0x00), bucket(0) {}

template <int KSIZE,class T,size_t (*hashfunc)(const char*,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter& ekmerhashmap<KSIZE,T,hashfunc>::iter::operator++()
{
  hitem=hitem->next;
  
  if (hitem) return(*this);

  ++bucket;
  while (hashmap->_hashitems[bucket]==0x00 && bucket < hashmap->_hashmask+1) ++bucket;

  if (bucket < hashmap->_hashmask+1)
    hitem=hashmap->_hashitems[bucket];
  return(*this); 
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
const char *ekmerhashmap<KSIZE,T,hashfunc>::iter::key() const
{
  lddieif(hitem==0x00,"trying to access end iterator");
  return(hitem->key);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
size_t ekmerhashmap<KSIZE,T,hashfunc>::iter::hash() const
{
  lddieif(hitem==0x00,"trying to access end iterator");
  return(hitem->hash);
}


template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T& ekmerhashmap<KSIZE,T,hashfunc>::iter::value() const
{
  lddieif(hitem==0x00,"trying to access end iterator");
  return(*hitem->value);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T& ekmerhashmap<KSIZE,T,hashfunc>::iter::operator*() const
{
  lddieif(hitem==0x00,"trying to access end iterator");
  return(*hitem->value);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T* ekmerhashmap<KSIZE,T,hashfunc>::iter::operator->() const
{
  lddieif(hitem==0x00,"trying to access end iterator");
  return(hitem->value);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
bool ekmerhashmap<KSIZE,T,hashfunc>::iter::operator==(const ekmerhashmap<KSIZE,T,hashfunc>::iter& i) const
{
  return(hashmap==i.hashmap && hitem==i.hitem);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
bool ekmerhashmap<KSIZE,T,hashfunc>::iter::operator!=(const ekmerhashmap<KSIZE,T,hashfunc>::iter& i) const
{
  return(hashmap!=i.hashmap || hitem!=i.hitem);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter& ekmerhashmap<KSIZE,T,hashfunc>::iter::operator=(const ekmerhashmap<KSIZE,T,hashfunc>::iter& i)
{
  hashmap=i.hashmap;
  bucket=i.bucket;
  hitem=i.hitem;
  return(*this);
}




template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
ekmerhashmap<KSIZE,T,hashfunc>::ekmerhashmap(): count(0)
{
  _hashmask = CSTRHASH_INIT_MASK;
  _hashitems=new ekmerhashitem<KSIZE,T>*[_hashmask+1];
  size_t i;
  for (i=0; i<_hashmask+1; ++i)
    _hashitems[i]=0x00;
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
ekmerhashmap<KSIZE,T,hashfunc>::ekmerhashmap(const ekmerhashmap<KSIZE,T,hashfunc>& oldhm): count(oldhm.count)
{
  _hashmask = oldhm._hashmask;
  _hashitems=new ekmerhashitem<KSIZE,T>*[_hashmask+1];

  ekmerhashitem<KSIZE,T> *oldhmitem;
  size_t i;
  for (i=0; i<_hashmask+1; ++i){
    _hashitems[i]=0x00;
    for (oldhmitem=oldhm._hashitems[i]; oldhmitem!=0x00; oldhmitem=oldhmitem->next)
      _hashitems[i]=new ekmerhashitem<KSIZE,T>(oldhmitem->key,oldhmitem->value,_hashitems[i],oldhmitem->hash);
  }
}




template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
ekmerhashmap<KSIZE,T,hashfunc>::~ekmerhashmap()
{
  clear();
  delete[] _hashitems;
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter ekmerhashmap<KSIZE,T,hashfunc>::get(const char * key) const
{
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter it;

  it.hashmap=this;
  it.bucket=hash(key) & _hashmask;
  it.hitem=_hashitems[it.bucket];

  while (it.hitem!=0x00){
    if (memcmp(key,it.hitem->key,KSIZE)==0)   // there is no collision
      return(it);
    it.hitem=it.hitem->next;
  }

  return(end());
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter ekmerhashmap<KSIZE,T,hashfunc>::gethash(size_t hash,const char * key) const
{
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter it;

  it.hashmap=this;
  it.bucket=hash & _hashmask;
  it.hitem=_hashitems[it.bucket];

  while (it.hitem!=0x00){
    if (memcmp(key,it.hitem->key,KSIZE)==0)   // there is no collision
      return(it);
    it.hitem=it.hitem->next;
  }

  return(end());
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter ekmerhashmap<KSIZE,T,hashfunc>::addhash(size_t hash,const char * key,const T& value)
{
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter it;

  it.hashmap=this;
  it.bucket=hash & _hashmask;
  it.hitem=_hashitems[it.bucket];

  while (it.hitem!=0x00){
    if (memcmp(key,it.hitem->key,KSIZE)==0){
      it.hitem->value=new T(value);
      return(it);
    }
    it.hitem=it.hitem->next;
  }

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();

  it.bucket=hash & _hashmask;

  // non existent value
  ++count;
  _hashitems[it.bucket]=new ekmerhashitem<KSIZE,T>(key,new T(value),_hashitems[it.bucket],hash);
  it.hitem=_hashitems[it.bucket];
  return(it);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter ekmerhashmap<KSIZE,T,hashfunc>::begin() const
{
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter i;
  i.hashmap=this;
  i.bucket=0;
  i.hitem=0x00;
  while (_hashitems[i.bucket]==0x00 && i.bucket<_hashmask+1) ++i.bucket;

  if (i.bucket==_hashmask+1)
    return(end());

  i.hitem=_hashitems[i.bucket];
  return(i);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
typename ekmerhashmap<KSIZE,T,hashfunc>::iter ekmerhashmap<KSIZE,T,hashfunc>::end() const
{
  typename ekmerhashmap<KSIZE,T,hashfunc>::iter i;
  i.hashmap=this;
  i.bucket=0;
  i.hitem=0x00;
  return(i);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
void ekmerhashmap<KSIZE,T,hashfunc>::reserve(size_t i)
{
  size_t a;
  size_t c=1;
  a=0x01;
  while (i>0){
    i=i>>1;
    a=(a<<1)|0x01;
    ++c;
  }

  if (c>=sizeof(size_t)*8){
    lwarn("reached limit of hash table index size");
    a=0x8000000000000000u-0x01u;
  }
  
  resizehash(a);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
void ekmerhashmap<KSIZE,T,hashfunc>::resizehash(size_t newmask) const
{
  size_t thashmask;
  ekmerhashitem<KSIZE,T> **thashitems;

  if (newmask < _hashmask) return;

  ldinfo("resizing hash table");

  if (newmask==0) 
    thashmask = (_hashmask << 1)|0x01;
  else
    thashmask = newmask;

  thashitems=new ekmerhashitem<KSIZE,T>*[thashmask+1u];
  ldieif(thashitems==0x00,"unable to allocate memory for hashmap");
  size_t i;
  for (i=0; i<thashmask+1u; ++i)
    thashitems[i]=0x00;

  size_t khash;
  ekmerhashitem<KSIZE,T>* hitem;

  typename ekmerhashmap<KSIZE,T,hashfunc>::iter it;
  for (it=begin(); it!=end(); ++it){
    khash = it.hash() & thashmask;
    hitem = it.hitem; //gethashitem(khash & _hashmask,it.key());

    if (hitem->prev)
      hitem->prev->next=hitem->next;
    else
      _hashitems[khash&_hashmask]=hitem->next;
    if (hitem->next)
      hitem->next->prev=hitem->prev;

    hitem->prev=0x00;
    hitem->next=thashitems[khash];
    thashitems[khash]=hitem;
    if (hitem->next)
      hitem->next->prev=hitem;
  } 

  delete[] _hashitems;
  _hashmask=thashmask;
  _hashitems=thashitems;
  ldinfo("finished resize");
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
ekmerhashitem<KSIZE,T>* ekmerhashmap<KSIZE,T,hashfunc>::gethashitem(size_t khash,const char *key) const
{
  ekmerhashitem<KSIZE,T>* hitem;
  hitem = _hashitems[khash];
  while (hitem != 0x00){
    if (memcmp(hitem->key,key,KSIZE)==0)
      return(hitem);
    hitem=hitem->next;
  }
  lderror("ekmerhashmap: did not find key");
  return(0x00); 
}


template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
void ekmerhashmap<KSIZE,T,hashfunc>::clear()
{
  size_t i;
  ekmerhashitem<KSIZE,T> *hitem;
  for (i=0; i<_hashmask+1; ++i){
    while (_hashitems[i]){
      hitem=_hashitems[i];
      _hashitems[i]=hitem->next;
      delete hitem->value;
      delete hitem;
    }
  }
  count=0;
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
bool ekmerhashmap<KSIZE,T,hashfunc>::exists(const char * key) const
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0)
      return(true);
    hitem=hitem->next;
  }
  // non existent value
  return(false);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
void ekmerhashmap<KSIZE,T,hashfunc>::erase(const typename ekmerhashmap<KSIZE,T,hashfunc>::iter& it)
{
  ldieif(it.hitem==0x00,"trying to delete empty iterator");
  if (it.hitem->prev) it.hitem->prev->next=it.hitem->next;
  else _hashitems[it.bucket]=it.hitem->next;
  if (it.hitem->next) it.hitem->next->prev=it.hitem->prev;
  delete it.hitem->value;
  delete it.hitem;
//  it.hitem=0x00;
  --count;
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
void ekmerhashmap<KSIZE,T,hashfunc>::erase(const char * key)
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  size_t j;
  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0){
      if (hitem->prev) hitem->prev->next=hitem->next;
      else _hashitems[i]=hitem->next;
      if (hitem->next) hitem->next->prev=hitem->prev;
      delete hitem->value;
      delete hitem;
      --count;
      return;
    }
    hitem=hitem->next;
  }
  lddie("tried to delete key from hashmap that does not exist");
  // non existent value
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
ekmerhashmap<KSIZE,T,hashfunc>& ekmerhashmap<KSIZE,T,hashfunc>::operator+=(const ekmerhashmap<KSIZE,T,hashfunc>& hm)
{
  size_t i;
  for (i=0; i<hm.size(); ++i)
    add(hm.keys(i),hm.values(hm.keys(i)));
  return(*this);
}


template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T& ekmerhashmap<KSIZE,T,hashfunc>::addref(const char * key,T* value)
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();


  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0){   // there is no collision
      hitem->value = value;
      return(*hitem->value);
    }
    hitem=hitem->next;
  }

  // non existent value
  ++count;
  _hashitems[i]=new ekmerhashitem<KSIZE,T>(key,value,_hashitems[i],um_hash);
  return(*_hashitems[i]->value);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T& ekmerhashmap<KSIZE,T,hashfunc>::add(const char * key,const T& value)
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0){   // there is no collision
      hitem->value=new T(value);
      return(*hitem->value);
    }
    hitem=hitem->next;
  }

  // non existent value
  ++count;
  _hashitems[i]=new ekmerhashitem<KSIZE,T>(key,new T(value),_hashitems[i],um_hash);
  return(*_hashitems[i]->value);
/*
  operator[](key) = value;
  return(operator[](key));
*/
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
const T& ekmerhashmap<KSIZE,T,hashfunc>::operator[](const char * key) const
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0)   // there is no collision
      return(*hitem->value);
    hitem=hitem->next;
  }

  // non existent value
//  ++count;
  _hashitems[i]=new ekmerhashitem<KSIZE,T>(key,new T,_hashitems[i],um_hash);
  return(*_hashitems[i]->value);
}

/*
template <int KSIZE,class T>
const T& ekmerhashmap<KSIZE,T,hashfunc>::operator[](int ind) const
{
  int i;
  ehashitem<T>* hitem;

  lddieif(ind > _keys.size(),"ekmerhashmap: index out of bounds");
  i=hash(_keys.at(ind));
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (_keys.at(ind) == hitem->key)   // there is no collision
      return(*hitem->value);
    hitem=hitem->next;
  }

  ldie("ekmerhashmap: index out of bounds: "+estr(ind));
//  return(*(T*)0x00);
}
*/

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T& ekmerhashmap<KSIZE,T,hashfunc>::operator[](const char * key)
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0)   // there is no collision
      return(*hitem->value);
    hitem=hitem->next;
  }

//  lerror("key not found");
//  throw "key not found";

  // non existent value
  ++count;
  _hashitems[i]=new ekmerhashitem<KSIZE,T>(key,new T,_hashitems[i],um_hash);
  return(*_hashitems[i]->value);
}

//#pragma GCC diagnostic ignored "-Wreturn-type"
template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
const T& ekmerhashmap<KSIZE,T,hashfunc>::values(const char * key) const
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0)   // there is no collision
      return(*hitem->value);
    hitem=hitem->next;
  }

  lerror("ekmerhashmap: key not found");
  throw "ekmerhashmap: key not found";
//  return(*(T*)0x00);
}

template <int KSIZE,class T,size_t (*hashfunc)(const char *,int)>
T& ekmerhashmap<KSIZE,T,hashfunc>::values(const char * key)
{
  size_t i;
  ekmerhashitem<KSIZE,T>* hitem;

  if (size() > (3u*(_hashmask+1))>>2u) resizehash();

  size_t um_hash=hash(key);
  i=um_hash & _hashmask;
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (memcmp(key,hitem->key,KSIZE)==0)   // there is no collision
      return(*hitem->value);
    hitem=hitem->next;
  }
  
  lerror("ekmerhashmap: key not found");
  throw "ekmerhashmap: key not found";
//  return(*(T*)0x00);
}

/*
template <int KSIZE,class T>
T& ekmerhashmap<KSIZE,T,hashfunc>::operator[](int ind)
{
  int i;
  ehashitem<T>* hitem;

  i=hash(_keys.at(ind));
  hitem=_hashitems[i];

  while (hitem!=0x00){
    if (_keys.at(ind) == hitem->key)   // there is no collision
      return(*hitem->value);
    hitem=hitem->next;
  }

  ldie("ekmerhashmap: index out of bounds: "+estr(ind));
  return(*(T*)0x00);
}
//#pragma GCC diagnostic warning "-Wreturn-type"
*/

/*
template <unsigned int (*hashfunc)(const evar&)>
void ekmerhashmap<evar,evar,hashfunc>::addvar(evar& key,evar& var)
{
  add(key,var);
}

template <unsigned int (*hashfunc)(const evar&)>
evar ekmerhashmap<evar,evar,hashfunc>::getvar(int i) const
{
  return(evar());
//  return(values(i));
}

template <unsigned int (*hashfunc)(const evar&)>
evar ekmerhashmap<evar,evar,hashfunc>::getvarkey(int i) const
{
  return(evar());
//  return(keys(i));
}
*/

#endif

