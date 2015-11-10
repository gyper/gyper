#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

// functions to split a string by a specific delimiter
#include <string>
#include <vector>
#include <sstream>
#include <seqan/basic.h>
#include <seqan/sequence.h>

// join a vector of elements by a delimiter object.  ostream<< must be defined
// for both class S and T and an ostream, as it is e.g. in the case of strings
// and character arrays
template<typename TSeqanVariable typename TContainer, typename TDelimiter>
seqan::String<TSeqanVariable> join (TContainer & elems, TDelimiter & delimiter)
{
  std::stringstream ss;
  typename TContainer::iterator e = elems.begin();
  ss << *e++;

  for (; e != elems.end(); ++e)
  {
    ss << delimiter << *e;
  }

  return ss.str();
}

#endif
