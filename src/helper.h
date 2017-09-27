//#include <limits.h>
//#include "DPPmcmc.h"
#include <Rcpp.h>

/*
template <typename T>
std::vector<T> rep(T a,int num_reps)
{
  std::vector<T> result(num_reps);
  for(int i=0;i<num_reps;i++)
  {
     result[i]=a;
  }
  return result;
}
*/

template <typename T>
std::vector<double> divideVectorByDouble(std::vector<T> vector1,double denominator){
  std::vector<double> returnVector(vector1.size());
  for(int i=0;i<vector1.size();i++){
    returnVector[i]=(((double)vector1[i])/denominator);
  }
  return returnVector;

}


template <typename T>
std::vector<T> whichAreEqual(std::vector<T> vector1,T val){

  //return is 0 based
  std::vector<T> returnVector(0);

  for(int i=0;i<vector1.size();i++){
    if (vector1[i]==val) returnVector.push_back(i);

  }
  return returnVector;
}


template <typename T>
std::vector<T> logVector(std::vector<T> vector1)
{
  std::vector<T> result(vector1.size());
  for(int i=0;i<vector1.size();i++)
  {
    result[i]=log(vector1[i]);
  }
  return result;
}

template <typename T>
std::vector<T> expVector(std::vector<T> vector1)
{
  std::vector<T> result(vector1.size());
  for(int i=0;i<vector1.size();i++)
  {
    result[i]=exp(vector1[i]);
  }
  return result;
}
template <typename T>
T sumVector(std::vector<T> vector1)
{
  T result=0;
  for(int i=0;i<vector1.size();i++)
  {
    result=result+vector1[i];
  }
  return result;
}
template <typename T>
std::vector<T> concatenateVectors(std::vector<T> vector1,std::vector<T> vector2){
  std::vector<T> concatenated(vector1.size()+vector2.size());

  int vector1size=vector1.size();
  for (int i=0;i<vector1size;i++){
    concatenated[i]=vector1[i];
  }

  for (int i=0;i<vector2.size();i++){
    concatenated[i+vector1size]=vector2[i];
  }
  return concatenated;
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
  // assert(a.size() == b.size());
  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(),
                 std::back_inserter(result), std::plus<T>());
  return result;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
  // assert(a.size() == b.size());
  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(),
                 std::back_inserter(result), std::minus<T>());
  return result;
}

template <typename T>
bool anyEqual(std::vector<T> vector1,double val){
  bool returnBool=FALSE;
  for(int i=0;i<vector1.size();i++){
    if(vector1[i]==val) returnBool=TRUE;
  }
  return returnBool;
}

template <typename T>
std::vector<T> removeElementAtPosition(std::vector<T> vector1,int element){
  //element should be 0 based.
  std::vector<T> resultVector(vector1.size()-1);

  int counter=0;
  for(int i=0;i<vector1.size();i++){
    if (i!=element)
     {
        resultVector[counter]=vector1[i];
        counter++;
     }
  }
  return resultVector;
}

template <typename T>
std::vector<T> elementsInRange(int rangeMin, int rangeMax,std::vector<T> vector1){
  //range should be 0 based.
  std::vector<T> resultVector(0);

  for(int i=rangeMin;i<=rangeMax;i++){
    resultVector.push_back(vector1[i]);
  }
  return resultVector;
}





