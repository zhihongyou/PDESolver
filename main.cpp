#include <iostream> 
#include <vector>
#include <cstdio>
#include <string>
#include "src/system/systemclass.h"
#include "src/system/systemclass.cpp"
#include "external/exprtk/exprtk.hpp"

using namespace std;

template <typename T>
struct myfunc : public exprtk::ifunction<T>
{
   myfunc()
   : exprtk::ifunction<T>(2)
   { exprtk::disable_has_side_effects(*this); }

   inline T operator()(const T& v1, const T& v2)
   {
      return T(1) + (v1 * v2) / T(3);
   }
};

template <typename T>
inline T myotherfunc(T v0, T v1, T v2)
{
   return std::abs(v0 - v1) * v2;
}

template <typename T>
void custom_function(std::string myExpr)
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>   expression_t;
   typedef exprtk::parser<T>       parser_t;

   std::string expression_string = myExpr;

   T x = T(1);
   T y = T(2);

   myfunc<T> mf;

   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);
   symbol_table.add_variable("y",y);
   symbol_table.add_function("myfunc",mf);
   symbol_table.add_function("otherfunc",myotherfunc);
   symbol_table.add_constants();

   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);

   T result = expression.value();
   printf("Result: %10.5f\n",result);
}
    

int main() {

    System mySys;
    
    // mySys.set();
    std::cout << mySys.Dimension <<endl;
    
    custom_function<double>("myfunc(sin(x / pi), otherfunc(3 * y, x / 2, x * y))");
    return 0;
};


