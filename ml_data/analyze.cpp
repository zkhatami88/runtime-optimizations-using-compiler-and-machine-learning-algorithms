namespace boost {
    struct irange
    {
      irange(unsigned a, unsigned b) {}
      
      int begin()
      {
          return 1;
      }
    };
}

namespace hpx { namespace parallel {
    
    template<typename P, typename B, typename E, typename F>
    void for_each(P p, B b, E e, F f)
    {
        
    }  
  
    template<typename P, typename B, typename E, typename F>
    void for_each_n(P p, B b, E e, F f)
    {
        
    } 


    double seq = 3;
    double par = 3;
}

namespace util
{
    struct high_resolution_timer
    {
        double elapsed()
        {
            return 0.;
        }
    };
}

    
    int finalize()
    {
        return 1;
    }

}

#include <iostream>


#include "functions/func.hpp"


int main(int argc, char* argv[])
{
    return 0;
}