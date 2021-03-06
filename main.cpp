#include <iostream>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

const ZZ n = conv<ZZ>(1018), N = conv<ZZ>(n + 1);  /* N = 1019 -- prime     */
 const ZZ alpha = conv<ZZ>(2);            /* generator             */
 const ZZ beta = conv<ZZ>(5);             /* 2^{10} = 1024 = 5*/


void new_xab( ZZ& x, ZZ& a, ZZ& b,ZZ p, ZZ n );
ZZ pollard_rho(ZZ G, ZZ g,ZZ p);
int main()
{
    ZZ x=conv<ZZ>(1), a=conv<ZZ>(0), b=conv<ZZ>(0);
   ZZ X=x, A=a, B=b;
   ZZ r;
   /*for(ZZ i = conv<ZZ>(1); i < N; ++i ) {
     new_xab( x, a, b );
     new_xab( X, A, B );
     new_xab( X, A, B );
     cout<<x<<" "<<a<<" "<<b<<" "<<X<<" "<<A<<" "<<B<<endl;
     if( x == X )break;

}
*/
    cout<<pollard_rho(N,alpha,n)<<endl;
return 0;
}




void new_xab( ZZ& x, ZZ& a, ZZ& b,ZZ p ,ZZ n ) {
   switch( x%3 ) {
   case 0: x = x*x     % p;  a =  a*2  % n;  b =  b*2  % n;  break;
   case 1: x = x*alpha % p;  a = (a+1) % n;                  break;
   case 2: x = x*beta  % p;                  b = (b+1) % n;  break;
   }
 }

 ZZ pollard_rho(ZZ G, ZZ g, ZZ p){
 ZZ x=conv<ZZ>(1), a=conv<ZZ>(0), b=conv<ZZ>(0);
   ZZ X=x, A=a, B=b;
   ZZ r;
   do{
     new_xab( x, a, b,N,n );
     new_xab( X, A, B,N,n );
     new_xab( X, A, B,N,n );
   }while(x!=X);
   return ((a-A)/(B-b));
 }


