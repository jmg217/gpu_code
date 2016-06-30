#include <iostream>


int function(int pointer){
/*int c=0;
for(int i=0; i<5; i++){
c+=pointer[i];
}*/
return pointer;
}

int main(){

int N=2;
int M=3;
int L=4;

double foo [N];

for(int j=0; j<N; j++){

foo[j]=2;
}

for(int i=0; i<N; i++){
std::cout<<foo[i]<<std::endl;
}

double Foo[N][M][L];

for(int k=0; (k < N); k++)
   {
      
      for(int l=0; (l < M); l++)
      {
       
        for(int m=0; m<L; m++){
                Foo[k][l][m] = k*l*m;
                }
        }
   }



/*
for(int x=0; (x < N); x++)
   {
std::cout<<"x"<<std::endl;
      for(int y=0; (y < M); y++)
      {
std::cout<<"y"<<std::endl;
        for(int z=0; z<L; z++){
                std::cout<<Foo[x][y][z]<<std::endl;
                }
        }
   }
*/

int* pointer;
pointer= new int[5];

for(int it=0; it<5; it++){
pointer[it]=it;
}

int c= function(pointer[1]);

std::cout<<c<<std::endl;

delete[] pointer;

return 0;
}
