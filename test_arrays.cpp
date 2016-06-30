#include <iostream>


int main(){

int N=100;

double foo [N];

for(int j=0; j<100; j++){

foo[j]=2;
}

for(int i=0; i<100; i++){
std::cout<<foo[i]<<std::endl;
}
return 0;
}
