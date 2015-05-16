// Doolittle uses unit diagonals for the lower triangle
void Doolittle(int d,double*S,double*D){
   for(int k=0;k<d;++k){
      for(int j=k;j<d;++j){
         double sum=0.;
         for(int p=0;p<k;++p)sum+=D[k*d+p]*D[p*d+j];
         D[k*d+j]=(S[k*d+j]-sum); // not dividing by diagonals
      }
      for(int i=k+1;i<d;++i){
         double sum=0.;
         for(int p=0;p<k;++p)sum+=D[i*d+p]*D[p*d+k];
         D[i*d+k]=(S[i*d+k]-sum)/D[k*d+k];
      }
   }
}
void solveDoolittle(int d,double*LU,double*b,double*x){
   double y[d];
   for(int i=0;i<d;++i){
      double sum=0.;
      for(int k=0;k<i;++k)sum+=LU[i*d+k]*y[k];
      y[i]=(b[i]-sum); // not dividing by diagonals
   }
   for(int i=d-1;i>=0;--i){
      double sum=0.;
      for(int k=i+1;k<d;++k)sum+=LU[i*d+k]*x[k];
      x[i]=(y[i]-sum)/LU[i*d+i];
   }
}
