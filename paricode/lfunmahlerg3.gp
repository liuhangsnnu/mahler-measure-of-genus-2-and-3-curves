\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ compare mahler measure of reciprocal genus 3 families and special value of L-functions
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

read("data3.gp");
read("mahler.gp");


{
getLengthG3(index) =
if(index==1,return(length(alist1)));
if(index==2,return(length(alist2)));
if(index==3,return(length(alist3)));
if(index==4,return(length(alist4)));
if(index==5 || index==6,return(length(alist5)));
if(index==7,return(length(alist7)));
}

{
getAByindexG3(num,index) =
if(index==1,return(alist1[num]));
if(index==2,return(alist2[num]));
if(index==3,return(alist3[num]));
if(index==4,return(alist4[num]));
if(index==5 || index==6,return(alist5[num]));
if(index==7,return(alist7[num]));
}

{
getCByindexG3(A,index) =
if(index==1 || index==2,return(x^(5-poldegree(A))*A));
if(index==3 || index==4 || index==5 || index==6,return(x^(6-poldegree(A))*A));
if(index==7,return(A));
}

{
getBG3ByindexG3(num,k,index) =
if(index==1,l=-k-1+llist1[num];return(x^5+1+k*x^4+l*x^3+l*x^2+k*x));
if(index==2,l=-k+llist2[num];return(k*x^4+l*x^3+l*x^2+k*x));
if(index==3,l=llist3[num];m=-2*k+mlist3[num];return(k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==4,l=-1+llist4[num];m=-2*k+mlist4[num];return(x^6+1+k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==5,l=-1+llist5[num];m=-2*k+mlist5[num];return(x^6+1+k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==6,l=-2+llist5[num];m=-2*k+mlist5[num];return(2*x^6+2+k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==7,l=llist7[num];m=-2*k+mlist7[num];return(2*x^8+2+klist7[num]*(x^7+x)+k*(x^6+x^2)+l*(x^5+x^3)+m*x^4));
}

{
getBG3(num,index) =
if(index==1,l=-k-1+llist1[num];return(x^5+1+k*x^4+l*x^3+l*x^2+k*x));
if(index==2,l=-k+llist2[num];return(k*x^4+l*x^3+l*x^2+k*x));
if(index==3,l=llist3[num];m=-2*k+mlist3[num];return(k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==4,l=-1+llist4[num];m=-2*k+mlist4[num];return(x^6+1+k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==5,l=-1+llist5[num];m=-2*k+mlist5[num];return(x^6+1+k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==6,l=-2+llist5[num];m=-2*k+mlist5[num];return(2*x^6+2+k*x^5+l*x^4+m*x^3+l*x^2+k*x));
if(index==7,l=llist7[num];m=-2*k+mlist7[num];return(2*x^8+2+klist7[num]*(x^7+x)+k*(x^6+x^2)+l*(x^5+x^3)+m*x^4));
}

\\ write results to the tex file
filepath="tableGenus3.tex";


{
calcmahlerG3(m,n,index) =
for(num=1,getLengthG3(index),
\\for(num=25,25,
print("===========================================================");
print("num = ", num);
A=getAByindexG3(num,index);
C=getCByindexG3(A,index);
cl=calcEllipticFactor(getBG3(num,index)^2-4*A*C);
for(k=m,n,
  B=getBG3ByindexG3(num,k,index);
  calclfunmahler(k,cl);
);
);
}

{
calcmahlerG3Tex(m,n,index) =
Y;X;
k=varlower("k");
for(num=1,getLengthG3(index),
print("===========================================================");
print("num = ", num);
A=getAByindexG3(num,index);
C=getCByindexG3(A,index);
cl=calcEllipticFactor(getBG3(num,index)^2-4*A*C);
my(strA="1");
my(strC=pol2Str(A));
if(A!=1,strA=pol2Str(A));
if(C/A!=1,strC=Str(Strtex(C/A),strC));
write(filepath,Str("$$","A=",strA,"$$"));
write(filepath,Str("$$","B=",Strtex(getBG3(num,index)),"$$"));
write(filepath,Str("$$","C=",strC,"$$"));
write(filepath,"\\begin{longtable}{|l|l|l|lllll|}");
write(filepath,"\\hline");
write(filepath,"$k$ & $s$ & $N_E$ & $a_1$ & $a_2$ & $a_3$ & $a_4$ & $a_5$\\\\");
write(filepath,"\\hline");
for(kk=m,n,
  B=getBG3ByindexG3(num,kk,index);
  print(kk);
  print(subst(eval(cl),k,kk));
  calclfunmahler(kk,subst(eval(cl),k,kk));
);
write(filepath,"\\hline");
write(filepath,"\\end{longtable}");
);
}

{
calcmahlerG3TexAll(m,n) =
write(filepath,"\\documentclass{amsart}");
write(filepath,"\\usepackage{longtable}");
write(filepath,"\\usepackage{amssymb}");
write(filepath,"\\begin{document}");
for(index=1,7,calcmahlerG3Tex(m,n,index));
write(filepath,"\\end{document}");
}

{
calcjInvariantG3(index) =
my(vec=vector(getLengthG3(index)));
for(num=1,getLengthG3(index),
my(A=getAByindexG3(num,index));
my(C=getCByindexG3(A,index));
my(B=getBG3(num,index));
my(cl=calcEllipticFactor(B^2-4*A*C));
my(E = ellinit(cl));
vec[num]=[E.j,A*y^2+B*y+C];
);
return(vec);
}
