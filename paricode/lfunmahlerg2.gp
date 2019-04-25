\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ compare mahler measure of reciprocal genus 2 families and special value of L-functions
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


read("data2.gp");
read("mahler.gp");


{
getLengthG2(index) =
if(index==1,return(length(alist1g2)));
if(index==2,return(length(alist2g2)));
if(index==3,return(length(alist3g2)));
if(index==4,return(length(alist4g2)));
if(index==5,return(length(alist5g2)));
if(index==6,return(length(alist6g2)));
if(index==7,return(length(alist7g2)));
if(index==8,return(length(alist8g2)));
if(index==9,return(length(alist9g2)));
}

{
getAByindexG2(num,index) =
if(index==1,return(alist1g2[num]));
if(index==2,return(alist2g2[num]));
if(index==3,return(alist3g2[num]));
if(index==4,return(alist4g2[num]));
if(index==5,return(alist5g2[num]));
if(index==6,return(alist6g2[num]));
if(index==7,return(alist7g2[num]));
if(index==8,return(alist8g2[num]));
if(index==9,return(alist9g2[num]));
}

{
getCByindexG2(A,index) =
if(index==1 || index==2,return(subst(A,x,1/x)*x^3));
if(index==3 || index==4 || index==5 || index==6 || index==7,return(subst(A,x,1/x)*x^4));
if(index==8 || index==9,return(subst(A,x,1/x)*x^6));
}

{
getBG2ByindexG2(num,k,index) =
if(index==1,return(x^3+1+k*x^2+k*x));
if(index==2,return(k*x^2+k*x));
if(index==3,l=-2*k+llist3g2[num];return(k*x^3+l*x^2+k*x));
if(index==4,l=2*k+llist4g2[num];return(x^4+1+k*x^3+l*x^2+k*x));
if(index==5,l=-2*k+llist5g2[num];return(x^4+1+k*x^3+l*x^2+k*x));
if(index==6,l=2*k+llist6g2[num];return(2*x^4+2+k*x^3+l*x^2+k*x));
if(index==7,l=-2*k+llist7g2[num];return(2*x^4+2+k*x^3+l*x^2+k*x));
if(index==8,m=2*k+mlist8g2[num];return(2*x^6+2+klist8g2[num]*(x^5+x)+k*(x^4+x^2)+m*x^3));
if(index==9,m=-2*k+mlist9g2[num];return(2*x^6+2+klist9g2[num]*(x^5+x)+k*(x^4+x^2)+m*x^3));
}

{
getBG2(num,index) =
if(index==1,return(x^3+1+k*x^2+k*x));
if(index==2,return(k*x^2+k*x));
if(index==3,return(k*x^3+(-2*k+llist3g2[num])*x^2+k*x));
if(index==4,return(x^4+1+k*x^3+(2*k+llist4g2[num])*x^2+k*x));
if(index==5,return(x^4+1+k*x^3+(-2*k+llist5g2[num])*x^2+k*x));
if(index==6,return(2*x^4+2+k*x^3+(2*k+llist6g2[num])*x^2+k*x));
if(index==7,return(2*x^4+2+k*x^3+(-2*k+llist7g2[num])*x^2+k*x));
if(index==8,return(2*x^6+2+klist8g2[num]*(x^5+x)+k*(x^4+x^2)+(2*k+mlist8g2[num])*x^3));
if(index==9,return(2*x^6+2+klist9g2[num]*(x^5+x)+k*(x^4+x^2)+(-2*k+mlist9g2[num])*x^3));

}

\\ write results to the tex file
filepath="tableGenus2.tex";

{
calcmahlerG2(m,n,index) =
for(num=1,getLengthG2(index),
print("===========================================================");
print("num = ", num);
A=getAByindexG2(num,index);
C=getCByindexG2(A,index);
cl=calcEllipticFactor(getBG2(num,index)^2-4*A*C);
for(k=m,n,
  B=getBG2ByindexG2(num,k,index);
  calclfunmahler(k,cl);
);
);
}


{
calcmahlerG2Tex(m,n,index) =
Y;X;
k=varlower("k");
for(num=1,getLengthG2(index),
print("===========================================================");
print("num = ", num);
A=getAByindexG2(num,index);
C=getCByindexG2(A,index);
cl=calcEllipticFactor(getBG2(num,index)^2-4*A*C);
my(strA="1");
my(strC=pol2Str(A));
if(A!=1,strA=pol2Str(A));
if(C/A!=1,strC=Str(Strtex(C/A),strC));
write(filepath,Str("$$","A=",strA,"$$"));
write(filepath,Str("$$","B=",Strtex(getBG2(num,index)),"$$"));
write(filepath,Str("$$","C=",strC,"$$"));
write(filepath,"\\begin{longtable}{|l|l|l|lllll|}");
write(filepath,"\\hline");
write(filepath,"$k$ & $s$ & $N_E$ & $a_1$ & $a_2$ & $a_3$ & $a_4$ & $a_5$\\\\");
write(filepath,"\\hline");
for(kk=m,n,
  B=getBG2ByindexG2(num,kk,index);
  print(kk);
  print(subst(eval(cl),k,kk));
  calclfunmahler(kk,subst(eval(cl),k,kk));
);
write(filepath,"\\hline");
write(filepath,"\\end{longtable}");
);
\\write(filepath,"\\end{document}");
}

{
calcmahlerG2TexAll(m,n) =
write(filepath,"\\documentclass{amsart}");
write(filepath,"\\usepackage{longtable}");
write(filepath,"\\usepackage{amssymb}");
write(filepath,"\\begin{document}");
for(index=1,9,calcmahlerG2Tex(m,n,index));
write(filepath,"\\end{document}");
}


{
calcjInvariantG2(index) =
my(vec=vector(getLengthG2(index)));
for(num=1,getLengthG2(index),
my(A=getAByindexG2(num,index));
my(C=getCByindexG2(A,index));
my(B=getBG2(num,index));
my(cl=calcEllipticFactor(B^2-4*A*C));
my(E = ellinit(cl));
vec[num]=[E.j,A*y^2+B*y+C];
);
return(vec);
}
