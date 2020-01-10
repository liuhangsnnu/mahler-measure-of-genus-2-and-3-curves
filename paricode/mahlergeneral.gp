\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ calculate the Mahler measure of general two variable polynomial
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

read("jGenus1.gp");
read("mahler.gp");

defaultprecision=default(realprecision);
\\ rootprecision is used to determine if the root is on the unit circle
rootprecision=10^(-ceil(defaultprecision/2));

P;
\\ all roots of P
{
logr(u)=
my(x=exp(2*Pi*I*u));
my(pol=correctpol(eval(P)));
my(roots);
if(poldegree(pol,y)==2,roots=abs(rootsd2(pol)),roots=abs(polroots(pol)));
\\roots=abs(polroots(pol));
\\this is to avoid log of 0
for(i=1,#roots,
 if(roots[i]==0,roots[i]=1);
);
\\if(#roots==2 && (roots[1]>1 && roots[2]>1), print(u));
return(vecsort(log(roots)));
}

{
rootsd2(pol)=
my(A=polcoeff(pol,2,y));
my(B=polcoeff(pol,1,y));
my(C=polcoeff(pol,0,y));
return([(-B+sqrt(B^2-4*A*C))/(2*A),(-B-sqrt(B^2-4*A*C))/(2*A)]);
}

{
\\ if the real part or imaginary part of coefficients of pol equals 0
\\ it should be 0, bug of Pari?
correctpol(pol) =
my(degree=poldegree(pol));
my(result=0);
for(i=0,degree,
  my(coeff=polcoeff(pol,i));
  if(real(coeff)==0,coeff=imag(coeff)*I);
  if(imag(coeff)==0,coeff=real(coeff));
  result=result+coeff*y^i
);
return(result);
}

{
\\calculate the sum of the integrand
integrandgeneral(u,allintervallist) =
my(result=0);
\\ for degenerating case
my(logr=logr(u));
for(i=1,poldegree(P,y),
  my(flag=0);
  intervallist=allintervallist[i];
  for(j=1,#intervallist,
    if(intervallist[j][1]<=u && u<=intervallist[j][2],
      \\interval contains the point u
      \\print(u);
      if(i<=#logr,result=result+logr[i]);
      break;
    )
  )
);
return(result);
}

{
plotintegrand() =
allintervallist=splitallinterval();
ploth(u=0,0.5,integrandgeneral(u,allintervallist));
}

{
splitallinterval() =
my(allintervallist=listcreate());
for(i=1,poldegree(P,y),
  listput(allintervallist,splitintervalgeneral(i));
);
return(allintervallist);
}

{
unionlist(allintervallist) =
my(ulist=listcreate());
for(i=1,#allintervallist,
  for(j=1,#allintervallist[i],
    if(!inList(ulist, allintervallist[i][j]),listput(ulist,allintervallist[i][j]));
  )
);
return(ulist);
}

{
\\split the interval by ramification points and change sign points
\\calculate the integration intervals of function in logr(x) according to ramification points and change sign points
splitintervalgeneral(index) =
my(criticalPointsList=calcCriticalPoints());
my(intervallist=listcreate());
for(i=1,#criticalPointsList-1,
  my(x=(criticalPointsList[i]+99*criticalPointsList[i+1])/100);
  if(index<=#logr(x) && logr(x)[index]>0 && criticalPointsList[i] != criticalPointsList[i+1],
    listput(intervallist,[criticalPointsList[i],criticalPointsList[i+1]]);
  );
);
return(intervallist);
}




{
calcCriticalPoints() =
my(l1=calcRamiPoints());
\\print("l1 = ", l1);
my(l2=calcChangeSignPoints());
\\print("l2 = ", l2);
my(jl=joinlist(l1,l2));
\\add 1/6 for calculate m(P_k(x-1,y)) more accurately
\\listput(jl,1/6.0);
listsort(jl);
\\print(jl);
if(#jl==0,listput(jl,0,1);listput(jl,1/2,2));
\\ the endpoints are 0 and 1/2
if(jl[#jl]!=1/2,listput(jl,1/2));
if(jl[1]!=0,listinsert(jl,0,1));
return(jl);
}



\\ find ramification points and change sign points on unit circle
{
\\ramification points
calcRamiPoints() =
my(DP=deriv(P,y));
my(resultant=polresultant(P,DP,y));
\\print(resultant);
my(rootuc=listcreate());
if(resultant==0, return(rootuc));
my(roots=polroots(resultant));
\\print(abs(roots));
for(i=1,#roots,
  if(isOnUnitCircle(roots[i])==true,
    my(root=imag(log(roots[i]))/2/Pi);
    if(root>0 && root<1/2 && !inList(rootuc, root), listput(rootuc,root));
  );
);
return(rootuc);
}

{
\\change sign points, P is non-reciprocal
calcChangeSignPoints() =
my(Q=x^poldegree(P,x)*y^poldegree(P,y)*subst(subst(P,x,1/x),y,1/y));
\\print(Q);
my(resultant=polresultant(P,Q,y));
\\print(resultant);
my(rootuc=listcreate());
if(resultant==0, return(rootuc));
my(roots=polroots(resultant));
\\print(abs(roots));
for(i=1,#roots,
  if(isOnUnitCircle(roots[i])==true,
    my(yroots=polroots(subst(P,x,roots[i])));
    for(j=1,#yroots,
      if(isOnUnitCircle(yroots[j])==true,
        my(root=imag(log(roots[i]))/2/Pi);
        if(root>0 && root<1/2 && !inList(rootuc, root), listput(rootuc,root));
        break;
      );
    );
  )
);
return(rootuc);
}

{
isOnUnitCircle(p) =
if(abs(abs(p)-1)<rootprecision,return(true),return(false));
}




\\find infinity points
\\split the interval again if it contains more than 2 infinity points of logr(x)
\\check if r equals 0 in the interval
\\ infinity points
{
calcInfinityPoints() =
my(P0=normalize(subst(P,y,0)));
l=imag(log(polroots(P0)))/2/Pi;
my(result=listcreate());
for(i=1,#l,
  if(l[i]>=0 && l[i]<=1/2,listput(result,l[i]));
);
listsort(result);
return(result);
}


{
calcmahlergeneral() =
my(degree=poldegree(P,y));
my(mahler=0);
my(allintervallist=splitallinterval());
\\print(allintervallist);
my(ulist=unionlist(allintervallist));
\\print(ulist);
for(i=1,#ulist,
  \\integrate if there is no singularity
  \\print(i);
  \\print(ulist[i]);
  mahler=mahler+intnum(x=ulist[i][1],ulist[i][2],integrandgeneral(x,allintervallist));
  \\integrate if there is a singularity, no need now
);
return(2*mahler);
}

{
checkreciprocal(P) =
my(P1=subst(subst(P,x,1/x),y,1/y)*x^poldegree(P,x)*y^poldegree(P,y));
if(P==P1 || P==-P1,
  A=polcoeff(P,2,y);
  B=polcoeff(P,1,y);
  C=polcoeff(P,0,y);
  return(1);
);
return(0);
}
