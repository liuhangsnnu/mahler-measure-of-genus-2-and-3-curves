\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ calculate the Mahler measure of reciprocal polynomial A*y^2+B*y+C
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

defaultprecision=default(realprecision);
\\ rootprecision is used to determine if the root is on the unit circle
rootprecision=10^(-ceil(defaultprecision/2));


A;B;C;
\\ tranform the functions a(x),b(x),c(x) to trigonometric polynomial
a(x)= transpol(A,x);
b(x)=transpol(B,x);
c(x)=transpol(C,x);
r(x)=(-b(x)+sqrt(b(x)^2-4*a(x)*c(x)))/(2*a(x));
n1(x)=log(abs((-b(x)-sqrt(b(x)^2-4*a(x)*c(x)))));
n2(x)=log(abs((-b(x)+sqrt(b(x)^2-4*a(x)*c(x)))));
\\ the integrand
integrand(x) =  log(abs(r(x)));
integranda(x) = log(2*abs(a(x)));
\\ use integration by parts to transform the function integranda
transinta(x,gamm,mult)=if(x!=gamm,(gamm-x)*dertranspol(A,x)/transpol(A,x),-mult);



\\ this function is used if the interval does not contain log singularity
mahlerint(x0,x1) = abs(intnum(x=x0,x1,integrand(x)));


{
\\ deal with log sigularity
\\ gamm is the singularity, i.e. root of a(x), mult is the multiplicity of gamm, there is only one singularity in [alph, beta]
\\ this function is used if the interval contains a log singularity
lsmahlerint(alph,beta,gamm,mult) =
\\ integration of integranda(x), integration by parts
my(inta,intn,intpart);
intpart=if(abs(gamm-(beta+alph)/2)>rootprecision,
  intnum(x=alph,beta,transinta(x,gamm,mult)),
  intnum(x=alph,(beta+gamm)/2,transinta(x,gamm,mult))+intnum(x=(beta+gamm)/2,beta,transinta(x,gamm,mult))
);
inta=(beta-gamm)*integranda(beta)+(gamm-alph)*integranda(alph)+intpart;
if(abs(gamm-(beta+alph)/2)>rootprecision,
  if(b(gamm)>0,
    intn=intnum(x=alph,beta,n1(x)),
    intn=intnum(x=alph,beta,n2(x))
  ),
  if(b(gamm)>0,
    intn=intnum(x=alph,(beta+alph)/2,n1(x))+intnum(x=(beta+alph)/2,beta,n1(x)),
    intn=intnum(x=alph,(beta+alph)/2,n2(x))+intnum(x=(beta+alph)/2,beta,n2(x))
  );
);
return(abs(intn-inta));
}

\\ elliminate the factor of x
normalize(pol) = pol/x^(valuation(pol,x));


{
\\ transform the reciprocal and anti-reciprocal polynomials to trigonometric polynomials
transpol(pol,x) =
if(pol==0,return(0));
pol=normalize(pol);
my(d=poldegree(pol));
my(tr=0);
\\ reciprocal
if(polcoeff(pol,d)==polcoeff(pol,0),
  for(i=0,ceil(d/2)-1,tr=tr+2*polcoeff(pol,i)*cos((d/2-i)*2*Pi*x));
  \\ add the middle term for d even
  if(d%2==0,tr=tr+polcoeff(pol,d/2));
);
\\ anti-reciprocal
if(polcoeff(pol,d)==-polcoeff(pol,0),
  for(i=0,ceil(d/2)-1,tr=tr-2*I*polcoeff(pol,i)*sin((d/2-i)*2*Pi*x));
);
return(tr);
}

{
\\ derivative of transpol(pol,x), consider reciprocal and anti-reciprocal polynomials
dertranspol(pol,x)=
my(d=poldegree(pol));
my(tr=0);
\\ reciprocal
if(polcoeff(pol,d)==polcoeff(pol,0),
  for(i=0,ceil(d/2)-1,tr=tr-2*(d/2-i)*2*Pi*polcoeff(pol,i)*sin((d/2-i)*2*Pi*x));
);
\\ anti-reciprocal
if(polcoeff(pol,d)==-polcoeff(pol,0),
  for(i=0,ceil(d/2)-1,tr=tr-2*I*(d/2-i)*2*Pi*polcoeff(pol,i)*cos((d/2-i)*2*Pi*x));
);
return(tr);
}

{
\\ degree of first nonzero term of a polynomial
polfirstdeg(pol)=
for(i=0,poldegree(pol),if(polcoeff(pol,i)!=0,return(i)));
}

{
\\ find the interval where abs(y) != 1 && 0<arg(x)<Pi and split the interval such that only the end points have log and square root singularity
calcinterval(delta) =
my(rod=polroots(delta));
\\ roots on unit circle
my(rouc=listcreate(length(rod)));
my(j=1);
\\ find all the roots on the unit circle exp(2*Pi*I*x)
for(i=1,length(rod),if(abs(abs(rod[i])-1)<rootprecision,listput(rouc,imag(log(rod[i]))/2/Pi,j);j=j+1));
listsort(rouc);
\\ x is in [0,1/2]
my(l=length(rouc));
for(i=1,l,if(rouc[l+1-i]<0,listpop(rouc,l+1-i)));
l=length(rouc);
if(l==0,listput(rouc,0,1);listput(rouc,1/2,2);l=2);
\\ the endpoints are 0 and 1/2
if(rouc[l]!=1,listput(rouc,1/2,l+1));
if(rouc[1]!=0,listinsert(rouc,0,1));
\\ all the intervals such that abs(y) != 1
l=length(rouc);
my(interval=listcreate(l));
for(i=1,l-1,
  my(x=exp(2*Pi*I*(rouc[i]+rouc[i+1])/2));
  if(real(eval(delta)/x^((poldegree(delta)+polfirstdeg(delta))/2))>0 && rouc[i] != rouc[i+1],listput(interval,[rouc[i],rouc[i+1]]));
);
\\ multiple roots, square root singularity
my(mr=multipleroots(rouc));
\\print("mr = " mr);
\\ roots of A, log singularity
my(ra=tolist(imag(log(polroots(A))/2/Pi)));
l=length(ra);
\\ we only need the roots with argument from 0 to Pi
for(i=1,l,if(ra[l+1-i]<0,listpop(ra,l+1-i)));
\\ split the interval if the interval contains more than one log sigularity, i.e. a(x)=0
return(splitintervallist(interval,ra));
}

{
\\ interval is a list, split the interval list
\\ futhermore split the interval when both the endpoints have log singularity
splitintervallist(interval, points) =
my(l=listcreate());
for(i=1,length(interval),l=joinlist(l, splitinterval(tolist(interval[i]), points)));
return(l);
}

{
\\ interval is an interval [a,b], split the interval when it contains more than one root of a(x)
splitinterval(interval, points) =
listsort(points);
\\ points in the interval
pii=listcreate();
if(length(points)!=0 && (interval[1]<points[1]) && (points[1]<interval[2]), listput(pii,points[1]));
for(i=2,length(points),if((interval[1]<points[i]) && (points[i]<interval[2]) && (points[i] != points[i-1]), listput(pii,points[i])));
my(l=listcreate());
if(length(pii)>1,
  listput(l,[interval[1],(pii[1]+pii[2])/2]);
  for(i=1,length(pii)-2,listput(l,[(pii[i]+pii[i+1])/2,(pii[i+1]+pii[i+2])/2]));
  listput(l,[(pii[length(pii)-1]+pii[length(pii)])/2,interval[2]]),
  listput(l,interval)
);
return(l);
}

{
tolist(l) =
l1=listcreate(length(l));
for(i=1,length(l),listput(l1,l[i]));
return(l1);
}


{
joinlist(l1,l2)=
my(l=listcreate());
for(i=1,length(l1),listput(l,l1[i]));
for(i=1,length(l2),listput(l,l2[i]));
return(l);
}

\\ check if the singlar points are in the interval
{
checksingularpoint(rouc, interval) =
\\ find multiple roots, i.e. singular points
my(mr=multipleroots(rouc));
checkintervalsing(mr,interval);
}


{
\\ check if the interval contains a sigular point
checkintervalsing(mr,interval) =
if(length(interval)==0 || length(mr)==0,return(0));
for(i=1,length(mr),
for(j=1,length(interval),
  if(mr[i]>=interval[j][1] && interval[j][2]>=mr[i],return(1));
);
);
}

{
\\ find multiple roots, rouc is sorted
multipleroots(rouc) =
my(mr=listcreate(length(rouc)));
for(i=1,length(rouc)-1,
  if(length(mr)==0,if(rouc[i]==rouc[i+1],listput(mr,rouc[i])),if(rouc[i]==rouc[i+1] && rouc[i]!=mr[length(mr)],listput(mr,rouc[i])));
);
return(mr);
}

{
\\ check if the interval contain a root of a(x) and the multiplicity of the root
\\ ro is 0 if there is no root in the interval, mult is the multiplicity and root is the root
checkintervala(interval) =
my(ra=tolist(imag(log(polroots(A))/2/Pi)));
\\ ro indicates if the interval contains a root of a(x), mult is the multiplicity of the number, root is the root in the interval
my(ro=0);
my(mult=0);
my(root=0);
for(i=1,length(ra),
  if((ra[i]-interval[1])>rootprecision && (ra[i]-interval[2])<-rootprecision,
  ro=1;mult=mult+1;root=ra[i])
);
result=[ro,mult,root];
return(result);
}


{
\\ calculate the elliptic factor Z_k/<sigma>, sigma:(x,y)->(1/x,1/y)
\\ Z_k is the curve defined by A*y^2+B*y+C
calcEllipticFactor(delta) =
if(delta%x^2==0,delta=delta/x^2);
my(flag=0);
if(delta%x^2==0,delta=delta/x^2);
\\ flag=-1 for genus 2, we need to use y^2 = x^3*h(1/x) instead
if(delta%(x-1)^2==0,flag=-1;delta=delta/(x-1)^2);
if(delta%(x+1)^2==0,delta=delta/(x+1)^2);
\\ when delta%(x+1)^4==0
if(delta%(x+1)^2==0,delta=delta/(x+1)^2);
if(delta%(x-1)^2==0,delta=delta/(x-1)^2);
my(deg=poldegree(delta)); \\ 5,6,7 or 8
\\ degenerating case, degree=5 is possible
if(deg<5, return([0,0,0,0,0]));
i=subst(subst(y^2-delta,x,(x+1)/(x-1)),y,y/(x-1)^(ceil(deg/2)))*(x-1)^(2*ceil(deg/2));
\\print("delta=",delta);
\\print("i=",i);
i=substpol(i,x^2,x);
D=i-y^2;
\\genus 2
if(poldegree(D)==3,
  c0=polcoeff(-D,0);
  c1=polcoeff(-D,1);
  c2=polcoeff(-D,2);
  c3=polcoeff(-D,3);
  if(flag==-1,c0=polcoeff(-D,3);c1=polcoeff(-D,2);c2=polcoeff(-D,1);c3=polcoeff(-D,0));
  \\print("c0 = ",c0);
  \\print("c1 = ",c1);
  \\print("c2 = ",c2);
  \\print("c3 = ",c3);
  \\ eliminate
  a6=c0*c3^2;
  a4=c1*c3;
  a2=c2;
  my(v=[0,a2,0,a4,a6]);
  return(v);
);
\\ genus 3
if(poldegree(D)==4,
  c0=polcoeff(-D,0);
  c1=polcoeff(-D,1);
  c2=polcoeff(-D,2);
  c3=polcoeff(-D,3);
  c4=polcoeff(-D,4);
  \\ Jacobian of genus 1 curve, see the paper of Boyd
  a2=c2;
  a4=c1*c3-4*c0*c4;
  a6=-(4*c0*c2*c4-c1^2*c4-c0*c3^2);
  my(v=[0,a2,0,a4,a6]);
  return(v);
);
return([0,0,0,0,0]);
}



{
\\ calculate the elliptic factor Z_k/<sigma2>
\\ Z_k is the curve defined by A*y^2+B*y+C
calcEllipticFactor2(delta) =
if(delta%x^2==0,delta=delta/x^2);
my(flag=0);
\\ flag=-1 for genus 2, we need to use y^2 = x^3*h(1/x) instead
if(delta%(x-1)^2==0,flag=-1;delta=delta/(x-1)^2);
if(delta%(x+1)^2==0,delta=delta/(x+1)^2);
my(deg=poldegree(delta)); \\ 5,6,7 or 8
\\ degenerating case, degree=5 is possible
if(deg<5, return([0,0,0,0,0]));
i=subst(subst(y^2-delta,x,(x+1)/(x-1)),y,y/(x-1)^(ceil(deg/2)))*(x-1)^(2*ceil(deg/2));
\\print("delta=",delta);
\\print("i=",i);
i=substpol(i,x^2,x);
D=i-y^2;
\\genus 2
if(poldegree(D)==3,
  c0=polcoeff(-D,0);
  c1=polcoeff(-D,1);
  c2=polcoeff(-D,2);
  c3=polcoeff(-D,3);
  if(flag==-1,c0=polcoeff(-D,3);c1=polcoeff(-D,2);c2=polcoeff(-D,1);c3=polcoeff(-D,0));
  \\print("c0 = ",c0);
  \\print("c1 = ",c1);
  \\print("c2 = ",c2);
  \\print("c3 = ",c3);
  \\ eliminate
  a2=c1;
  a4=c2*c0;
  a6=c3*c0^2;
  my(v=[0,a2,0,a4,a6]);
  return(v);
);
}




{
\\modify the interval if it contains 0 and A%(x-1)==0 or it contains 1/2 and A%(x+1)==0
modifyinterval(interval) =
if(length(interval)!=0,
  my(modifyinterval=listcreate());
  \\a(0)=0
  if(interval[1][1]==0 && A%(x-1)==0,
    ci=checkintervala(interval[1]);
    \\there is no singularity
    if(ci[1]==0,
      listput(modifyinterval,[-interval[1][2],interval[1][2]]);
    );
    \\there is a singularity
    if(ci[1]!=0,
      point=ci[3];
      listput(modifyinterval,[-point/2,point/2]);
      listput(modifyinterval,[point/2,interval[1][2]]);
    );
    for(i=2,length(interval),
      listput(modifyinterval,interval[i]);
    );
    return(modifyinterval);
  );
  \\a(1/2)=0
  len=length(interval);
  if(interval[len][2]==1/2 && A%(x+1)==0,
    for(i=1,len-1,
      listput(modifyinterval,interval[i]);
    );
    ci=checkintervala(interval[len]);
    \\there is no singularity
    if(ci[1]==0,
      listput(modifyinterval,[interval[len][1],1-interval[len][1]]);
    );
    \\there is a singularity
    if(ci[1]!=0,
      point=ci[3];
      listput(modifyinterval,[interval[len][1],(point+1/2)/2]);
      listput(modifyinterval,[(point+1/2)/2,1-(point+1/2)/2]);
    );
    return(modifyinterval);
  );
);
return(interval);
}


{
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ this function calculate the Mahler measure of A*y^2+B*y+C
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
calcmahler() =
my(delta=B^2-4*A*C);
my(interval = calcinterval(delta));
interval = modifyinterval(interval);
\\print(interval);
\\ check whether the singular point is in the interval, todo
my(mm=0);
for(i=1,length(interval),
  ci=checkintervala(interval[i]);
  \\ if the interval is in [0,1/2], we need to multiply the result by 2
  m=2;
  \\ if the interval is not in [0,1/2], m=1
  if((i==1&&interval[i][1]<0)||(i==length(interval)&&interval[i][2]>1/2),m=1);
  \\ if the interval does not contain a root of a(x), use mahlerint, otherwise use lsmahlerint
  if(ci[1]==0,mm=mm+m*mahlerint(interval[i][1],interval[i][2]),mm=mm+m*lsmahlerint(interval[i][1],interval[i][2],ci[3],ci[2]));
);
return(mm);
}


{
\\ compare l-value and mahler measure and write the results to tex file
calclfunmahler(k,cl) =
my(d=gcd(A,B));
if(d!=1,A=A/d;B=B/d;C=C/d);
\\ calculate the Mahler measure
mm=0;
if(!(poldegree(A)==0 && poldegree(B)==0 && poldegree(C)==0),mm=calcmahler());
print("k = ", k, ", mahler measure is ", mm);
\\ get A,B,C back
if(d!=1,A=A*d;B=B*d;C=C*d);
\\ calculate the elliptic curve factor in the Jacobian of the curve
a2=eval(cl[2]);
a4=eval(cl[4]);
a6=eval(cl[5]);
disc = eval(18*a2*a4*a6 - 4*a2^3*a6 + a2^2*a4^2 - 4*a4^3 - 27*a6^2);
\\ if discriminant != 0
if(disc !=0,
  my(E = ellinit(eval(cl)));
  \\ calculate the conductor and minimal Weierstrass equation
  egr=ellglobalred(E);
  N=egr[1];
  ecc=ellchangecurve(E,egr[2]);
  \\ calculate the L-value of corresponding elliptic curve
  my(s=lfun(E, 0, 1)/mm);
  print("s = ", s, ", rational number close to s = ", bestappr(s));
  if(denominator(bestappr(s))<100,
  \\ one should calculate a1,a2,a3,a4,a6 of minimal model
    write(filepath,Str(k,"&",bestappr(s),"&",N,"&",ecc.a1,"&",ecc.a2,"&",ecc.a3,"&",ecc.a4,"&",ecc.a6,"\\\\"));
  );
);
\\degenerating cases
delta=B^2-4*A*C;
if(disc==0 && delta!=0,f=factor(delta);
  co = polcoeff(delta,poldegree(delta))/calcC(f);
  dlvalueList=listcreate();
  dlist=listcreate();
  listput(dlist,-3);listput(dlist,-4);listput(dlist,-24);
  my(d);
  listput(dlvalueList,mm);
  listput(dlvalueList,lfund(-3));
  listput(dlvalueList,lfund(-4));
  listput(dlvalueList,lfund(-24));
  for(i=1,length(f[,1]),
    \\ find discriminant of imaginary quadratic field where the boundary of Deninger cycle lies in
    if(f[i,2]==2 && poldegree(f[i,1])==1,dd=co*calcD(f,i);d=coredisc(dd*denominator(dd)^2));
    if(poldegree(f[i,1])==2,d=coredisc(polcoeff(f[i,1],1)^2-4*polcoeff(f[i,1],0)*polcoeff(f[i,1],2)));
    \\if(f[i,2]==2 && poldegree(f[i,1])==2,d=coredisc(polcoeff(f[i,1],1)^2-4*polcoeff(f[i,1],0)*polcoeff(f[i,1],2)));
    if(!inList(dlist,d),
      if(d<0,listput(dlist,d);listput(dlvalueList,lfund(d)));
    );
  );
  my(knowdlist=[-3,-4,-7,-8,-15,-19,-24,-39,-40,-55,-120]);
  rv=lindep(vector(length(dlvalueList), i, dlvalueList[i]));
  \\ length of rv except 0
  lrv=0;
  \\ discriminant
  finald=0;
  for(i=1,#rv,
    if(rv[i]!=0,
      lrv=lrv+1;
      if(i!=1,finald=dlist[i-1]);
    );
  );
  \\ Test for new example of Chinburg's conjecture
  if(lrv==2 && inList(knowdlist, finald)==0,
    print(A,"    ",B,"    ",C);
    print(dlvalueList);
    print(dlist);
    print(rv);
  );

  coeffm=rv[1];
  relationStr="$m=";
  my(signrelation="");
  if(coeffm!=0 && abs(coeffm)<1000,
    my(num=0);
    my(flag=0);
    for(i=2,#rv,
      if(rv[i]!=0,
        if(-rv[i]/coeffm>0 && flag,signrelation="+");
        relationStr=Str(relationStr,signrelation,-rv[i]/coeffm,"d_{",abs(dlist[i-1]),"}");
        num=num+1;flag=1;
      );
    );
    if(num!=0,
      relationStr=Str(relationStr,"$"),
      relationStr=Str(relationStr,0,"$")
    );
    my(curveEqu="");
    if(num!=0,
      write(filepath,k,"&",relationStr,"&","&","\\multicolumn{5}{c|}{",curveEqu,"}","\\\\");
    );
  );
);
}



\\ some functions

{
\\ calculate d_{-f}
lfund(f)=(-f)^(1.5)*3/2/Pi^3*lfun(x^2-f, 2);
}

{
\\ integration of functions with vertical tangent
intvtfunction(a,b,f) =
my(u=a+t^2);
my(v=b-t^2);
my(c=(a+b)/2);
left=intnum(t=0,sqrt(c-a),f(a+t^2)*2*t);
right=intnum(t=-sqrt(b-c),0,-f(b-t^2)*2*t);
\\print(left,"   ",right);
return(left+right);
}

{
  inList(list, value)=for(i=1,#list, if(list[i]==value, return(i))); 0;
}



{
calcD(f,index) =
  result=1;
  p=-polcoeff(f[index,1],0)/polcoeff(f[index,1],1);
  for(i=1,length(f[,1]),if(index!=i,result=result*subst(f[i,1],x,p)^f[i,2]));
  return(result);
}

{
  \\ product of leading coefficient
  calcC(f) =
  result=1;
  for(i=1,length(f[,1]),result=result*polcoeff(f[i,1],poldegree(f[i,1]))^f[i,2]);
  return(result);
}


{
\\ polynomial A to tex
pol2Str(pol) =
my(f=factor(pol));
my(polStr="");
for(i=1,#f[,1],
  my(lp="(");
  my(rp=")");
  if(f[i,1]==x,lp="";rp="");
  if(f[i,2]!=1,
    polStr=Str(polStr,Str(lp,Strtex(f[i,1]),rp,"^{",f[i,2],"}")),
    polStr=Str(polStr,Str(lp,Strtex(f[i,1]),rp))
  );
);
return(polStr);
}

{
\\ polynomial A*y^2+B*y+C to tex
pol2Str2Var(P) =
if(checkreciprocal(P)==1,
  my(A=polcoeff(P,2,y));
  my(B=polcoeff(P,1,y));
  my(C=polcoeff(P,0,y));
  return(Str(pol2Str(A),"(y^2+",Strtex(C/A),")+(",Strtex(B),")y"));
);
my(coef,coefstr,d);
my(str="");
for(i=0,poldegree(P,y),
  coef=polcoeff(P,poldegree(P,y)-i,y);
  if(coef!=0,
      if(coef==x^poldegree(coef,x) && coef!=1,
          coefstr=Str(Strtex(coef)),
          if(coef==1,coefstr="",coefstr=Str("(",Strtex(coef),")"));
      );
      str=Str(str,coefstr);
      if(i!=poldegree(P,y),
        d=poldegree(P,y)-i;
        if(d==1,d="");
        str=Str(str,"y^{",d,"}+");
      );
  );
);
return(str);
}
