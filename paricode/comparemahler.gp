\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ compare the mahler measure of different families
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



read("jGenus1.gp");
read("lfunmahlerg2.gp");
read("lfunmahlerg3.gp");
read("data2.gp");
read("data3.gp");
read("mahler.gp");
read("mahlergeneral.gp");


\\calculate the Mahler measure for integers in [-bound,bound]
bound=200;
\\ write results to the file
filepath="comparemahler.txt";

{
getcurve(i) =
lg1=#polylist;
lg21=lg1+#alist1g2;
lg22=lg21+#alist2g2;
lg23=lg22+#alist3g2;
lg24=lg23+#alist4g2;
lg25=lg24+#alist5g2;
lg26=lg25+#alist6g2;
lg27=lg26+#alist7g2;
lg28=lg27+#alist8g2;
lg29=lg28+#alist9g2;
lg31=lg29+#alist1;
lg32=lg31+#alist2;
lg33=lg32+#alist3;
lg34=lg33+#alist4;
lg35=lg34+#alist5;
lg36=lg35+#alist5;
lg37=lg36+#alist7;
\\print(lg1," ",lg21," ",lg22," ",lg23," ",lg24," ",lg25," ",lg26," ",lg27," ",lg28," ",lg29," ",lg31," ",lg32," ",lg33," ",lg34," ",lg35," ",lg36," ",lg37);
if(i<=lg1,return(i));
my(num);
my(index);
my(genus);
\\genus2
if(lg1<i && i<=lg21, num=i-lg1; index=1; genus=2);
if(lg21<i && i<=lg22, num=i-lg21; index=2; genus=2);
if(lg22<i && i<=lg23, num=i-lg22; index=3; genus=2);
if(lg23<i && i<=lg24, num=i-lg23; index=4; genus=2);
if(lg24<i && i<=lg25, num=i-lg24; index=5; genus=2);
if(lg25<i && i<=lg26, num=i-lg25; index=6; genus=2);
if(lg26<i && i<=lg27, num=i-lg26; index=7; genus=2);
if(lg27<i && i<=lg28, num=i-lg27; index=8; genus=2);
if(lg28<i && i<=lg29, num=i-lg28; index=9; genus=2);
\\genus3
if(lg29<i && i<=lg31, num=i-lg29; index=1; genus=3);
if(lg31<i && i<=lg32, num=i-lg31; index=2; genus=3);
if(lg32<i && i<=lg33, num=i-lg32; index=3; genus=3);
if(lg33<i && i<=lg34, num=i-lg33; index=4; genus=3);
if(lg34<i && i<=lg35, num=i-lg34; index=5; genus=3);
if(lg35<i && i<=lg36, num=i-lg35; index=6; genus=3);
if(lg36<i && i<=lg37, num=i-lg36; index=7; genus=3);
return([num,index,genus]);
}

{
comparemahler(relationequalvector) =
relationlist=relationequalvector[1];
equallist=relationequalvector[2];
my(mlist=listcreate());
\\get all the curves
for(i=1,#equallist,
  listput(mlist,calcmahlermeasure(getcurve(equallist[i])));
);
\\print(mlist);
my(a);my(b);
my(relation1);my(relation2);
my(resultlist=listcreate());
for(i=1,#equallist,
  for(j=i+1,#equallist,
    \\print("=========",i,"    ",j,"==========");
    if(i==1,relation1=relationlist[j-1];
      \\k2=a*k1+b
      a=-polcoeff(relation1,1,k1)/polcoeff(relation1,1,k2);
      b=-subst(subst(relation1,k1,0),k2,0)/polcoeff(relation1,1,k2);
    );
    \\calculate the relation
    if(i!=1,
      relation1=relationlist[i-1];
      relation2=relationlist[j-1];
      a1=polcoeff(relation1,1,k1);
      a2=polcoeff(relation2,1,k1);
      b1=polcoeff(relation1,1,k2);
      b2=polcoeff(relation2,1,k2);
      c1=subst(subst(relation1,k1,0),k2,0);
      c2=subst(subst(relation2,k1,0),k2,0);
      \\k2=a*k1+b
      a=(a2*b1)/(a1*b2);
      b=(a2*c1-a1*c2)/(a1*b2);
      \\print(relation1);
      \\print(relation2);
      \\print(a1,"   ",a2,"   ",b1,"   ",b2,"   ",c1,"   ",c2,"   ",a,"   ",b);
    );
    \\print(a,"    ",b);
    my(m1list=mlist[i]);
    my(m2list=mlist[j]);
    my(mahlerrelationlist=listcreate());
    my(mink1=0);
    my(maxk1=0);
    for(k=1,#m1list,
      \\k1 and k2 are parameters in the equation
      my(k1=k-bound-1);
      my(k2=a*k1+b);
      kk=k2+1+bound;
      \\print(k,"   ",kk);
      \\find the range where the quotient is rational
      if(kk>0 && kk<=#m2list,
        \\ first k1
        if(mink1==0,mink1=k1);
        \\ last k1
        maxk1=k1;
        if(m1list[k]!=0 && abs(m2list[kk])>10^(-30),
          relation=m1list[k]/m2list[kk];
          \\print(k1,"    ",k2,"    ",m1list[k],"     ",m2list[kk],"       ",relation);
          my(appr=bestappr(relation,10^30));
          if(denominator(appr)<5,
            listput(mahlerrelationlist,[k1,appr])
          );
        );
      );
    );
    \\print(mink1,"      ",maxk1);
    \\print(splitrelationinterval(mahlerrelationlist,mink1,maxk1));
    if(#mahlerrelationlist!=0,
      listput(resultlist,[[i,j],[a,b],splitrelationinterval(mahlerrelationlist,mink1,maxk1)]);
    );
  );
);
return(resultlist);
}




{
splitrelationinterval(relationlist,mink1,maxk1) =
intervallist = listcreate();
my(begin=1);
my(a,b);
for(i=2,#relationlist,
  if(relationlist[i][1]!=relationlist[i-1][1]+1 || relationlist[i][2]!=relationlist[i-1][2],
    a=relationlist[begin][1];
    b=relationlist[i-1][1];
    if(a==mink1,a=-oo);
    if(b==maxk1,b=+oo);
    listput(intervallist,[a,b,relationlist[begin][2]]);
    begin=i;
  );
);
if(#relationlist!=0,
  a=relationlist[begin][1];
  b=relationlist[#relationlist][1];
  if(a==mink1,a=-oo);
  if(b==maxk1,b=+oo);
  listput(intervallist,[a,b,relationlist[begin][2]]);
);
return(intervallist);
}


{
toTexCompare(relationequalvector) =
Y;X;
k=varlower("k");
my(comparelist=comparemahler(relationequalvector));
if(#comparelist==0,return(false));
if(checkcomparelistgenus1(comparelist,relationequalvector)==true,return(false));
my(pollist=getPolList(relationequalvector[2]));
for(i=1,#pollist,
  if(checkpolrelation(#pollist+1-i,comparelist)==true,
    write(filepath,Str("\\(\\displaystyle ",getPolSymbol(i),"_k=",pol2Str2Var(pollist[#pollist+1-i]),"\\)\\\\"));
  );
);
for(i=1,#comparelist,
  my(compare=comparelist[#comparelist+1-i]);
  symbol1=getPolSymbol(#pollist+1-compare[1][1]);
  symbol2=getPolSymbol(#pollist+1-compare[1][2]);
  subscript=compare[2][1]*k+compare[2][2];
  relationstr=getRelationString(compare[3],symbol1,symbol2,subscript);
  write(filepath,relationstr);
);
return(true);
}

{
\\ check if mahler measure of a family has relation with others
checkpolrelation(index,comparelist)=
my(flag=false);
for(i=1,#comparelist,
  if((index==comparelist[i][1][1])||(index==comparelist[i][1][2])
  && (checkallpointinterval(comparelist[i][3])==false),
    flag=true;
  );
);
return(flag);
}

{
checkallpointinterval(compare)=
my(flag=true);
for(i=1,#compare,
  if(compare[i][1]!=compare[i][2],
    flag=false;
  );
);
return(flag);
}


{
\\ check if all the relations are between genus 1 families
checkcomparelistgenus1(comparelist,relationequalvector)=
my(flag=true);
for(i=1,#comparelist,
  if(relationequalvector[2][comparelist[i][1][1]]>#polylist ||
      relationequalvector[2][comparelist[i][1][2]]>#polylist,
      flag=false
  );
);
return(flag);
}

{
\\ check if all the polynomials define curves of genus 1
checkgenus1(relationequalvector) =
my(flag=true);
for(i=1,#relationequalvector[2],
  if(relationequalvector[2][i]>#polylist,flag=false);
);
return(flag);
}

{
\\check if there is a genus 3 family
checkgenus3(relationequalvector) =
my(flag=false);
for(i=1,#relationequalvector[2],
  if(relationequalvector[2][i]>#polylist+#alist1g2+#alist2g2+#alist3g2+#alist4g2+#alist5g2+#alist6g2+#alist7g2,
    flag=true;
  );
);
return(flag);
}


{
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ this function writes the all the relation between
\\ the Mahler measure of different families to a tex file
\\ the input comparealljlist is the result of compareallj()
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
toTexAll(comparealljlist) =
write(filepath,"\\begin{longtable}{|l|}");
write(filepath,"\\hline");
for(i=1,#comparealljlist,
  print("==========================================");
  print(i);
  \\ genus 2 and genus 3
  \\if(checkgenus1(comparealljlist[i])==false,
  \\ genus 3 case
  if(checkgenus1(comparealljlist[i])==false && checkgenus3(comparealljlist[i])==true,
  \\ genus 2 case
  \\if(checkgenus1(comparealljlist[i])==false && checkgenus3(comparealljlist[i])==false,
      result=toTexCompare(comparealljlist[i]);
      if(result!=false,
        write(filepath,"\\hline");
      );
  );
);
write(filepath,"\\end{longtable}");
}

{
getRelationString(intervallist,symbol1,symbol2,subscript) =
my(resultstr="\\(\\displaystyle ");
for(i=1,#intervallist,
  if(intervallist[i][1]!=intervallist[i][2],
    my(ratio=intervallist[i][3]);
    if(ratio==1,ratio="");
    my(relationstr=Str("m(",symbol1,"_k) = ", ratio, "m(", symbol2, "_{",Strtex(subscript),"}),"));
    my(interval=Str(Strtex(intervallist[i][1])," \\leqslant k \\leqslant ",Strtex(intervallist[i][2])));
    resultstr=Str(resultstr,Str(relationstr,interval));
    if(i!=#intervallist,resultstr=Str(resultstr,",\\quad "));
  );
);
resultstr=Str(resultstr,"\\)\\\\");
\\print(resultstr);
return(resultstr);
}

{
getPolSymbol(i) =
if(i==1, return("P"));
if(i==2, return("Q"));
if(i==3, return("R"));
if(i==4, return("S"));
if(i==5, return("T"));
if(i==6, return("U"));
if(i==7, return("V"));
if(i==8, return("W"));
if(i==9, return("F"));
if(i==10, return("G"));
if(i==11, return("H"));
if(i==12, return("I"));
if(i==13, return("J"));
if(i==14, return("L"));
if(i==15, return("M"));
if(i==16, return("N"));
}

{
getPolList(equallist) =
my(resultlist=listcreate());
for(i=1,#equallist,
  P=getPol(getcurve(equallist[i]));
  \\listput(resultlist,subst(subst(P,y,Y),x,X));
  listput(resultlist,P);
);
return(resultlist);
}

{
getPol(curve) =
my(pol);
if(#curve==1,pol=polylist[curve]+k*x*y);
if(#curve!=1, genus=curve[3];
  my(num=curve[1]);
  my(index=curve[2]);
  my(A,B,C);
  if(genus==2,
    A=getAByindexG2(num,index);
    C=getCByindexG2(A,index);
    B=getBG2(num,index);
  );
  if(genus==3,
    A=getAByindexG3(num,index);
    C=getCByindexG3(A,index);
    B=getBG3(num,index);
  );
  pol=A*y^2+B*y+C;
);
return(pol);
}


{
calcmahlermeasure(curve) =
m=-bound;n=bound;
my(l=listcreate());
\\genus1
if(#curve==1,
  my(pol=polylist[curve]);
  for(k=m,n,
    P=pol+k*x*y;
    \\print(P);
    \\plotintegrand();
    my(mm=0);
    \\if(numfactor(P)==1,
    \\print(reci);
      reci=checkreciprocal(P);
      if(reci!=0,
        \\print(A,"    ",B);
        my(d=gcd(A,B));
        if(d!=1,A=A/d;B=B/d;C=C/d);
        if(!(poldegree(A)==0 && poldegree(B)==0 && poldegree(C)==0),
          \\print(A,"    ",B);
          mm=calcmahler();
          if(d!=1,A=A*d;B=B*d;C=C*d);
        );
      );
      if(reci==0,mm=calcmahlergeneral());
    \\);
    listput(l,mm);
  );
);
if(#curve!=1, genus=curve[3];
  my(num=curve[1]);
  my(index=curve[2]);
  if(genus==2,
    A=getAByindexG2(num,index);
    C=getCByindexG2(A,index);
    \\print(A,"   ",getBG2(num,index));
    for(k=m,n,
      B=getBG2ByindexG2(num,k,index);
      \\print(A*y^2+B*y+C);
      \\print(A,"   ",B,"   ",C);
      my(d=gcd(A,B));
      if(d!=1,A=A/d;B=B/d;C=C/d);
      my(mm=0);
      if(!(poldegree(A)==0 && poldegree(B)==0 && poldegree(C)==0),
        \\print(A,"    ",B);
        mm=calcmahler();
        \\print(k,"       ",mm);
      );
      if(d!=1,A=A*d;B=B*d;C=C*d);
      listput(l,mm);
    );
  );
  if(genus==3,
    A=getAByindexG3(num,index);
    C=getCByindexG3(A,index);
    \\print(A,"   ",getBG3(num,index));
    for(k=m,n,
      B=getBG3ByindexG3(num,k,index);
      my(d=gcd(A,B));
      if(d!=1,A=A/d;B=B/d;C=C/d);
      my(mm=0);
      if(!(poldegree(A)==0 && poldegree(B)==0 && poldegree(C)==0),
        mm=calcmahler();
      );
      if(d!=1,A=A*d;B=B*d;C=C*d);
      listput(l,mm);
    );
  );
);
return(l);
}
