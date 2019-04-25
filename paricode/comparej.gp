\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ compare the j-invariant of genus 1 factors
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


read("jGenus1.gp");
read("lfunmahlerg2.gp");
read("lfunmahlerg3.gp");


\\ j-invariant of elliptic curve factor of genus2 curve
{
\\ result [[j,defining equation],[j,defining equation]...]
getj(genus) =
my(l=listcreate());
if(genus==1, for(i=1,#jGenus1List,listput(l,[jGenus1List[i],polylist[i]+k*x*y])));
if(genus==2,
  for(i=1,9,
    my(jlist=calcjInvariantG2(i));
    for(j=1,#jlist,
      listput(l,jlist[j]);
    );
  );
);
if(genus==3,
  for(i=1,7,
    my(jlist=calcjInvariantG3(i));
    for(j=1,#jlist,
      listput(l,jlist[j]);
    );
  );
);
return(l);
}




{
comparej(genus1,genus2) =
my(jlist1=getj(genus1));
if(genus1==genus2,
  for(i=1,#jlist1,
    for(j=1,#jlist1,
      my(ji=subst(jlist1[i][1],k,k));
      my(jj=subst(jlist1[j][1],k,k));
      degni=poldegree(numerator(ji));
      degnj=poldegree(numerator(jj));
      degdi=poldegree(denominator(ji));
      degdj=poldegree(denominator(jj));
      if(i<j && degni==degnj && degdi==degdj,
        compfunc=compareRationalFunction(ji,jj);
        if(compfunc!=0,print(i);print(j);print(compfunc););
      );
    );
  );
);
if(genus1!=genus2,
  my(jlist2=getj(genus2));
  for(i=1,#jlist1,
    for(j=1,#jlist2,
      my(ji=subst(jlist1[i][1],k,k));
      my(jj=subst(jlist2[j][1],k,k));
      degni=poldegree(numerator(ji));
      degnj=poldegree(numerator(jj));
      degdi=poldegree(denominator(ji));
      degdj=poldegree(denominator(jj));
      if(degni==degnj && degdi==degdj,
        compfunc=compareRationalFunction(ji,jj);
        if(compfunc!=0,print(i);print(j);print(compfunc););
      );
    );
  );
);
}

{
getallj() =
my(list1=getj(1));
my(list2=getj(2));
my(list3=getj(3));
return(joinlist(joinlist(list1,list2),list3));
}

{
compareallj() =
my(jlist1=getallj());
my(equalnumberlist=listcreate());
my(resultlist=listcreate());
for(i=1,#jlist1,
  \\ number in all the equallists
  if(inList(equalnumberlist, i),next;);
  my(equallist=listcreate());
  my(relationlist=listcreate());
  for(j=1,#jlist1,
    my(ji=subst(jlist1[i][1],k,k));
    my(jj=subst(jlist1[j][1],k,k));
    degni=poldegree(numerator(ji));
    degnj=poldegree(numerator(jj));
    degdi=poldegree(denominator(ji));
    degdj=poldegree(denominator(jj));
    if(i<j && degni==degnj && degdi==degdj,
      compfunc=compareRationalFunction(ji,jj);
      if(compfunc!=0,
        listput(relationlist,compfunc);
        if(#equallist==0,listput(equallist,i));
        listput(equallist,j);
        print(compfunc);
      );
    );
  );
  if(#equallist>1,listput(resultlist,[relationlist,equallist]));
  equalnumberlist=joinlist(equalnumberlist,equallist);
  if(#equallist!=0,print(equallist));
);
return(resultlist);
}



{
compareRationalFunction(a,b) =
na=subst(numerator(a),k,k1);
da=subst(denominator(a),k,k1);
nb=subst(numerator(b),k,k2);
db=subst(denominator(b),k,k2);
my(c=na*db-nb*da);
my(c0=polcoeff(polcoeff(c,0),0));
my(d=[-30..30]);
if(c0!=0,d=divisors(c0));
for(i=1,#d,
  linear=k1-k2-d[i];
  if(c%(linear)==0,return(linear);break);
  linear=k1-k2+d[i];
  if(c%(linear)==0,return(linear);break);
  linear=k1+k2-d[i];
  if(c%(linear)==0,return(linear);break);
  linear=k1+k2+d[i];
  if(c%(linear)==0,return(linear);break);
);
return(0);
}



{
compareRationalFunctiontest(a,b) =
na=subst(numerator(a),k,k1);
da=subst(denominator(a),k,k1);
nb=subst(numerator(b),k,k2);
db=subst(denominator(b),k,k2);
my(c=na*db-nb*da);
my(d=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]);
my(linearlist=listcreate());
for(i=1,#d,
  linear=k1-k2-d[i];
  if(c%(linear)==0,listput(linearlist,linear));
  if(d[i]!=0,
    linear=k1-k2+d[i];
    if(c%(linear)==0,listput(linearlist,linear));
  );
  linear=k1+k2-d[i];
  if(c%(linear)==0,listput(linearlist,linear));
  if(d[i]!=0,
    linear=k1+k2+d[i];
    if(c%(linear)==0,listput(linearlist,linear));
  );
);
return(linearlist);
}


{
comparealljtest() =
my(jlist1=getallj());
my(equalnumberlist=listcreate());
my(resultlist=listcreate());
for(i=1,#jlist1,
  \\ number in all the equallists
  if(inList(equalnumberlist, i),next;);
  my(equallist=listcreate());
  my(relationlist=listcreate());
  for(j=i+1,#jlist1,
    my(ji=subst(jlist1[i][1],k,k));
    my(jj=subst(jlist1[j][1],k,k));
    degni=poldegree(numerator(ji));
    degnj=poldegree(numerator(jj));
    degdi=poldegree(denominator(ji));
    degdj=poldegree(denominator(jj));
    if(degni==degnj && degdi==degdj,
      compfunc=compareRationalFunctiontest(ji,jj);
      if(#compfunc>=2,
        listput(relationlist,compfunc);
        if(#equallist==0,listput(equallist,i));
        listput(equallist,j);
        print(i,"     ",j,"      ",compfunc);
      );
    );
  );
  if(#equallist>1,listput(resultlist,[relationlist,equallist]));
  equalnumberlist=joinlist(equalnumberlist,equallist);
);
return(resultlist);
}
