
#include "diag.h"


    
void transpose(int n, ccomplex *a)
{
  int i, j;
  ccomplex zt, *aii, *aij, *aji;

  for(i=0,aii=a;i<n-1;i++,aii+=n+1)
    for(j=i+1,aij=aii+1,aji=aii+n;j<n;j++,aij++,aji+=n){
      zt   = *aij;
      *aij = *aji;
      *aji = zt;
    }
}

void diag(int n, ccomplex *a, double *e)
{
  int  i,j,k,l,m;
  double c,d,f,g,h,r,s,t,u,v,w,x,y,*ei,*ej,*ek,*el;
  ccomplex zs,zt,zu,*akkp1,*akj,*aki,*aji,*aij,*b,*bi,*bj,*bk,*bl;

  b = malloc(n * sizeof(ccomplex));

  for(k=0,akkp1=a+1,bk=b,ek=e;k<n-2;k++,akkp1+=n+1,bk++,ek++){
    (*ek) = (akkp1-1)->re;
    for(j=k+1,s=0.0,akj=akkp1;j<n;j++,akj++)
      s += (akj->re)*(akj->re) + (akj->im)*(akj->im) ;
    s = sqrt(s);
    zt.re = s; zt.im = 0.0;
    r = sqrt((akkp1->re)*(akkp1->re) + (akkp1->im)*(akkp1->im));
    if(r>0.0){ zt.re = (akkp1->re)*s/r ; zt.im = (akkp1->im)*s/r; }
    akkp1->re += zt.re; akkp1->im += zt.im;
    bk->re = -zt.re; bk->im = zt.im;
    h = (r+s)*s;
    if(h>(2.e-16)){
      for(i=k+1,t=0.0,aki=akkp1,bi=bk+1;i<n;i++,aki++,bi++){
	zt.re = 0.0; zt.im = 0.0;
	for(j=k+1,aji=aki+n,akj=akkp1;j<=i;j++,aji+=n,akj++){
	  zt.re += (aji->re)*(akj->re) - (aji->im)*(akj->im);
	  zt.im += (aji->re)*(akj->im) + (aji->im)*(akj->re); }
	for(aij=aji-n+1;j<n;j++,aij++,akj++){
	  zt.re += (aij->re)*(akj->re) + (aij->im)*(akj->im);
	  zt.im += (aij->re)*(akj->im) - (aij->im)*(akj->re); }
	bi->re = zt.re/h; bi->im = zt.im/h;
	t += (bi->re)*(aki->re) + (bi->im)*(aki->im);
      }
      u = t*0.5/h;
      for(j=k+1,akj=akkp1,bj=bk+1;j<n;j++,akj++,bj++){
	(bj->re) = (akj->re)*u - (bj->re); (bj->im) = (akj->im)*u - (bj->im);
	for(i=k+1,aij=akj+n,aki=akkp1,bi=bk+1;i<=j;i++,aij+=n,aki++,bi++){
	  (aij->re) += (bi->re)*(akj->re) + (aki->re)*(bj->re)
	            + (bi->im)*(akj->im) + (aki->im)*(bj->im);
	  (aij->im) += (bi->re)*(akj->im) + (aki->re)*(bj->im)
	            - (bi->im)*(akj->re) - (aki->im)*(bj->re);
	}
      }
    }
    (akkp1-1)->re = h; (akkp1-1)->im = 0.0;
  }
  e[n-1] = (akkp1+n)->re;
  e[n-2] = (akkp1-1)->re;
  (b[n-2]).re = (akkp1->re); (b[n-2]).im = -(akkp1->im);
  (b[n-1]).re = 0.0;        (b[n-1]).im = 0.0;
  g = fabs(e[0]);
  for(j=1,bj=b+1,ej=e+1;j<n;j++,bj++,ej++){
    v = fabs((bj-1)->re)+fabs((bj-1)->im)+fabs(*ej);
    g = (v>g) ? v : g;
  }
  if(g<(1.e-16)) exit(0);
  d = (2.e-16)*g;
  k = n-1; bk = b+k; ek = e+k;
  while(k>0){
    for(l=k,bl=bk,el=ek;l>0;l--,bl--,el--) 
      if((fabs((bl-1)->re)+fabs((bl-1)->im))<d) break;
    if(l!=k){
      w  = (*(ek-1)+*ek)*0.5;
      r  = *ek - w;
      y  = sqrt((bk-1)->re*(bk-1)->re+(bk-1)->im*(bk-1)->im+r*r);
      x  = (r>=0.0) ? fabs(y)+w : -fabs(y)+w;
      f  = (*el-=x);
      zu = *bl;
      r  = sqrt(zu.re*zu.re+zu.im*zu.im+f*f);
      for(j=l,bj=bl,ej=el;j<k;j++,bj++,ej++){
	if(j>l){
	  r = sqrt((bj->re)*(bj->re)+(bj->im)*(bj->im)+(*ej)*(*ej));
	  (bj-1)->re = zs.re*r; (bj-1)->im = zs.im*r;
	  f = (*ej)*c;
	  zu.re = (bj->re)*c; zu.im = (bj->im)*c;
	}
	c = (*ej)/r;
	zs.re = (bj->re)/r; zs.im = (bj->im)/r;
	w = (*(ej+1))-x;
	(*ej) = (c*zu.re+w*zs.re)*zs.re + (c*zu.im+w*zs.im)*zs.im + f + x;
	(*(ej+1)) = c*w - zu.re*zs.re - zu.im*zs.im;
      }
      (bk-1)->re = zs.re*(*ek); (bk-1)->im = zs.im*(*ek);
      (*ek) = (*ek)*c + x;
    } else k--, bk--, ek--;
//    if (sigIntFlag) return;
  }
  free(b);
  j = n-1;
  while(j>0){
    k = 0; l = 0; m = 0;
    for(i=1;i<=j;i++){
      if(e[i]<e[k]){ m = k; l = i; }
      else k = i;
    } 
    if(l==m) j = l-1;
    else{ w = e[l]; e[l] = e[m]; e[m] = w; }
    j = l-1;
  }
}

