#include "random.h"

long setseed(void) {
	time_t lt;
	lt=time(NULL);
	return (long)lt;
}


/*C routine random number generator from
Numerical Recipes in C: Press et al*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(void) {
	int j;
	long k;
	extern long *idum;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7; j>=0; j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
			
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j]=*idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double ran2(double min, double max)
{
	double r = ran2();
	r *= (max - min);
	r += min;
	return r;
}

int ran2(int min, int max)
{
	double r = ran2();
	r *= (double)(max - min);
	r += min;
	return (int)r;
}

int ran2(int max)
{
	double r = ran2();
	r *= (double)max;
	return (int)r;
}

// Returns LnGam(xx)
// From Numerical Recipes in C
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// AdamA: Normal Distrubted random number generator
// From Numerical Recipes in C, Press et al. 1988
// Returns a normally distributed random number with zero mean and unit variance
double normrnd(void)
{
	extern long *idum;
	static int iset=0;
	static double gset;
	double fac, r, v1, v2;

	if (iset == 0)
	{
		do
		{
			v1 = 2.0 * ran2() - 1.0;
			v2 = 2.0 * ran2() - 1.0;
			r = v1 * v1 + v2 * v2;
		} while (r >= 1.0);
		
		fac = sqrt(-2.0*log(r) / r);
		gset = v1 * fac;
		iset = 1;
		return v2*fac;
	} 
	else
	{
		iset = 0;
		return gset;
	}
}

// Returns random numbers distributed as a gamma distribution
// From Devroye 1986
double gamrnd(double ia, double ib)
{
	double b, c, d, u, v, w, x, y, z;
	int accept;
	double ret = 0;

	if (ia == 1.0)
	{
		// gamma is exponetial (Devroye pg 405)
		ret = -ib * log(ran2());
	}
	else if ((ia < 1) && (ia > 0))
	{
		c = 1 / ia;
		d = 1 / (1 - ia);
		accept = 0;
		do 
		{
			u = ran2();
			v = ran2();
			x = pow(u, c);
			y = pow(v, d);
			z = x + y;
			if (z <= 1.0)
				accept = 1;
		} while (accept != 1);

		ret = -ib * log(ran2()) * x / z;
	}
	else if (ia > 1)
	{
		b = ia - 1;
		c = 3.0 * ia - 0.75;
		accept = 0;
		do
		{
			u = ran2();
			v = ran2();
			w = u * (1 - u);
			y = (u - 0.5) * sqrt(c / w);
			x = b + y;
			if (x >= 0.0)
			{
				z = 64.0 * pow(w, 3) * pow(v, 2);
				if (z <= (1 - 2 * y * y / x))
				{ accept = 1; }
				else
				{	
					if (log(z) <= (2 * (b * log(x / b) - y)))
					{ accept = 1; }
				}
			}
		} while (accept != 1);

			ret = ib * x;
	}
	return ret;
}


// AdamA: Poisson Distrubted random number generator
// For Numerical Recipes in C, Press et al. 1988
// Returns a poisson distributed random number with mean xm

int poissrnd(double xm)
{
	//float gammln(float xx);
	static double sq, alxm, g, oldm=(-1.0);
	double em,t,y;
	
	if (xm < 12.0) { 
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm); 
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran2();
		} while (t > g);
	} else { 
		if (xm != oldm) { 
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g = xm*alxm - gammln(xm + 1.0);
		}
		do {
			do { 
				y=tan(PI * ran2());
				em=sq*y+xm; 
			} while (em < 0.0); 
			em=floor(em); 
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran2() > t);
	}
	return (int)em;
}

double exprnd(double mu)
{
	return -mu * log(ran2());
}

int bnlrnd(double pp, int n)
{
	int j;
	static int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;
	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;		//This is the mean of the deviate to be produced.
	if (n < 25) { 
		bnl = 0.0;
		for (j=1;j<=n;j++)
			if (ran2() < p) ++bnl;
	} else if (am < 1.0) { 
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran2();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) { 
			en=n; 
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) { 
			pc=1.0-p; 
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc); 
		do {
			do {
				angle=PI*ran2();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0)); 
			em=floor(em); 
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran2() > t); bnl=em; 
	}
	if (p != pp) bnl=n-bnl; 
	return (int)bnl;
}

