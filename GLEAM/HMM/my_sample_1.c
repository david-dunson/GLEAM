#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //Declarations
mxArray *pData;
double *pValues, *z;
int i;
double cumsum[3],u;

pData = prhs[0];

pValues = mxGetPr(pData);

plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type

z = mxGetPr(plhs[0]); //out z

cumsum[0]=pValues[0];
for(i=1;i<3;i++)
{
    cumsum[i]=cumsum[i-1]+pValues[i];
}

u = rand()/(RAND_MAX+1.0);

if(u<cumsum[0]/cumsum[2])
{
    *z=1;
}
else if(u<cumsum[1]/cumsum[2])
{
    *z=2;
}
else
{
    *z=3;
}


return;    
}