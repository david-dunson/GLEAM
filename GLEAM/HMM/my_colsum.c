#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//Declarations
mxArray *mmData ;
double *mmValues, *outArray;
int i,j;


//Copy input pointer mm 
mmData = prhs[0];
//Get matrix x
mmValues = mxGetPr(mmData);

//Allocate memory and assign output pointer
plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL); //mxReal is our data-type
//Get a pointer to the data space in our newly allocated memory
outArray = mxGetPr(plhs[0]);

for(j=0;j<3;j++) // column
{
  outArray[j] = mmValues[j]+mmValues[3+j]+mmValues[6+j];
}

return;
}