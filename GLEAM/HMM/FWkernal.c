#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//Declarations
mxArray *laData, *mmData, *raData;
double *laValues,*mmValues, *raValues, *outArray;
int i,j;


//Copy input pointer la mm ra
laData = prhs[0];
mmData = prhs[1];
raData = prhs[2];

//Get matrix x
laValues = mxGetPr(laData);
mmValues = mxGetPr(mmData);
raValues = mxGetPr(raData);


//Allocate memory and assign output pointer
plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL); //mxReal is our data-type

//Get a pointer to the data space in our newly allocated memory
outArray = mxGetPr(plhs[0]);

//
for(i=0;i<3;i++) // row
{
    for(j=0;j<3;j++) // column
    {
        outArray[(i*3)+j] = laValues[j]*mmValues[(i*3)+j]*raValues[i];
    }
}
            
    
    return;
}
