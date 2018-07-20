#include "./matrices.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <sstream>
#include <ctime>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <string>
#include<cstring>

using namespace std;
using namespace matrices;


double point(double istar, int maxsize, int minsize)
{

    int intistar = istar;
    double value = max(minsize,min(intistar,maxsize));
    return value;
}

double interpolate1d(double hS,  double ds , int smax, array &vec)
{

    int s1;
    int s2;
    

    double ds1;
 

    ds1 = hS / ds;


    s1 = point(ds1,smax-1,0);
    s2 = s1+1;


    double step1 = ((s2 * ds - hS) / double(s2 * ds - s1 * ds)) * vec[s1] + ((hS - s1 * ds) / double(s2 * ds - s1 * ds)) * vec[s2];
   
    return step1;
}


array triSolver(matrix &mat)
{

    double a1 = 0;
    double b1 = 0; // these are current row coeffs
    double c1 = 0;
    double d1 = 0;

    double a2 = 0; // these are next row coeffs
    double b2 = 0;
    double c2 = 0;
    double d2 = 0;

    array solution(mat.size()); // since the matrix is from i = 1 to n-1
    // static int count;
    // count++;


    for (int i = 0; i != mat.size() - 1; i++)
    {

        a1 = mat[i][0];
        b1 = mat[i][1];
        c1 = mat[i][2];
        d1 = mat[i][3];
        a2 = mat[i+1][0];
        b2 = mat[i+1][1];
        c2 = mat[i+1][2];
        d2 = mat[i+1][3];

        if (i == 0)
        {

            mat[i+1][0] = 0;
            mat[i+1][1] = b2; //b1 = 0
            mat[i+1][2] = c2; //c1 = 0
            mat[i+1][3] = d2 - (a2*d1/a1);

        }
        else if (i > 0 && i < (mat.size() - 2))
        {

            mat[i + 1][0] = 0;
            mat[i + 1][1] = (b2) - (a2*c1)/b1;
            mat[i + 1][2] = (c2);
            mat[i + 1][3] = d2 - (a2*d1/b1);

        }


    }
    solution[mat.size() - 1] = mat[mat.size() - 1][3] / mat[mat.size() - 1][2];

    for (int m = mat.size() - 2; m >= 1; m--)
    {
        solution[m] = (1. / mat[m][1])*(mat[m][3] - mat[m][2] * solution[m + 1]);
     //  cout << " solution " << solution[100] << endl ;
    }

    solution[0] = (mat[0][3]) / (mat[0][0]);
   // cout << solution[0] << endl;
    return (solution);
}
array hjbDiscretisation(int i, double dS, double g, double f ,double dtow, double vol) //
{

    array coeff(3);
    double alpha; // i-1 pivot location
    double beta; // i
    double gamma; //i+1
    double alphacentral;
    double betacentral;
    double gammacentral;
    double alphaforward;
    double betaforward;
    double gammaforward;
    double gammaback;
    double alphaback;

    alphacentral = -(vol * vol)*(i * i) / (2.) - f / ((2.)/dS);
    gammacentral = -(vol*vol)*((i)*(i)) / (2.) + f/ (2./dS);

    alphaforward = -(vol*vol)*((i)*(i))/2.;
    gammaforward = -(vol*vol) * i * i / (2.) + f/dS;

    alphaback = -((vol*vol)*i*i)/2. - f/dS;
    gammaback = -(vol * vol)*((i)*(i)) / 2.;

    if(alphacentral <=0 && gammacentral <=0 )
    {
        alpha = alphacentral;
        gamma = gammacentral;
    }
    else if(alphaforward<=0 && gammaforward <=0)
    {
        alpha = alphaforward;
        gamma = gammaforward;
    }
    else
    {
        alpha = alphaback;
        gamma = gammaback;
    }

        beta = -(gamma + alpha) + g ;

        coeff[0] = alpha*dtow;

        coeff[1] = 1 + dtow*beta;
        coeff[2] = gamma*dtow;

        return coeff;
}

int semiLagrangian(int gS, int gT, double SMAX, double T, double vol, double r , double K, array &result, array &initresult)
{

    double ds = SMAX/(double)gS;

    double dt = T/(double)gT;

    array initResult(gS+1);
    array preResult(gS+1);
   

    matrix tri(gS+1,array(4));

    for(int i = 0; i != gS+1; i++)
    {
        initResult[i] = i*ds - K > 0 ? i*ds - K :  0;
    }

    preResult = initResult;
    
    for(int t = 0;t != gT; t++)
    {
        
        for (int i = 0; i != gS+1; i++)
        {
            array hjb(3);
            hjb= hjbDiscretisation(i, ds, r, -r*i*ds ,dt, vol);
            if(i == 0)
            {

                tri[i][0] = 1;
                tri[i][1] = 0;
                tri[i][2] = 0;
                tri[i][3] = 0;
                continue;
            }


            if(i!=0 && i!=gS)
            {
                
                tri[i][0] = hjb[0];
                tri[i][1] = hjb[1];
                tri[i][2] = hjb[2];
                double hS = i*ds*exp(r*dt);
       
                tri[i][3] = interpolate1d(hS, ds, gS, preResult);
             
            }


            if(i==gS)
            {

                tri[i][0] = 0;
                tri[i][1] = 0;
                tri[i][2] = 1.;
                tri[i][3] = i*ds - K*exp(-r*((t+1)*dt));

                continue;
            }

        }
        preResult = triSolver(tri);
    }
//   for(int i = 0; i!=gS+1; i++)
//    {
//        cout << i*ds << "  : " <<preResult[i] << endl;
//    }
    result = preResult;
    initresult = initResult;
    return 1;
}

void printFile(const array &results, int gS,double SMAX)
{

    double ds = SMAX/(double)gS;
    char rname [300];
    sprintf(rname, "TestFileEuroOption");
    ofstream dataFile(rname, ios::out);
    dataFile.precision(10);
    cout.precision(10);

    if (!dataFile) {
        cout << " Error " << endl;
    }

  //  dataFile << "# " << setw(5)<<  "Alpha" << setw(15) << "S" << setw(15) << "V" << '\n' << endl;
    //cout <<"# " << setw(4) << "Alpha"<< setw(15) << "S" << setw(15) <<"V"  <<'\n' <<endl;

  
        for (int i = 0; i <= results.size() - 1; i++) {

            dataFile << showpoint << setw(6)  << ds * (i) << setw(25) << results[i] << endl;
            // cout <<showpoint <<setw(6) << k*da <<setw(15) << ds*(i)<< setw(15)  << results[k][i] << endl;
        }

        dataFile << endl;
        //cout << endl;
    
   
}


int main()
{
    int gS = 1000;

    int gT = 100;

    double SMAX = 10.;

    double T = 1.;

    double vol = 0.4;

    double r = 0.3;

    
    double K = 1.;

    array results(gS+1);
    array initResult(gS+1);
    //    for(int i = 0; i!=gS+1; i++)
//    {
//        cout << i << "  1 : " <<results[gS+1] << endl;
//    }
    semiLagrangian(gS,gT,SMAX,T, vol, r ,K, results, initResult);
   cout << "result at 2  =  " << interpolate1d(2,  10/(double)gS , gS, results)<< endl;
    printFile(results,gS,SMAX);
   // printFile(initResult,gS,SMAX);
    return (0);
}
