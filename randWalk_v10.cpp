//Neha Khetan : SOCM Lab, November 2013 //
// Random walk with drift
// modified on 13 December 2014: Switch case for different gradient shapes incorporated, 0: homograd (set by homowd), 1: exponential grad (set by rHalf), 2: sigmoid (set by rHals/sigSlope)
// modified on 12 feb 2015: adding cases for initiation: 0 = at center , 1: periphery of circle set by oocyte radius,  2: random nucleation (set by oocyte and nucleus radius).
// 18 feb 2015: adding box- muller transformation fucntion modified from Numerical recipes in C
// 20 Jan 2016: Added switch case for step choice: 0 = existing model ; 1 = vt * ( 1 + ( 1 -w) * RN ) ; 1 =  vt  + ( sqrt( vt )* ( 1 -w) * RN )
// 25 Jan 2016: Modified rwc.c as randWalk_v1.cpp  and using the MT 64 bit the inbuilt function in C++ : http://www.cplusplus.com/reference/random/mersenne_twister_engine/ , http://stackoverflow.com/questions/22923551/generating-number-0-1-using-mersenne-twister-c
// output file : Instead of dynamics static ==> randomwalkCoordinates.txt , wr_values.txt , wd_values.txt

//26 Jan: NOTE CHANGED THE ORDER OF values in the input parameter file. r1, s1, r2, s2 as line NO. 17 , 18, 19 , 20.

// 11 Feb 2016 : Modifying the Model and gradient Shape!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// hill added...which did not work earlier in 2014- however added now ------------ hill attract n reflect
// 



#include <stdio.h>
#include <stdlib.h>
#include<cmath>
#include<string.h>
#include<time.h>
#define pi 3.14
#include <iostream> 
#include "mt_rand.h"
#include <sstream>
#include <fstream>


using namespace std;

FILE *c;
FILE *status;
FILE *reflectMap;
FILE *grad;
FILE *step;
FILE *gradMap;



// Function to calculate radial distance when x and y coordinates are the input parameters //
double radial(double a, double b){
return(sqrt((a*a)+(b*b)));
}


// box muller transformation for converting uniform random number to gaussian random number
float bmt( float mu, float sd )
{				       
	static int iset =0;
	static float gset;
	float fac,w,v1,v2,rsq;
	
	if (iset==0)		        // use value from previous call 
	{
	float rn1,rn2;
		do{ 
			rn1= getUniRand();
			rn2= getUniRand();

			v1 = 2.0 * rn1 - 1.0;
			v2 = 2.0 * rn2 - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while ( rsq >= 1.0 || rsq ==0);
		fac = sqrt( (-2.0 * log( rsq ) ) / rsq );
		gset= v1*fac;
		iset =1;
		// return v2*fac;   // returns mean =0, sd =+/- 1
		return (v2*fac*sd +mu);
	}
	else{
		iset =0;
		//return gset;
		return (gset*sd +mu);
		}
}

// =========================================================================================================
// 
//=========  	 	Functions defining:  Gradient for Asymmetry ========================================

// Hill: Attractive Function
double hill_attract( double r_val , double km , double n ){
	double w =0.00;
	w = 1 / ( 1 + pow( ( r_val / km ) , n ))  ;
	   
	return w;
	}


// Hill: Reflective Function
double hill_reflect( double r_val , double km , double n ){
	double w =0.00;
	w = 1 / ( 1 + pow( ( ( 30 - r_val )/ ( km) ) , n ))  ;   
	return w;
	}


// exponential gradient : parameters (1): rhalf, i.e. half maximal range of the exp decay
	double expGradient(double r_val, double rhalf){
		double w= 0.00;
		//w = exp(-r_val* 0.693/r_half);
		//w = exp(-(r_val-nucleus_radius)* 0.693/rhalf);  // 13 dec             
		w = exp(-(log(2)*r_val/rhalf));  //ca: 13 dec             
		return w;
	}

// sigmoid gradient, parameters(2): rhalf: inflection point in case of steep gradients, s: slope
	double sig_attract( double r_val , double rhalf , double s ){
		double w=0.00;
		// w = 1/(1+exp(-(rc-r_val)/s));      13 dec
		// w = 1 / (1+exp((r_val-nucleus_radius)-rhalf)/s));     // 13 dec 
		w = 1 / (1+exp((r_val-rhalf)/s));     // ca: 13 dec 
		return w;
	}


// Reflecting Sigmoid boundary condition
	double sig_reflect(double r_val,  double refHalf ,double refSlope, double domain ){
		double w =0.00;
		w = 1 / (1+ exp(-( r_val - ( domain - refHalf)  )/ refSlope ));    		
		return w;
	}



// =======================================================================================================================================================
// *******************************************************************************************************************************************************

// stores the coordinates in x and y: for any given t: values for (t-1) and t are stored
double oldx[200]={}, oldy[200]={};
 
int main()
{

// files to write into:
std::stringstream ss;               // output file for randomwalk coordinates
std::stringstream wrsim;         // output file for wmap
std::stringstream wdsim;



//=========================================================================================================================================================
//      							            Reading the parameters from the inputParameter.txt file
//=========================================================================================================================================================

	    //Variable Parameters:  Read from the file 
		FILE *in;
		float number, totalTime;								
		int gradient , initialization , stepChoice , boundaryChoice , ref_gradient;

		float array[21];
		double Np,dt,wd=0.00,pstop, velocityInput, rHalf,wMin,wMax, D, homoWd=0, sigSlope=0 ,std ,wr;
		double oocyte_radius, nucleus_radius,pmembraneattach = 0.00 ,xx,yy,capture , refSlope, refHalf;
		double stepDiffusive, stepBallistic, step_net;




		in= fopen("inputParameter.txt","r");
		int count = 0;
		while( fscanf( in, "%f", &number) ==1 )							// reads the file till the end, stores the content in number
		{
			array[count] = number;								        // the content is added to the array
			count = count+1;
		}


		//setting the parameters as read from the input file {inputparameter.txt}		
		dt 					    	 = array[0];    //time interval betw sim steps
		refHalf                      = array[1];   // half value of the reflection function
		refSlope			         = array[2];   // Slope value for the reflection function
		wMin					 	 = array[3];    //scaling factor (not used)
		wMax 				    	 = array[4];    //scaling factor (not used)
		oocyte_radius 			 	 = array[5];    //whole cell radius (including the nucleus, chromatin...etc)
		nucleus_radius   		 	 = array[6];    //nuclear radius (inner region)
		totalTime			    	 = array[7];    //simulation time (sec)
		Np 				 	    	 = array[8];    //no. of particles (assumed to be constant for now)
		velocityInput       		 = array[9];    //particle velocity (ballistic!)
		pstop				    	 = array[10];    //prob. of stopping (default value = 0)
		pmembraneattach	 			 = array[11];   //cortial/cellwall attachment prob. (default=0)
		gradient			         = array[12];   //boolean to turn it on/off (wd gradient 0=off, 1=on)
		boundaryChoice		 		 = array[13];   //boolean boundary condition (bc) 0: no bc, 1:reflecting bc, with memb attach
		stepChoice 			     	 = array[14];   //step size distributions (not implemented, or used)
		homoWd 				     	 = array[15];   //Wd in case of no gradient. 0=Wd is 0, If (0 < homWd <= 1) will be used to set the Wd value
		initialization		     	 = array[16];   //location of particles at initial time, 0: center ; 1 =  on periphery, 2:random uniform (both r and theta are randomly set
		std					 		 = array[17];   // for std for gaussian steps.		
		rHalf				         = array[18];    // with respect to  chrom edge.
		sigSlope			         = array[19];   // for sigmoid gradient, this would be the slope value
		ref_gradient 				 = array[20];   // 0 => no reflection gradient, 1 => sigmoid gradient , 2 => hill gradient
		
	
         
		// noSteps and step size calculation
		int noSteps;
		noSteps = round(totalTime/dt);
	
									
//---------------------------------------------------------------------------------------------------------------------------------------------------------
//      							            End Of Parameter File
//=========================================================================================================================================================



//=========================================================================================================================================================
// 										Variables and Output Files used During the Code run
// ---------------------------------------------------------------------------------------------------------------------------------------------------------

	//variables used to store
	double rn1, rn2,rn3, rn4;			// RNs generated : initialization, stopping prob, random angles, membrane attachment
	double theta_rand_ini,theta_rand;		// Angles computed in radians for initialization and random
	double theta_drift, drift,phi, theta_net;
	double x,y;					// calculates the step distance covered: x & y respectively where, x=  stepDiffusive*cosine(theta_rand)& y=stepDiffusive *si(theta_rand)
	int statusParticle;					      // to store status of particle
	double r_previous,r_present,r_reverse,x_rev, y_rev,newx,newy;		// for boundary condition
		
		
	// Output files	created and written into
	step		= fopen("stepSize.txt","w");
	status 		= fopen("statusParticles.txt","w");
    gradMap 	= fopen("wd_map.out", "w");
	reflectMap  = fopen("wr_map.out", "w");
	


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
//                                                 End of variable declarations and Output files
//======================================================================================================================== //



// ======================================================================================================================
//                                    WRITE the Gradient MAPs:  13 December 2014
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ***** CENTRIPETAL Gradient map


double domain = oocyte_radius - nucleus_radius;
double k =0;
double wMap =0.00;
double wMapNorm =0.00;
double wrMap = 0.00;

	switch(gradient){
		case 0:{
			fprintf(gradMap,"%0.2lf\n",homoWd);   
			}break;
				   
		case 1:{
			// exponential grad
			double ymin = expGradient((oocyte_radius-nucleus_radius), rHalf);     
			double ymax = expGradient((nucleus_radius-nucleus_radius),rHalf);
					
			for(k= 0; k <= domain+0.1;k=k+0.1){	
				wMap = expGradient(k,rHalf);
				wMapNorm = ((wMax-wMin) *((wMap-ymin)/(ymax-ymin))) + wMin;
				fprintf(gradMap,"%0.2lf\t%0.5lf\n",k,wMapNorm);   
			}
			}break;
		case 2:{
		    	// sigmoid grad
		    	double ymin = sig_attract((oocyte_radius-nucleus_radius) , rHalf , sigSlope );
      			double ymax = sig_attract((nucleus_radius-nucleus_radius), rHalf , sigSlope );
			for(k= 0; k <= domain+0.1;k=k+0.1){
				wMap = sig_attract( k ,rHalf ,sigSlope );
				wMapNorm = ((wMax-wMin) *((wMap-ymin)/(ymax-ymin))) + wMin;
				fprintf(gradMap,"%0.2lf\t%0.5lf\n",k,wMapNorm);   
			}
			}break;
					
					
		case 3: {
			// Attractive hill function
			for(k= 0; k <= domain+0.1;k=k+0.1){
				wMap = hill_attract( k , rHalf , sigSlope );
				fprintf( gradMap,"%0.2lf\t%0.5lf\n",k,wMap ); 
			}
			} break;
			default:{
			printf("Hi! No gradient map implemented \n");
			}
		}


// ****  Reflection Gradient map


	switch(ref_gradient){
			
			case 0:{
					wrMap = 0.00;
					fprintf( reflectMap,"%0.2lf\t%0.2lf\n" , k , wrMap );
					}break;
			
			
			case 1:{
				double rbmax = sig_reflect((oocyte_radius  - nucleus_radius) , refHalf , refSlope  , domain );
				double rbmin = sig_reflect((nucleus_radius - nucleus_radius) , refHalf , refSlope  , domain );
				for(k= 0; k <= domain+0.1; k=k+0.1){
					wrMap = sig_reflect( k,  refHalf , refSlope  , domain );
					wrMap = ((( wrMap - rbmin)/( rbmax - rbmin )));
					fprintf(reflectMap,"%0.2lf\t%0.2lf\n" , k , wrMap );

					}
				
				} break;
				
			case 2: {
				// Reflect hill function
				for(k= 0; k <= domain+0.1;k=k+0.1){
					wrMap = hill_reflect( k , refHalf , refSlope );
					fprintf( reflectMap,"%0.2lf\t%0.2lf\n" , k , wrMap  );
					printf( "%0.2lf\t%0.2lf\n" , k , wrMap  );
					};
					}break;
					
					
			   default:{
			printf("Hi! No REFELCTION gradient map implemented \n");
			}
		}
					
					

//                                     End of Printing Gradient Maps into a file
//======================================================================================================================== //





// ========================================================================================================================
//                                    INITIALIZATION OF PARTICLES 
// .........  ( 0 = at  centre , 1 = periphery , 2 = uniform random in space between nucleus and oocyte)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// last modified: removing redundancy: 23 Jan 2016


// outer loop is for time i.e. At time = 0, inner loop is for all the particles,so run through each particle.
switch(initialization){

	// at center
	case 0: {
	
			int i,j,part;
			for (i=0;i<1;i++)
			{
			double capture =0.00;
			//fprintf(out,"%0.3lf\t",(i*dt));
			ss << ( i * dt )  << "\t" ;
			
			
			for ( j=0 ; j< Np ; j++)
				{	
				oldx[j] = 0.00;
				oldy[j] = 0.00;
				// fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", oldx[j],oldy[j],capture);
				ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;
				}
			
			}

		}break;				
			
	case 1: {
			
			// at periphery					
			int i,j,part;
			for (i=0;i<1;i++)
				{
					double capture =0.00;
					//fprintf(out,"%0.3lf\t",(i*dt));
					ss << ( i * dt )  << "\t" ;
			
					double theta_rand_ini , rn1 ;
	
						for (j=0; j<Np ; j++)
						{	
							rn1 = getUniRand();
							theta_rand_ini = rn1*360*pi/180;
			
							oldx[j]= oocyte_radius * cos (theta_rand_ini);
							oldy[j]= oocyte_radius * sin (theta_rand_ini);
							// fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", oldx[j],oldy[j],capture);
							ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;
						}
				}
			}break;
			
	
	 case 2: {
			
			 // random nucleation
			int i,j,part;
			for (i=0;i<1;i++)
			{
			double capture =0.00;
			//fprintf(out,"%0.3lf\t",(i*dt));
			ss << ( i * dt )  << "\t" ;
			
			double initialStep= 0.00 , rr , theta_rand_ini, rn1 ;
	
				for (j=0; j<Np; j++)
				{	
					rr = getUniRand();
					
					rn1 = getUniRand();
					theta_rand_ini = rn1*360*pi/180;
					oldx[j]= (((oocyte_radius-nucleus_radius) * rr) + nucleus_radius) * cos (theta_rand_ini);
					oldy[j]= (((oocyte_radius-nucleus_radius) * rr) + nucleus_radius) * sin (theta_rand_ini);
					ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;				 
				 }
			}		
		}break;	
		
		default:
			printf("\n Please enter a valid initialization condition !! \n");	
	}
				
		
	
//-------------------------------------------------------------------------------------------------------------------------------
//                                                  End of INITIALIZATION
//======================================================================================================================== //




//======================================================================================================================== //
//                                    BEGIN TO DO THE RANDOM WALK!!! :):) !! 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//Random Walk

	int i = 1 ;
	int j = 0 ;

	for ( i=1; i<= noSteps; i++ )
	{
	double newx[200]={}, newy[200]={}; 
	
	//fprintf(out,"\n%0.3lf\t",(i*dt));
	ss << "\n" << ( i * dt )  << "\t" ;
	
	//fprintf(status,"\n%0.3lf\t",(i*dt));
	//fprintf(step,"\n%0.3lf\t",(i*dt));
		for (j=0; j<Np; j++)
	    	{
		r_previous, r_present, xx=0,yy=0;

	    // Determine the State of the particle in delta time step //
		       
		rn2 =  getUniRand();
		theta_rand = rn2*360*pi/180;


		// Stop particle State
		if (rn2 < pstop)
		    {
		    capture =0.00;
		        statusParticle =1;                                  		// if true , status = Stop = 1
		        newx[j]= oldx[j];
		        newy[j]= oldy[j];
		        //fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", newx[j],newy[j],capture);
		        //fprintf(status,"%d\t",statusParticle);                         // prints the status of the particle
		        ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;	
		    }
		// Go particle State		 
		else
		    {
		        statusParticle=0;                                    	      // iffalse, status = Move = 0
		        //fprintf(status,"%d\t",statusParticle);
		       
		        //............ Random walk................. //
		                             //Return a pseudorandom double in [0,1) with 32 bits of randomness : n3 for Random Walk Angles
			    rn3        =  getUniRand();
		        theta_rand = rn3*360*pi/180;
			               
		        
		        
		        //----------------------------------------------------------------------------- //
		        //........... End of random walk........... //

		       
		      //----------------------------------------------------------------------------- //
		        //............Drift Calculation ............ //
		            phi =(atan2 (oldy[j],oldx[j]))*180/pi;          //(in degrees) partice at its intial position makes an angle phi with + x axis
		            // Pushing the particle towards the center
		            if (phi >0)
		            {
		            drift= phi - 180;
		            }
		            else if ( phi < 0 ){
		            drift = phi + 180;
		            }
		            else if ( phi ==0 ){
		            drift = 180;
		            }
		            theta_drift = drift *pi/180;  // angles in radians
		        //........... drift angle calculation ........// 




		// To calculate the position of the particle at (t-1) time, which will determine the Wd experienced in the move during the period (t-1) to t	
		
		r_previous = radial(oldx[j],oldy[j]);   // ensures that for each particle it's own previous position is considered for wd calculation


//....................................................................................................
// To implement a switch case scenario fro gradient implementation....13 dec 2014
switch(gradient){
	case 0: {
			wd = homoWd;
			wdsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wd <<  "\n" ;
			} break;

	
	case 1:{
	// Exponential Gradient: calls the function expGradient: expGradient(r_value, r_half) where both are with respect to the chromatin edge.
	// r_values are subtracted by nucleus radius: value as r_value ranges from [nucleus radius to oocyte radius, subtracting nuc. radius from each to 		shift origin from center to nuclear radius
			double w=0;
			double ymin = expGradient((oocyte_radius-nucleus_radius), rHalf);     
			double ymax = expGradient((nucleus_radius-nucleus_radius),rHalf);
			// for min-max normalization
			//w = expGradient(r_previous,rHalf);
			// on 13 dec. commenting out
			// wd = ((wMax-wMin) *((w-ymin)/(ymax-ymin))) + wMin;		//normalized weight of drift for all sorts of functions  13 dec
			
			
			
			
			w = expGradient((r_previous-nucleus_radius),rHalf); 								// 13 dec 2014
			wd = ((wMax-wMin) *((w-ymin)/(ymax-ymin))) + wMin;		//normalized weight of drift for all sorts of functions  13 dec		 			
			wdsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wd <<  "\n" ;
			} break;

	// Sigmoid gradient
	case 2: {
	 
			double w=0.00;
			double ymin = sig_attract( (oocyte_radius-nucleus_radius)  ,rHalf  ,sigSlope );
			double ymax = sig_attract( (nucleus_radius-nucleus_radius) ,rHalf  ,sigSlope);
			
			//w = sig_attract(r_previous,s,rHalf); // 13 dec....
			//wd = ((wMax-wMin) *((w-ymin)/(ymax-ymin))) + wMin;		//normalized weight of drift for all sorts of functions
			//fprintf(wdvalues,"%0.2lf\t%0.2lf\t%0.2lf\t%0.5lf\n",oldx[j],oldy[j],r_previous,wd);
		
			
			w = sig_attract( (r_previous-nucleus_radius) ,rHalf ,sigSlope );  // 13 dev 2014
			wd = ((wMax-wMin) *((w-ymin)/(ymax-ymin))) + wMin;		//normalized weight of drift for all sorts of functions
			wdsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wd <<  "\n" ;

			}break;
			
			
			
	case 3: {
            wd = 0;
			wd = hill_attract( (r_previous - nucleus_radius) ,rHalf ,sigSlope ); 

			wdsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wd <<  "\n" ;
			// printf( " i am in hill attract\n");
			}break;
	default :
			
			printf("\n Invalid gradient input entered! \n");
			break;	
			}
		
// -----------------------------

// FOR REFLECTION
switch(ref_gradient){
			case 0:{
					double wr = 0.00;	
					wrsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wr <<  "\n" ;
				} break;
				
			case 1:{
					
					wr = 0.00;	
					 double rbmax = sig_reflect((oocyte_radius  - nucleus_radius) , refHalf , refSlope  , domain );
					 double rbmin = sig_reflect((nucleus_radius - nucleus_radius) , refHalf , refSlope  , domain );
					 wr = sig_reflect((r_previous-nucleus_radius) , refHalf , refSlope  , domain);      // 13 dev 2014
					wr = (( wr -rbmin )/( rbmax-rbmin )) ;		//normalized weight of drift for all sorts of functions
					wrsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wr <<  "\n" ;
				} break;
			case 2:{
					 wr = 0.00;
					wr = hill_reflect(  ( r_previous - nucleus_radius)   , refHalf  , refSlope);
					wrsim << oldx[j] <<  "\t" << oldy[j] <<  "\t" <<  r_previous <<  "\t" <<  wr <<  "\n" ;
					} break;
				
			default :
				printf("\n Invalid Ref gradient input entered...during simulation! \n");
				break;	
			}
				


// Net Angle Calculation :
// angle in radians   
// nk: 21 jan 2016: theta_net = atan2(((wd* sin (theta_drift)) +((( 1 -wd)*sin(theta_rand))), ((wd*cos(theta_drift)) + ( 1 - wd)*cos(theta_rand)));  

// // change on 22 Jan 2016 : Adding reflection angle

// case 2b
double numY =     ( wd*sin(theta_drift) ) +  ( ( 1 -wd )*( ( 1 -wr )*sin(theta_rand) + ( wr*sin(theta_drift)) ));
double denX =     ( wd*cos(theta_drift) ) +  ( ( 1 -wd )*( ( 1 -wr )*cos(theta_rand) + ( wr*cos(theta_drift)) ));

// // CASE 3
// double numY =     ( ( wd + wr ) *sin(theta_drift) )  +   ( ( 1 - wd - wr )*sin(theta_rand) ) ;
// double denX =     ( ( wd + wr ) *cos(theta_drift) )  +   ( ( 1 - wd - wr )*cos(theta_rand) ) ;


theta_net = atan2(  numY , denX );      


// ======================================================================================================================
// ==============================   STEP CHOICE CALCULATION ====================================	
						 		
						 	switch(stepChoice){
							case 0: {		
				
				
									// Diffusion Coefficient: set by scaling the input ballistic motion velocity. That val. is itself derived from motor velocities...
					
								    D =	0.006 ;
									stepDiffusive = sqrt( 4  * D * dt);
									stepBallistic = velocityInput *dt ;
									// step_net = ( wd*stepBallistic ) +  ( ( wMax+wMin-wd ) * stepDiffusive ) ; 
									
									//case 2									
									step_net = ( wd*stepBallistic ) +  ( ( 1 - wd ) * (   ( wr * stepBallistic ) +  ( ( 1 -wr) *  stepDiffusive ) ) ) ; 
									//printf("%0.3lf \t %0.3lf \t %0.3lf \t%0.3lf\n " ,  r_previous, wd ,  wd*stepBallistic , ( 1 - wd) *stepDiffusive);

									// case 3:
									//step_net = ( ( wd + wr ) * stepBallistic  )  +  ( ( 1 - wd - wr ) * stepDiffusive )   ; 
									
									//fprintf(step,"%f\n", step_net);
									}break;
						
							case 1:
							 		{
							 		// assumption is that diffusive and drif tstep are same i.e.: mean ( drift movement ) and SD ( spread of diffusive movement ) is SAME.
							 		// or the diffusive magnitude is comparable to the drift magnitude within 1 SD
							 		//   v*t  * ( 1 + { ( 1 - w)* RN( u =0, std = 1 )}
							 		
							 		
							 		float rng_gauss;
									rng_gauss = bmt( 0,std );
							 		step_net = (velocityInput*dt) * ( 1 +  ((1 - wd) * rng_gauss) );
							 		//printf("%f\n",step_net);
							 		//fprintf(step,"%f\n", step_net);
									}break;
						
							case 2:
							 		{
							 		// assumption is that diffusive and drift are not really same
							 		// :the spread due to the diffusion is in magnitude equated to the sqrt ( mean ) or implying that: mean is comparable to the variance of the diffusive
							 		//   v*t   + { sqrt(v*t) * ( 1 - w)* RN( u =0, std = 1 )}
							 		
							 		float rng_gauss;	
							 		rng_gauss = bmt( 0, std );
							 		step_net = (velocityInput*dt) + ( sqrt(velocityInput*dt) * ((1 - wd) * rng_gauss) );
							 		// printf("%f\n",step_net);
							 		//fprintf(step,"%f\n", step_net);
									}break;
							case 3:
							 		{
							 		// assumption is that diffusive and drift are not really same
							 		// :the spread due to the diffusion is in magnitude equated to the sqrt ( mean ) or implying that: mean is comparable to the variance of the diffusive
							 		//   v*t   + { sqrt(v*t) * ( 1 - w)* RN( u =0, std = 1 )}
							 		
							 		float rng_gauss;	
							 		rng_gauss = bmt( 0, std );
							 		step_net =     ( wd * (velocityInput*dt) ) + ( (velocityInput*dt) * ((1 - wd) * rng_gauss) );
							 		// printf("%f\n",step_net);
							 		//fprintf(step,"%f\n", step_net);
									}break; 	
							 	
							 	
							 default :
						
									printf("\n PLease enter a valid Step choice !!! \n");
									break;	
									}	
				 	
			 		
			 		 

		       
		        // Calculate delta step size    
		            x = step_net * cos(theta_net);
		            y = step_net * sin(theta_net);
		       

		        /* Calculate New position */
		            xx = oldx[j] + x;
		            yy = oldy[j] + y;


///////////////////////////////////////////// Boundary Condition ///////////////

			



double capture =0.00;

			if (boundaryChoice==1.00)
			{
			// To calculate the position of the particle at (t-1) time, which will determine the Wd experienced in the move during the period (t-1) to t	
			r_previous = radial(oldx[j],oldy[j]);   // ensures that for each particle it's own previous position is considered for wd calculation
			
			r_present  = radial(xx,yy);  
 					
			if(r_present<= nucleus_radius)
				
				{
				
					if(r_previous =nucleus_radius)
						{
						capture=1.00;	
				            	newx[j]= oldx[j];
				            	newy[j]= oldy[j];
					    	// fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", newx[j],newy[j],capture);
				            	ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;	
	 					}
	 				else
	 					{
	 			
	 					capture=1.00;	
				            	newx[j]= nucleus_radius*cos(theta_net);;
				            	newy[j]= nucleus_radius*sin(theta_net);
					    	//fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", newx[j],newy[j],capture);
					    	ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;	
					    	}
	 				
	 			
	 			} 			

			else
				{

				    if(r_previous =oocyte_radius)
				    {
				    	   
					     
					    rn4=  getUniRand();
						if (rn4 < pmembraneattach)
						{
						    newx[j]= oldx[j];
						    newy[j]=oldy[j];

						}
						else
						{
						    if (r_present <= oocyte_radius)
						    {
						        newx[j]=oldx[j]+x;
						        newy[j]=oldy[j]+y;
	
						    }
						    else
						    {
						  
						     x_rev=oldx[j]-x;
						     y_rev=oldy[j]-y;
		            		   
		            		 
		            		
						    r_reverse= radial(x_rev,y_rev);
						        if (r_reverse <=oocyte_radius)
						        {
						         newx[j]= x_rev;
						         newy[j]= y_rev;
						        }
						        else
						        {
						        newx[j]= oldx[j];
						        newy[j]= oldy[j];
						        }
						    }
												
						}

						// fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", newx[j],newy[j],capture);
						ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;	
					
					    }

				    else
				    {
					
				        if (r_present <=oocyte_radius)
				        {
				        newx[j]=oldx[j]+x;
				        newy[j]=oldy[j]+y;
				        }
				        else
				        {
				        x_rev=oldx[j]-x;
				        y_rev=oldy[j]-y;
				        r_reverse= radial(x_rev,y_rev);
				            if (r_reverse <=oocyte_radius)
				            {
				            newx[j]= x_rev;
				            newy[j]= y_rev;
				            }
				            else
				            {
				            newx[j]= oldx[j];
				            newy[j]=oldy[j];
				            }
				        }

					//fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", newx[j],newy[j],capture);
				    ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;	
				    }
			  }     
			/////////////////////// end of boundary condition
			 
			}
			else
			{
			newx[j]=oldx[j]+x;
			newy[j]=oldy[j]+y;
			//fprintf(out,"%0.2lf\t%0.2lf\t%0.2lf\t", newx[j],newy[j],capture);
			ss << oldx[j]  << "\t" << oldy[j] << "\t" << capture << "\t" ;	
			//fprintf(step,"%0.2lf\t%0.2lf\t", r_present,(r_present-r_previous));
			}
		
			oldx[j]	=newx[j];
			oldy[j]	=newy[j];
			
			
	    	} 
	   
	    }  // particle loop ends
}  // time loop end





// ====================================== FILE OUTPUTS ========
std::string s  = ss.str();
std::string w1 = wdsim.str();
std::string w2 = wrsim.str();


// Output the FilerandomwalkCoordinates
ofstream outfile;
outfile.open("randomwalkCoordinates.out");
outfile << s ;
outfile.close();

// Output for the file Wr values in simulation
ofstream wrval;
wrval.open("wr_sim.out");
wrval << w2 ;
wrval.close();


// Output for the file Wr values in simulation
ofstream wdval;
wdval.open("wd_sim.out");
wdval << w1 ;
wdval.close();


return 0 ;
};
