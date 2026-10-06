// ==============================================================================
// File Name: sis_star_ni.c
// Author: Juliane T. de Moraes
// Updated: 2026-09
// License: MIT
// Usage: ./sis_star_ni beta taumed
// Description: This code computes the lifespan in a star graph of intermittent nodes
// See Ref.: Moraes, J. T. and and Ferreira, S. C. (2026). Intermittent quarantine 
// suppresses epidemic spreading beyond simple contact reduction [Preprint].
// arXiv. https://doi.org/10.48550/arXiv.2610.00429
// ==============================================================================

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

// ---------------------------------------------------------
// KISS pseudo-random number generator

#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836

#define Nmax 1000000000

unsigned int z, jsr, w, jcong;

double rand01_kiss(){
   z = 69069*z+1327217885;
   jsr^= (jsr<<13); jsr^=(jsr>>17); jsr^=(jsr<<5);
   w = 18000 * (w & 65535) + (w >> 16);
   jcong = 30903 * (jcong & 65535) + (jcong >> 16);
   return( (z + jsr + (w << 16) +jcong) * 2.328306436538696e-10);
}

void kissinit(long idum){
   int k;
   idum= abs(1099087573 * idum);
   if (idum == 0) idum = 1;
   if (idum >= IM) idum = IM-1;

   k=(idum)/IQ;
   idum=IA*(idum-k*IQ)-IR*k;
   if (idum < 0) idum += IM;
   if (idum < 1) z=idum+1; else z=idum;
   k=(idum)/IQ;
   idum=IA*(idum-k*IQ)-IR*k;
   if (idum < 0) idum += IM;
   if (idum < 1) jsr=idum+1; else jsr=idum;
   k=(idum)/IQ;
   idum=IA*(idum-k*IQ)-IR*k;
   if (idum < 0) idum += IM;
   if (idum < 1) w=idum+1; else w=idum;
   k=(idum)/IQ;
   idum=IA*(idum-k*IQ)-IR*k;
   if (idum < 0) idum += IM;
   if (idum < 1) jcong=idum+1; else jcong=idum;

 } 
 // ------------------------------------------------

int main(int argc,char *argv[]){

  if (argc != 3){
    printf("\n Usage: %s beta taumed \n\n",argv[0]);
    return 1;
  }
  
  int N,L,e,*adj,*beg,*dgr,g,m,i,semente;
  char filename[4*64];
  
  semente = 10;
  kissinit(semente);  

  dgr = (int *) malloc(Nmax * sizeof(int));
  beg = (int *) malloc(Nmax * sizeof(int));
  adj = (int *) malloc(Nmax * sizeof(int));
		
  double lambda,deltat,t,p,Z,tau,freqtau,somatau,tls,beta,tau0,tau_med,r,*t_ch,prob;
  int Ni,*I,viz,Q,qq,Nn,*sigma,b,q,rr,runs,out,*EST,f;     

  I = (int *) malloc(Nmax * sizeof(int));
  sigma = (int *) malloc(Nmax * sizeof(int));
	
  EST = (int *) malloc(Nmax * sizeof(int));
  
  t_ch = (double *) malloc(Nmax * sizeof(double));       
  
  // constants for the inter-event time distribution
    
  beta = atof(argv[1]);         // exponent for the power-law distribution 
    
  tau0 = 0.05;         // constant of the power-law distribution 
    
  tau_med = atof(argv[2]);
    
  r = beta*(tau_med - tau0) + tau0 - 2*tau_med;
  
  lambda = 0.1;
  

  time_t seconds;
   
  runs = pow(10,4);
  
  int nn = 25;          // Number of points
  double m0 = 30;      // Initial value
  double mf = 5000;     // Final value
  double R = pow(mf/m0, 1.0/(nn - 1));

  double t_rlx = pow(10,2);

  m = m0;
  int n = 1;
  
  do{         //hub degree loop
  
  
    N = m + 1;
    
    for(i=1;i<=N;i++){
	
      dgr[i] = 0;
      beg[i] = 0;
	
    }
    
    dgr[1] = m; //hub 
	
    for(i=2;i<=m+1;i++){
	
      dgr[i]++;
		
    }

    L=0;

    for(g=1;g<=N;g++){

      beg[g]=L;
      L=L+dgr[g];
				
    }
    


    for(e=2;e<=m+1;e++){
        
      adj[e-1]=e;
          
    }
    
    for(e=m+1;e<=L;e++){
    
      adj[e]=1;
    
    }
    
    semente++;
    kissinit(semente);  
	
	
    sprintf(filename, "sis_star_ni_lifespanxm_lambda%.1f_beta%.1f_taumed%.1f_tau0_%.2f_1e4runs.dat",lambda,beta,tau_med,tau0);
    FILE *arq1;
    arq1 = fopen( filename, "a+");
	      
    somatau = 0.0;
    freqtau = 0.0;
    
    for(rr=1;rr<=runs;rr++){		
					
      seconds = time(NULL);
      semente=seconds+rr;
      kissinit(semente);
		
      t=t_rlx;	
      tls = 0.0;
		
    
      for(f=1;f<=N;f++){
      
        prob = rand01_kiss(semente);
      	
        if(prob>0.5){EST[f]=1;}else{EST[f]=0;}
      	
        t_ch[f] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r;
		
      
      }
      
      for(q=1;q<=N;q++){
      
        if(t_rlx>t_ch[q]){         // it means the node q should already change its state
      	  
          do{
      	  
            if(EST[q]==1){
          
              EST[q]=0;                
                                               
            }else{EST[q]=1;} 
   
            t_ch[q]+= (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r;
               
          }while(t_ch[q]<t_rlx); 		  
        }
      }
      
      // EPIDEMIC PROCESS
      
      q=0;
      b=0;

      lambda=0.1;
      
      Ni=1;				
      Nn=m;				
      
      I[1] = 1;			        
      sigma[1] = 1;
		
      for(qq=2;qq<=N;qq++){			
	
        sigma[qq]=0;			
	
      }
      
      out = 0;          
      
      do{
      
        Z=rand01_kiss(semente);
    
        p=(Ni)*pow((Ni+(lambda*Nn)),-1);	
               
        
        if(p>Z){		
			
	  Q=(Ni)*rand01_kiss(semente)+1;		

	  q=I[Q];
			
	  if(t>t_ch[q]){		
                    
            do{
                                           
              if(EST[q]==1){	
                            
                EST[q]=0;
                            
              }else{		
                            
                EST[q]=1;
                                
              }
                        
              t_ch[q] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[q];
                        
            }while(t_ch[q]<t); 		    
          }	
	
	  sigma[q]=0;
	  I[Q]=I[Ni];

	  Ni--;
	  Nn=Nn-dgr[q];
			
	  if(Ni==0){		
				
	    tls = t;			
				
	    somatau = somatau + tls;		
				
	    freqtau = freqtau + 1.0;
		
	    out = 1;			
	
	  }

        }else{			
					
          if(sigma[1]==0){		
			
	    if(t>t_ch[1]){
                    
              do{          
                                        
                if(EST[1]==1){	
                            
                  EST[1]=0;
                                                                  
                }else{		
                            
                  EST[1]=1;          
               
                }
                        
              
                t_ch[1] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[1];
                        
              }while(t_ch[1]<t); 		
                	 
            }
				
	    if(EST[1]==1){	
				
	      Ni++;          
              I[Ni]=1;          
              Nn=Nn+dgr[1];          
              sigma[1]=1;
			
	    }
			
          }else{				//hub
			
	    b=(dgr[1])*rand01_kiss(semente)+1; 	
      
	    viz = adj[beg[1]+b];
					
	    if(t>t_ch[viz]){
                        
              do{
                                                     
                if(EST[viz]==1){	
                                
                  EST[viz]=0;                                           
       
                }else{	
                                
                  EST[viz]=1;
			
		}

                t_ch[viz] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[viz];
         
              }while(t_ch[viz]<t); 		
            }
                    
            if((sigma[viz]==0) && (EST[viz]==1)){ 
                      
              Ni++;
              I[Ni]=viz;          
              Nn=Nn+dgr[viz];          
              sigma[viz]=1;
             
            }	
          }
	}  
	
	deltat=-log(rand01_kiss(semente)+pow(10,-8))*pow((Ni+(lambda*Nn)),-1);
	
        t+=deltat;		
          
      }while(out!=1);

      
      
    }   

    tau = somatau*pow(freqtau,-1) - t_rlx;

    fprintf(arq1,"%d %e \n",m,tau);	
    
    fclose(arq1);
  
    n++;
    
    m=m0*pow(R,n);
  
  }while(m<=mf); 


  free(I);
  free(sigma);

  free(EST);

  free(dgr);
  free(beg);
  free(adj);
  
  free(t_ch);

  return (0);
}

