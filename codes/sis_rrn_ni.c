// ==============================================================================
// File Name: sis_rrn_ni.c
// Author: Juliane T. de Moraes
// Updated: 2026-09
// License: MIT
// Usage: ./sir_rrn_ni beta taumed lambda_ini lambda_fin delta_lamb expon
// Description: This code computes the susceptibility and the quasi-stationary 
// average of prevalence as function of infection rate for the SIS model in a
// node intermittent random regular network.
// See Ref.: Moraes, J. T. and and Ferreira, S. C. (2026). Intermittent quarantine 
// suppresses epidemic spreading beyond simple contact reduction [Preprint].
// arXiv. https://doi.org/10.48550/arXiv.2610.00429
// ==============================================================================

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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

  int N,L,f,e,*adj,a,somadgr,n,g,*S,q,o,stop,*beg,*dgr,*DGR,cont2,alerta,aux,*sigma,s,u,links,k;
  int semente;
  char filename[4*64];
	
  int expon,smp;
  
  if (argc != 7){
    printf("\n Usage: %s  beta taumed lambda_ini lambda_fin delta_lamb expon \n\n",argv[0]);
    return 1;
  }
  
  expon = atoi(argv[6]);	// N = 10^expon size of original net
  smp = 1;                      // net sample
	
  N=pow(10,expon);
  
  semente = 10*smp;

  kissinit(semente);  

  dgr = (int *) malloc(Nmax * sizeof(int));
  beg = (int *) malloc(Nmax * sizeof(int));
  adj = (int *) malloc(Nmax * sizeof(int));
  S = (int *) malloc(Nmax * sizeof(int));
  DGR = (int *) malloc(Nmax * sizeof(int));
  sigma = (int *) malloc(Nmax * sizeof(int));

  q=0;
  o=0;
  aux=0;
  a=0;
  s=0;
  u=0;

  
  /***************************************************************************/
  // Building the network using UCM 
  // Ref. Catanzaro et al. PRE vol. 71, 2, pp. 1-4, 2025
  /***************************************************************************/
	
  do{

    a=0;
    somadgr=0;
    semente++;
    kissinit(semente);  
	
    do{
            
      k=4;
      a++;
      dgr[a]=k;
      somadgr+=dgr[a];
     

    }while(a<N);
    
    aux = somadgr % 2 ;

  }while(aux != 0);

  L=0;

  for(g=1;g<=N;g++){

    beg[g]=L;
    L=L+dgr[g];
	
  }


  do{

    cont2=0;

    for(s=1;s<=N;s++){

      for(o=beg[s]+1;o<=beg[s]+dgr[s];o++){

	S[o]=s;

      }	
    }


    for(e=1;e<=N;e++){

      DGR[e]=dgr[e]; 
	
    }

	
    for(e=1;e<=L;e++){

      adj[e]=0;
	
    } 

	
    n=L;

    alerta=0;

    links=0;

    do{

      stop=0;

      q=(n)*rand01_kiss(semente)+1;
      u=(n)*rand01_kiss(semente)+1;
      
      if(q!=u && S[q]!=S[u] && DGR[S[u]]!=0 && DGR[S[q]]!=0){

	for(f=beg[S[u]]+1;f<=beg[S[u]]+dgr[S[u]];f++){

	  if(S[q]==adj[f]){

	    stop=1;

	  }

	}


      if(stop!=1){
    		
	adj[beg[S[u]]+DGR[S[u]]]=S[q];  
        DGR[S[u]]--;

	adj[beg[S[q]]+DGR[S[q]]]=S[u]; 
	DGR[S[q]]--;

	cont2=0;
		      
        links++;

	if(q==n){

	  n--;
	  S[u]=S[n];
          n--;

	}else if(u==n){

	  n--;
	  S[q]=S[n];
	  n--;	

	}else if(q!=n && u!=n){
				
	  S[q]=S[n];
          n--;
          S[u]=S[n];
	  n--;

	}

      }

      if((n==2 && stop==1) || (n==0) ){

	alerta=1;		 
			
      }
	
    }else{

      cont2++;
			
      if((n==2 && DGR[S[q]]==0)||(n==2 && DGR[S[u]]==0)||(n==1)||(cont2==1000000 )||(n==2 && S[q]==S[u])||n==0){

	alerta=1;	
	    
      }
    }

    if((n==0) && links<(L/2)){

      n=L;

      alerta=1;

        }
		      
    }while(alerta!=1);		
	
  }while(n>0);

  
  // finished the building of the network
  /******************************************************************************/

  /******************************************************************************/
  // SIS model through the OGA in an intermittent net
  // Ref.: Cota and Ferreira, Comp. Phys. Comm. vol 19, pp. 303-312, 2017.
  // Ref.: Moraes, J. T. and and Ferreira, S. C. (2026). Intermittent quarantine 
  // suppresses epidemic spreading beyond simple contact reduction [Preprint].
  // arXiv. https://doi.org/10.48550/arXiv.2610.00429
  /******************************************************************************/


  double lambda,deltat,t,p,tmax,Z,tmed,mediarho,mediarho2,X,normp,somap;
  double beta,*t_ch,tau_med,tau0,r,lambda_ini,lambda_fin,delta_lambda,*P,prob,t_rlx;
  int Ni,*I,viz,Q,qq,t0,Nn,*EST,N_act,b,i;    
    
  // variables of the epidemic process 

  I = (int *) malloc(Nmax * sizeof(int));
   
  // variables of the intermittency
    
  EST = (int *) malloc(Nmax * sizeof(int));   // EST def. the state of the node EST = 1 the node is active, EST = 0 inactive 
  
  P = (double *) malloc(Nmax * sizeof(double));
  
  t_ch = (double *) malloc(Nmax * sizeof(double));       
  // constants for the inter-event time distribution
    
  beta = atof(argv[1]);         // exponent for the power-law distribution 
    
  tau0 = 0.05;         // constant of the power-law distribution 
    
  tau_med = atof(argv[2]);
    
  r = beta*(tau_med - tau0) + tau0 - 2*tau_med ; 
    
  q=0;
  o=0;
  u=0;
  a=1;
    
  lambda_ini = atof(argv[3]);
  lambda_fin = atof(argv[4]);
  delta_lambda = atof(argv[5]);
  
  lambda = lambda_ini;    // initial value of lambda
    
  do{               //  loop in lambda 
      
    for(i=1;i<=N;i++){
    
      P[i]=0.0;     
    
    }
    
    t_rlx = pow(10,2);
    t0 = t_rlx+5*pow(10,5);     // transient time 
    tmed = 5*pow(10,5);    // time used for averaging the quantities 
    tmax = tmed+t0;
    
    t = 0.0;
    
    //Setting the initial conditions of the intermittent network
    N_act=0;
    
    for(f=1;f<=N;f++){
      
      prob = rand01_kiss(semente);
      	
      if(prob>0.5){EST[f]=1;}else{EST[f]=0;}
      	
      t_ch[f] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r;
		
      N_act+=EST[f];
      
    }
    
    

    s = 0;
    
    sprintf(filename, "sis_rrn_ni_lambdaxrho_beta%.1f_taumed%.1f_N10e%d.dat",beta,tau_med,expon);
    FILE *arq3;  
    arq3 = fopen( filename, "a+");
        
    sprintf(filename, "sis_rrn_ni_lambdaxsuscet_beta%.1f_taumed%.1f_N10e%d.dat",beta,tau_med,expon);
    FILE *arq4;
    arq4 = fopen( filename, "a+");
      
        
    Ni=N;
    Nn=somadgr;			        //  Nn is the sum of the degrees of the infected nodes
        
    for(qq=1;qq<=Ni;qq++){		//  Initial condition for the epidemic model
            
        I[qq]=qq;			//I[N+1] is the list of infected node 
        sigma[qq]=1;		        //sigma is a state list, sigma[i]=1 i is infected, sigma[i]=0 i is susceptible       
      
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

      
    //printf("lambda = %lf \n",lambda);         // if you want to see what is the lambda being computed
    
    
      
    do{     //loop in time
            
      Z=rand01_kiss(semente);
        
      p=(lambda*Nn)*pow((Ni+(lambda*Nn)),-1);
        
      deltat=-log(rand01_kiss(semente)+pow(10,-8))*pow((Ni+(lambda*Nn)),-1);
            
        
      if(t>t0){
        
        P[Ni]+=deltat;     
          
      }

                
      if(Z<p){		// Infection event
      
        Q=(Ni)*rand01_kiss(semente)+1;		

        q=I[Q];          	      
   	            
        if(t>t_ch[q]){         // it means the node q should already change its state
		        
          do{
                
            if(EST[q]==1){    // EST[q]==1 : q is active 	
                            
              EST[q]=0;                
                             
            }else{          // q is inactive 
                            
              EST[q]=1;
              
            } 
   
            t_ch[q] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[q];
               
          }while(t_ch[q]<t); 		  
        }
            
        if(EST[q]==1){		// the node is active, so the infection occurs 
  
          b=(dgr[q])*rand01_kiss(semente)+1;    //randomly choosing between neighbors of q
      	  viz = adj[beg[q]+b];
      	  
   	  if(t>t_ch[viz]){         // it means the node viz should already change its state
		        
            do{
                
              if(EST[viz]==1){    // EST[viz]==1 : viz is active 	
                            
                EST[viz]=0;                  
                                          
              }else{          // viz is inactive 
                            
                EST[viz]=1;

              } 
   
              t_ch[viz] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[viz];
               
            }while(t_ch[viz]<t); 		  
          }
      		
          if((sigma[viz]==0) && (EST[viz]==1)){   // if the neighbor is active and susceptible the infection occurs
 
            Ni++;
            I[Ni]=viz;
            Nn=Nn+dgr[viz];
            sigma[viz]=1;
             
          }
                            
        }
                
      }else{ 		//q becomes susceptible
            
        Q=(Ni)*rand01_kiss(semente)+1;		
	q=I[Q];	
	
        if(t>t_ch[q]){         // it means the node q should already change its state
		        
          do{
                
            if(EST[q]==1){    // EST[q]==1 : q is active 	
                            
              EST[q]=0;                
                                  
            }else{          // q is inactive 
                            
              EST[q]=1;
              
            } 
   
            t_ch[q] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[q];
               
          }while(t_ch[q]<t); 		  
        }
                
        sigma[q]=0;
        I[Q]=I[Ni];  
        Ni--;
        Nn=Nn-dgr[q];

        if(Ni==0){		//absorbing state
	
	  Ni++;	
	  Nn=Nn+dgr[q];
	  I[Ni]=q;
	  sigma[q]=1;
	       	  
        } 
      }
  
      t+=deltat;

    }while(t<tmax);
        
    mediarho=0.0;
    mediarho2=0.0;  
    somap=0.0;
      
    for(i=1;i<=N;i++){
      
      somap+=P[i];          
    
    }
      
      
    normp = pow(somap,-1);   
        
    for(i=1;i<=N;i++){

      mediarho+=i*normp*P[i];
      mediarho2+=pow(i,2)*normp*P[i];

    }
     
    X=(mediarho2-pow(mediarho,2))*pow(mediarho,-1);
        
        
    fprintf(arq3,"%lf %e %e \n",lambda,mediarho,mediarho2);
    fprintf(arq4,"%lf %e \n",lambda,X);
        
    a++;
        
    lambda+=delta_lambda;
         
    fclose(arq3);
    fclose(arq4);

  }while(lambda<=lambda_fin);
    

    free(I);
    free(sigma);
    free(EST);
    free(dgr);
    free(beg);
    free(adj);
    free(P);
    free(t_ch);
    free(S);
    free(DGR);
  
    return (0);
}

