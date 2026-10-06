// ==============================================================================
// File Name: sir_ucm_ni.c
// Author: Juliane T. de Moraes
// Updated: 2026-09
// License: MIT
// Usage: ./sir_ucm_ni gamma beta lambda_ini lambda_fin delta_lamb expon smp taumed
// Description: This code computes the variability as function of infection rate for 
// the SIR model in a node intermittent network with power-law degree distribution.
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

  int N,L,f,e,*adj,a,somadgr,n,g,*S,q,o,stop,*beg,*dgr,*DGR,cont2,alerta,aux,*sigma,s,u,*hub,rr,num,kmax,links,v,ra;
  double gamma,epsilon,*p_k,somap_k,C,k,AA,st_moment,nd_moment,k0,kc,pp,aa;
  int semente;
  char filename[4*64];
	
  int expon,smp;
  
  if (argc != 9){
    printf("\n Usage: %s  gamma beta lambda_ini lambda_fin delta_lamb expon smp taumed \n\n",argv[0]);
    return 1;
  }
  
  expon = atoi(argv[6]);	// N = 10^expon: size of original net
  smp = atoi(argv[7]);          // net sample
  N=pow(10,expon);
  semente = 10*smp;

  kissinit(semente);  

  dgr = (int *) malloc(Nmax * sizeof(int));
  beg = (int *) malloc(Nmax * sizeof(int));
  adj = (int *) malloc(Nmax * sizeof(int));
  S = (int *) malloc(Nmax * sizeof(int));
  DGR = (int *) malloc(Nmax * sizeof(int));
  sigma = (int *) malloc(Nmax * sizeof(int));
  hub = (int *) malloc(Nmax * sizeof(int));
  p_k = (double *) malloc(Nmax * sizeof(double));

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
        	      
  gamma = atof(argv[1]);
  k0 = 3.0;
  kc=2*sqrt(N);                             			    
  C = 1.0/(1.0-gamma);
  
  // this two files will give info about the network backbone:
  
  sprintf(filename, "sis_ucm%.2f_parameters_N10e%d_net%d.dat",gamma,expon,smp);
  FILE *arq1;  
  arq1 = fopen( filename, "a+");	
  
  sprintf(filename, "sis_ucm%.2f_degreedistrib_N10e%d_net%d.dat",gamma,expon, smp);
  FILE *arq2;  
  arq2 = fopen( filename, "a+");
	
  do{

    a=0;
    somadgr=0;
    semente++;
    kissinit(semente);  
	
    do{
	
      epsilon=rand01_kiss(semente);      
      k = k0*pow((1-(1-pow((k0*pow(kc,-1)),(gamma-1)))*epsilon),C);    
	
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

  kmax=1;
  num=0;

  for(rr=1;rr<=N;rr++){

    if(dgr[rr]>kmax){ 

      kmax=dgr[rr];

    }
  }


  for(rr=1;rr<=N;rr++){

    if(dgr[rr]==kmax){   

      num++;

      hub[num]=rr;
      
      fprintf(arq1,"hub: %d \n",rr);
	
    }
  }
  
  for(f=1;f<=N;f++){    

    p_k[f]=0;    
	
  }

  for(e=1;e<=N;e++){  

    p_k[dgr[e]]++;					//creating the distribution

  }

  somap_k=0.0;

  for(e=1;e<=N;e++){

    somap_k+=p_k[e];

  }

  AA=pow(somap_k,-1);
  

  for(e=1;e<=N;e++){  

    if(p_k[e]!=0){		
		
      fprintf(arq2,"%d %e \n",e,AA*p_k[e]);  // normalized degree distribution

    }
  }

  st_moment=0.0;
  nd_moment=0.0;

  for(v=k0;v<=kmax;v++){

    st_moment+=v*AA*p_k[v];
    nd_moment+=(v*v)*AA*p_k[v];

  }

  fprintf(arq1,"N = %d \nst moment = %lf \nnd moment = %lf \nkmax = %d \n# hubs = %d \n",N,st_moment,nd_moment,kmax,num);

  fclose(arq1);
  fclose(arq2);
  
  // finished the building of the network
  /******************************************************************************/

  
  /******************************************************************************/
  // SIR model through the OGA in an intermittent net
  // Ref.: Cota and Ferreira, Comp. Phys. Comm. vol 19, pp. 303-312, 2017.
  // Ref.: Moraes, J. T. and and Ferreira, S. C. (2026). Intermittent quarantine 
  // suppresses epidemic spreading beyond simple contact reduction [Preprint].
  // arXiv. https://doi.org/10.48550/arXiv.2610.00429
  /******************************************************************************/
	
  double lambda,deltat,t,p,rho_rec,Z,somarho_rec,somarho2_rec,*mediarho_rec,*mediarho2_rec,*Delta,*mediat,somat;
  double beta,tau0,*t_ch,lambda_ini,lambda_fin,delta_lambda,tau_med,r,prob;
  int Ni,*I,Q,b,Nn,ams,ams_max,ini,Nr,*R,i,*EST,viz;     
	
  // variables of the intermittency
    
  EST = (int *) malloc(Nmax * sizeof(int));   // if EST = 1 the node is active, EST = 0 inactive 
  t_ch = (double *) malloc(Nmax * sizeof(double));       
  
  
  // constants for the inter-event time distribution
    
  beta = atof(argv[2]);         // exponent for the power-law distribution 
    
  tau0 = 0.05;         // constant of the power-law distribution 
    
  tau_med = atof(argv[8]);
    
  r = beta*(tau_med - tau0) + tau0 - 2*tau_med ; 
    
  o=0;
  u=0;
  b=0;
  a=1;
  q=0;
  aux=0;
  s=0;

  I = (int *) malloc(Nmax * sizeof(int));   // infected nodes' list
  R = (int *) malloc(Nmax * sizeof(int));   // recovered nodes' list
  mediarho_rec = (double *) malloc(Nmax * sizeof(double));
  mediarho2_rec = (double *) malloc(Nmax * sizeof(double));	
  Delta = (double *) malloc(Nmax * sizeof(double));		// Delta is the epidemic variability
  mediat = (double *) malloc(Nmax * sizeof(double));

  lambda_ini = atof(argv[3]);
  lambda_fin = atof(argv[4]);
  delta_lambda = atof(argv[5]);
  
  lambda = lambda_ini;    // initial value of lambda
    
			
  ams_max = pow(10,5);		// how many times the absorbing state is reached 


  for(i=1;i<=N;i++){ 		// all the nodes, except for one, begin susceptible

    sigma[i] = 0;

  }

  do{			//loop in lambda
  
	
    sprintf(filename, "sir_ucm_ni_2-75_beta%.1f_lambdaxrho-rec_N10e%d_net%d_taumed%.1f_trlx0.dat",beta,expon, smp,tau_med);
    FILE *arq2;         
    arq2 = fopen( filename, "a+");
        
    sprintf(filename, "sir_ucm_ni_2-75_beta%.1f_lambdaxDelta_N10e%d_net%d_taumed%.1f_trlx0.dat",beta,expon, smp,tau_med);
    FILE *arq3;         
    arq3 = fopen( filename, "a+");
        
    sprintf(filename, "sir_ucm_ni_2-75_beta%.1f_lambdaxtmax-med_N10e%d_net%d_taumed%.1f_trlx0.dat",beta,expon, smp,tau_med);
    FILE *arq4;         
    arq4 = fopen( filename, "a+");


    ams = 1; 		
  
    somarho_rec = 0.0;	
    somarho2_rec = 0.0;
    somat = 0.0;
    double t_rlx = pow(10,2);

    t=0.0;
	
	
    do{             //loop samples
    
    
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
    
      ini = (N)*rand01_kiss(semente)+1;	// randomly choose a node to infect

		
      Ni = 1;		        
      I[1] = ini;		
      sigma[ini] = 1; 
      Nn = dgr[ini];		

      Nr=0;				

      do{                 //loop sir: Ni!=0

        Z = rand01_kiss(semente);
        p = lambda*Nn*pow((Ni+lambda*Nn),-1);

        if(Z<p){		// An infection will occur

	  do{			// choosing a node proportionally to its degree
	    
	    ra=0;

	    Q=(Ni)*rand01_kiss(semente)+1;		
	    q=I[Q];	

	    pp = rand01_kiss(semente);
	    aa = dgr[q]*pow(kmax,-1);	

	    if(pp<aa){ra=1;}
	
	  }while(ra!=1);
	    
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
        
            
          if(EST[q]==1){		         
  
            b=(dgr[q])*rand01_kiss(semente)+1;    
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
	    
            if((sigma[viz]==0) && (EST[viz]==1)){   

              Ni++;
              I[Ni]=viz;
              Nn=Nn+dgr[viz];
              sigma[viz]=1;
	   
	    }
          }			
        
        }else{		                          // A recovery will occur

          Q=(Ni)*rand01_kiss(semente)+1;		//randomly choose between the infected nodes
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

	  sigma[q]=2;			// recover q
	  I[Q]=I[Ni];
	  Ni--;
	  Nn=Nn-dgr[q];

	  Nr++;				// increase the number of recovered   
	  R[Nr]=q;			// place it in a recovered nodes' list
	
        } 
	
        if(Ni!=0){

          deltat=-log(rand01_kiss(semente)+pow(10,-8))*pow((Ni+(lambda*Nn)),-1);
	  t+=deltat;
		
        }
      }while(Ni!=0);  	
	
	
      rho_rec = Nr*pow(N,-1);			
      somarho_rec+=rho_rec; 	
      somarho2_rec+=rho_rec*rho_rec;
      somat+=t;
	
      // resetting the initial conditions      

      for(i=1;i<=Nr;i++){		

        sigma[R[i]]=0;			

      }	

      ams++;
    
    }while(ams<=ams_max);		

    mediat[a] = somat*pow(ams_max,-1); 
    mediarho_rec[a] = somarho_rec*pow(ams_max,-1); 		
    mediarho2_rec[a] = somarho2_rec*pow(ams_max,-1); 	

    Delta[a] = pow((mediarho2_rec[a] - (mediarho_rec[a]*mediarho_rec[a])),0.5)*pow(mediarho_rec[a],-1);

    fprintf(arq2,"%lf %e \n",lambda,mediarho_rec[a]);	
    fprintf(arq3,"%lf %e \n",lambda,Delta[a]);	
    fprintf(arq4,"%lf %e \n",lambda,mediat[a]);	

    a++;						

    lambda+=delta_lambda;


    fclose(arq2);
    fclose(arq3);
    fclose(arq4);

  }while(lambda<=lambda_fin);


  free(I);
  free(R);
  free(sigma);	
  free(mediarho_rec);
  free(mediarho2_rec);
  free(mediat);
  free(Delta);
  free(EST);
  free(t_ch);
  free(hub);
  free(DGR);
  free(S);
  free(dgr);
  free(beg);
  free(p_k);
  free(adj);

  return (0);
}


