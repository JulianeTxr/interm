// ==============================================================================
// File Name: sis_ucm_li.c
// Author: Juliane T. de Moraes
// Updated: 2026-09
// License: MIT
// Usage: ./sis_ucm_li gamma beta lambda_ini lambda_fin delta_lamb expon smp taumed
// Description: This code computes the susceptibility and the quasi-stationary 
// average of prevalence as function of infection rate for the SIS model in a link
// intermittent network with power-law degree distribution.
// See Ref.: referencia do artigo no arxiv
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

  int N,L,f,e,*adj,a,somadgr,n,g,*S,q,o,stop,*beg,*dgr,*DGR,cont2,alerta,aux,*sigma,s,u,*hub,rr,num,kmax,links,v,h;
  double gamma,epsilon,*p_k,somap_k,C,k,AA,st_moment,nd_moment,k0,kc,change,prob;
  int semente;
  char filename[4*64];
	
  int expon,smp;
  
  if (argc != 9){
    printf("\n Usage: %s  gamma beta lambda_ini lambda_fin delta_lamb expon smp taumed \n\n",argv[0]);
    return 1;
  }
  
  expon = atoi(argv[6]);	// N = 10^expon size of original net
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
      
      if(k>kc){printf("k = %lf \n",k);}
	
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
  // SIS model through the OGA in an intermittent net
  // Ref.: Cota and Ferreira, Comp. Phys. Comm. vol 19, pp. 303-312, 2017.
  // Ref.: incluir ref arxiv paper interm. 
  /******************************************************************************/
  
  
  double lambda,deltat,t,p,tmax,Z,tmed,mediarho,mediarho2,X,normp,somap,t_rlx;
  double beta,*t_ch,tau_med,tau0,r,lambda_ini,lambda_fin,delta_lambda,*P,aa,pp;
  int Ni,*I,viz,Q,qq,t0,Nn,*EST,L_act,b,i,aresta,ra;    
    
  // variables of the epidemic process 

  I = (int *) malloc(Nmax * sizeof(int));
   
  // variables of the intermittency
    
  EST = (int *) malloc(Nmax * sizeof(int));   // EST def. the state of the node EST = 1 the node is active, EST = 0 inactive 
  
  P = (double *) malloc(Nmax * sizeof(double));
  
  t_ch = (double *) malloc(Nmax * sizeof(double));  
       
  // constants for the inter-event time distribution
    
  beta = atof(argv[2]);         // exponent for the power-law distribution 
    
  tau0 = 0.05;         // constant of the power-law distribution 
    
  tau_med = atof(argv[8]);
    
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
    t0 = t_rlx+5*(pow(10,5));     // transient time 
    tmed = 5*pow(10,5);           // time used for averaging the quantities 
    tmax = tmed+t0;
    
    t = 0.0;
    
    //distributing the time of the first change of each link, according to the inter-event time distribution
    
    for(i=1;i<=N;i++){
    
      for(f=beg[i]+1;f<=beg[i]+dgr[i];f++){
      	
      	prob = rand01_kiss(semente);
      	
      	if(prob>0.5){EST[f]=1;}else{EST[f]=0;}
      	
	change = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r;
	t_ch[f]=change;
        
        for(g=beg[adj[f]]+1;g<=beg[adj[f]]+dgr[adj[f]];g++){
          
          if(adj[g]==i){ 
            
            EST[g]=EST[f];
            t_ch[g]=change;
          
          } 
        }	
      }
    }
    
    L_act=0;
        
    for(qq=1;qq<=L;qq++){
    
    L_act+=EST[qq];
    
    }
    
    // Iterating the network before the epidemics for a time = t_rlx
    
    for(q=1;q<=L;q++){
      
      if(t_rlx>t_ch[q]){         // it means the link q should already change its state
      	  
        do{
      	  
          if(EST[q]==1){
          
            EST[q]=0;                
                                               
          }else{EST[q]=1;} 
   
          t_ch[q]+= (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r;
               
        }while(t_ch[q]<t_rlx); 		  
      }
    }
      
    s = 0;
    
    aresta = 0;
    
    sprintf(filename, "sis_ucm_li_gamma%.2f_lambdaxrho_beta%.1f_taumed%.1f_N1e%d.dat",gamma,beta,tau_med,expon);
    FILE *arq3;  
    arq3 = fopen( filename, "a+");
        
    sprintf(filename, "sis_ucm_li_ucm_gamma%.2f_lambdaxsuscet_beta%.1f_taumed%.1f_N1e%d.dat",gamma,beta,tau_med,expon);
    FILE *arq4;
    arq4 = fopen( filename, "a+");
      
        
    Ni=N;
    Nn=somadgr;			        //  Nn is the sum of the degrees of the infected nodes
        
    for(qq=1;qq<=Ni;qq++){		//  Initial condition for the epidemic model
            
        I[qq]=qq;			//I[N+1] is the list of infected node 
        sigma[qq]=1;		        //sigma is a state list, sigma[i]=1 i is infected, sigma[i]=0 i is susceptible       
      
    }
      
    //printf("lambda = %lf \n",lambda);         // if you want to see what is the lambda being computed
      
    do{                                 //loop in time
            
      Z=rand01_kiss(semente);
        
      p=(lambda*Nn)*pow((Ni+(lambda*Nn)),-1);
        
      deltat=-log(rand01_kiss(semente)+pow(10,-8))*pow((Ni+(lambda*Nn)),-1);
            
        
      if(t>t0){
        
        P[Ni]+=deltat;     
          
      }

                
      if(Z<p){		// Infection event
                	      
        do{			//choosing proportionally to degree

	  ra=0;

	  Q=(Ni)*rand01_kiss(semente)+1;		

	  q=I[Q];	

	  pp = rand01_kiss(semente);

	  aa = dgr[q]*pow(kmax,-1);	

	  if(pp<aa){

	    ra=1;

	  }
	
	}while(ra!=1);
				
	b=(dgr[q])*rand01_kiss(semente)+1; 	
        viz = adj[beg[q]+b];          
   	                  	  
        if(sigma[viz]==0){ 
                    
        aresta = beg[q]+b;
                       
        if(t>t_ch[aresta]){      
                    	
          do{
                                                           
            if(EST[aresta]==1){	
                            
              EST[aresta]=0;
                                                                                                                      
            }else{		
                            
              EST[aresta]=1;
                  
            }
                        
            t_ch[aresta] = (tau0+r)*pow((1-rand01_kiss(semente)),(pow((-beta+1),-1))) - r + t_ch[aresta];
                      
          }while(t_ch[aresta]<t);      
            
          for(g=beg[adj[aresta]]+1;g<=beg[adj[aresta]]+dgr[adj[aresta]];g++){         // the ends of the link are found in the adjacency list
          
            if(adj[g]==q){ 
              
              EST[g] = EST[aresta];           
              t_ch[g] = t_ch[aresta];
            
              
            }     
          }                    
        }
      	  
        if(EST[aresta]==1){         //if the link is active the infection occurs
        
          Ni++;
                        
          I[Ni]=viz;
                        
          Nn=Nn+dgr[viz];
                        
          sigma[viz]=1;
                        
        }
                            
      }
                
    }else{ 		//q becomes susceptible
            
        Q=(Ni)*rand01_kiss(semente)+1;		
	q=I[Q];	
	         
        sigma[q]=0;
        I[Q]=I[Ni];  
        Ni--;
        Nn=Nn-dgr[q];

        if(Ni==0){		//absorbing state - hub reactivation
	
	  h=(num)*rand01_kiss(semente)+1; 	  
  	  Ni++;	
	  Nn=Nn+dgr[hub[h]];
	  I[Ni]=hub[h];
	  sigma[hub[h]]=1;
	       	  
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
      mediarho2+=(i*i)*normp*P[i];

    }
     
    X=(mediarho2-(mediarho*mediarho))/(mediarho);
        
        
    fprintf(arq3,"%lf %e \n",lambda,mediarho);
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

    free(p_k);
    free(S);
    free(DGR);
    free(hub);

    return (0);
}

