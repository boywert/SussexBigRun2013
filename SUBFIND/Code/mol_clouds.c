#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "proto.h"
#include "allvars.h"

#ifdef MOL_CLOUDS


/* lifetime/scatter of MCs in internal time units */
#define MOL_CLOUDS_LIFETIME_SCATTER 0.0*MOL_CLOUDS_LIFETIME   // SCATTER OPTION UNTESTED 
/* trace this MC particle (index is internal MC ID) */
#define MCTRACER_INDEX 0

/* correction factor for correct mean mass fraction */
//#define MOL_CLOUDS_MASS_CORRECTION                          // ONLY DELTA_MASS TESTED
/* do on/off mass instead of Gaussian growth */
#define DELTA_MASS

#ifdef MOL_CLOUDS_MASS_CORRECTION
 #define KAPPA 0.855624
#else
 #define KAPPA 1.0
#endif


/* sort structures */
static struct sort_data
{
  MyIDType ID;  
  int ind;
}
 *sortdata, *sortdatasmall;

static struct sort_dataD
{
  MyIDType ID1;
  MyIDType ID2;  
  int ind;
}
 *sortdata1, *sortdata2;


/* sort functions */
int compare_sort_data(const void *a, const void *b)
{
  if(((struct sort_data *) a)->ID < ((struct sort_data *) b)->ID)
    return -1;

  if(((struct sort_data *) a)->ID > ((struct sort_data *) b)->ID)
    return +1;

  return 0;
}

/* compare functions */
int compare_sort_data1(const void *a, const void *b)
{
  if(((struct sort_dataD *) a)->ID1 < ((struct sort_dataD *) b)->ID1)
    return -1;

  if(((struct sort_dataD *) a)->ID1 > ((struct sort_dataD *) b)->ID1)
    return +1;

  return 0;
}

/* compare functions */
int compare_sort_data2(const void *a, const void *b)
{
  if(((struct sort_dataD *) a)->ID2 < ((struct sort_dataD *) b)->ID2)
    return -1;

  if(((struct sort_dataD *) a)->ID2 > ((struct sort_dataD *) b)->ID2)
    return +1;

  return 0;
}

/* mass of MC as a function of time */
double mass_function(double t, double t_born, double t_lifetime)
{
#ifdef DELTA_MASS
  return All.MOL_CLOUDS_MassMCs;
#else
 double t_m = t_born + t_lifetime/2.0;
 
 return All.MOL_CLOUDS_m1_tilde + All.MOL_CLOUDS_m2_tilde * exp(-(t-t_m)*(t-t_m)/((t_lifetime/sqrt(2.0))*(t_lifetime/sqrt(2.0))));
#endif 
}

/* generate uniformly distributed random lifetimes for the MCs */
double generate_lifetime(MyIDType id)
{
 return 2.0*MOL_CLOUDS_LIFETIME_SCATTER*(drand48()-0.5)  +  MOL_CLOUDS_LIFETIME;
 //return 2.0*MOL_CLOUDS_LIFETIME_SCATTER*(get_random_number(id)-0.5)  +  MOL_CLOUDS_LIFETIME;
}


/* INIT CLOUD DATA */
void init_mol_clouds()
{
 int i;
 double mean_lifetime_local=0.0, mean_lifetime_global=0.0, variance_lifetime_local=0.0, variance_lifetime_global=0.0;
 double total_mc_mass_local=0.0, total_mc_mass_global=0.0, total_star_mass_local=0.0, total_star_mass_global=0.0;
 double initial_timeborn[header.npartTotal[3]], initial_lifetime[header.npartTotal[3]];
 
 All.MOL_CLOUDS_NumMCs = header.npartTotal[3];
 All.MOL_CLOUDS_NumStars = header.npartTotal[2];

 if (header.mass[3]!=0.0)
  All.MOL_CLOUDS_MassMCs = header.mass[3];
 else
  {
   double mc_mass_local=0.0, mc_mass_global=0.0;
   for (i=0; i<NumPart; i++)
    {
     if (P[i].Type==3)
      mc_mass_local=P[i].Mass;
    } 
   MPI_Allreduce(&mc_mass_local, &mc_mass_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   All.MOL_CLOUDS_MassMCs = mc_mass_global;
  }
 All.MOL_CLOUDS_MassStars = header.mass[2];
 
 All.MOL_CLOUDS_m1_tilde = (KAPPA*All.MOL_CLOUDS_MassStars-All.MOL_CLOUDS_MassMCs*exp(-0.5))/(KAPPA-exp(-0.5));
 All.MOL_CLOUDS_m2_tilde = (All.MOL_CLOUDS_MassMCs-All.MOL_CLOUDS_MassStars)/(KAPPA-exp(-0.5));

 
 srand48(230381);


 /* let every task generate the same random numbers for MC startup */
 for (i=0; i<All.MOL_CLOUDS_NumMCs; i++)
   {
    initial_lifetime[i]=generate_lifetime(i);                         /* assign initial lifetime */
    //initial_timeborn[i]=-get_random_number(i)*initial_lifetime[i];  /* cloud was born some time in the past within its lifetime */
    initial_timeborn[i]=-drand48()*initial_lifetime[i];               /* cloud was born some time in the past within its lifetime */
   }

 for (i=0; i<NumPart; i++)
   { 
    /* shift IDs to start at 1 (depends on ICs)*/ 
    P[i].ID++; 

    if (P[i].Type==2)    
     total_star_mass_local+=P[i].Mass;
    
    if (P[i].Type==3)
     {
      /* MC */
      P[i].MOL_CLOUDS_index=P[i].ID-All.MOL_CLOUDS_NumStars-1; 
      P[i].MOL_CLOUDS_LifeTime = initial_lifetime[P[i].MOL_CLOUDS_index];
      P[i].MOL_CLOUDS_TimeBorn = initial_timeborn[P[i].MOL_CLOUDS_index];
      P[i].Mass=mass_function(All.Time, P[i].MOL_CLOUDS_TimeBorn, P[i].MOL_CLOUDS_LifeTime); /* calculate cloud mass based on time of birth and age */
      mean_lifetime_local+=P[i].MOL_CLOUDS_LifeTime/All.MOL_CLOUDS_NumMCs;
      total_mc_mass_local+=P[i].Mass;
     }
    else
     {
      /* star */
      P[i].MOL_CLOUDS_LifeTime = 0.0;     
      P[i].MOL_CLOUDS_TimeBorn = 0.0;           
      P[i].MOL_CLOUDS_index=-1;  
     } 
   }

  MPI_Allreduce(&mean_lifetime_local, &mean_lifetime_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (i=0; i<NumPart; i++)
   {
    if (P[i].Type!=3)
     continue;
    variance_lifetime_local+=(P[i].MOL_CLOUDS_LifeTime-mean_lifetime_global)*(P[i].MOL_CLOUDS_LifeTime-mean_lifetime_global)/All.MOL_CLOUDS_NumMCs;
   } 
 
  MPI_Allreduce(&total_mc_mass_local, &total_mc_mass_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_star_mass_local, &total_star_mass_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
  MPI_Allreduce(&variance_lifetime_local, &variance_lifetime_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 
  if (ThisTask==0)
   {
    printf(" MOL_CLOUDS_NumMCs=%d\n", All.MOL_CLOUDS_NumMCs);
    printf(" MOL_CLOUDS_NumStars=%d\n", All.MOL_CLOUDS_NumStars);    
    printf(" MOL_CLOUDS_MassMCs=%g\n", All.MOL_CLOUDS_MassMCs);        
    printf(" MOL_CLOUDS_MassStars=%g\n", All.MOL_CLOUDS_MassStars);            
    printf(" MOL_CLOUDS_m1_tilde=%g\n", All.MOL_CLOUDS_m1_tilde);        
    printf(" MOL_CLOUDS_m2_tilde=%g\n", All.MOL_CLOUDS_m2_tilde);            
    printf(" KAPPA=%g\n", KAPPA);                
    printf(" mean_lifetime=%g (%g)\n", mean_lifetime_global, MOL_CLOUDS_LIFETIME);     
    printf(" stdev_lifetime=%g (%g)\n", sqrt(variance_lifetime_global), MOL_CLOUDS_LIFETIME_SCATTER);         
    printf(" mass_function(0, 0, MOL_CLOUDS_LIFETIME)=%g\n", mass_function(0, 0, MOL_CLOUDS_LIFETIME));            
    printf(" mass_function(MOL_CLOUDS_LIFETIME/2, 0, MOL_CLOUDS_LIFETIME)=%g\n", mass_function(MOL_CLOUDS_LIFETIME/2, 0, MOL_CLOUDS_LIFETIME));                
    printf(" mass_function(MOL_CLOUDS_LIFETIME, 0, MOL_CLOUDS_LIFETIME)=%g\n", mass_function(MOL_CLOUDS_LIFETIME, 0, MOL_CLOUDS_LIFETIME));                    
    printf(" time = %g  total mass in MCs = %g ( fraction = %g )   total mass in stars = %g ( fraction = %g )\n", 
             All.Time, 
	     total_mc_mass_global, total_mc_mass_global/(total_mc_mass_global + total_star_mass_global), 
 	     total_star_mass_global, total_star_mass_global/(total_mc_mass_global + total_star_mass_global));

   }    
}




/* CREATE / DESTROY CLOUDS */
int do_mol_clouds(void)
{
 int i, j;
 int ID1_List_local[All.MOL_CLOUDS_NumMCs], ID1_List_global[All.MOL_CLOUDS_NumMCs];
 int ID2_List_local[All.MOL_CLOUDS_NumMCs], ID2_List_global[All.MOL_CLOUDS_NumMCs]; 
 double lifetimes_local[All.MOL_CLOUDS_NumMCs], lifetimes_global[All.MOL_CLOUDS_NumMCs];
 double dyingtimes_local[All.MOL_CLOUDS_NumMCs], dyingtimes_global[All.MOL_CLOUDS_NumMCs];

 int  minID_local=All.MOL_CLOUDS_NumMCs+All.MOL_CLOUDS_NumStars, minID_global=All.MOL_CLOUDS_NumMCs+All.MOL_CLOUDS_NumStars, maxID_local=0, maxID_global=0;
 int converted_stars_to_MCs_local=0, converted_stars_to_MCs_global=0, converted_MCs_to_stars_local=0, converted_MCs_to_stars_global=0;  
 int MCcount_local=0, MCcount_global=0, MCs_to_convert_local=0, MCs_to_convert_global=0;
 double total_mc_mass_local=0.0, total_mc_mass_global=0.0 ,total_star_mass_local=0.0, total_star_mass_global=0.0;
 double tracer_mass_local=0.0, tracer_mass_global=0.0; 
 double tracer_pos_local[]={0.0,0.0,0.0}, tracer_pos_global[]={0.0,0.0,0.0};
 double t_born_min_local=All.TimeMax, t_born_min_global=All.TimeMax, t_born_max_local=All.TimeBegin, t_born_max_global=All.TimeBegin; 
 double t_lifetime_min_local=All.TimeMax, t_lifetime_min_global=All.TimeMax;
 double t_lifetime_max_local=All.TimeBegin, t_lifetime_max_global=All.TimeBegin; 
 int found_tracer_local=0, found_tracer_global=0;
 int flag;
 int select = All.MOL_CLOUDS_NumMCs+1, app=0;

 

 /* update mass of clouds */ 
 for (i=0; i<NumPart; i++)
   {
    if (P[i].ID>maxID_local)
     maxID_local=P[i].ID;

    if (P[i].ID<minID_local)
     minID_local=P[i].ID;
 
    if (P[i].Type==2)    
     total_star_mass_local+=P[i].Mass;
     
    if (P[i].Type!=3)
     continue;

    P[i].Mass=mass_function(All.Time, P[i].MOL_CLOUDS_TimeBorn, P[i].MOL_CLOUDS_LifeTime);
    total_mc_mass_local+=P[i].Mass;
    
    if (P[i].MOL_CLOUDS_TimeBorn>t_born_max_local)
     t_born_max_local=P[i].MOL_CLOUDS_TimeBorn;

    if (P[i].MOL_CLOUDS_TimeBorn<t_born_min_local)
     t_born_min_local=P[i].MOL_CLOUDS_TimeBorn;

    if (P[i].MOL_CLOUDS_LifeTime>t_lifetime_max_local)
     t_lifetime_max_local=P[i].MOL_CLOUDS_LifeTime;

    if (P[i].MOL_CLOUDS_LifeTime<t_lifetime_min_local)
     t_lifetime_min_local=P[i].MOL_CLOUDS_LifeTime;

    if (P[i].MOL_CLOUDS_index==MCTRACER_INDEX)
     {
      tracer_mass_local=P[i].Mass;
      tracer_pos_local[0]=P[i].Pos[0];
      tracer_pos_local[1]=P[i].Pos[1];
      tracer_pos_local[2]=P[i].Pos[2];
      found_tracer_local++;
     }

   } 

 /* get statistics */
 MPI_Allreduce(&total_mc_mass_local, &total_mc_mass_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(&total_star_mass_local, &total_star_mass_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
 MPI_Allreduce(&tracer_mass_local, &tracer_mass_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
 MPI_Allreduce(&tracer_pos_local[0], &tracer_pos_global[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(&found_tracer_local, &found_tracer_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(&t_born_max_local, &t_born_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);   
 MPI_Allreduce(&t_born_min_local, &t_born_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);    
 MPI_Allreduce(&t_lifetime_max_local, &t_lifetime_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);   
 MPI_Allreduce(&t_lifetime_min_local, &t_lifetime_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);    
 MPI_Allreduce(&maxID_local, &maxID_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 MPI_Allreduce(&minID_local, &minID_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

 if (ThisTask==0)
  {
   printf(" time = %g  total mass in MCs = %g ( fraction = %g )  total mass in stars = %g ( fraction = %g )  t_born_diff = %g ( %g %g )  t_lifetime_diff = %g ( %g %g )\n", 
            All.Time, 
	    total_mc_mass_global, total_mc_mass_global/(total_mc_mass_global + total_star_mass_global), 
	    total_star_mass_global, total_star_mass_global/(total_mc_mass_global + total_star_mass_global), 
	    t_born_max_global-t_born_min_global, t_born_min_global, t_born_max_global, 
	    t_lifetime_max_global-t_lifetime_min_global, t_lifetime_min_global, t_lifetime_max_global);	    
   printf(" minID_global=%d  maxID_global=%d\n", minID_global, maxID_global);
   printf(" mctracer %g %g  %g  %g  %g %d\n", All.Time, tracer_mass_global, tracer_pos_global[0], tracer_pos_global[1], tracer_pos_global[2], found_tracer_global);
  }

 /* zero arrays */
 for (j=0; j<All.MOL_CLOUDS_NumMCs; j++)
   {
    ID1_List_local[j]=0; ID1_List_global[j]=0;
    ID2_List_local[j]=0; ID2_List_global[j]=0;
    lifetimes_local[j]=0; lifetimes_global[j]=0;
    dyingtimes_local[j]=0; dyingtimes_global[j]=0;
   }

  /* search for MCs that need to be removed */ 
  for (i=0; i<NumPart; i++)
    {
     if (P[i].ID>All.MOL_CLOUDS_NumStars)
       MCcount_local++;
       
     if (P[i].Type!=3)
      continue;

     dyingtimes_local[P[i].MOL_CLOUDS_index]=P[i].MOL_CLOUDS_TimeBorn + P[i].MOL_CLOUDS_LifeTime;     

     if (All.Time - P[i].MOL_CLOUDS_TimeBorn > P[i].MOL_CLOUDS_LifeTime)
      {
       ID2_List_local[P[i].MOL_CLOUDS_index]=P[i].ID;        
       MCs_to_convert_local++;
      }
    } 

 MPI_Allreduce(&MCs_to_convert_local, &MCs_to_convert_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(&MCcount_local, &MCcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(ID2_List_local, ID2_List_global, All.MOL_CLOUDS_NumMCs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(dyingtimes_local, dyingtimes_global, All.MOL_CLOUDS_NumMCs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 /* no MCs to convert, do nothing and return */
 if (MCs_to_convert_global==0)
  return 0;

 if (ThisTask==0)
  {
   printf(" identified clouds to destroy\n");
   fflush(stdout);
  }


 sortdata = (struct sort_data *) mymalloc("sortdata", sizeof(struct sort_data) * NumPart);
 sortdatasmall = (struct sort_data *) mymalloc("sortdata", sizeof(struct sort_data) * All.MOL_CLOUDS_NumMCs); 
 sortdata1 = (struct sort_dataD *) mymalloc("sortdata1", sizeof(struct sort_dataD) * All.MOL_CLOUDS_NumMCs); 
 sortdata2 = (struct sort_dataD *) mymalloc("sortdata2", sizeof(struct sort_dataD) * All.MOL_CLOUDS_NumMCs); 


 /* let task 0 generate a unique list of new stellar (distinct) IDs and lifetimes for all MCs */
 if (ThisTask == 0)
  {
   for (j=0; j<All.MOL_CLOUDS_NumStars; j++)
     //if ((MyIDType) floor((All.MOL_CLOUDS_NumStars - j) * get_random_number(j))+1 < select)
     if ((MyIDType) floor((All.MOL_CLOUDS_NumStars - j) * drand48())+1 < select)
      {
       ID1_List_local[app]=j+1; //IDs start at 1!
       select--;
       app++;
      }

   for (j=0; j<All.MOL_CLOUDS_NumMCs; j++)
    lifetimes_local[j]=generate_lifetime(MCs_to_convert_global+j);

   printf(" SAMPLE OF NEW IDS: %d %d %d\n", ID1_List_local[0], ID1_List_local[All.MOL_CLOUDS_NumMCs/2], ID1_List_local[All.MOL_CLOUDS_NumMCs-1]);
  } 

 /* spread the list of IDs and lifetimes to other tasks */
 MPI_Allreduce(ID1_List_local, ID1_List_global, All.MOL_CLOUDS_NumMCs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);   
 MPI_Allreduce(lifetimes_local, lifetimes_global, All.MOL_CLOUDS_NumMCs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 


 /* filter out those MCs that need to be removed */
 for (j=0; j<All.MOL_CLOUDS_NumMCs; j++)
   if (ID2_List_global[j]==0)
     ID1_List_global[j]=0;
 
 if (ThisTask==0)
  {
   printf(" generated new IDs\n");
   fflush(stdout);
  } 

 /* now switch the IDs between MC and star particles (we distinguish MCs from stars by ID) */
 for (i=0; i<NumPart; i++)
   {
    sortdata[i].ID=P[i].ID;
    sortdata[i].ind=i;    
   }


 for (j=0; j<All.MOL_CLOUDS_NumMCs; j++)
   {
    sortdata1[j].ID1=ID1_List_global[j];
    sortdata1[j].ID2=ID2_List_global[j];    
    sortdata1[j].ind=j;    
   }


 for (j=0; j<All.MOL_CLOUDS_NumMCs; j++)
   {
    sortdata2[j].ID1=ID1_List_global[j];
    sortdata2[j].ID2=ID2_List_global[j];    
    sortdata2[j].ind=j;    
   }


 qsort(sortdata, NumPart, sizeof(struct sort_data), compare_sort_data);
 qsort(sortdata1, All.MOL_CLOUDS_NumMCs, sizeof(struct sort_dataD), compare_sort_data1); 
 qsort(sortdata2, All.MOL_CLOUDS_NumMCs, sizeof(struct sort_dataD), compare_sort_data2); 
 
 
 /* check that sorting was done correctly */ 
 if (ThisTask==0)
  for (j=0; j<All.MOL_CLOUDS_NumMCs-1; j++)
  {  
   if (((sortdata1[j].ID1==sortdata1[j+1].ID1) && (sortdata1[j].ID1!=0)) || ((sortdata2[j].ID2==sortdata2[j+1].ID2) && (sortdata2[j].ID2!=0)) || (sortdata1[j].ID1>sortdata1[j+1].ID1) || (sortdata2[j].ID2>sortdata2[j+1].ID2))
    printf(" WARNING\n");
  }
 
 /* EXCHANGE IDs and set unique MOL_CLOUDS_index */
 flag=1; j=0; i=0;
 while(flag)
 {
  if (sortdata2[j].ID2 == sortdata[i].ID)
   {
    P[sortdata[i].ind].ID=sortdata2[j].ID1;       
    P[sortdata[i].ind].MOL_CLOUDS_index=-1;
    j++; i++;
    converted_MCs_to_stars_local++;   
   }
  else
   {
     if (sortdata2[j].ID2 < sortdata[i].ID)
      j++;
    else
      i++;
   } 
  if ((j>=All.MOL_CLOUDS_NumMCs)  || (i>=NumPart))
   flag=0;
 }

 flag=1; j=0; i=0; 
 while(flag)
 { 
  if (sortdata1[j].ID1 == sortdata[i].ID)
   {
    P[sortdata[i].ind].ID=sortdata1[j].ID2;       
    P[sortdata[i].ind].MOL_CLOUDS_index=sortdata1[j].ind;
    j++; i++;
    converted_stars_to_MCs_local++;
   }
  else
   {
    if (sortdata1[j].ID1 < sortdata[i].ID)
      j++;
    else
      i++;
   } 
  if ((j>=All.MOL_CLOUDS_NumMCs)  || (i>=NumPart))
   flag=0;
 }

 if (ThisTask==0)
  {
   printf(" switched IDs\n");
   fflush(stdout);
  } 


 MPI_Allreduce(&converted_stars_to_MCs_local, &converted_stars_to_MCs_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
 MPI_Allreduce(&converted_MCs_to_stars_local, &converted_MCs_to_stars_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  

 if (ThisTask==0)
  {
   printf(" MOL_CLOUDS: time=%g stars_to_MCs(ID)=%d MCs_to_stars(ID)=%d\n", All.Time, converted_stars_to_MCs_global, converted_MCs_to_stars_global);
   printf(" MOL_CLOUDS: MCs died=%d total number of MCs=%d  \n", MCs_to_convert_global, MCcount_global);  
  } 


 /* based on the ID we can now change the particle properties accordingly */
 converted_stars_to_MCs_local=0; converted_stars_to_MCs_global=0; converted_MCs_to_stars_local=0; converted_MCs_to_stars_global=0;  


 /* EXCHANGE particle properties */
 for (i=0; i<NumPart; i++)
   {
    /* MC to star (prop) */
    if ((P[i].ID<=All.MOL_CLOUDS_NumStars) && (P[i].Type==3))
     {
      P[i].Type=2;
      P[i].Mass=All.MOL_CLOUDS_MassStars; /* stars have the mass of stellar particles */
      P[i].MOL_CLOUDS_LifeTime = 0.0;
      P[i].MOL_CLOUDS_TimeBorn = 0.0;            
      converted_MCs_to_stars_local++;
     }
    /* star to MC (prop) */
    if ((P[i].ID>All.MOL_CLOUDS_NumStars) && (P[i].Type==2))
     {
      P[i].Type=3;
      P[i].Mass=All.MOL_CLOUDS_MassStars; /* at the beginning MCs have the mass of stars and then grow in time */
      P[i].MOL_CLOUDS_LifeTime = lifetimes_global[P[i].MOL_CLOUDS_index];
      P[i].MOL_CLOUDS_TimeBorn = dyingtimes_global[P[i].MOL_CLOUDS_index];      
      converted_stars_to_MCs_local++;
     }
   }

 MPI_Allreduce(&converted_stars_to_MCs_local, &converted_stars_to_MCs_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
 MPI_Allreduce(&converted_MCs_to_stars_local, &converted_MCs_to_stars_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  

 if (ThisTask==0)
  {
   printf(" MOL_CLOUDS: stars_to_MCs(prop)=%d MCs_to_stars(prop)=%d\n", converted_stars_to_MCs_global, converted_MCs_to_stars_global); 
  } 
 
 myfree(sortdata2); 
 myfree(sortdata1);
 myfree(sortdatasmall); 
 myfree(sortdata); 
  
 /* to be sure: check IDs (can be removed to make code faster) */
 test_id_uniqueness();
 
 return 1;
}

#endif
