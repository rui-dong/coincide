#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include "aims.h"

#define PI 3.1415926
#define N_VECTOR_MAX 400000        // max number of matched vector pairs
#define N_SOL_FACTOR 1           // number of solution per rotation angle
#define BIG_NUM 65535
#define EPS 1e-4
#define DEBUG_MODE 0

/*
 * N_SOL_FACTOR=1
 * output only ONE solution per rotation angle. (finding the smallest area)
 * SORT the solution_final, output the smallest 5 configuration
 */

int brav_lattice ( double , double , double , double , double ) ;

class coincide {
    //private
    int ** solution_tmp;
    double ** norm_tmp;
    int n_sol;
    int brav_1, brav_2; // 1 SQR 2 REC 3 HEX 4 CRT 5 OBL
    double cos_angle,angle;

    void eliminate_linear ( void );
    int find_min_area ( int );
    void sort_solution ( void );
    
  public:
    int fix_const; // FIX constraint 
	double angle_const;
    double angle_start, angle_end, angle_step;
    int num_angle,num_solution;
    double norm_a1, norm_a2, norm_b1, norm_b2, angle1, angle2;
    int N;
    int N_atom_limit;
    double tolerance, angle_tolerance, max_tilt;
    double a[2][2],b[2][2];
    int nat1,nat2;
    int ** solution;
    int ** solution_final;
    double ** strain;

    coincide ( aims, aims, double, double, double );
    void loop ( void );
    void solution_output ( char * );
    void strain_calculator ( int );

    FILE * debugFile, * outputFile;
 
};

coincide::coincide ( aims layer1, aims layer2, double angle_start_in, double angle_end_in, double angle_step_in) {
  int i,j; 

  double angle_limit[5][5] = { 
  	{ 45, 45, 15, 45, 90},
  	{ 45, 90, 30, 90, 180},
  	{ 15, 30, 30, 30, 60},
  	{ 45, 90, 30, 90, 180},
  	{ 90, 180, 60, 180, 180},
	};

  N=10;
  N_atom_limit=BIG_NUM;
  angle_start = angle_start_in ;
  angle_end = angle_end_in ;
  angle_step = angle_step_in ;

  tolerance=0.02;
  angle_tolerance=0.05;
  max_tilt=0;

  brav_1=0;
  brav_2=0;
  fix_const=0;
  angle_const=0;
  //layer1.output("tp1.in");
  //layer2.output("tp2.in");

  if ( DEBUG_MODE ) {
    printf("===============\n");
    for ( i=0;i<4;i++) {
      for ( j=0;j<4;j++) {
        printf("%f ",layer1.alatt[i][j]);
      }
      printf("\n");
    }
    for ( i=0;i<3;i++) {
      for ( j=0;j<3;j++) {
        printf("%f ",layer2.alatt[i][j]);
      }
      printf("\n");
    }
  }

  a[0][0]=layer1.alatt[0][0]; 
  a[0][1]=layer1.alatt[0][1];
  a[1][0]=layer1.alatt[1][0];
  a[1][1]=layer1.alatt[1][1];
  b[0][0]=layer2.alatt[0][0];
  b[0][1]=layer2.alatt[0][1];
  b[1][0]=layer2.alatt[1][0];
  b[1][1]=layer2.alatt[1][1];

  if ( DEBUG_MODE ) {
    printf("%f\n",a[0][0]);
    printf("%f\n",a[0][1]);
    printf("%f\n",a[1][0]);
    printf("%f\n",a[1][1]);
    printf("%f\n",b[0][0]);
    printf("%f\n",b[0][1]);
    printf("%f\n",b[1][0]);
    printf("%f\n",b[1][1]);
  }

  nat1=layer1.nat;
  nat2=layer2.nat;
  
  norm_a1=sqrt( a[0][0] * a[0][0] + a[0][1] * a[0][1] );
  norm_a2=sqrt( a[1][0] * a[1][0] + a[1][1] * a[1][1] );
  cos_angle=( a[0][0] * a[1][0] + a[0][1] * a[1][1] ) / ( norm_a1 * norm_a2 );
  if ( cos_angle > 1 ) 
    cos_angle = 1;
  else if ( cos_angle < -1 )
    cos_angle = -1;
  angle1=acos(cos_angle)*180/PI;

  norm_b1=sqrt( b[0][0] * b[0][0] + b[0][1] * b[0][1] );
  norm_b2=sqrt( b[1][0] * b[1][0] + b[1][1] * b[1][1] );
  cos_angle=( b[0][0] * b[1][0] + b[0][1] * b[1][1] ) / (  norm_b1 * norm_b2 );
  if ( cos_angle > 1 )
    cos_angle = 1;
  else if ( cos_angle < -1 )
    cos_angle = -1;
  angle2=acos(cos_angle)*180/PI;

  brav_1 = brav_lattice ( norm_a1, norm_a2, angle1, tolerance, angle_tolerance );
  brav_2 = brav_lattice ( norm_b1, norm_b2, angle2, tolerance, angle_tolerance );
  
  if ( angle_end == angle_start ) {
	angle_end = angle_start + angle_limit[brav_1-1][brav_2-1];
    printf("Automatic angle range for lattice %d & %d: %f to %f\n",brav_1,brav_2,angle_start,angle_end);
  }


  if ( DEBUG_MODE ) {
    printf("cell 1: %f-%f-%f;cell 2: %f-%f-%f\n",angle1,norm_a1,norm_a2,angle2,norm_b1,norm_b2);
    printf("%d %d\n",brav_1,brav_2);
  }

  if ( brav_1 == 1 && brav_2 == 1 ) { // 
    fix_const=1;
	angle_const=90;
    printf("!!Warning: Both of the lattices are SQR.\n");
    printf("!!Warning: Only output square matched lattices!\n");
  }
  
  if ( brav_1 == 3 && brav_2 == 3 ) { // 
    fix_const=2;
	angle_const=120;
    printf("!!Warning: Both of the lattices are HEX.\n");
    printf("!!Warning: Only output hexagonal matched lattices!\n");
  }


  // initialization of the system
  // n1 n2 m1 m2 norm1 norm2 
  solution_tmp = (int **) malloc ( N_VECTOR_MAX * sizeof (int*) );
  norm_tmp = (double **) malloc ( N_VECTOR_MAX * sizeof (double*) );
  for (i=0;i<N_VECTOR_MAX;i++) {
    solution_tmp[i]= (int *) malloc (4 * sizeof (int));
    norm_tmp[i]= (double *) malloc (2 * sizeof (double));
  }

  solution = (int **) malloc ( N_VECTOR_MAX * sizeof (int*) );
  for (i=0;i<N_VECTOR_MAX;i++) {
    solution[i]= (int *) malloc (8 * sizeof (int));
  }

  num_angle=int(angle_end-angle_start)/angle_step+1;
  num_solution = 0; //N_SOL_FACTOR * num_angle ;
  solution_final= (int **) malloc ( N_SOL_FACTOR * num_angle * sizeof (int *));
  strain= (double **) malloc ( N_SOL_FACTOR * num_angle * sizeof (double *));
  for(i=0;i< num_angle * N_SOL_FACTOR; i++) {
    solution_final[i]= (int *) malloc (11 * sizeof(int));
    strain[i]= (double *) malloc (3 * sizeof(double));
    for (j=0;j<11;j++) {
      solution_final[i][j]=0;
    }
  }
  
  //printf("number number numer: %f %f %f %d\n",angle_start, angle_end, angle_step, num_angle);

}


void coincide::eliminate_linear ( void ) {
  double f1,f2;
  double min_norm;
  int flag,min_t,min_i;
  int i,j,k;
  int t;
  double norm_t;

  if ( DEBUG_MODE ) { // matched pairs before sorting
    fprintf (outputFile, "matched pairs before sorting\n");
    for (i=0;i<n_sol;i++) {
      for (k=0;k<4;k++) {
        fprintf (outputFile, "%d ", solution_tmp[i][k]);
      }
      for (k=0;k<2;k++) {
        fprintf (outputFile, "%f ", norm_tmp[i][k]);
      }
      fprintf (outputFile, "\n");
    }
  }



  // sorting ascending normA
  // sorting the array (insert)
  
  for (i=0;i<n_sol-1;i++) {
    min_norm=norm_tmp[i][0];
    //min_i=solution_tmp[i][3];
    min_t=i;
    for (j=i+1;j<n_sol;j++) {
      if ( norm_tmp[j][0] < min_norm ) { 
        min_norm = norm_tmp[j][0];
      //if ( solution_tmp[j][3] < min_i ) { 
        //min_i = solution_tmp[j][3];
        min_t = j;
      }
    }
    fprintf (outputFile, "%f %d\n",min_norm,min_t);
    
    for (k=0;k<4;k++) {
      t=solution_tmp[i][k];
      solution_tmp[i][k]=solution_tmp[min_t][k];
      solution_tmp[min_t][k]=t;
    }
    for (k=0;k<2;k++) {
      norm_t=norm_tmp[i][k];
      norm_tmp[i][k]=norm_tmp[min_t][k];
      norm_tmp[min_t][k]=norm_t;
    }
  } 

  
  if ( DEBUG_MODE ) { // matched pairs before processing
    fprintf (outputFile, "matched pairs before processing\n");
    for (i=0;i<n_sol;i++) {
      for (k=0;k<4;k++) {
        fprintf (outputFile, "%d ", solution_tmp[i][k]);
      }
      for (k=0;k<2;k++) {
        fprintf (outputFile, "%f ", norm_tmp[i][k]);
      }
      fprintf (outputFile, "\n");
    }
  }






  /*
    for the pair of (n1,m1) and (n2,m2)
    if "at least one of them is 0"
      if "both of them are 0" (because the other pair cannot be both 0, they are liean dependent, set flat to 1, and eliminate j_th solution later)
      else: "only one of them is 0" (the pair cannot be depedent)
    else: "neither of them is 0" (calculate the ratio, use for later comparison)

    if "flag == 1" || ratio1=ratio2 depedent: (elimilate j_th solution)

  */

  //
  
  
  for (i=0;i<n_sol;i++) {
    for (j=i+1;j<n_sol;j++) {
      if ( solution_tmp[j][0] == 0 && solution_tmp[j][1] == 0 ) {
        continue ; // j_th vector is already eliminated
      }
      flag = 0;
      // pair 1
      if ( solution_tmp[i][0] * solution_tmp[j][0] == 0 ) { // at least one is 0
        if ( solution_tmp[i][0] + solution_tmp[j][0] == 0 ) { // both are 0
	  flag=1;
	} else { // only one is 0
	  continue; //not linear dependent
	}
      } else { // neither is 0
	f1 = double(solution_tmp[i][0]) / double(solution_tmp[j][0]) ;
      }
      // pair 2
      if ( solution_tmp[i][1] * solution_tmp[j][1] == 0 ) { // at least one is 0
        if ( solution_tmp[i][1] + solution_tmp[j][1] == 0 ) { // both are 0
	  flag=1;
	} else { // only one is 0
	  continue; //not linear dependent
	}
      } else { // neither is 0
	f2 = double(solution_tmp[i][1]) / double(solution_tmp[j][1]) ;
      }
      // check the flag and ratio
      if ( flag == 1 || fabs(f1-f2) < EPS ) {
        solution_tmp[j][0]=0;
        solution_tmp[j][1]=0;
        solution_tmp[j][2]=0;
        solution_tmp[j][3]=0;
      }
    }
  }

  
  if ( DEBUG_MODE ) { // matched pairs after processing
    fprintf (outputFile, "matched pairs after processing\n");
    for (i=0;i<n_sol;i++) {
      for (k=0;k<4;k++) {
        fprintf (outputFile, "%d ", solution_tmp[i][k]);
      }
      for (k=0;k<2;k++) {
        fprintf (outputFile, "%f ", norm_tmp[i][k]);
      }
      fprintf (outputFile, "\n");
    }
  }

  // re-organise the solution_tmp
  
  flag=0;
  for (i=0;i<n_sol;i++) {
    if ( solution_tmp[i][0] == 0 && solution_tmp[i][1] == 0 ) { 
      continue; // empty slot: i++
    } else { // non-empty slot:
      if ( flag == i ) { // every slot is occupied so far, don't copy, flag++ and i++
        flag++;
	continue;
      } else { // need to copy, and flag++
        solution_tmp[flag][0]=solution_tmp[i][0];
        solution_tmp[flag][1]=solution_tmp[i][1];
        solution_tmp[flag][2]=solution_tmp[i][2];
        solution_tmp[flag][3]=solution_tmp[i][3];
        norm_tmp[flag][0]=norm_tmp[i][0];
        norm_tmp[flag][1]=norm_tmp[i][1];
	flag++;
      }
    }
  }
  n_sol=flag;

  

  
  if ( DEBUG_MODE ) { // matched pairs after re-organization
    fprintf (outputFile, "matched pairs after re-organization\n");
    for (i=0;i<n_sol;i++) {
      for (k=0;k<4;k++) {
        fprintf (outputFile, "%d ", solution_tmp[i][k]);
      }
      for (k=0;k<2;k++) {
        fprintf (outputFile, "%f ", norm_tmp[i][k]);
      }
      fprintf (outputFile, "\n");
    }
  }

}

int coincide::find_min_area ( int i_angle ) {
  double Am[2],Am_p[2];
  int count,k;
  int vec_i=-1,vec_j=-1;
  int i,j,t,m1,m2,m1_p,m2_p,min_area,area;
  int n1, n2;
  double f1,f2;
  double min_angle=90;
  double normA, normB;
  int nat;
  int flag;
  min_area=BIG_NUM;
  int * sorting_tmp;

  for (i=0;i<n_sol;i++) {
    for (j=i+1;j<n_sol;j++) {
      m1=solution_tmp[i][0];
      m2=solution_tmp[i][1];
      m1_p=solution_tmp[j][0];
      m2_p=solution_tmp[j][1];

      Am[0] = m1*a[0][0] + m2*a[1][0];
      Am[1] = m1*a[0][1] + m2*a[1][1];
      Am_p[0] = m1_p * a[0][0] + m2_p * a[1][0];
      Am_p[1] = m1_p * a[0][1] + m2_p * a[1][1];

      area = abs ( m1 * m2_p - m1_p * m2 ) ;
      normA = sqrt(Am[0]*Am[0]+Am[1]*Am[1]);
      normB = sqrt(Am_p[0]*Am_p[0]+Am_p[1]*Am_p[1]);
      cos_angle=( Am[0] * Am_p[0] + Am[1] * Am_p[1]) / ( normA * normB );

      if ( cos_angle > 1 )
        cos_angle = 1;
      else if ( cos_angle < -1 )
        cos_angle = -1;
      
      angle=acos(cos_angle)*180/PI;

      
      flag=0;
      if ( area >=1 && area < min_area ) {
        flag=1;
      }
      if ( fix_const > 0 ) {
	    if ( fabs( fabs(angle)-angle_const ) > angle_tolerance || fabs(normA-normB) > tolerance ) {
          flag=0;
		}
      }
	  if ( fabs(angle) < max_tilt || fabs(angle) > 180-max_tilt ) {
	    flag=0;
	  }
      if ( flag == 1 ) {
        min_area = area;
		vec_i=i; vec_j=j;
      }
    }
  }
  
  //check for the total of atoms
  if ( vec_i!=-1 ) {

    //k= i_angle * N_SOL_FACTOR +i ;
    // get the solution
    //while (sorting_tmp[j] < 0) {
    //  j++;
    //}

    //min_area = solution[sorting_tmp[i]][0];
    //vec_i = solution[sorting_tmp[i]][1];
    //vec_j = solution[sorting_tmp[i]][2]; 
    //j++;
         
    nat = min_area * nat1;
    area = abs ( solution_tmp[vec_i][2]*solution_tmp[vec_j][3]-solution_tmp[vec_j][2]*solution_tmp[vec_i][3] );
    nat += area * nat2;

    if ( nat >= N_atom_limit ) { 
      return 0;
    } 

    //printf ("nat = %d\n",nat);
    //printf ("i_angle: %d  config: %d  total solution: %d\n",i_angle,i,num_solution);

    solution_final[ num_solution ][0] = solution_tmp[vec_i][0] ;
    solution_final[ num_solution ][1] = solution_tmp[vec_i][1] ;
    solution_final[ num_solution ][2] = solution_tmp[vec_j][0] ;
    solution_final[ num_solution ][3] = solution_tmp[vec_j][1] ;
    solution_final[ num_solution ][4] = solution_tmp[vec_i][2] ;
    solution_final[ num_solution ][5] = solution_tmp[vec_i][3] ;
    solution_final[ num_solution ][6] = solution_tmp[vec_j][2] ;
    solution_final[ num_solution ][7] = solution_tmp[vec_j][3] ;
    solution_final[ num_solution ][8] = i_angle;
    solution_final[ num_solution ][9] = nat;
    solution_final[ num_solution ][10] = area;
    // calculate some properties, strain
    strain_calculator(num_solution);

    num_solution++;
    printf ("find solution at i_angle=%d, total solution is now %d\n",i_angle,num_solution);

  
    //printf("surface 1: A1=%d*a1+%d*a2, A2=%d*a1+%d*a2.\n",solution[i_angle][0],solution[i_angle][1],solution[i_angle][2],solution[i_angle][3]);
    //printf("surface 2: B1=%d*b1+%d*b2, B2=%d*b1+%d*b2.\n",solution[i_angle][4],solution[i_angle][5],solution[i_angle][6],solution[i_angle][7]);
  }
  //free(sorting_tmp);
  return 1; 
}


void coincide::solution_output ( char * file_out ) {
  int i,k,j;
  double angle;
  FILE * fid_log;

  /* // normal version
  for (i=0;i<num_angle;i++) {
    printf("angle=%f: surface 1: ",angle_start+angle_step*i);
    for ( k=0;k<4;k++ ) {
      printf("%d ",solution[i][k]);
    }
    printf("; surface 2: ");
    for ( ;k<8;k++ ) {
      printf("%d ",solution[i][k]);
    }
    printf("\n");
  }
  */

  // reuquested by Tom Archer
  /*
  double a1[2],b1[2];
  double norm_a, norm_b, c,err;
  for (i=0;i<num_angle;i++) {
    
    a1[0]=solution_final[i][0]*a[0][0]+solution_final[i][1]*a[1][0];
    a1[1]=solution_final[i][0]*a[0][1]+solution_final[i][1]*a[1][1];

    b1[0]=solution_final[i][4]*b[0][0]+solution_final[i][5]*b[1][0];
    b1[1]=solution_final[i][4]*b[0][1]+solution_final[i][5]*b[1][1];

    norm_a=vec_norm(a1,2);
    norm_b=vec_norm(b1,2);

    if ( norm_a < EPS ) {
      continue;
    }

    c=(norm_a+norm_b)/2;
    err=(norm_a-norm_b)/c;
    printf("%f  %f  %f\n",angle_start+angle_step*i,c*c*sqrt(3),err);
  } */
  
  fid_log=fopen(file_out,"w");
  fprintf(fid_log,"# number angle nat area strain_x strain_y strain_xy vector index x8 \n");
  for (i=0;i<num_solution;i++) {
    angle=angle_start+angle_step*double(solution_final[i][8]);
    fprintf(fid_log,"%d  %f  %d  %d  %f  %f  %f ",i,angle,solution_final[i][9],solution_final[i][10],strain[i][0],strain[i][1],strain[i][2]);
    for (j=0;j<8;j++) {
      fprintf(fid_log,"%d ",solution_final[i][j]);
    }
    fprintf(fid_log,"\n");
  }
  fclose(fid_log);

}


void coincide::loop ( void ) {
  //FILE *outputFile, *inputFile;
  if ( DEBUG_MODE ) {
    debugFile = fopen ("debug.tmp", "w");
  }
  
  double xA_1, yA_1, xA_2, yA_2;
  double xB_1, yB_1, xB_2, yB_2;

  xA_1=a[0][0]; yA_1=a[0][1];
  xA_2=a[1][0]; yA_2=a[1][1];
  xB_1=b[0][0]; yB_1=b[0][1];
  xB_2=b[1][0]; yB_2=b[1][1];
  
  int Nmax; Nmax=N;


  //Loops through the angles and vectors in order to seek coincidences
  //and prints the output file with the solutions for each angle
   
  if ( DEBUG_MODE ) { outputFile = fopen ("coincidences.tmp", "w"); }
    
  int m1, m2, n1, n2;
  double angle, angleRad, angle_Am_MBn, cosAngle;
  double xAm, yAm, xMBn, yMBn;
  double norm, normA, normB;
  int i_angle; 


  /* //system check
  if ( (if_hex_1 != 1 || if_hex_2 != 1 ) && fix_hex == 1 ) {
    printf("warning: at least one of the cells are not hexagonal, output should NOT be fixed to hexagonal!!\n");
  } */


  i_angle=0;
  //printf("before executing: %f %f %f\n",angle_start,angle_end,angle_step);
  for (angle = angle_start; angle < angle_end + angle_step; angle = angle + angle_step) {
    angleRad = angle*PI/180;
    
    n_sol=0;
    if ( DEBUG_MODE ) { fprintf (outputFile, "\n%.2f\n", angle); }
    
    for (m1 = -Nmax; m1 <= Nmax; m1++) {
      for (m2 = -Nmax; m2 <= Nmax; m2++) {
        for (n1 = -Nmax; n1 <= Nmax; n1++) {
          for (n2 = -Nmax; n2 < Nmax; n2++) {
            // |Am - MBn| < tolerance
            xAm = m1*xA_1 + m2*xA_2;
            xMBn = n1*(xB_1*cos(angleRad) + yB_1*sin(angleRad)) + n2*(xB_2*cos(angleRad) + yB_2*sin(angleRad));
            
            yAm = m1*yA_1 + m2*yA_2;
            yMBn = n1*(-xB_1*sin(angleRad) + yB_1*cos(angleRad)) + n2*(-xB_2*sin(angleRad) + yB_2*cos(angleRad));
            
            normA = sqrtf (pow(xAm, 2) + pow(yAm, 2));
            normB = sqrtf (pow(xMBn, 2) + pow(yMBn, 2));
            
            if (normA >= normB)
              norm = normB;
            else
              norm = normA;
            
            // Angle between the two vectors
            cosAngle = (xAm*xMBn + yAm*yMBn)/(normA*normB);
            
            // To prevent rounding errors, which lead to cos > 1
            if (cosAngle > 1.0)
              cosAngle = 1.0;
              
            angle_Am_MBn = 180*acos(cosAngle)/PI;
           
            if (sqrtf (pow(xAm - xMBn, 2) + pow(yAm - yMBn, 2))/norm < tolerance && fabs(angle_Am_MBn) < angle_tolerance) {
              if ( DEBUG_MODE ) { fprintf (outputFile, "%d %d %d %d\n", m1, m2, n1, n2); }
              solution_tmp[n_sol][0]= m1;
              solution_tmp[n_sol][1]= m2;
              solution_tmp[n_sol][2]= n1;
              solution_tmp[n_sol][3]= n2;
              norm_tmp[n_sol][0]= normA;
              norm_tmp[n_sol][1]= normB;
	      ++n_sol;
	    }
	    if ( n_sol >= N_VECTOR_MAX ) {
	      printf("error: n_sol overflow!!!!\n");
	    }
          }
        }
      }
    }
    // eliminate_linear(); 
    find_min_area(i_angle);
    i_angle++;
  }
  
  // sort solution_final
  int min_nat,min_j,tmp,k,i,j;
  for ( i=0;i<num_solution-1;i++ ) {
    min_nat=solution_final[i][9];
    min_j=i;
    for ( j=i+1;j<num_solution;j++ ) {
      if ( solution_final[j][9] < min_nat ) {
        min_j=j;
	min_nat=solution_final[j][9];
      }
    }
    if (min_j!=i) {
      for (k=0;k<11;k++) {
        tmp=solution_final[i][k];
        solution_final[i][k]=solution_final[min_j][k];
        solution_final[min_j][k]=tmp;
      }
    }
  }
  
  if ( DEBUG_MODE ) {
    fclose (outputFile);  
    fclose (debugFile);
  }
 
}

void coincide::strain_calculator ( int k ) {
  double ref[3][3],target[3][3];
  double ref_p[6], target_p[6];
  double n,n_p;
  double y2norm,y2x,y2y;
  int i,j;
    
    ref[0][0]=solution_final[k][0]*a[0][0]+solution_final[k][1]*a[1][0];
    ref[0][1]=solution_final[k][0]*a[0][1]+solution_final[k][1]*a[1][1];
    ref[1][0]=solution_final[k][2]*a[0][0]+solution_final[k][3]*a[1][0];
    ref[1][1]=solution_final[k][2]*a[0][1]+solution_final[k][3]*a[1][1];
    ref[0][2]=0; ref[1][2]=0;
    ref[2][0]=0; ref[2][1]=0; ref[2][2]=1;

    target[0][0]=solution_final[k][4]*b[0][0]+solution_final[k][5]*b[1][0];
    target[0][1]=solution_final[k][4]*b[0][1]+solution_final[k][5]*b[1][1];
    target[1][0]=solution_final[k][6]*b[0][0]+solution_final[k][7]*b[1][0];
    target[1][1]=solution_final[k][6]*b[0][1]+solution_final[k][7]*b[1][1];
    target[0][2]=0; target[1][2]=0;
    target[2][0]=0; target[2][1]=0; target[2][2]=1;


    if ( DEBUG_MODE ) {
      for ( i=0;i<3;i++) {
        for ( j=0;j<3;j++) {
          printf("%f ",ref[i][j]);
        }
        printf("\n");
      }
      for ( i=0;i<3;i++) {
        for ( j=0;j<3;j++) {
          printf("%f ",target[i][j]);
        }
        printf("\n");
      }
    }


    alatt2cellp(ref,ref_p);
    alatt2cellp(target,target_p);
    
    if ( DEBUG_MODE ) {
      for ( j=0;j<6;j++) {
        printf("%f ",ref_p[j]);
      }
      printf("\n");
      for ( j=0;j<6;j++) {
        printf("%f ",target_p[j]);
      }
      printf("\n");
    }

    strain[k][0]=(ref_p[0]-target_p[0])/target_p[0]; // e_xx

    cellp2alatt(ref_p,ref);
    cellp2alatt(target_p,target);

    if ( DEBUG_MODE ) {
      for ( i=0;i<3;i++) {
        for ( j=0;j<3;j++) {
          printf("%f ",ref[i][j]);
        }
        printf("\n");
      }
      for ( i=0;i<3;i++) {
        for ( j=0;j<3;j++) {
          printf("%f ",target[i][j]);
        }
        printf("\n");
      }
    }


    n=-ref[1][0]/ref[0][0]; 
    n_p=-target[1][0]/target[0][0]; 

    y2x= n_p * ref[0][0] + ref[1][0]; 
    y2y= n_p * ref[0][1] + ref[1][1]; 
    y2norm= sqrt( y2x*y2x + y2y*y2y );

    if (DEBUG_MODE) {
      printf("n = %f\n",n);
      printf("n_p = %f\n",n_p);
      printf(" y2x = %f\n",y2x);
      printf(" y2y = %f\n",y2y);
      printf(" y2norm= %f\n",y2norm);
    }

    //strain[k][1]= ( target[1][1] - y2norm ) / y2norm;
    strain[k][1]= ( y2norm - target[1][1] ) / target[1][1];
    strain[k][2]=atan(y2x/y2y)/2;

    
    if ( DEBUG_MODE ) {
      printf("exx=%f eyy=%f exy=%f\n",strain[k][0],strain[k][1],strain[k][2]);
    }

    //n=ref_p[1] * sin( (ref_p[3]-90)/180*PI ) / ref_p[0] ; 
    //n_p=target_p[1] * sin( (target_p[3]-90)/180*PI ) / target_p[0] ; 

    //y1y = target_p[1] * sin ( target_p[3]/180*PI );
    //y2y = ref_p[1] * sin ( ref_p[3]/180*PI );
    //y2x = n_p * ref_p[0] + ref_p[1] * cos ( ref_p[3]/180*PI );

    //strain[k][1]=(sqrt( y2x*y2x + y2y*y2y )-y1y)/y1y; // e_yy
    //strain[k][2]=atan(y2x/y2y)/2;


    //c=(norm_a+norm_b)/2;
    //err=(norm_a-norm_b)/c;
    //printf("%f  %f  %f\n",angle_start+angle_step*i,c*c*sqrt(3),err);

}


int brav_lattice ( double a1, double a2, double angle, double tol, double angle_tol ) {
	int brav=5;
	if ( fabs(angle-120) < angle_tol || fabs(angle-60) < angle_tol ) {
		if ( fabs( a1 - a2 ) < tol ) {
			brav=3;
			return brav;
		}
	}
	if ( fabs(angle-90) < angle_tol ) {
		if ( fabs( a1 - a2 ) < tol ) {
			brav=1;
		} else {
			brav=2;
		}
		return brav;
	}
	if ( fabs( a1 - a2 ) < tol ) {
		brav=4;
	}
	return brav;
}





