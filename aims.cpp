#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <cstdlib>
#include <math.h>

#define NAT_MAX 500000
#define EPS 1e-4
#define PI 3.1415926
#define DEBUG_AIMS 0

using namespace std;

void alatt2cellp (double [3][3] , double [6] ) ;
void cellp2alatt (double [6] , double [3][3] ) ;
double vec_norm( double * , int ) ;
double vec_dot( double * , double * , int ) ;
void average_cell( double [3][3] , double [3][3] , double [3][3] ) ;

class aims {
    int i,j,k,l,m,n;
  public:
    int N, M, sys_type;
    char prefix [64];
    char prefix_out [64];

    int nat,nlatt;
    double xyz[NAT_MAX][3], xyz_frac[NAT_MAX][3];
    char name[NAT_MAX][10];
    int frac_flag[NAT_MAX]; // mixed "atom" and "atom_frac" are allowed
    double alatt[3][3];
    double cell_p[6];
  
    int if_frac; // 1: output fractional coordinates
    int fix_axis[3]; // 1: fix cell dimension in this direction
    
    aims (char *);
    aims ();
    virtual ~aims() {};
    int output (char *);
    void xyz2frac ( void );
    void frac2xyz ( void );
    void multiple( int, int, int, int );
    void simple_multiple( int, int, int);
    void supercell ( int [4] );
    void displace ( double [3] );
};


int aims_merge ( aims * , aims , aims , int , int , int) ; 
int aims_merge_pointer ( aims * , aims *, aims *, int , int , int) ; 


aims::aims () { // to create an empty aims object
  int i,j;

  nat=0; nlatt=0; if_frac=0;
  fix_axis[0]=0; fix_axis[1]=0; fix_axis[2]=0;
  
  //for (i=0;i<NAT_MAX;i++) {
  //  for (j=0;j<3;j++) {
  //    xyz[i][j]=0;
  //    xyz_frac[i][j]=0;
  //    name[i][j]=32;
  //  }
  //}

  sprintf(prefix_out,"geometry");
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      alatt[i][j]=0;
    }
  }

}

aims::aims (char * prefix) {
  // read in AIMS geometry file named PREFIX.in
  int i;
  double x, y, z;
  char filein[128];
  char tmp[256];
  std::string line;

  strcpy(filein,prefix);
  strcat(filein,".in");
  

  ifstream fidin ( filein );
  if (!fidin.is_open()) {
    printf("file cannot be open!"); abort();
  }

  nat=0; nlatt=0; if_frac=0;
  fix_axis[0]=0; fix_axis[1]=0; fix_axis[2]=0;

  sprintf(prefix_out,"geometry");


  while ( getline(fidin,line) ) {
    strcpy(tmp,line.c_str());
    n=line.size();
    if (n==0) { //empty line
      continue;
    } else if (tmp[0]==35) { //comment line
      continue;
    }
    //printf ("%s\n",tmp);

    // atomic positions
    n=line.find("atom_frac");
    if (n!=-1) {
      //printf("%d\n",n);
      frac_flag[nat]=1;
      sscanf(&tmp[n+10],"%lf %lf %lf %s",&xyz[nat][0],&xyz[nat][1],&xyz[nat][2],name[nat]);
      nat+=1;
    } else {
      n=line.find("atom");
      if (n!=-1) {
        //printf("%d\n",n);
        sscanf(&tmp[n+5],"%lf %lf %lf %s",&xyz[nat][0],&xyz[nat][1],&xyz[nat][2],name[nat]);
        nat+=1;
      } else {
        // lattice constants
        n=line.find("lattice_vector");
        if (n!=-1) {
          //printf("%d\n",n);
          sscanf(&tmp[n+14],"%lf %lf %lf",&alatt[nlatt][0],&alatt[nlatt][1],&alatt[nlatt][2]);
          nlatt+=1;
        }
      }
    }   
  }

  // check before proceed
  if ( nat>NAT_MAX ) {
    printf("number of atom overflow!"); abort();
  } else if ( nat==0 ) {
    printf("no atom found"); abort();
  } else if ( nlatt>3 ) {
    printf("more than three lattice vectors!"); abort();
  }




  // convert all atomic positions to xyz
  for(i=0;i<nat;i++) {
    if (frac_flag[i]!=0) {
      xyz_frac[i][0]=xyz[i][0]; 
      xyz_frac[i][1]=xyz[i][1]; 
      xyz_frac[i][2]=xyz[i][2];

      
      xyz[i][0]=xyz_frac[i][0]*alatt[0][0]+xyz_frac[i][1]*alatt[1][0];
      xyz[i][1]=xyz_frac[i][0]*alatt[0][1]+xyz_frac[i][1]*alatt[1][1];
      xyz[i][2]=xyz_frac[i][2]*alatt[2][2];

      /*//full version
      xyz[i][0]=xyz_frac[0]*alatt[0][0]+xyz_frac[1]*alatt[1][0]+xyz_frac[2]*alatt[2][0];
      xyz[i][1]=xyz_frac[0]*alatt[0][1]+xyz_frac[1]*alatt[1][1]+xyz_frac[2]*alatt[2][1];
      xyz[i][2]=xyz_frac[0]*alatt[0][2]+xyz_frac[1]*alatt[1][2]+xyz_frac[2]*alatt[2][2];
      */
    }    
  }
  // create fractional position for the rest
  for(i=0;i<nat;i++) {
    if (frac_flag[i]==0) {
      // special for 2D materials
      xyz_frac[i][0]=(alatt[1][1]*xyz[i][0]-alatt[1][0]*xyz[i][1])/(alatt[0][0]*alatt[1][1]-alatt[1][0]*alatt[0][1]); 
      xyz_frac[i][1]=(alatt[0][1]*xyz[i][0]-alatt[0][0]*xyz[i][1])/(alatt[1][0]*alatt[0][1]-alatt[0][0]*alatt[1][1]); 
      xyz_frac[i][2]=xyz[i][2]/alatt[2][2];

      // full version (......)
    }
  }

  printf("---- finish reading from FHI-AIMS geometry file %s.in\n",prefix);
  
  alatt2cellp(alatt,cell_p);
  printf("-- Number of atoms: %d\n",nat);
  printf("-- cell paramters: a=%4.2f; b=%4.2f; c=%4.2f; bc=%4.1f; ac=%4.1f; ab=%4.1f\n",cell_p[0],cell_p[1],cell_p[2],cell_p[4],cell_p[5],cell_p[3]);

  // print to log file
  /*
  printf("nat=%d\n",nat);
  for (i=0;i<nat;i++) {
    printf("%f %f %f %s %d\n",xyz[i][0],xyz[i][1],xyz[i][2],name[i],frac_flag[i]);
  }
  printf("nlatt=%d\n",nlatt);
  for (i=0;i<3;i++) {
    printf("%f %f %f\n",alatt[i][0],alatt[i][1],alatt[i][2]);
  } */

}

void aims::simple_multiple( int nx, int ny, int nz ) {
  int flag=0;
  int count,i,j,ia;

  count=nat;
  for ( k=0;k<nz;k++) {
    for ( j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        if ( i==0 && j==0 && k==0 ) {
          continue;
        }
        for ( ia=0;ia<nat;ia++ ) {
            xyz_frac[count][0]=xyz_frac[ia][0]+i;
            xyz_frac[count][1]=xyz_frac[ia][1]+j;
            xyz_frac[count][2]=xyz_frac[ia][2]+k;
            sprintf(name[count],"%s",name[ia]);
            ++count;
        }
      }
    }
  }

  nat=count;
  frac2xyz();
  for ( i=0;i<3;i++ ) {
    alatt[0][i]*=nx;
    alatt[1][i]*=ny;
    alatt[2][i]*=nz;
  }
  xyz2frac();

}


void aims::multiple( int x_min, int x_max, int y_min, int y_max ) {
  int flag=0;
  int count,i,j,ia;

  count=nat;
  for (i=x_min;i<=x_max;i++) {
    for ( j=y_min;j<=y_max;j++) {
      if ( i==0 && j==0 ) {
        flag=1;
	continue;
      }
      for ( ia=0;ia<nat;ia++ ) {
          xyz_frac[count][0]=xyz_frac[ia][0]+i;
          xyz_frac[count][1]=xyz_frac[ia][1]+j;
          xyz_frac[count][2]=xyz_frac[ia][2];
          sprintf(name[count],"%s",name[ia]);
	  ++count;
      }
    }
  }

  if ( flag==0 ) {
    printf("warning: original cell is not supposed to be included in the large cell, however it is in!\n");
  }

  nat=count;

  frac2xyz();

}

void aims::supercell( int m[4] ) {
  //int m1,m2,n1,n2
  int i,j,count,area,flag;
  int x_min,x_max,y_min,y_max,nat_prim,nat_area;

  nat_prim=nat;

  // get the "rectangular" that enclose the target super cell
  if ( m[0]*m[2] < 0) { // m1 m3 opposite
    x_min=min(m[0],m[2]);
    x_max=max(m[0],m[2]);
  } else if ( m[0]+m[2] < 0 ) { // m1 m3 both negative (including 0/-)
    x_min=m[0]+m[2];
    x_max=0;
  } else { // m1 m3 both positive (including 0/+)
    x_min=0;
    x_max=m[0]+m[2];
  }

  if ( m[1]*m[3] < 0) { // m2 m4 opposite
    y_min=min(m[1],m[3]);
    y_max=max(m[1],m[3]);
  } else if ( m[1]+m[3] < 0 ) { // m2 m4 both negative (including 0/-)
    y_min=m[1]+m[3];
    y_max=0;
  } else { // m2 m4 both positive (including 0/+)
    y_min=0;
    y_max=m[1]+m[3];
  }
  
  //printf("basis 1 from %d to %d\nbasis 2 from %d to %d\n",x_min,x_max,y_min,y_max);
  multiple(x_min,x_max,y_min,y_max);
  //output("tmp1.in");

  // cut the supercell from the "rectangular"
  // new lattice constant
  double a1,a2,b1,b2;
  a1=m[0]*alatt[0][0]+m[1]*alatt[1][0];
  a2=m[0]*alatt[0][1]+m[1]*alatt[1][1];
  b1=m[2]*alatt[0][0]+m[3]*alatt[1][0];
  b2=m[2]*alatt[0][1]+m[3]*alatt[1][1];

  area = m[0]*m[3]-m[1]*m[2];
  if ( area > 0 ) {
    alatt[0][0]=a1; alatt[0][1]=a2;
    alatt[1][0]=b1; alatt[1][1]=b2;
  } else {
    alatt[0][0]=b1; alatt[0][1]=b2;
    alatt[1][0]=a1; alatt[1][1]=a2;
    area=-area;
  }
  
  xyz2frac();  //convert to fractional using new lattice
  
  nat_area=area*nat_prim;

  double eps=EPS*10;
  double dx,dy,dz;
  count=0; //
  for (i=0;i<nat;i++) {
    if ( -eps <= xyz_frac[i][0] && xyz_frac[i][0] < 1+eps && -eps <= xyz_frac[i][1] && xyz_frac[i][1] < 1+eps )  {
      frac_flag[count]=i;
      ++count;
    }
  }
  
  for (i=0;i<count;i++) {
    if ( frac_flag[i]<0 ) { continue; } // deleted
    for (j=i+1;j<count;j++) {
      if ( frac_flag[j]<0 ) { continue; } //deleted

      dx = fabs( xyz_frac[frac_flag[i]][0] - xyz_frac[frac_flag[j]][0] );
      dy = fabs( xyz_frac[frac_flag[i]][1] - xyz_frac[frac_flag[j]][1] );
      dz = fabs( xyz_frac[frac_flag[i]][2] - xyz_frac[frac_flag[j]][2] );
      if ( dy<EPS || ( dy>1-EPS && dy<1+EPS ) ) {
        if ( dx<EPS || ( dx>1-EPS && dx<1+EPS ) ) {
	  if ( dz<EPS ) {
          frac_flag[j]=-1; 
	  }
        }
      }
    }
  } // duplicated deleted
  for ( i=0,j=0;j<count;j++ ) {
    if ( frac_flag[j]<0 ) { continue; }
    if ( i!=j ) {
      frac_flag[i]=frac_flag[j]; // move the last one here
    }
    i++;
  }
  count=i;

  nat=count;
  for (i=0;i<count;i++) {
    j=frac_flag[i];
    if (i==j) { continue; }
    xyz_frac[i][0]=xyz_frac[j][0];
    xyz_frac[i][1]=xyz_frac[j][1];
    xyz_frac[i][2]=xyz_frac[j][2];
    sprintf(name[i],"%s",name[j]);
  }

  frac2xyz();

  // possibly a function of sorting

}

void aims::displace ( double vector[3] ) {
  int i,k;
  for ( k=0;k<3;k++ ) {
    if ( fabs(vector[k]) < EPS ) { continue; }
    for ( i=0;i<nat;i++ ) {
      xyz[i][k]+=vector[k];
    }
  }
  xyz2frac();
}


void aims::frac2xyz ( void ) {
  int i;
  for(i=0;i<nat;i++) {
    xyz[i][0]=xyz_frac[i][0]*alatt[0][0]+xyz_frac[i][1]*alatt[1][0];
    xyz[i][1]=xyz_frac[i][0]*alatt[0][1]+xyz_frac[i][1]*alatt[1][1];
    xyz[i][2]=xyz_frac[i][2]*alatt[2][2];

    /*//full version
    xyz[i][0]=xyz_frac[0]*alatt[0][0]+xyz_frac[1]*alatt[1][0]+xyz_frac[2]*alatt[2][0];
    xyz[i][1]=xyz_frac[0]*alatt[0][1]+xyz_frac[1]*alatt[1][1]+xyz_frac[2]*alatt[2][1];
    xyz[i][2]=xyz_frac[0]*alatt[0][2]+xyz_frac[1]*alatt[1][2]+xyz_frac[2]*alatt[2][2];
    */
  }
}

void aims::xyz2frac (void) {
  int i;
  for(i=0;i<nat;i++) {
    // special for 2D materials
    xyz_frac[i][0]=(alatt[1][1]*xyz[i][0]-alatt[1][0]*xyz[i][1])/(alatt[0][0]*alatt[1][1]-alatt[1][0]*alatt[0][1]); 
    xyz_frac[i][1]=(alatt[0][1]*xyz[i][0]-alatt[0][0]*xyz[i][1])/(alatt[1][0]*alatt[0][1]-alatt[0][0]*alatt[1][1]); 
    xyz_frac[i][2]=xyz[i][2]/alatt[2][2];

    // full version (......)
  }

}

int aims::output( char * fileout ){
  int i;
  //char * fileout;
  FILE * fidout;

  //fileout=sprintf("%s.in",prefix_out);
  fidout=fopen(fileout,"w");
  //fidout=fopen("geometry.in","w");

  fprintf(fidout,"#lattice constants\n");
  for(i=0;i<3;i++) {
    fprintf(fidout,"lattice_vector %14.9f %14.9f %14.9f\n",alatt[i][0],alatt[i][1],alatt[i][2]);
    if ( fix_axis[i]==1 ) {
      fprintf(fidout,"constrain_relaxation .true.\n");
    }
  }

  fprintf(fidout,"#atomic positions\n");
  if ( if_frac ) {
    for(i=0;i<nat;i++) {
      fprintf(fidout,"atom_frac %14.9f %14.9f %14.9f %s\n",xyz_frac[i][0],xyz_frac[i][1],xyz_frac[i][2],name[i]);
    }
  } else {
    for(i=0;i<nat;i++) {
      fprintf(fidout,"atom %14.9f %14.9f %14.9f %s\n",xyz[i][0],xyz[i][1],xyz[i][2],name[i]);
    }
  }

  alatt2cellp(alatt,cell_p);
  fprintf(fidout,"#cell_parameter %d %14.9f %14.9f %14.9f %8.4f %8.4f %8.4f\n",nat,cell_p[0],cell_p[1],cell_p[2],cell_p[3],cell_p[4],cell_p[5]);

  fclose(fidout);
  return 1;

}


int aims_merge ( aims * target , aims source_1, aims source_2, int coord_type_1, int coord_type_2 , int cell_type ) { 
  /* 
	merge two aims structures "source_1" and "source_2" to "target"
	coord_type: 0=use xyz directly; 1=use xyz_frac
	cell_type: 0= use "average" cell (deform)
	          1/2 = use the cell from source_1/2
	
  */

int i,j,count;

  // in deformed cell, use fractional coord by force
  if ( cell_type==0 ) { 
    coord_type_1=1; coord_type_2=1;
  }

  if ( cell_type == 1) {
    alatt2cellp(source_1.alatt,target->cell_p);
    cellp2alatt(target->cell_p,target->alatt);
  } else if (cell_type ==2 ){
    alatt2cellp(source_2.alatt,target->cell_p);
    cellp2alatt(target->cell_p,target->alatt);
  } else {
    average_cell(target->alatt,source_1.alatt,source_2.alatt);
  }

  count=0;
  if (coord_type_1==0) {
    for(i=0;i<source_1.nat;i++) {
      target->xyz[count][0]=source_1.xyz[i][0];
      target->xyz[count][1]=source_1.xyz[i][1];
      target->xyz[count][2]=source_1.xyz[i][2];
      sprintf(target->name[count],"%s",source_1.name[i]);
      ++count;
    }
  } else {
    for(i=0;i<source_1.nat;i++) {
      target->xyz_frac[count][0]=source_1.xyz_frac[i][0];
      target->xyz_frac[count][1]=source_1.xyz_frac[i][1];
      target->xyz_frac[count][2]=source_1.xyz_frac[i][2];
      sprintf(target->name[count],"%s",source_1.name[i]);
      ++count;
    }
  }

  if (coord_type_2==0) {
    for(i=0;i<source_2.nat;i++) {
      target->xyz[count][0]=source_2.xyz[i][0];
      target->xyz[count][1]=source_2.xyz[i][1];
      target->xyz[count][2]=source_2.xyz[i][2];
      sprintf(target->name[count],"%s",source_2.name[i]);
      ++count;
    }
  } else {
     for(i=0;i<source_2.nat;i++) {
      target->xyz_frac[count][0]=source_2.xyz_frac[i][0];
      target->xyz_frac[count][1]=source_2.xyz_frac[i][1];
      target->xyz_frac[count][2]=source_2.xyz_frac[i][2];
      sprintf(target->name[count],"%s",source_2.name[i]);
      ++count;
    }
  }
  
  target->nat=source_1.nat+source_2.nat;
  if (target->nat != count)
    printf("warning: number of atom mismatch in at least one of the structures");

  // convert to coord using new lattice constants
  if ( coord_type_1 == 0 && coord_type_2 == 0 ) {
    target->xyz2frac();
  } else {
    target->frac2xyz();
  }
}



int aims_merge_pointer ( aims * target , aims * source_1, aims * source_2, int coord_type_1, int coord_type_2 , int cell_type ) { 
  /* 
	merge two aims structures "source_1" and "source_2" to "target"
	coord_type: 0=use xyz directly; 1=use xyz_frac
	cell_type: 0= use "average" cell (deform)
	          1/2 = use the cell from source_1/2
	
  */

int i,j,count;

  // in deformed cell, use fractional coord by force
  if ( cell_type==0 ) { 
    coord_type_1=1; coord_type_2=1;
  }

  if ( cell_type == 1) {
    alatt2cellp(source_1->alatt,target->cell_p);
    cellp2alatt(target->cell_p,target->alatt);
  } else if (cell_type ==2 ){
    for (i=0;i<3;i++) {
      for (j=0;j<3;j++) {
        target->alatt[i][j]=source_2->alatt[i][j];
      }
    }
  } else {
    average_cell(target->alatt,source_1->alatt,source_2->alatt);
  }

  count=0;
  if (coord_type_1==0) {
    for(i=0;i<source_1->nat;i++) {
      target->xyz[count][0]=source_1->xyz[i][0];
      target->xyz[count][1]=source_1->xyz[i][1];
      target->xyz[count][2]=source_1->xyz[i][2];
      sprintf(target->name[count],"%s",source_1->name[i]);
      ++count;
    }
  } else {
    for(i=0;i<source_1->nat;i++) {
      target->xyz_frac[count][0]=source_1->xyz_frac[i][0];
      target->xyz_frac[count][1]=source_1->xyz_frac[i][1];
      target->xyz_frac[count][2]=source_1->xyz_frac[i][2];
      sprintf(target->name[count],"%s",source_1->name[i]);
      ++count;
    }
  }

  if (coord_type_2==0) {
    for(i=0;i<source_2->nat;i++) {
      target->xyz[count][0]=source_2->xyz[i][0];
      target->xyz[count][1]=source_2->xyz[i][1];
      target->xyz[count][2]=source_2->xyz[i][2];
      sprintf(target->name[count],"%s",source_2->name[i]);
      ++count;
    }
  } else {
     for(i=0;i<source_2->nat;i++) {
      target->xyz_frac[count][0]=source_2->xyz_frac[i][0];
      target->xyz_frac[count][1]=source_2->xyz_frac[i][1];
      target->xyz_frac[count][2]=source_2->xyz_frac[i][2];
      sprintf(target->name[count],"%s",source_2->name[i]);
      ++count;
    }
  }
  
  target->nat=source_1->nat+source_2->nat;
  if (target->nat != count)
    printf("warning: number of atom mismatch in at least one of the structures");

  // convert to coord using new lattice constants
  if ( coord_type_1 == 0 && coord_type_2 == 0 ) {
    target->xyz2frac();
  } else {
    target->frac2xyz();
  }
} 








void average_cell( double alatt_0[3][3], double alatt_1[3][3], double alatt_2[3][3]) {
  /*
	function to calculate the average of two sets of lattices
	alatt_? are 3x3 matrices
	subscript: 0=target, 1/2=source
  */
  int i;
  double cell_para_0[6];
  double cell_para_1[6];
  double cell_para_2[6];
  // convert 3x3 matrix to 6 parameters
  alatt2cellp(alatt_1,cell_para_1); 
  alatt2cellp(alatt_2,cell_para_2);

  //for(i=1;i<6;i++) {
  //  cell_para_0[i]=(cell_para_1[i]+cell_para_2[i])/2;
  //}

  // 2D special cases
  cell_para_0[0]=(cell_para_1[0]+cell_para_2[0])/2;
  cell_para_0[1]=(cell_para_1[1]+cell_para_2[1])/2;
  cell_para_0[3]=(cell_para_1[3]+cell_para_2[3])/2;
  cell_para_0[2]=max(cell_para_1[2],cell_para_2[2]);

  cellp2alatt(cell_para_0,alatt_0);

}

void alatt2cellp (double alatt[3][3], double cellp[6]) {
  int i,j;
  /*
  for (i=0;i<3;i++) {
    cellp[i]=vec_norm(alatt[i],3);
  }
  cellp[3]=vec_dot(alatt[0],alatt[1],3)/(vec_norm(alatt[0],3)*vec_norm(alatt[1],3));
  cellp[4]=vec_dot(alatt[0],alatt[2],3)/(vec_norm(alatt[0],3)*vec_norm(alatt[2],3));
  cellp[5]=vec_dot(alatt[1],alatt[2],3)/(vec_norm(alatt[1],3)*vec_norm(alatt[2],3));
  for (i=3;i<6;i++) {
    cellp[i]=acos(cellp[i])*180/PI;
  }
  */

  //2D special case
  for (i=0;i<2;i++) {
    cellp[i]=vec_norm(alatt[i],3);
  }
  cellp[2]=alatt[2][2];
  cellp[3]=vec_dot(alatt[0],alatt[1],3)/(vec_norm(alatt[0],3)*vec_norm(alatt[1],3));
  cellp[3]=acos(cellp[3])*180/PI;
  cellp[4]=90; cellp[5]=90;
  
  
  if ( DEBUG_AIMS ) {
    for (i=0;i<6;i++) {
      printf("%f ",cellp[i]);
    }
    printf("\n");
  }

}


void cellp2alatt (double cellp[6], double alatt[3][3]) {
  int i,j;

  for (i=3;i<6;i++) {
    cellp[i]=cellp[i]/180*PI;
  }

  alatt[0][0]=cellp[0];
  alatt[0][1]=0;
  alatt[0][2]=0;

  alatt[1][0]=cellp[1]*cos(cellp[3]);
  alatt[1][1]=cellp[1]*sin(cellp[3]);
  alatt[1][2]=0;

  // 2D special case
  alatt[2][2]=cellp[2];
  alatt[2][1]=0;
  alatt[2][0]=0;

  for (i=3;i<6;i++) {
    cellp[i]=cellp[i]*180/PI;
  }
 

  if ( DEBUG_AIMS ) {
    for (i=0;i<2;i++) {
      for (j=0;j<2;j++) {
        printf("%f ",alatt[i][j]);
      }
      printf("\n");
    }
  }
}


double vec_norm( double * vec, int n) {
  int i;
  double total=0;
  for (i=0;i<n;i++) {
    total += vec[i]*vec[i];
  }
  return sqrt(total);
}


double vec_dot( double * vec1, double * vec2, int n) {
  int i;
  double total=0;
  for (i=0;i<n;i++) {
    total += vec1[i]*vec2[i];
  }
  return total;
}













