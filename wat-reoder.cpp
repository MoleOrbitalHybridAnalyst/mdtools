#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*-----

step 1: read in control variables from file "control"
step 2: read in OW list from file "OList" and build HW list
step 3: read in coordinate from argv[1]
step 4: reorder hydrogens to build water
step 5: build hydronium
step 6: put excess proton at the end
step 7: check the coo- donor
step 8: output xyz

THIS IS ONLY FOR ONE PROTON!
NO MULTIPRO!

-----*/

int natom, nOW, nHW, OC1, OC2, OT, HT;
int *OList, *HList, *numH;
double *x, *xsave;
double pbcx, pbcy, pbcz, hpbcx, hpbcy, hpbcz;

// copy j-th data of x to i-th data of xsave
#define SAVEX(i,j) { int _t = (i)*3, _tt = (j)*3; xsave[_t]=x[_tt]; xsave[_t+1]=x[_tt+1]; xsave[_t+2]=x[_tt+2]; }
// copy j-th data of xsave to i-th data of x
#define LOADX(i,j) { int _t = (i)*3, _tt = (j)*3; x[_t]=xsave[_tt]; x[_t+1]=xsave[_tt+1]; x[_t+2]=xsave[_tt+2]; }

// All atoms must be wrapped into a single PBC cell!
double r2ij(double *x1, double *x2)
{
  double dx = x1[0] - x2[0];
  double dy = x1[1] - x2[1];
  double dz = x1[2] - x2[2];
  
  while(dx<-hpbcx) dx+= pbcx; while(dx>hpbcx) dx-= pbcx;
  while(dy<-hpbcy) dy+= pbcy; while(dy>hpbcy) dy-= pbcy;
  while(dz<-hpbcz) dz+= pbcz; while(dz>hpbcz) dz-= pbcz;
  
  return dx*dx+dy*dy+dz*dz;
}

int main(int argc, char **argv)
{
  // step 1: read control
  
  FILE *fp;
  
  fp = fopen("control","r");
  fscanf(fp,"%d %d %d %d %d %d\n", &natom, &nOW, &OC1, &OC2, &OT, &HT);
  OC1--; OC2--; OT--; HT--;
  fscanf(fp,"%lf %lf %lf\n", &pbcx, &pbcy, &pbcz);
  hpbcx = pbcx*0.5; hpbcy = pbcy*0.5; hpbcz = pbcz*0.5;
  fclose(fp);
  
  x = new double [natom*3];
  xsave = new double [natom*3];
  OList = new int [nOW];
  numH  = new int [nOW];
  
  nHW = nOW*2+1;
  HList = new int [nHW]; 
  HList[0] = HT;
  
  // step 2: read OW list
  
  fp = fopen("OList","r");
  for(int i=0; i<nOW; i++)
  {
    fscanf(fp,"%d\n", OList+i);
    OList[i]--;
    numH[i]=0;
    
    HList[i*2+1] = OList[i]+1;
    HList[i*2+2] = OList[i]+2;
  }
  fclose(fp);
  
  // step 3: read coord
  
  fp = fopen(argv[1],"r");
  char line[100];
  // discard first two lines
  fgets(line,99,fp);
  fgets(line,99,fp);
    
  for(int i=0; i<natom; i++)
  {
    fscanf(fp, "%s %lf %lf %lf\n", line, x+i*3, x+i*3+1, x+i*3+2);
    SAVEX(i,i);
  }
  fclose(fp);
  
  // step 4: build water
  
  int hh = -1;
  double cut = 1.2 * 1.2;
  
  // For each hydrogen, find oxygens within "cut"
  for(int i=0; i<nHW; i++)
  {
    double* xH = xsave + HList[i]*3;
    int target = -1;
    
    for(int j=0; j<nOW; j++) if(numH[j]<2)
    {
      int iO = OList[j];
      double r2 = r2ij(xsave+iO*3, xH);     
      if(r2<cut) { target = j; break; }
    }
    
    if(target==-1) 
    { 
      printf("EXCESS PROTON %d\n", HList[i]+1);
      if(hh == -1) hh = HList[i];
      else {
		printf("WRONG\n");
		return 255;
	  }
    }
    else
    {
      LOADX(OList[target]+numH[target]+1, HList[i]);
      numH[target]++;
    }
  } 
  
  LOADX(HT, hh);
  
  // step 5: build hydronium
  
  int oo = -1;
  double r2min = 100.0;
  
  // For the excess proton, find the closest oxygen
  for(int i=0; i<nOW; i++)
  {
    int iO = OList[i];
    double r2 = r2ij(xsave+iO*3, xsave+hh*3);
    if(r2<r2min)
    {
      r2min = r2;
      oo = iO;
    }
  }
  
  for(int i=0; i<natom; i++) SAVEX(i,i);
  
  printf("EXCESS OXYGEN %d\n", oo+1);
  // SWAP
  SAVEX(oo, OT);  SAVEX(oo+1, OT+1);  SAVEX(oo+2, OT+2);
  SAVEX(OT, oo);  SAVEX(OT+1, oo+1);  SAVEX(OT+2, oo+2);
  
  // step 6: sort again

  for(int i=0; i<natom; i++) LOADX(i,i);
//  for(int i=0; i<HT; i++) LOADX(i, i);
//  for(int i=HT; i<natom-1; i++) LOADX(i, i+1);
//  // Put excess proton to the last one?
//  LOADX(natom-1, HT);

//  // step 7: remap oc1 oc2
//  
//  int list1[3]; list1[0]=natom-3; list1[1]=natom-2; list1[2]=natom-1; // hydronium hydrogens
//  int list2[2]; list2[0]=OC1; list2[1]=OC2;
//  
//  int ox, hx;
//  double rmin=1000;
//  
//  for(int i=0; i<3; i++) for(int j=0; j<2; j++)
//  {
//    double rr = r2ij(x+list1[i]*3, x+list2[j]*3);
//    if(rr<rmin)
//    {
//      rmin = rr;
//      hx = list1[i];
//      ox = list2[j];
//    }
//  }
//  
//#define SWAP(i,j) { double _t, *xx=x+(i)*3, *yy=x+(j)*3; _t=xx[0]; xx[0]=yy[0]; yy[0]=_t; _t=xx[1]; xx[1]=yy[1]; yy[1]=_t; _t=xx[2]; xx[2]=yy[2]; yy[2]=_t; }
//  
//  if(natom-1!=hx) SWAP(natom-1, hx);
//  if(ox!=OC1) SWAP(ox,OC1);
//  double roo = r2ij(x+OC1*3, x+(natom-4)*3);
//  printf("%s %lf %lf\n", argv[1], sqrt(roo), sqrt(rmin));  
  
  // step 8: output xyz
  
  fp = fopen("output.xyz","w");
  fprintf(fp,"%d\n\n", natom);
  for(int i=0; i<natom; i++) fprintf(fp, "X %10lf %10lf %10lf\n", x[i*3], x[i*3+1], x[i*3+2]);
  fclose(fp);
  
  return 0;
}
