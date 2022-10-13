#!/usr/bin/gawk -f
#
# usage:
# gmxstrain.awk -v id=ID -v namx=NMAX -v dd=DD -v to=T0
#
# this expects a pdb file id.pdb and the corresponding gmx topology file named id.top to be present
# other variables to set on the command line:
# nmax ... number of strain values
# dd ..... strain increment, stress tensor will be calculated at strain = -nmax*dd, ..., -dd, 0, dd, ..., nmax*dd
#          for each of 3 compressions and 3 shears.
# t0 ..... time (in ps) to discard from trajectory when calculating average stress tensors
#
# number of time steps to run each MD simulation and temperature need to be set in the mdp files
#
BEGIN {
  pi=3.141592653589793238462643383279502884197;
}{
#
  print "read", $0;
#
} END {
#
# do NPT run to establish lattice parameters at P=0
#
  system("gmx grompp -f mdnpt.mdp -c "id".pdb -p "id".top -o "id"-npt.tpr");
  system("gmx mdrun -v -deffnm "id"-npt -c tmp-"id"-npt.pdb > er-"id"-npt 2>&1");
#
  system("cat 0 | gmx trjconv -f tmp-"id"-npt.pdb -s "id"-npt.tpr -pbc nojump -o "id"-npt.pdb");
#
  system("gmx trjconv -f "id"-npt.trr -o lat-"id"-npt.pdb -s "id"-npt.tpr -n 1.ndx");
  system("grep CRYS lat-"id"-npt.pdb > crys-"id"-npt");
#
  afn="crys-"id"-npt";
  while ((getline line < afn) > 0) {
    split(line,aa);
    if(length(aa)>6) {
      nread++;
      aeq[nread]=aa[2];
      beq[nread]=aa[3];
      ceq[nread]=aa[4];
      alfeq[nread]=aa[5];
      beteq[nread]=aa[6];
      gameq[nread]=aa[7];
    }
  }

  close("crys-"id"-npt");
#
#
# discard first 20% of simulation time and calc average lattice parameters
# a,b,c,ala,beta,gamma from the rest
#
  for(i=int(nread/5);i<=nread;i++) {
    nnread++;
    sa+=aeq[i];
    sb+=beq[i];
    sc+=ceq[i];
    salf+=alfeq[i];
    sbet+=beteq[i];
    sgam+=gameq[i];
  }
  a=sa/nnread;
  b=sb/nnread;
  c=sc/nnread;
  alf=salf/nnread;
  bet=sbet/nnread;
  gam=sgam/nnread;
#
  lx=a;
  xy=b*cos(gam/180.0*pi);
  xz=c*cos(bet/180.0*pi);
  ly=sqrt(b*b-xy*xy);
  yz=(b*c*cos(alf/180.0*pi)-xy*xz)/ly;
  lz=sqrt(c*c-xz*xz-yz*yz);
#
# do NVT run for stress tensor in box with avg lattice parameters at P=0
#
  printf "CRYST1%9.3f%9.3f%9.3f%8.2f%8.2f%8.2f P 1       1\n",a,b,c,alf,bet,gam > "in-"id"-nvt.pdb";
  system("grep ATOM "id"-npt.pdb >> in-"id"-nvt.pdb");
#
  system("gmx grompp -f mdnvt.mdp -c in-"id"-nvt.pdb -p "id".top -o "id"-nvt.tpr");
  system("gmx mdrun -v -deffnm "id"-nvt -c tmp-"id"-nvt.pdb > er-"id"-nvt 2>&1");
#
  system("cat 0 | gmx trjconv -f tmp-"id"-nvt.pdb -s "id"-nvt.tpr -pbc nojump -o "id"-nvt.pdb");
  system("grep ATOM "id"-nvt.pdb > 0-"id".pdb");
  fncoor = "0-"id".pdb";
#
  system("cat 4stress | gmx energy -f "id"-nvt.edr -xvg none -o stress-"id"-nvt");
  nread=0;
  afn="stress-"id"-nvt.xvg";
  while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sss0[1]+=aa[2];
        sss0[2]+=aa[3];
        sss0[3]+=aa[4];
        sss0[4]+=aa[5];
        sss0[5]+=aa[6];
        sss0[6]+=aa[7];
        sss0[7]+=aa[8];
        sss0[8]+=aa[9];
        sss0[9]+=aa[10];
      }
    }
    sss0[1]/=nread;
    sss0[2]/=nread;
    sss0[3]/=nread;
    sss0[4]/=nread;
    sss0[5]/=nread;
    sss0[6]/=nread;
    sss0[7]/=nread;
    sss0[8]/=nread;
    sss0[9]/=nread;
#
#
#
  for(n=1;n<=nmax;n++) {
#
    delta = n*dd;
    lastdelta=(n-1)*dd;
#
# compress x
#
    new_lx = lx-delta*lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy-delta*xy;
    new_xz = xz-delta*xz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("sq1-%s-%06d", id, delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("sq1-%s-%06d.pdb", id, lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        ssq1[1]+=aa[2];
        ssq1[2]+=aa[3];
        ssq1[3]+=aa[4];
        ssq1[4]+=aa[5];
        ssq1[5]+=aa[6];
        ssq1[6]+=aa[7];
        ssq1[7]+=aa[8];
        ssq1[8]+=aa[9];
        ssq1[9]+=aa[10];
      }
    }
    ssq1[1]/=nread;
    ssq1[2]/=nread;
    ssq1[3]/=nread;
    ssq1[4]/=nread;
    ssq1[5]/=nread;
    ssq1[6]/=nread;
    ssq1[7]/=nread;
    ssq1[8]/=nread;
    ssq1[9]/=nread;
#
# compress y
#
    new_lx = lx;
    new_ly = ly-delta*ly;
    new_lz = lz;
    new_xy = xy;
    new_xz = xz;
    new_yz = yz-delta*yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("sq2-%s-%06d", id, delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("sq2-%s-%06d.pdb", id, lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        ssq2[1]+=aa[2];
        ssq2[2]+=aa[3];
        ssq2[3]+=aa[4];
        ssq2[4]+=aa[5];
        ssq2[5]+=aa[6];
        ssq2[6]+=aa[7];
        ssq2[7]+=aa[8];
        ssq2[8]+=aa[9];
        ssq2[9]+=aa[10];
      }
    }
    ssq2[1]/=nread;
    ssq2[2]/=nread;
    ssq2[3]/=nread;
    ssq2[4]/=nread;
    ssq2[5]/=nread;
    ssq2[6]/=nread;
    ssq2[7]/=nread;
    ssq2[8]/=nread;
    ssq2[9]/=nread;

#
# compress z
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz-delta*lz;
    new_xy = xy;
    new_xz = xz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("sq3-%s-%06d", id, delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("sq3-%s-%06d.pdb", id, lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        ssq3[1]+=aa[2];
        ssq3[2]+=aa[3];
        ssq3[3]+=aa[4];
        ssq3[4]+=aa[5];
        ssq3[5]+=aa[6];
        ssq3[6]+=aa[7];
        ssq3[7]+=aa[8];
        ssq3[8]+=aa[9];
        ssq3[9]+=aa[10];
      }
    }
    ssq3[1]/=nread;
    ssq3[2]/=nread;
    ssq3[3]/=nread;
    ssq3[4]/=nread;
    ssq3[5]/=nread;
    ssq3[6]/=nread;
    ssq3[7]/=nread;
    ssq3[8]/=nread;
    ssq3[9]/=nread;
#
# shear 1 (4) rotate around x towards negative y
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy;
    new_xz = xz;
    new_yz = yz-delta*lz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("sq4-%s-%06d", id, delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("sq4-%s-%06d.pdb", id, lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        ssq4[1]+=aa[2];
        ssq4[2]+=aa[3];
        ssq4[3]+=aa[4];
        ssq4[4]+=aa[5];
        ssq4[5]+=aa[6];
        ssq4[6]+=aa[7];
        ssq4[7]+=aa[8];
        ssq4[8]+=aa[9];
        ssq4[9]+=aa[10];
      }
    }
    ssq4[1]/=nread;
    ssq4[2]/=nread;
    ssq4[3]/=nread;
    ssq4[4]/=nread;
    ssq4[5]/=nread;
    ssq4[6]/=nread;
    ssq4[7]/=nread;
    ssq4[8]/=nread;
    ssq4[9]/=nread;
#
# shear 2 (5) rotate around y towards negative x
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy;
    new_xz = xz-delta*lz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("sq5-%s-%06d", id, delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("sq5-%s-%06d.pdb", id, lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        ssq5[1]+=aa[2];
        ssq5[2]+=aa[3];
        ssq5[3]+=aa[4];
        ssq5[4]+=aa[5];
        ssq5[5]+=aa[6];
        ssq5[6]+=aa[7];
        ssq5[7]+=aa[8];
        ssq5[8]+=aa[9];
        ssq5[9]+=aa[10];
      }
    }
    ssq5[1]/=nread;
    ssq5[2]/=nread;
    ssq5[3]/=nread;
    ssq5[4]/=nread;
    ssq5[5]/=nread;
    ssq5[6]/=nread;
    ssq5[7]/=nread;
    ssq5[8]/=nread;
    ssq5[9]/=nread;
#
# shear 3 (6) rotate around z towards negative x
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy-delta*ly;
    new_xz = xz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("sq6-%s-%06d", id, delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("sq6-%s-%06d.pdb", id, lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        ssq6[1]+=aa[2];
        ssq6[2]+=aa[3];
        ssq6[3]+=aa[4];
        ssq6[4]+=aa[5];
        ssq6[5]+=aa[6];
        ssq6[6]+=aa[7];
        ssq6[7]+=aa[8];
        ssq6[8]+=aa[9];
        ssq6[9]+=aa[10];
      }
    }
    ssq6[1]/=nread;
    ssq6[2]/=nread;
    ssq6[3]/=nread;
    ssq6[4]/=nread;
    ssq6[5]/=nread;
    ssq6[6]/=nread;
    ssq6[7]/=nread;
    ssq6[8]/=nread;
    ssq6[9]/=nread;
#
# positive deformations
#
    delta = -n*dd;
    lastdelta=-(n-1)*dd;
#
# expand x
#
    new_lx = lx-delta*lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy-delta*xy;
    new_xz = xz-delta*xz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("ex1-%s-%06d", id, -delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("ex1-%s-%06d.pdb", id, -lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sex1[1]+=aa[2];
        sex1[2]+=aa[3];
        sex1[3]+=aa[4];
        sex1[4]+=aa[5];
        sex1[5]+=aa[6];
        sex1[6]+=aa[7];
        sex1[7]+=aa[8];
        sex1[8]+=aa[9];
        sex1[9]+=aa[10];
      }
    }
    sex1[1]/=nread;
    sex1[2]/=nread;
    sex1[3]/=nread;
    sex1[4]/=nread;
    sex1[5]/=nread;
    sex1[6]/=nread;
    sex1[7]/=nread;
    sex1[8]/=nread;
    sex1[9]/=nread;
#
# expand y
#
    new_lx = lx;
    new_ly = ly-delta*ly;
    new_lz = lz;
    new_xy = xy;
    new_xz = xz;
    new_yz = yz-delta*yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("ex2-%s-%06d", id, -delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("ex2-%s-%06d.pdb", id, -lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sex2[1]+=aa[2];
        sex2[2]+=aa[3];
        sex2[3]+=aa[4];
        sex2[4]+=aa[5];
        sex2[5]+=aa[6];
        sex2[6]+=aa[7];
        sex2[7]+=aa[8];
        sex2[8]+=aa[9];
        sex2[9]+=aa[10];
      }
    }
    sex2[1]/=nread;
    sex2[2]/=nread;
    sex2[3]/=nread;
    sex2[4]/=nread;
    sex2[5]/=nread;
    sex2[6]/=nread;
    sex2[7]/=nread;
    sex2[8]/=nread;
    sex2[9]/=nread;
#
# expand z
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz-delta*lz;
    new_xy = xy;
    new_xz = xz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("ex3-%s-%06d", id, -delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("ex3-%s-%06d.pdb", id, -lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sex3[1]+=aa[2];
        sex3[2]+=aa[3];
        sex3[3]+=aa[4];
        sex3[4]+=aa[5];
        sex3[5]+=aa[6];
        sex3[6]+=aa[7];
        sex3[7]+=aa[8];
        sex3[8]+=aa[9];
        sex3[9]+=aa[10];
      }
    }
    sex3[1]/=nread;
    sex3[2]/=nread;
    sex3[3]/=nread;
    sex3[4]/=nread;
    sex3[5]/=nread;
    sex3[6]/=nread;
    sex3[7]/=nread;
    sex3[8]/=nread;
    sex3[9]/=nread;
#
# shear 1 (4) rotate around x towards positive y
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy;
    new_xz = xz;
    new_yz = yz-delta*lz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("ex4-%s-%06d", id, -delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("ex4-%s-%06d.pdb", id, -lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sex4[1]+=aa[2];
        sex4[2]+=aa[3];
        sex4[3]+=aa[4];
        sex4[4]+=aa[5];
        sex4[5]+=aa[6];
        sex4[6]+=aa[7];
        sex4[7]+=aa[8];
        sex4[8]+=aa[9];
        sex4[9]+=aa[10];
      }
    }
    sex4[1]/=nread;
    sex4[2]/=nread;
    sex4[3]/=nread;
    sex4[4]/=nread;
    sex4[5]/=nread;
    sex4[6]/=nread;
    sex4[7]/=nread;
    sex4[8]/=nread;
    sex4[9]/=nread;
#
# shear 2 (5) rotate around y towards positive x
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy;
    new_xz = xz-delta*lz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("ex5-%s-%06d", id, -delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("ex5-%s-%06d.pdb", id, -lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sex5[1]+=aa[2];
        sex5[2]+=aa[3];
        sex5[3]+=aa[4];
        sex5[4]+=aa[5];
        sex5[5]+=aa[6];
        sex5[6]+=aa[7];
        sex5[7]+=aa[8];
        sex5[8]+=aa[9];
        sex5[9]+=aa[10];
      }
    }
    sex5[1]/=nread;
    sex5[2]/=nread;
    sex5[3]/=nread;
    sex5[4]/=nread;
    sex5[5]/=nread;
    sex5[6]/=nread;
    sex5[7]/=nread;
    sex5[8]/=nread;
    sex5[9]/=nread;
#
#
# shear 3 (6) rotate around z towards positive x
#
    new_lx = lx;
    new_ly = ly;
    new_lz = lz;
    new_xy = xy-delta*ly;
    new_xz = xz;
    new_yz = yz;
#
    new_a = new_lx;
    new_b = sqrt(new_ly*new_ly+new_xy*new_xy);
    new_c = sqrt(new_lz*new_lz+new_xz*new_xz+new_yz*new_yz);
    new_alfa = 180.0/pi*acos((new_xy*new_xz+new_ly*new_yz)/(new_b*new_c));
    new_beta = 180.0/pi*acos(new_xz/new_c);
    new_gama = 180.0/pi*acos(new_xy/new_b);
#
    cid=sprintf("ex6-%s-%06d", id, -delta*1000000);
    if(n==1) pdbcoords=fncoor; else pdbcoords=sprintf("ex6-%s-%06d.pdb", id, -lastdelta*1000000);
    printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",new_a,new_b,new_c,new_alfa,new_beta,new_gama > "in-"cid".pdb";
    system("grep -v CRYS "pdbcoords" >> in-"cid".pdb");
    system("gmx grompp -f mdnvt.mdp -c in-"cid".pdb -p "id".top -o "cid".tpr");
    system("gmx mdrun -v -deffnm "cid" -c "cid".pdb > er-"cid" 2>&1");
#
    system("cat 4stress | gmx energy -f "cid".edr -xvg none -o stress-"cid);
    nread=0;
    afn="stress-"cid".xvg";
    while ((getline line < afn) > 0) {
      split(line,aa);
      if(length(aa)==10 && aa[1]>t0) {
        nread++;
        sex6[1]+=aa[2];
        sex6[2]+=aa[3];
        sex6[3]+=aa[4];
        sex6[4]+=aa[5];
        sex6[5]+=aa[6];
        sex6[6]+=aa[7];
        sex6[7]+=aa[8];
        sex6[8]+=aa[9];
        sex6[9]+=aa[10];
      }
    }
    sex6[1]/=nread;
    sex6[2]/=nread;
    sex6[3]/=nread;
    sex6[4]/=nread;
    sex6[5]/=nread;
    sex6[6]/=nread;
    sex6[7]/=nread;
    sex6[8]/=nread;
    sex6[9]/=nread;
#
# calcualte elastic constants for strain = n*dd
#
    lx0=a;
    ly0=b;
    lz0=c;
#
    len0[1]=lx0;
    len0[2]=ly0;
    len0[3]=lz0;
    len0[4]=lz0;
    len0[5]=lz0;
    len0[6]=ly0;

    pxx0=sss0[1];
    pyy0=sss0[5];
    pzz0=sss0[9];
    pxy0=0.5*(sss0[2]+sss0[4]);
    pxz0=0.5*(sss0[3]+sss0[7]);
    pyz0=0.5*(sss0[6]+sss0[8]);

    pmxx[1]=ssq1[1];
    pmyy[1]=ssq1[5];
    pmzz[1]=ssq1[9];
    pmxy[1]=0.5*(ssq1[2]+ssq1[4]);
    pmxz[1]=0.5*(ssq1[3]+ssq1[7]);
    pmyz[1]=0.5*(ssq1[6]+ssq1[8]);
    ppxx[1]=sex1[1];
    ppyy[1]=sex1[5];
    ppzz[1]=sex1[9];
    ppxy[1]=0.5*(sex1[2]+sex1[4]);
    ppxz[1]=0.5*(sex1[3]+sex1[7]);
    ppyz[1]=0.5*(sex1[6]+sex1[8]);

    pmxx[2]=ssq2[1];
    pmyy[2]=ssq2[5];
    pmzz[2]=ssq2[9];
    pmxy[2]=0.5*(ssq2[2]+ssq2[4]);
    pmxz[2]=0.5*(ssq2[3]+ssq2[7]);
    pmyz[2]=0.5*(ssq2[6]+ssq2[8]);
    ppxx[2]=sex2[1];
    ppyy[2]=sex2[5];
    ppzz[2]=sex2[9];
    ppxy[2]=0.5*(sex2[2]+sex2[4]);
    ppxz[2]=0.5*(sex2[3]+sex2[7]);
    ppyz[2]=0.5*(sex2[6]+sex2[8]);

    pmxx[3]=ssq3[1];
    pmyy[3]=ssq3[5];
    pmzz[3]=ssq3[9];
    pmxy[3]=0.5*(ssq3[2]+ssq3[4]);
    pmxz[3]=0.5*(ssq3[3]+ssq3[7]);
    pmyz[3]=0.5*(ssq3[6]+ssq3[8]);
    ppxx[3]=sex3[1];
    ppyy[3]=sex3[5];
    ppzz[3]=sex3[9];
    ppxy[3]=0.5*(sex3[2]+sex3[4]);
    ppxz[3]=0.5*(sex3[3]+sex3[7]);
    ppyz[3]=0.5*(sex3[6]+sex3[8]);

    pmxx[4]=ssq4[1];
    pmyy[4]=ssq4[5];
    pmzz[4]=ssq4[9];
    pmxy[4]=0.5*(ssq4[2]+ssq4[4]);
    pmxz[4]=0.5*(ssq4[3]+ssq4[7]);
    pmyz[4]=0.5*(ssq4[6]+ssq4[8]);
    ppxx[4]=sex4[1];
    ppyy[4]=sex4[5];
    ppzz[4]=sex4[9];
    ppxy[4]=0.5*(sex4[2]+sex4[4]);
    ppxz[4]=0.5*(sex4[3]+sex4[7]);
    ppyz[4]=0.5*(sex4[6]+sex4[8]);

    pmxx[5]=ssq5[1];
    pmyy[5]=ssq5[5];
    pmzz[5]=ssq5[9];
    pmxy[5]=0.5*(ssq5[2]+ssq5[4]);
    pmxz[5]=0.5*(ssq5[3]+ssq5[7]);
    pmyz[5]=0.5*(ssq5[6]+ssq5[8]);
    ppxx[5]=sex5[1];
    ppyy[5]=sex5[5];
    ppzz[5]=sex5[9];
    ppxy[5]=0.5*(sex5[2]+sex5[4]);
    ppxz[5]=0.5*(sex5[3]+sex5[7]);
    ppyz[5]=0.5*(sex5[6]+sex5[8]);

    pmxx[6]=ssq6[1];
    pmyy[6]=ssq6[5];
    pmzz[6]=ssq6[9];
    pmxy[6]=0.5*(ssq6[2]+ssq6[4]);
    pmxz[6]=0.5*(ssq6[3]+ssq6[7]);
    pmyz[6]=0.5*(ssq6[6]+ssq6[8]);
    ppxx[6]=sex6[1];
    ppyy[6]=sex6[5];
    ppzz[6]=sex6[9];
    ppxy[6]=0.5*(sex6[2]+sex6[4]);
    ppxz[6]=0.5*(sex6[3]+sex6[7]);
    ppyz[6]=0.5*(sex6[6]+sex6[8]);
#
    cfac=0.0001;
    up=n*dd;
    for(dir=1;dir<=6;dir++) {
#
      delta=-up*len0[dir];
      deltaxy=-up*xy;
      deltaxz=-up*xz;
      deltayz=-up*yz;

      C1neg=-(pmxx[dir]-pxx0)/(delta/len0[dir])*cfac;
      C2neg=-(pmyy[dir]-pyy0)/(delta/len0[dir])*cfac;
      C3neg=-(pmzz[dir]-pzz0)/(delta/len0[dir])*cfac;
      C4neg=-(pmyz[dir]-pyz0)/(delta/len0[dir])*cfac;
      C5neg=-(pmxz[dir]-pxz0)/(delta/len0[dir])*cfac;
      C6neg=-(pmxy[dir]-pxy0)/(delta/len0[dir])*cfac;

      delta=up*len0[dir];
      deltaxy=up*xy;
      deltaxz=up*xz;
      deltayz=up*yz;

      C1pos=-(ppxx[dir]-pxx0)/(delta/len0[dir])*cfac;
      C2pos=-(ppyy[dir]-pyy0)/(delta/len0[dir])*cfac;
      C3pos=-(ppzz[dir]-pzz0)/(delta/len0[dir])*cfac;
      C4pos=-(ppyz[dir]-pyz0)/(delta/len0[dir])*cfac;
      C5pos=-(ppxz[dir]-pxz0)/(delta/len0[dir])*cfac;
      C6pos=-(ppxy[dir]-pxy0)/(delta/len0[dir])*cfac;

      C1[dir]=0.5*(C1neg+C1pos);
      C2[dir]=0.5*(C2neg+C2pos);
      C3[dir]=0.5*(C3neg+C3pos);
      C4[dir]=0.5*(C4neg+C4pos);
      C5[dir]=0.5*(C5neg+C5pos);
      C6[dir]=0.5*(C6neg+C6pos);
    }

    C11all = C1[1];
    C22all = C2[2];
    C33all = C3[3];
    C12all = 0.5*(C1[2]+C2[1]);
    C13all = 0.5*(C1[3]+C3[1]);
    C23all = 0.5*(C2[3]+C3[2]);
    C44all = C4[4];
    C55all = C5[5];
    C66all = C6[6];
    C14all = 0.5*(C1[4]+C4[1]);
    C15all = 0.5*(C1[5]+C5[1]);
    C16all = 0.5*(C1[6]+C6[1]);
    C24all = 0.5*(C2[4]+C4[2]);
    C25all = 0.5*(C2[5]+C5[2]);
    C26all = 0.5*(C2[6]+C6[2]);
    C34all = 0.5*(C3[4]+C4[3]);
    C35all = 0.5*(C3[5]+C5[3]);
    C36all = 0.5*(C3[6]+C6[3]);
    C45all = 0.5*(C4[5]+C5[4]);
    C46all = 0.5*(C4[6]+C6[4]);
    C56all = 0.5*(C5[6]+C6[5]);

    ofn=sprintf("ec-%s-%06d", id, up*1000000);

    printf "C11 = %12.4f\n", C11all  > ofn;
    printf "C22 = %12.4f\n", C22all >> ofn;
    printf "C33 = %12.4f\n", C33all >> ofn;
    printf "C12 = %12.4f\n", C12all >> ofn;
    printf "C13 = %12.4f\n", C13all >> ofn;
    printf "C23 = %12.4f\n", C23all >> ofn;
    printf "C44 = %12.4f\n", C44all >> ofn;
    printf "C55 = %12.4f\n", C55all >> ofn;
    printf "C66 = %12.4f\n", C66all >> ofn;
    printf "C14 = %12.4f\n", C14all >> ofn;
    printf "C15 = %12.4f\n", C15all >> ofn;
    printf "C16 = %12.4f\n", C16all >> ofn;
    printf "C24 = %12.4f\n", C24all >> ofn;
    printf "C25 = %12.4f\n", C25all >> ofn;
    printf "C26 = %12.4f\n", C26all >> ofn;
    printf "C34 = %12.4f\n", C34all >> ofn;
    printf "C35 = %12.4f\n", C35all >> ofn;
    printf "C36 = %12.4f\n", C36all >> ofn;
    printf "C45 = %12.4f\n", C45all >> ofn;
    printf "C46 = %12.4f\n", C46all >> ofn;
    printf "C56 = %12.4f\n", C56all >> ofn;
  }
}

function acos(x) { return atan2(sqrt(1-x*x), x) }
