/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "vec.h"
#include "copyrite.h"
#include "smalloc.h"
#include "macros.h"
#include "symtab.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "pdb2top.h"
#include "gpp_nextnb.h"
#include "topdirs.h"
#include "toputil.h"
#include "h_db.h"
#include "pgutil.h"
#include "resall.h"
#include "topio.h"
#include "string2.h"
#include "physics.h"
#include "pdbio.h"
#include "gen_ad.h"
#include "filenm.h"
#include "index.h"
#include "gen_vsite.h"
#include "add_par.h"
#include "toputil.h"
#include "fflibutil.h"
#include "strdb.h"

/* this must correspond to enum in pdb2top.h */
const char *hh[ehisNR]   = { "HISA", "HISB", "HISH", "HIS1" };

static int missing_atoms(t_restp *rp, int resind,
			 t_atoms *at, int i0, int i, bool bCTer)
{
  int  j,k,nmiss;
  char *name;
  bool bFound, bRet;

  nmiss = 0;
  for (j=0; j<rp->natom; j++) {
    name=*(rp->atomname[j]);
    /* if ((name[0]!='H') && (name[0]!='h') && (!bCTer || (name[0]!='O'))) { */
    bFound=FALSE;
    for (k=i0; k<i; k++) 
      bFound=(bFound || !strcasecmp(*(at->atomname[k]),name));
    if (!bFound) {
      nmiss++;
      fprintf(stderr,"\nWARNING: "
	      "atom %s is missing in residue %s %d in the pdb file\n",
	      name,*(at->resinfo[resind].name),at->resinfo[resind].nr);
      if (name[0]=='H' || name[0]=='h')
	fprintf(stderr,"         You might need to add atom %s to the hydrogen database of residue %s\n"
		       "         in the file ff???.hdb (see the manual)\n",
		name,*(at->resinfo[resind].name));
      fprintf(stderr,"\n");
      }
      /* } */
  }
  
  return nmiss;
}

bool is_int(double x)
{
  const double tol = 1e-4;
  int   ix;
  
  if (x < 0)
    x=-x;
  ix=gmx_nint(x);
  
  return (fabs(x-ix) < tol);
}

void
choose_ff(const char *ffsel,
          char *forcefield, int ff_maxlen,
          char *ffdir, int ffdir_maxlen)
{
    int  nff;
    char **ffdirs,**ffs,*ptr;
    int  i,sel;
    char buf[STRLEN],*doc_dir;
    FILE *fp;
    char *pret;

    nff = fflib_search_file_in_dirend(fflib_forcefield_itp(),
                                      fflib_forcefield_dir_ext(),
                                      &ffdirs);

    if (nff == 0)
    {
        gmx_fatal(FARGS,"No force fields found (files with name '%s' in subdirectories ending on '%s')",
                  fflib_forcefield_itp(),fflib_forcefield_dir_ext());
    }

    /* Store the force field names in ffs */
    snew(ffs,nff);
    for(i=0; i<nff; i++)
    {
        /* Remove the path from the ffdir name */
        ptr = strrchr(ffdirs[i],DIR_SEPARATOR);
        if (ptr == 0)
        {
            ffs[i] = strdup(ffdirs[i]);
        }
        else
        {
            ffs[i] = strdup(ptr+1);
        }
        /* Remove the extension from the ffdir name */
        ffs[i][strlen(ffs[i])-strlen(fflib_forcefield_dir_ext())] = '\0';
    }

    if (ffsel != NULL)
    {
        sel = -1;
        for(i=0; i<nff; i++)
        {
            if (strcmp(ffs[i],ffsel) == 0)
            {
                sel = i;
            }
        }
        if (sel == -1)
        {
            gmx_fatal(FARGS,"Could not find force field '%s'",ffsel);
        }
    }
    else if (nff > 1)
    {
        printf("\nSelect the Force Field:\n");
        for(i=0; (i<nff); i++)
        {
            sprintf(buf,"%s%c%s",
                    ffdirs[i],DIR_SEPARATOR,fflib_forcefield_doc());
            doc_dir = low_gmxlibfn(buf,FALSE);
            if (doc_dir != NULL)
            {
                /* We don't use fflib_open, because we don't want printf's */
                fp = ffopen(doc_dir,"r");
                get_a_line(fp,buf,STRLEN);
                ffclose(fp);
                sfree(doc_dir);
                printf("%2d: %s\n",i,buf);
            }
            else
            {
                printf("%2d: %s\n",i,ffs[i]);
            }
        }
        do
        {
            pret = fgets(buf,STRLEN,stdin);
            
            if(pret != NULL)
            {
                sscanf(buf,"%d",&sel);
            }
        }
        while ( pret==NULL || (sel < 0) || (sel >= nff));
    }
    else
    {
        sel = 0;
    }

    if (strlen(ffs[sel]) >= ff_maxlen)
    {
        gmx_fatal(FARGS,"Length of force field name (%d) >= maxlen (%d)",
                  strlen(ffs[sel]),ff_maxlen);
    }
    strcpy(forcefield,ffs[sel]);

    if (strlen(ffdirs[sel]) >= ffdir_maxlen)
    {
        gmx_fatal(FARGS,"Length of force field dir (%d) >= maxlen (%d)",
                  strlen(ffdirs[sel]),ffdir_maxlen);
    }
    strcpy(ffdir,ffdirs[sel]);

    for(i=0; (i<nff); i++)
    {
        sfree(ffdirs[i]);
        sfree(ffs[i]);
    }
    sfree(ffdirs);
    sfree(ffs);
}

static int name2type(t_atoms *at, int **cgnr, gpp_atomtype_t atype, 
		     t_restp restp[])
{
  int     i,j,prevresind,resind,i0,prevcg,cg,curcg;
  char    *name;
  bool    bProt, bNterm;
  double  qt;
  int     nmissat;
  t_aa_names *aan;
  
  nmissat = 0;

  resind=-1;
  bProt=FALSE;
  bNterm=FALSE;
  i0=0;
  snew(*cgnr,at->nr);
  qt=0;
  prevcg=NOTSET;
  curcg=0;
  cg=-1;
  j=NOTSET;
  aan = get_aa_names();
  
  for(i=0; (i<at->nr); i++) {
    prevresind=resind;
    if (at->atom[i].resind != resind) {
      resind = at->atom[i].resind;
      bProt = is_protein(aan,*(at->resinfo[resind].name));
      bNterm=bProt && (resind == 0);
      if (resind > 0)
	nmissat += 
	  missing_atoms(&restp[prevresind],prevresind,at,i0,i,
			(!bProt && is_protein(aan,restp[prevresind].resname)));
      i0=i;
    }
    if (at->atom[i].m == 0) {
      if (debug)
	fprintf(debug,"atom %d%s: curcg=%d, prevcg=%d, cg=%d\n",
		i+1,*(at->atomname[i]),curcg,prevcg,
		j==NOTSET ? NOTSET : restp[resind].cgnr[j]);
      qt=0;
      prevcg=cg;
      name=*(at->atomname[i]);
      j=search_jtype(&restp[resind],name,bNterm);
      at->atom[i].type = restp[resind].atom[j].type;
      at->atom[i].q    = restp[resind].atom[j].q;
      at->atom[i].m    = get_atomtype_massA(restp[resind].atom[j].type,
					    atype);
      cg = restp[resind].cgnr[j];
      if ( (cg != prevcg) || (resind != prevresind) )
	curcg++;
    } else {
      if (debug)
	fprintf(debug,"atom %d%s: curcg=%d, qt=%g, is_int=%d\n",
		i+1,*(at->atomname[i]),curcg,qt,is_int(qt));
      cg=-1;
      if (is_int(qt)) {
	qt=0;
	curcg++;
      }
      qt+=at->atom[i].q;
    }
    (*cgnr)[i]=curcg;
    at->atom[i].typeB = at->atom[i].type;
    at->atom[i].qB    = at->atom[i].q;
    at->atom[i].mB    = at->atom[i].m;
  }
  nmissat += missing_atoms(&restp[resind],resind,at,i0,i,
			   (!bProt || is_protein(aan,restp[resind].resname)));
  done_aa_names(&aan);
			   
  return nmissat;
}

static void print_top_heavy_H(FILE *out, real mHmult)
{
  if (mHmult == 2.0) 
    fprintf(out,"; Using deuterium instead of hydrogen\n\n");
  else if (mHmult == 4.0)
    fprintf(out,"#define HEAVY_H\n\n");
  else if (mHmult != 1.0)
    fprintf(stderr,"WARNING: unsupported proton mass multiplier (%g) "
	    "in pdb2top\n",mHmult);
}

void print_top_comment(FILE *out,const char *filename,const char *title,bool bITP)
{
  char tmp[256]; 
  
  nice_header(out,filename);
  fprintf(out,";\tThis is your %stopology file\n",bITP ? "include " : "");
  cool_quote(tmp,255,NULL);
  fprintf(out,";\t%s\n",title[0]?title:tmp);
  fprintf(out,";\n");
}

void print_top_header(FILE *out,const char *filename, 
                      const char *title,bool bITP,const char *ffdir,real mHmult)
{
    print_top_comment(out,filename,title,bITP);
    
    print_top_heavy_H(out, mHmult);
    fprintf(out,"; Include forcefield parameters\n");
    fprintf(out,"#include \"%s%c%s\"\n\n",
            ffdir,DIR_SEPARATOR,fflib_forcefield_itp());
}

static void print_top_posre(FILE *out,const char *pr)
{
  fprintf(out,"; Include Position restraint file\n");
  fprintf(out,"#ifdef POSRES\n");
  fprintf(out,"#include \"%s\"\n",pr);
  fprintf(out,"#endif\n\n");
}
  
static void print_top_water(FILE *out,const char *ffdir,const char *water)
{
    char buf[STRLEN];

  fprintf(out,"; Include water topology\n");
  fprintf(out,"#include \"%s%c%s.itp\"\n",ffdir,DIR_SEPARATOR,water);
  fprintf(out,"\n");
  fprintf(out,"#ifdef POSRES_WATER\n");
  fprintf(out,"; Position restraint for each water oxygen\n");
  fprintf(out,"[ position_restraints ]\n");
  fprintf(out,";%3s %5s %9s %10s %10s\n","i","funct","fcx","fcy","fcz");
  fprintf(out,"%4d %4d %10g %10g %10g\n",1,1,1000.0,1000.0,1000.0);
  fprintf(out,"#endif\n");
  fprintf(out,"\n");

    sprintf(buf,"%s%c%s",ffdir,DIR_SEPARATOR,"ions.itp");
    if (fflib_fexist(buf))
    {
        fprintf(out,"; Include topology for ions\n");
        fprintf(out,"#include \"%s%cions.itp\"\n",ffdir,DIR_SEPARATOR);
        fprintf(out,"\n");
    }
}

static void print_top_system(FILE *out, const char *title)
{
  fprintf(out,"[ %s ]\n",dir2str(d_system));
  fprintf(out,"; Name\n");
  fprintf(out,"%s\n\n",title[0]?title:"Protein");
}

void print_top_mols(FILE *out,
                    const char *title, const char *ffdir, const char *water,
                    int nincl, char **incls, int nmol, t_mols *mols)
{
  int  i;
  char *incl;

  if (nincl>0) {
    fprintf(out,"; Include chain topologies\n");
    for (i=0; (i<nincl); i++) {
        incl = strrchr(incls[i],DIR_SEPARATOR);
        if (incl == NULL) {
            incl = incls[i];
        } else {
            /* Remove the path from the include name */
            incl = incl + 1;
        }
      fprintf(out,"#include \"%s\"\n",incl);
    }
    fprintf(out,"\n");
  }

    if (water)
    {
      print_top_water(out,ffdir,water);
    }
    print_top_system(out, title);
  
  if (nmol) {
    fprintf(out,"[ %s ]\n",dir2str(d_molecules));
    fprintf(out,"; %-15s %5s\n","Compound","#mols");
    for (i=0; (i<nmol); i++)
      fprintf(out,"%-15s %5d\n",mols[i].name,mols[i].nr);
  }
}

void write_top(FILE *out, char *pr,char *molname,
	       t_atoms *at,int bts[],t_params plist[],t_excls excls[],
	       gpp_atomtype_t atype,int *cgnr, int nrexcl)
     /* NOTE: nrexcl is not the size of *excl! */
{
  if (at && atype && cgnr) {
    fprintf(out,"[ %s ]\n",dir2str(d_moleculetype));
    fprintf(out,"; %-15s %5s\n","Name","nrexcl");
    fprintf(out,"%-15s %5d\n\n",molname?molname:"Protein",nrexcl);
    
    print_atoms(out, atype, at, cgnr);
    print_bondeds(out,at->nr,d_bonds,      F_BONDS,    bts[ebtsBONDS], plist);
    print_bondeds(out,at->nr,d_constraints,F_CONSTR,   0,              plist);
    print_bondeds(out,at->nr,d_constraints,F_CONSTRNC, 0,              plist);
    print_bondeds(out,at->nr,d_pairs,      F_LJ14,     0,              plist);
    print_excl(out,at->nr,excls);
    print_bondeds(out,at->nr,d_angles,     F_ANGLES,   bts[ebtsANGLES],plist);
    print_bondeds(out,at->nr,d_dihedrals,  F_PDIHS,    bts[ebtsPDIHS], plist);
    print_bondeds(out,at->nr,d_dihedrals,  F_IDIHS,    bts[ebtsIDIHS], plist);
    print_bondeds(out,at->nr,d_cmap,       F_CMAP,     bts[ebtsCMAP],  plist);
    print_bondeds(out,at->nr,d_polarization,F_POLARIZATION,   0,       plist);
    print_bondeds(out,at->nr,d_thole_polarization,F_THOLE_POL,0,       plist);
    print_bondeds(out,at->nr,d_vsites2,    F_VSITE2,   0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3,   0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3FD, 0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3FAD,0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3OUT,0,              plist);
    print_bondeds(out,at->nr,d_vsites4,    F_VSITE4FD, 0,              plist);
    print_bondeds(out,at->nr,d_vsites4,    F_VSITE4FDN, 0,             plist);
    
    if (pr)
      print_top_posre(out,pr);
  }
}

static atom_id search_res_atom(const char *type,int resind,
			       int natom,t_atom at[],
			       char ** const *aname,
			       const char *bondtype,bool bAllowMissing)
{
  int i;

  for(i=0; (i<natom); i++)
    if (at[i].resind == resind)
      return search_atom(type,i,natom,at,aname,bondtype,bAllowMissing);
  
  return NO_ATID;
}

static void do_ssbonds(t_params *ps,int natoms,t_atom atom[],char **aname[],
		       int nssbonds,t_ssbond *ssbonds,bool bAllowMissing)
{
  int     i,ri,rj;
  atom_id ai,aj;
  
  for(i=0; (i<nssbonds); i++) {
    ri = ssbonds[i].res1;
    rj = ssbonds[i].res2;
    ai = search_res_atom(ssbonds[i].a1,ri,natoms,atom,aname,
			 "special bond",bAllowMissing);
    aj = search_res_atom(ssbonds[i].a2,rj,natoms,atom,aname,
			 "special bond",bAllowMissing);
    if ((ai == NO_ATID) || (aj == NO_ATID))
      gmx_fatal(FARGS,"Trying to make impossible special bond (%s-%s)!",
		  ssbonds[i].a1,ssbonds[i].a2);
    add_param(ps,ai,aj,NULL,NULL);
  }
}

static bool inter_res_bond(const t_rbonded *b)
{
    return (b->AI[0] == '-' || b->AI[0] == '+' ||
            b->AJ[0] == '-' || b->AJ[0] == '+');
}

static void at2bonds(t_params *psb, t_hackblock *hb,
                     int natoms, t_atom atom[], char **aname[], 
                     int nres, rvec x[], 
                     real long_bond_dist, real short_bond_dist,
                     bool bAllowMissing)
{
  int     resind,i,j,k;
  atom_id ai,aj;
  real    dist2, long_bond_dist2, short_bond_dist2;
  const char *ptr;

  long_bond_dist2  = sqr(long_bond_dist);
  short_bond_dist2 = sqr(short_bond_dist);

  if (debug)
    ptr = "bond";
  else
    ptr = "check";

  fprintf(stderr,"Making bonds...\n");
  i=0;
  for(resind=0; (resind < nres) && (i<natoms); resind++) {
    /* add bonds from list of bonded interactions */
    for(j=0; j < hb[resind].rb[ebtsBONDS].nb; j++) {
      /* Unfortunately we can not issue errors or warnings
       * for missing atoms in bonds, as the hydrogens and terminal atoms
       * have not been added yet.
       */
      ai=search_atom(hb[resind].rb[ebtsBONDS].b[j].AI,i,natoms,atom,aname,
		     ptr,TRUE);
      aj=search_atom(hb[resind].rb[ebtsBONDS].b[j].AJ,i,natoms,atom,aname,
		     ptr,TRUE);
      if (ai != NO_ATID && aj != NO_ATID) {
          dist2 = distance2(x[ai],x[aj]);
          if (dist2 > long_bond_dist2 )
          {
              fprintf(stderr,"Warning: Long Bond (%d-%d = %g nm)\n",
                      ai+1,aj+1,sqrt(dist2));
          }
          else if (dist2 < short_bond_dist2 )
          {
              fprintf(stderr,"Warning: Short Bond (%d-%d = %g nm)\n",
                      ai+1,aj+1,sqrt(dist2));
          }
          add_param(psb,ai,aj,NULL,hb[resind].rb[ebtsBONDS].b[j].s);
      }
    }
    /* add bonds from list of hacks (each added atom gets a bond) */
    while( (i<natoms) && (atom[i].resind == resind) ) {
      for(j=0; j < hb[resind].nhack; j++)
	if ( ( hb[resind].hack[j].tp > 0 ||
	       hb[resind].hack[j].oname==NULL ) &&
	     strcmp(hb[resind].hack[j].AI,*(aname[i])) == 0 ) {
	  switch(hb[resind].hack[j].tp) {
	  case 9:          /* COOH terminus */
	    add_param(psb,i,i+1,NULL,NULL);     /* C-O  */
	    add_param(psb,i,i+2,NULL,NULL);     /* C-OA */
	    add_param(psb,i+2,i+3,NULL,NULL);   /* OA-H */
	    break;
	  default:
	    for(k=0; (k<hb[resind].hack[j].nr); k++)
	      add_param(psb,i,i+k+1,NULL,NULL);
	  }
	}
      i++;
    }
    /* we're now at the start of the next residue */
  }
}

static int pcompar(const void *a, const void *b)
{
  t_param *pa,*pb;
  int     d;
  pa=(t_param *)a;
  pb=(t_param *)b;
  
  d = pa->AI - pb->AI;
  if (d == 0) 
    d = pa->AJ - pb->AJ;
  if (d == 0)
    return strlen(pb->s) - strlen(pa->s);
  else
    return d;
}

static void clean_bonds(t_params *ps)
{
  int     i,j;
  atom_id a;
  
  if (ps->nr > 0) {
    /* swap atomnumbers in bond if first larger than second: */
    for(i=0; (i<ps->nr); i++)
      if ( ps->param[i].AJ < ps->param[i].AI ) {
	a = ps->param[i].AI;
	ps->param[i].AI = ps->param[i].AJ;
	ps->param[i].AJ = a;
      }
    
    /* Sort bonds */
    qsort(ps->param,ps->nr,(size_t)sizeof(ps->param[0]),pcompar);
    
    /* remove doubles, keep the first one always. */
    j = 1;
    for(i=1; (i<ps->nr); i++) {
      if ((ps->param[i].AI != ps->param[j-1].AI) ||
	  (ps->param[i].AJ != ps->param[j-1].AJ) ) {
	cp_param(&(ps->param[j]),&(ps->param[i]));
	j++;
      } 
    }
    fprintf(stderr,"Number of bonds was %d, now %d\n",ps->nr,j);
    ps->nr=j;
  }
  else
    fprintf(stderr,"No bonds\n");
}

void print_sums(t_atoms *atoms, bool bSystem)
{
  double m,qtot;
  int    i;
  const char *where;
  
  if (bSystem)
    where=" in system";
  else
    where="";
  
  m=0;
  qtot=0;
  for(i=0; (i<atoms->nr); i++) {
    m+=atoms->atom[i].m;
    qtot+=atoms->atom[i].q;
  }
  
  fprintf(stderr,"Total mass%s %.3f a.m.u.\n",where,m);
  fprintf(stderr,"Total charge%s %.3f e\n",where,qtot);
}

static void check_restp_type(const char *name,int t1,int t2)
{
    if (t1 != t2)
    {
        gmx_fatal(FARGS,"Residues in one molecule have a different '%s' type: %d and %d",name,t1,t2);
    }
}

static void check_restp_types(t_restp *r0,t_restp *r1)
{
    int i;

    check_restp_type("all dihedrals",r0->bAlldih,r1->bAlldih);
    check_restp_type("nrexcl",r0->nrexcl,r1->nrexcl);
    check_restp_type("HH14",r0->HH14,r1->HH14);
    check_restp_type("remove dihedrals",r0->bRemoveDih,r1->bRemoveDih);

    for(i=0; i<ebtsNR; i++)
    {
        check_restp_type(btsNames[i],r0->rb[i].type,r1->rb[i].type);
    }
}

void add_atom_to_restp(t_restp *restp,int resnr,int at_start,const t_hack *hack)
{
    char buf[STRLEN];
    int  k;
    const char *Hnum="123456";

    /*if (debug)
    {
        fprintf(debug,"adding atom(s) %s to atom %s in res %d%s in rtp\n",
                hack->nname,
                *restp->atomname[at_start], resnr, restp->resname);
                }*/
    strcpy(buf, hack->nname);
    buf[strlen(buf)+1]='\0';
    if ( hack->nr > 1 )
    {
        buf[strlen(buf)]='-';
    }
    /* make space */
    restp->natom += hack->nr;
    srenew(restp->atom,     restp->natom);
    srenew(restp->atomname, restp->natom);
    srenew(restp->cgnr,     restp->natom);
    /* shift the rest */
    for(k=restp->natom-1; k > at_start+hack->nr; k--)
    {
        restp->atom[k] =
            restp->atom    [k - hack->nr];
        restp->atomname[k] =
            restp->atomname[k - hack->nr];
        restp->cgnr[k] =
            restp->cgnr    [k - hack->nr];
    }
    /* now add them */
    for(k=0; k < hack->nr; k++)
    {
        /* set counter in atomname */
        if ( hack->nr > 1 )
        {
            buf[strlen(buf)-1] = Hnum[k];
        }
        snew( restp->atomname[at_start+1+k], 1);
        restp->atom     [at_start+1+k] = *hack->atom;
        *restp->atomname[at_start+1+k] = strdup(buf);
        if ( hack->cgnr != NOTSET )
        {
            restp->cgnr   [at_start+1+k] = hack->cgnr;
        }
        else
        {
            restp->cgnr   [at_start+1+k] = restp->cgnr[at_start];
        }
    }
}

void get_hackblocks_rtp(t_hackblock **hb, t_restp **restp, 
			       int nrtp, t_restp rtp[],
			       int nres, t_resinfo *resinfo, 
			       int nterpairs,
			       t_hackblock **ntdb, t_hackblock **ctdb,
			       int *rn, int *rc)
{
  int i, j, k, l;
  t_restp *res;
  char buf[STRLEN];
  const char *Hnum="123456";
  int tern,terc;
  bool bN,bC,bRM;

  snew(*hb,nres);
  snew(*restp,nres);
  /* first the termini */
  for(i=0; i<nterpairs; i++) {
      if (rn[i] >= 0 && ntdb[i] != NULL) {
          copy_t_hackblock(ntdb[i], &(*hb)[rn[i]]);
      }
      if (rc[i] >= 0 && ctdb[i] != NULL) {
          merge_t_hackblock(ctdb[i], &(*hb)[rc[i]]);
      }
  }  

  /* then the whole rtp */
  for(i=0; i < nres; i++) {
    res = search_rtp(*resinfo[i].name,nrtp,rtp);
    copy_t_restp(res, &(*restp)[i]);

    /* Check that we do not have different bonded types in one molecule */
    check_restp_types(&(*restp)[0],&(*restp)[i]);

    tern = -1;
    for(j=0; j<nterpairs && tern==-1; j++) {
        if (i == rn[j]) {
            tern = j;
        }
    }
    terc = -1;
    for(j=0; j<nterpairs && terc == -1; j++) {
        if (i == rc[j]) {
            terc = j;
        }
    }
    bRM = merge_t_bondeds(res->rb, (*hb)[i].rb,tern>=0,terc>=0);

    if (bRM && ((tern >= 0 && ntdb[tern] == NULL) ||
                (terc >= 0 && ctdb[terc] == NULL))) {
        gmx_fatal(FARGS,"There is a dangling bond at at least one of the terminal ends and the force field does not provide terminal entries or files. Edit a .n.tdb and/or .c.tdb file.");
    }
    if (bRM && ((tern >= 0 && ntdb[tern]->nhack == 0) ||
                (terc >= 0 && ctdb[terc]->nhack == 0))) {
        gmx_fatal(FARGS,"There is a dangling bond at at least one of the terminal ends. Select a proper terminal entry.");
    }
  }
  
  /* now perform t_hack's on t_restp's,
     i.e. add's and deletes from termini database will be 
     added to/removed from residue topology 
     we'll do this on one big dirty loop, so it won't make easy reading! */
    for(i=0; i < nres; i++)
    {
        for(j=0; j < (*hb)[i].nhack; j++)
        {
            if ( (*hb)[i].hack[j].nr )
            {
                /* find atom in restp */
                for(l=0; l < (*restp)[i].natom; l++)
                    if ( ( (*hb)[i].hack[j].oname==NULL && 
                           strcmp((*hb)[i].hack[j].AI, *(*restp)[i].atomname[l])==0 ) ||
                         ( (*hb)[i].hack[j].oname!=NULL &&
                           strcmp((*hb)[i].hack[j].oname,*(*restp)[i].atomname[l])==0 ) )
                        break;
                if (l == (*restp)[i].natom)
                {
                    /* If we are doing an atom rename only, we don't need
                     * to generate a fatal error if the old name is not found
                     * in the rtp.
                     */
                    /* Deleting can happen also only on the input atoms,
                     * not necessarily always on the rtp entry.
                     */
                    if (!((*hb)[i].hack[j].oname != NULL &&
                          (*hb)[i].hack[j].nname != NULL) &&
                        !((*hb)[i].hack[j].oname != NULL &&
                          (*hb)[i].hack[j].nname == NULL))
                    {
                        gmx_fatal(FARGS,"atom %s not found in residue %d%s "
                                  "while combining tdb and rtp",
                                  (*hb)[i].hack[j].oname!=NULL ? 
                                  (*hb)[i].hack[j].oname : (*hb)[i].hack[j].AI, 
                                  i+1,*resinfo[i].name);
                    }
                }
                else
                {
                    if ( (*hb)[i].hack[j].oname==NULL ) {
                        /* we're adding: */
                        add_atom_to_restp(&(*restp)[i],resinfo[i].nr,l,
                                          &(*hb)[i].hack[j]);
                    }
                    else
                    {
                        /* oname != NULL */
                        if ( (*hb)[i].hack[j].nname==NULL ) {
                            /* we're deleting */
                            if (debug) 
                                fprintf(debug, "deleting atom %s from res %d%s in rtp\n",
                                        *(*restp)[i].atomname[l], 
                                        i+1,(*restp)[i].resname);
                            /* shift the rest */
                            (*restp)[i].natom--;
                            for(k=l; k < (*restp)[i].natom; k++) {
                                (*restp)[i].atom    [k] = (*restp)[i].atom    [k+1];
                                (*restp)[i].atomname[k] = (*restp)[i].atomname[k+1];
                                (*restp)[i].cgnr    [k] = (*restp)[i].cgnr    [k+1];
                            }
                            /* give back space */
                            srenew((*restp)[i].atom,     (*restp)[i].natom);
                            srenew((*restp)[i].atomname, (*restp)[i].natom);
                            srenew((*restp)[i].cgnr,     (*restp)[i].natom);
                        } else { /* nname != NULL */
                            /* we're replacing */
                            if (debug) 
                                fprintf(debug, "replacing atom %s by %s in res %d%s in rtp\n",
                                        *(*restp)[i].atomname[l], (*hb)[i].hack[j].nname, 
                                        i+1,(*restp)[i].resname);
                            snew( (*restp)[i].atomname[l], 1);
                            (*restp)[i].atom[l]      =       *(*hb)[i].hack[j].atom;
                            *(*restp)[i].atomname[l] = strdup((*hb)[i].hack[j].nname);
                            if ( (*hb)[i].hack[j].cgnr != NOTSET )
                                (*restp)[i].cgnr   [l] = (*hb)[i].hack[j].cgnr;
                        }
                    }
                }
            }
        }
    }
}

static bool atomname_cmp_nr(const char *anm,t_hack *hack,int *nr)
{

  if (hack->nr == 1) {
    *nr = 0;
    
    return (strcasecmp(anm,hack->nname) != 0);
  } else {
    if (isdigit(anm[strlen(anm)-1])) {
      *nr = anm[strlen(anm)-1] - '0';
    } else {
      *nr = 0;
    }
    if (*nr <= 0 || *nr > hack->nr) {
      return FALSE;
    } else {
      return (strlen(anm) == strlen(hack->nname) + 1 &&
	      strncasecmp(anm,hack->nname,strlen(hack->nname)) == 0);
    }
  }
}

static bool match_atomnames_with_rtp_atom(t_atoms *pdba,int atind,
                                          t_restp *rptr,t_hackblock *hbr,
                                          bool bVerbose)
{
    int  resnr;
    int  i,j,k;
    char *oldnm,*resnm,*newnm;
    int  anmnr;
    char *start_at,buf[STRLEN];
    int  start_nr;
    bool bFoundInAdd;
    bool bDeleted;

    oldnm = *pdba->atomname[atind];
    resnr = pdba->resinfo[pdba->atom[atind].resind].nr;

    bDeleted = FALSE;
    for(j=0; j<hbr->nhack; j++)
    {
        if (hbr->hack[j].oname != NULL && hbr->hack[j].nname != NULL &&
            strcasecmp(oldnm,hbr->hack[j].oname) == 0)
        {
            /* This atom still has the old name, rename it */
            /* We need to find the add hack that can add this atom
             * to find out after which atom it should be added.
             */
            newnm = hbr->hack[j].nname;
            bFoundInAdd = FALSE;
            for(k=0; k<hbr->nhack; k++) {
                if (hbr->hack[k].oname == NULL && hbr->hack[k].nname != NULL) {
                    if (atomname_cmp_nr(newnm,&hbr->hack[k],&anmnr))
                    {
                        if (anmnr <= 1)
                        {
                            start_at = hbr->hack[k].a[0];
                        }
                        else
                        {
                            sprintf(buf,"%s%d",hbr->hack[k].nname,anmnr-1);
                            start_at = buf;
                        }
                        for(start_nr=0; start_nr<rptr->natom; start_nr++)
                        {
                            if (strcasecmp(start_at,(*rptr->atomname[start_nr])) == 0)
                            {
                                break;
                            }
                        }
                        if (start_nr == rptr->natom)
                        {
                            gmx_fatal(FARGS,"Could not find atom '%s' in residue building block '%s' to add atom '%s' to",
                                      start_at,rptr->resname,newnm);
                        }
                        /* We can add the atom after atom start_nr */
                        add_atom_to_restp(rptr,resnr,start_nr,
                                          &hbr->hack[j]);
                    }
                    else
                    {
                        gmx_fatal(FARGS,"Could not find an 'add' entry for atom named '%s' corresponding to the 'replace' entry from atom name '%s' to '%s' for tdb or hdb database of residue type '%s'",
                                  newnm,
                                  hbr->hack[j].oname,hbr->hack[j].nname,
                                  rptr->resname);
                    }
                }
            }
            if (!bFoundInAdd)
            {
                for(k=0; k<rptr->natom; k++)
                {
                    if (strcasecmp(newnm,*rptr->atomname[k]) == 0)
                    {
                        break;
                    }
                }
                if (k == rptr->natom) {
                    gmx_fatal(FARGS,"Renamed atom '%s' (old name '%s') not found in rtp entry '%s' for residue %d",
                              oldnm,newnm,rptr->resname,resnr);
                }
            }
            
            if (bVerbose)
            {
                printf("Renaming atom '%s' in residue '%s' %d to '%s'\n",
                       oldnm,rptr->resname,resnr,newnm);
            }
            /* Rename the atom in pdba */
            snew(pdba->atomname[atind],1);
            *pdba->atomname[atind] = strdup(newnm);
        }
        else if (hbr->hack[j].oname != NULL && hbr->hack[j].nname == NULL &&
                 strcasecmp(oldnm,hbr->hack[j].oname) == 0)
        {
            /* This is a delete entry, check if this atom is present
             * in the rtp entry of this residue.
             */
            for(k=0; k<rptr->natom; k++)
            {
                if (strcasecmp(oldnm,*rptr->atomname[k]) == 0)
                {
                    break;
                }
            }
            if (k == rptr->natom)
            {
                /* This atom is not present in the rtp entry,
                 * delete is now from the input pdba.
                 */
                if (bVerbose)
                {
                    printf("Deleting atom '%s' in residue '%s'\n",
                           oldnm,rptr->resname,resnr);
                }
                sfree(pdba->atomname[atind]);
                for(k=atind+1; k<pdba->nr; k++)
                {
                    pdba->atom[k-1]     = pdba->atom[k];
                    pdba->atomname[k-1] = pdba->atomname[k];
                }
                pdba->nr--;
                bDeleted = TRUE;
            }
        }
    }

    return bDeleted;
}
    
void match_atomnames_with_rtp(t_restp restp[],t_hackblock hb[],
                              t_atoms *pdba,
                              bool bVerbose)
{
    int  i,j,k;
    char *oldnm,*resnm,*newnm;
    int  resnr;
    t_restp *rptr;
    t_hackblock *hbr;
    int  anmnr;
    char *start_at,buf[STRLEN];
    int  start_nr;
    bool bFoundInAdd;
    
    for(i=0; i<pdba->nr; i++)
    {
        oldnm = *pdba->atomname[i];
        resnm = *pdba->resinfo[pdba->atom[i].resind].name;
        resnr = pdba->resinfo[pdba->atom[i].resind].nr;
        rptr  = &restp[pdba->atom[i].resind];
        for(j=0; (j<rptr->natom); j++)
        {
            if (strcasecmp(oldnm,*(rptr->atomname[j])) == 0)
            {
                break;
            }
        }
        if (j == rptr->natom)
        {
            /* Not found yet, check if we have to rename this atom */
            if (match_atomnames_with_rtp_atom(pdba,i,
                                              rptr,&(hb[pdba->atom[i].resind]),
                                              bVerbose))
            {
                /* We deleted this atom, decrease the atom counter by 1. */
                i--;
            }
        }
    }
}

void gen_cmap(t_params *psb, t_restp *restp, int natoms, t_atom atom[], char **aname[], int nres)
{
	int     residx,i,ii,j,k;
	atom_id ai,aj,ak,al,am;
	const char *ptr;
	
	if (debug)
		ptr = "cmap";
	else
		ptr = "check";
	
	fprintf(stderr,"Making cmap torsions...");
	i=0;
	/* End loop at nres-1, since the very last residue does not have a +N atom, and
	 * therefore we get a valgrind invalid 4 byte read error with atom am */
	for(residx=0; residx<nres-1; residx++)
	{
		/* Add CMAP terms from the list of CMAP interactions */
		for(j=0;j<restp[residx].rb[ebtsCMAP].nb; j++)
		{
			ai=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[0],i,natoms,atom,aname,
						   ptr,TRUE);
			aj=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[1],i,natoms,atom,aname,
						   ptr,TRUE);
			ak=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[2],i,natoms,atom,aname,
						   ptr,TRUE);
			al=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[3],i,natoms,atom,aname,
						   ptr,TRUE);
			am=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[4],i,natoms,atom,aname,
						   ptr,TRUE);
			
			/* The first and last residues no not have cmap torsions */
			if(ai!=NO_ATID && aj!=NO_ATID && ak!=NO_ATID && al!=NO_ATID && am!=NO_ATID)
			{
				add_cmap_param(psb,ai,aj,ak,al,am,restp[residx].rb[ebtsCMAP].b[j].s);
			}
		}
		
		if(residx<nres-1)
		{
			while(atom[i].resind<residx+1)
			{
				i++;
			}
		}
	}
	
	/* Start the next residue */
}

static void 
scrub_charge_groups(int *cgnr, int natoms)
{
	int i;
	
	for(i=0;i<natoms;i++)
	{
		cgnr[i]=i+1;
	}
}


void pdb2top(FILE *top_file, char *posre_fn, char *molname,
             t_atoms *atoms, rvec **x, gpp_atomtype_t atype, t_symtab *tab,
             int nrtp, t_restp rtp[],
             t_restp *restp, t_hackblock *hb,
             int nterpairs,t_hackblock **ntdb, t_hackblock **ctdb,
             int *rn, int *rc, bool bAllowMissing,
             bool bVsites, bool bVsiteAromatics,
             const char *ff, const char *ffdir, real mHmult,
             int nssbonds, t_ssbond *ssbonds,
             real long_bond_dist, real short_bond_dist,
             bool bDeuterate, bool bChargeGroups, bool bCmap,
             bool bRenumRes)
{
    /*
  t_hackblock *hb;
  t_restp  *restp;
    */
  t_params plist[F_NRE];
  t_excls  *excls;
  t_nextnb nnb;
  int      *cgnr;
  int      *vsite_type;
  int      i,nmissat;
  int      bts[ebtsNR];
  
  init_plist(plist);

  /* lookup hackblocks and rtp for all residues */
  //get_hackblocks_rtp(&hb, &restp, nrtp, rtp, atoms->nres, atoms->resinfo, 
  //	     nterpairs, ntdb, ctdb, rn, rc);
  /* ideally, now we would not need the rtp itself anymore, but do 
     everything using the hb and restp arrays. Unfortunately, that 
     requires some re-thinking of code in gen_vsite.c, which I won't 
     do now :( AF 26-7-99 */
  
  if (debug) {
    print_resall(debug, atoms->nres, restp, atype);
    dump_hb(debug, atoms->nres, hb);
  }
  
  /* Make bonds */
  at2bonds(&(plist[F_BONDS]), hb, 
           atoms->nr, atoms->atom, atoms->atomname, atoms->nres, *x, 
           long_bond_dist, short_bond_dist, bAllowMissing);
  
  /* specbonds: disulphide bonds & heme-his */
  do_ssbonds(&(plist[F_BONDS]),
	     atoms->nr, atoms->atom, atoms->atomname, nssbonds, ssbonds,
	     bAllowMissing);
  
  nmissat = name2type(atoms, &cgnr, atype, restp);
  if (nmissat) {
    if (bAllowMissing)
      fprintf(stderr,"There were %d missing atoms in molecule %s\n",
	      nmissat,molname);
    else
      gmx_fatal(FARGS,"There were %d missing atoms in molecule %s, if you want to use this incomplete topology anyhow, use the option -missing",
		  nmissat,molname);
  }
  
  /* Cleanup bonds (sort and rm doubles) */ 
  clean_bonds(&(plist[F_BONDS]));
  
  snew(vsite_type,atoms->nr);
  for(i=0; i<atoms->nr; i++)
    vsite_type[i]=NOTSET;
  if (bVsites) {
    /* determine which atoms will be vsites and add dummy masses 
       also renumber atom numbers in plist[0..F_NRE]! */
    do_vsites(nrtp, rtp, atype, atoms, tab, x, plist, 
              &vsite_type, &cgnr, mHmult, bVsiteAromatics, ffdir);
  }
  
  /* Make Angles and Dihedrals */
  fprintf(stderr,"Generating angles, dihedrals and pairs...\n");
  snew(excls,atoms->nr);
  init_nnb(&nnb,atoms->nr,4);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,restp[0].nrexcl,restp[0].HH14,
          plist,excls,hb,restp[0].bAlldih,restp[0].bRemoveDih,
          bAllowMissing);
  done_nnb(&nnb);
  
    /* Make CMAP */
    if (TRUE == bCmap)
    {
		gen_cmap(&(plist[F_CMAP]), restp, atoms->nr, atoms->atom, atoms->atomname, atoms->nres);
        if (plist[F_CMAP].nr > 0)
        {
            fprintf(stderr, "There are %4d cmap torsion pairs\n",
                    plist[F_CMAP].nr);
        }
    }

  /* set mass of all remaining hydrogen atoms */
  if (mHmult != 1.0)
    do_h_mass(&(plist[F_BONDS]),vsite_type,atoms,mHmult,bDeuterate);
  sfree(vsite_type);
  
  /* Cleanup bonds (sort and rm doubles) */ 
  /* clean_bonds(&(plist[F_BONDS]));*/
   
  fprintf(stderr,
	  "There are %4d dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d virtual sites\n",
	  plist[F_PDIHS].nr, plist[F_IDIHS].nr, plist[F_ANGLES].nr,
	  plist[F_LJ14].nr, plist[F_BONDS].nr,
	  plist[F_VSITE2].nr +
	  plist[F_VSITE3].nr +
	  plist[F_VSITE3FD].nr +
	  plist[F_VSITE3FAD].nr +
	  plist[F_VSITE3OUT].nr +
      plist[F_VSITE4FD].nr +
      plist[F_VSITE4FDN].nr );
  
  print_sums(atoms, FALSE);
  
  if (FALSE == bChargeGroups)
  {
	  scrub_charge_groups(cgnr, atoms->nr);
  }

    if (bRenumRes)
    {
        for(i=0; i<atoms->nres; i++) 
        {
            atoms->resinfo[i].nr = i + 1;
            atoms->resinfo[i].ic = ' ';
        }
    }
	
  if (top_file) {
    fprintf(stderr,"Writing topology\n");
    /* We can copy the bonded types from the first restp,
     * since the types have to be identical for all residues in one molecule.
     */
    for(i=0; i<ebtsNR; i++) {
        bts[i] = restp[0].rb[i].type;
    }
    write_top(top_file, posre_fn, molname,
	      atoms, bts, plist, excls, atype, cgnr, restp[0].nrexcl);
  }
  
  /* cleaning up */
  free_t_hackblock(atoms->nres, &hb);
  free_t_restp(atoms->nres, &restp);
    
  /* we should clean up hb and restp here, but that is a *L*O*T* of work! */
  sfree(cgnr);
  for (i=0; i<F_NRE; i++)
    sfree(plist[i].param);
  for (i=0; i<atoms->nr; i++)
    sfree(excls[i].e);
  sfree(excls);
}
