from pathlib import Path
from copy import deepcopy
from zipfile import ZipFile
import os, csv, json, hashlib, shutil, re
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[key]='1'
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from docx import Document
from docx.shared import Inches, Pt
from lxml import etree
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'phase4c'
for name in ['figures', 'supplementary']:
 (OUT/name).mkdir(parents=True, exist_ok=True)
ORDER=['BRCA','CCRCC','COAD','GBM','HNSCC','LSCC','LUAD','OV','PDAC','UCEC']
COL={'A':'#2878B5','B':'#C45A5A'}
PAL=dict(zip(ORDER,['#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#56B4E9','#332288','#882255','#44AA99','#999933']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
records=[]; legends={}; tables=[]
def read(p): return pd.read_csv(ROOT/p)
def record(kind,id,row,field,path,value,display=None,t='BOTH'):
 records.append(dict(kind=kind,id=id,panel_or_row=row,variable=field,source_file=path,source_field=field,transformation=t,exact_value=str(value),displayed_value=str(value if display is None else display),status='PASS'))
def register_frame(kind,id,df,path):
 for i,row in df.iterrows():
  for col,value in row.items():
   if pd.notna(value): record(kind,id,str(i),col,path,value,t=str(row.get('transformation','BOTH')))
def save(fig,name,supp=False):
 folder=OUT/('supplementary' if supp else 'figures')
 for ext,dpi in [('png',180),('pdf',180),('tiff',600)]:
  fig.savefig(folder/f'{name}.{ext}',dpi=dpi,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
 plt.close(fig)
def panel(ax,label): ax.set_title(label,loc='left',fontsize=10,fontweight='bold',pad=10)
def threshold(ax,value=.02): ax.axvline(value,color='#555555',ls='--',lw=.9,zorder=0)
def intervals(ax,rows,fields=('mu','ci_low','ci_high'),offset=0,color=None):
 for i,(_,r) in enumerate(rows.iterrows()):
  m,l,h=(float(r[f]) for f in fields)
  ax.errorbar(m,i+offset,xerr=[[m-l],[h-m]],fmt='o',ms=4,color=color or COL[r.transformation],capsize=2,lw=1)
cs=read('phase1a/results/cancer_level_summary.csv'); meta=read('phase1a/results/random_effects_summary.csv'); loo=read('phase1a/results/leave_one_cancer_out.csv'); common=read('phase1a/results/common_support_summary.csv')
pooled={t:meta[(meta.transformation==t)&(meta.analysis=='all_ten')].iloc[0] for t in 'AB'}
# Population and exact source-backed coverage.
pair=read('phase0_5/manifests/reconstructed_pairing.csv').query('in_phase0_triple == 1')
cov=read('phase1b/results/protein_coverage_comparison.csv').set_index('cancer')
counts=[len(pair[pair.cancer==c]) for c in ORDER]
assert sum(counts)==1027
fig,axes=plt.subplots(2,2,figsize=(7.2,5.7)); fig.subplots_adjust(hspace=.45,wspace=.22)
texts=[('A  Study population','10 CPTAC cancer cohorts\n1,027 case-level matched tumors\n\n'+'\n'.join('   '.join(f'{ORDER[i]}  {counts[i]}' for i in range(j,min(j+2,10))) for j in range(0,10,2))),('B  Measurement layers','Gene-level RNA\nPublished RNA covariate sensitivity\n\nTranscript-level RNA\nReconstructed total G and principal share P\n\nGene-level protein abundance\nPublished protein outcome'),('C  Principal-transcript definition','GENCODE v34\n↓\nUnique MANE Select\nUnique APPRIS principal-1 fallback\nfor genes without MANE annotations\n↓\n10,795 candidate genes'),('D  Statistical analysis','M0: protein ~ log₂(G) + WES purity + CNV\nM1: M0 + principal-transcript proportion P\n↓\nGene-level partial R², within each cancer\nDelete-one-case jackknife confidence interval\n↓\nREML synthesis across cancer means\nHartung–Knapp CI and prediction interval')]
for ax,(title,txt) in zip(axes.flat,texts):
 ax.axis('off'); panel(ax,title); ax.text(0,.96,txt,va='top',linespacing=1.55,transform=ax.transAxes,fontsize=8)
fig.text(.02,.015,'Case-level matching; aliquot-level correspondence was not established.',fontsize=8)
save(fig,'Figure1_study_design')
for c,n in zip(ORDER,counts): record('figure','Figure1','A',c,'phase0_5/manifests/reconstructed_pairing.csv',n)
for f,v,p in [('total_cases',1027,'phase1a/reports/cohort_reconstruction_report.md'),('candidate_genes',10795,'phase1a/results/cancer_level_summary.csv'),('GENCODE_version',34,'phase0_6/reports/statistical_model_freeze.md')]: record('figure','Figure1','A/C',f,p,v)
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,3.7),gridspec_kw={'width_ratios':[1.4,1]});fig.subplots_adjust(wspace=.38,bottom=.22)
a.axis('off');panel(a,'A  Nested models');a.text(0,.95,'Outcome: gene-level protein abundance\n\nM0: protein ~ log₂(G) + WES purity + CNV\nM1: M0 + principal-transcript proportion P\n\nG = sum of reconstructed linear transcript values\nP = principal-transcript value / G\n\npartial R² = (RSS_M0 − RSS_M1) / RSS_M0',va='top',linespacing=1.55,fontsize=8)
panel(b,'B  Pooled estimates')
for i,t in enumerate('AB'):
 r=pooled[t]; b.errorbar(r.mu,i,xerr=[[r.mu-r.ci_low],[r.ci_high-r.mu]],fmt='o',color=COL[t],capsize=4); b.text(.034,i-.30,f'{r.mu:.4f}\n[{r.ci_low:.4f}, {r.ci_high:.4f}]',ha='right',va='center',fontsize=7)
threshold(b);b.set(yticks=[0,1],yticklabels=['Literal A','Literal B'],xlim=(.015,.035),ylim=(1.6,-.6),xlabel='Pooled mean partial R²');b.set_xticks([.02,.025,.03])
fig.text(.02,.03,'Two co-primary measurement representations. Error bars show 95% Hartung–Knapp CIs.',fontsize=8)
save(fig,'Figure2_incremental_explanation')
register_frame('figure','Figure2',meta.query("analysis == 'all_ten'")[['transformation','mu','ci_low','ci_high']],'phase1a/results/random_effects_summary.csv')
fits={t:pd.read_csv(ROOT/f'phase1a/results/primary_literal_{t}_gene_fits.tsv.gz',sep='\t') for t in 'AB'}
fig,axes=plt.subplots(2,2,figsize=(7.2,6),gridspec_kw={'width_ratios':[1.2,1]});fig.subplots_adjust(wspace=.42,hspace=.38,bottom=.18)
for j,t in enumerate('AB'):
 for k,limit in [(0,.10)]:
  ax=axes[j,k]
  for c in ORDER:
   d=fits[t];v=np.sort(d.loc[(d.cancer==c)&(d.fit_status=='OK'),'partial_r2'].to_numpy());ax.plot(v,np.arange(1,len(v)+1)/len(v),color=PAL[c],lw=1,label=c)
  ax.set(xlim=(0,limit),ylim=(0,1),xlabel='Gene-level partial R²',ylabel=f'Literal {t} · cumulative fraction');threshold(ax)
  panel(ax,('A  Main range' if k==0 else 'B  Full range') if j==0 else '')
 ax=axes[j,1];d=cs[cs.transformation==t].set_index('cancer').loc[ORDER]
 ax.plot(d.mean_partial_r2,range(10),'o',color=COL[t],label='Mean',ms=4);ax.plot(d.median_partial_r2,range(10),'s',color=COL[t],mfc='white',label='Median',ms=4)
 ax.set(yticks=range(10),yticklabels=ORDER,xlim=(0,.038),ylim=(9.7,-.7),xlabel='Partial R²');threshold(ax);panel(ax,'B  Mean and median' if j==0 else '');
fig.legend([Line2D([],[],color=PAL[c]) for c in ORDER],ORDER,ncol=5,loc='lower center',frameon=False,fontsize=8)
save(fig,'Figure3_partial_R2_distributions')
# Full range supplement retains every successful fit without selecting by effect size.
fig,axes=plt.subplots(1,2,figsize=(7.2,3.5))
for ax,t in zip(axes,'AB'):
 for c in ORDER:
  d=fits[t];v=np.sort(d.loc[(d.cancer==c)&(d.fit_status=='OK'),'partial_r2']);ax.plot(v,np.arange(1,len(v)+1)/len(v),color=PAL[c],lw=1)
 ax.set(xlim=(0,1),ylim=(0,1),xlabel='Gene-level partial R²',ylabel='Cumulative fraction');panel(ax,'Literal '+t)
fig.legend([Line2D([],[],color=PAL[c]) for c in ORDER],ORDER,ncol=5,loc='lower center',frameon=False,fontsize=7);fig.subplots_adjust(bottom=.25,wspace=.3)
save(fig,'FigureS2_full_range_distributions',True)
for t,d in fits.items():
 ok=d[d.fit_status=='OK']; ok[['cancer','partial_r2']].to_csv(OUT/f'manifests/ECDF_points_literal_{t}.csv.gz',index=False)
 for c in ORDER:
  n=int((ok.cancer==c).sum()); assert n==int(cs.query('cancer == @c and transformation == @t').iloc[0].n_success)
  record('figure','Figure3/FigureS2',c,'successful_gene_fits',f'phase1a/results/primary_literal_{t}_gene_fits.tsv.gz',n,t=t)
register_frame('figure','Figure3',cs[['cancer','transformation','mean_partial_r2','median_partial_r2']],'phase1a/results/cancer_level_summary.csv')
fig,axes=plt.subplots(1,2,figsize=(8,5.7));fig.subplots_adjust(wspace=.55,bottom=.20)
for ax,t,label in zip(axes,'AB',['A','B']):
 d=cs[cs.transformation==t].set_index('cancer').loc[ORDER].reset_index();intervals(ax,d,('mean_partial_r2','ci_low','ci_high'),color=COL[t])
 r=pooled[t]; ax.errorbar(r.mu,11,xerr=[[r.mu-r.ci_low],[r.ci_high-r.mu]],fmt='D',ms=6,lw=3,color=COL[t],capsize=3)
 ax.hlines(12,r.pi_low,r.pi_high,color='#888888',lw=1.3);ax.plot(r.mu,12,'|',color=COL[t],ms=8);threshold(ax)
 ax.set(yticks=list(range(10))+[11,12],yticklabels=[f'{c}  (n={int(n)})' for c,n in zip(ORDER,d.n_success)]+['REML pooled','95% PI'],ylim=(13,-.8),xlim=(.012,.036),xlabel='Mean partial R²');panel(ax,f'{label}  Literal {t}')
 ax.text(0,-.19,f'Pooled {r.mu:.4f} · 95% CI {r.ci_low:.4f}–{r.ci_high:.4f}\n95% PI {r.pi_low:.4f}–{r.pi_high:.4f} · I² {r.I2:.1%}',transform=ax.transAxes,fontsize=8,linespacing=1.5)
fig.text(.04,.015,'n = successful gene fits. Dashed line = 0.02. Thick bar = pooled CI. Thin grey line = PI.',fontsize=8)
save(fig,'Figure4_cancer_heterogeneity')
register_frame('figure','Figure4',cs[['cancer','transformation','n_success','mean_partial_r2','ci_low','ci_high']],'phase1a/results/cancer_level_summary.csv');register_frame('figure','Figure4',meta.query("analysis == 'all_ten'"),'phase1a/results/random_effects_summary.csv')
fig=plt.figure(figsize=(7.2,7));gs=fig.add_gridspec(2,2,height_ratios=[1,1.4]);axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])];fig.subplots_adjust(wspace=.65,hspace=.65,bottom=.13)
ax=axes[0]
for r in common.itertuples():
 ax.scatter(r.mean_A,r.mean_B,color=PAL[r.cancer],s=25)
ax.plot([.017,.032],[.017,.032],ls='--',color='#555555');ax.set(xlim=(.017,.032),ylim=(.017,.032),xlabel='Literal A mean partial R²',ylabel='Literal B mean partial R²');ax.set_aspect('equal');panel(ax,'A  Paired support')
for ax,kind in [(axes[1],'summary'),(axes[2],'loo')]:
 for t,offset in [('A',-.14),('B',.14)]:
  rows=meta[(meta.transformation==t)].set_index('analysis').loc[['all_ten','paired_support','exclude_OV']].reset_index() if kind=='summary' else loo[loo.transformation==t].set_index('left_out').loc[ORDER].reset_index()
  intervals(ax,rows,offset=offset,color=COL[t])
 labels=['All ten','Paired support','Excluding OV'] if kind=='summary' else ['Omit '+c for c in ORDER]
 ax.set(yticks=range(len(labels)),yticklabels=labels,ylim=(len(labels)-.4,-.6),xlim=(.019,.030),xlabel='Pooled mean partial R²');threshold(ax);panel(ax,'B  Primary sensitivities' if kind=='summary' else 'C  Leave one cancer out')
fig.legend([Line2D([],[],marker='o',color=COL[t]) for t in 'AB'],['Literal A','Literal B'],loc='lower center',ncol=2,frameon=False)
save(fig,'Figure5_transformation_robustness')
for df,path in [(common,'common_support_summary.csv'),(meta,'random_effects_summary.csv'),(loo,'leave_one_cancer_out.csv')]:register_frame('figure','Figure5',df,'phase1a/results/'+path)
struct=[]
for t in 'AB':
 r=read(f'phase2b/results/primary_S2_vs_S1_literal_{t}.csv').query("quantity == 'adjusted_S2_minus_S1'").copy();struct.append(r);register_frame('figure','FigureS1',r,f'phase2b/results/primary_S2_vs_S1_literal_{t}.csv')
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,3.3),gridspec_kw={'width_ratios':[1.2,1]});fig.subplots_adjust(wspace=.35)
a.axis('off');panel(a,'A  Predicted categories');a.text(0,.92,'S2\nDifferent predicted protein sequence\nand supported Pfam, signal-peptide\nor transmembrane annotation difference\n\nS1\nDifferent predicted sequence without\na supported difference in those features',va='top',linespacing=1.5)
d=pd.concat(struct);intervals(b,d,('estimate','ci_low','ci_high'));threshold(b,.005);b.axvline(0,color='#BBBBBB',lw=.8);b.set(xlim=(-.0002,.007),yticks=[0,1],yticklabels=['Literal A','Literal B'],ylim=(1.5,-.5),xlabel='Adjusted S2 − S1 difference');panel(b,'B  Secondary association');save(fig,'FigureS1_structural_annotation',True)
sens=read('phase1a/results/sensitivity_summary.csv')
labels=list(dict.fromkeys(sens.sensitivity));fig,ax=plt.subplots(figsize=(9,4.5))
for j,label in enumerate(labels):
 for t,offset in [('A',-.15),('B',.15),('C',.30)]:
  rows=sens[(sens.sensitivity==label)&(sens.transformation==t)]
  vals=rows.mean_partial_r2.dropna()
  if len(vals):ax.scatter(vals,np.full(len(vals),j+offset),s=13,color=COL.get(t,'#777777'),marker='o' if t!='C' else '^')
ax.legend([Line2D([],[],marker='o',color=COL['A']),Line2D([],[],marker='o',color=COL['B']),Line2D([],[],marker='^',color='#777777')],['Literal A','Literal B','Candidate C'],loc='lower right',fontsize=7,frameon=False);ax.text(.001,labels.index('S4_clr'),'Not run: undefined contrast',fontsize=7,va='center');ax.set(yticks=range(len(labels)),yticklabels=labels,xlabel='Stored mean partial R²',ylim=(len(labels)-.4,-.6));threshold(ax);ax.tick_params(axis='y',labelsize=7);fig.subplots_adjust(left=.34,bottom=.18);save(fig,'FigureS3_expanded_sensitivities',True);register_frame('figure','FigureS3',sens,'phase1a/results/sensitivity_summary.csv')
