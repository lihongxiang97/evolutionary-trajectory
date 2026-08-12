#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data" / "main_figures" / "figure2" / "trajectory_classifications_three_methods.tsv"
OUT = ROOT / "output" / "figure2"
OUT.mkdir(parents=True, exist_ok=True)

CLASSES = ["Conservation", "Neofunctionalization(Child)",
           "Neofunctionalization(Parent)", "Specialization"]
SHORT = {"Conservation":"Con", "Neofunctionalization(Child)":"Neo(C)",
         "Neofunctionalization(Parent)":"Neo(P)", "Specialization":"Spe"}
MCOL = {"Conservation":"#4C78A8",
        "Neofunctionalization(Child)":"#F2A541",
        "Neofunctionalization(Parent)":"#59A14F",
        "Specialization":"#D65F5F"}
RCOL = {"P":"#4472C4", "C":"#E07A5F", "A":"#595959"}

mpl.rcParams.update({"font.family":"Arial", "font.size":8, "axes.linewidth":.7,
                     "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none"})

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]
    q=np.minimum.accumulate((r*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    z=np.empty_like(q); z[o]=np.minimum(q,1); return z

def star(p):
    return "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<.05 else "ns"

def clean(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

def panel_a(ax, d):
    tab=pd.crosstab(d.duplication_type,d.two_method_agreement).reindex(
        index=["TD","PD","TRD"],columns=CLASSES,fill_value=0)
    pct=tab.div(tab.sum(1),axis=0)*100
    bottom=np.zeros(3)
    for m in CLASSES:
        ax.bar(["TD","PD","TRD"],pct[m],bottom=bottom,width=.62,
               color=MCOL[m],edgecolor="white",linewidth=.7,label=SHORT[m])
        bottom += pct[m].to_numpy()
    ax.set_ylim(0,100); ax.set_ylabel("Percentage (%)")
    ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.16),
              handlelength=1.3,columnspacing=1.0)
    clean(ax)
    return tab,pct

def mini_profile(ax, x, y, values, color, width=.075, height=.14):
    n=len(values); bw=width/n
    for i,v in enumerate(values):
        ax.add_patch(Rectangle((x-width/2+i*bw,y-height/2),bw*.78,height*v,
                               facecolor=color,edgecolor="none"))

def panel_b(ax):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    names=["Con","Neo(C)","Neo(P)","Spe"]
    xs=np.linspace(.13,.87,4)
    anc=np.array([.75,.35,.85,.55,.65,.40])
    patterns={
        "Con":(anc*.94,anc*1.03),
        "Neo(C)":(anc,np.array([.18,.85,.20,.18,.78,.22])),
        "Neo(P)":(np.array([.18,.85,.20,.18,.78,.22]),anc),
        "Spe":(np.array([.72,.12,.18,.76,.15,.12]),np.array([.10,.66,.75,.12,.72,.65]))
    }
    for x,name in zip(xs,names):
        # Explicit ancestor-to-duplicates genealogy.
        ax.plot([x,x],[.72,.60],color="#666",lw=1)
        ax.plot([x-.045,x+.045],[.60,.60],color="#666",lw=1)
        ax.plot([x-.045,x-.045],[.60,.43],color="#666",lw=1)
        ax.plot([x+.045,x+.045],[.60,.43],color="#666",lw=1)
        ax.add_patch(Circle((x,.76),.018,facecolor=RCOL["A"],edgecolor="none"))
        ax.add_patch(Circle((x-.045,.40),.018,facecolor=RCOL["P"],edgecolor="none"))
        ax.add_patch(Circle((x+.045,.40),.018,facecolor=RCOL["C"],edgecolor="none"))
        mini_profile(ax,x,.86,anc,RCOL["A"])
        mini_profile(ax,x-.045,.27,patterns[name][0],RCOL["P"])
        mini_profile(ax,x+.045,.27,patterns[name][1],RCOL["C"])
        ax.text(x,.08,name,ha="center",va="center",fontsize=8)
        if x==xs[0]:
            ax.text(x-.078,.40,"P",ha="center",va="center",color=RCOL["P"],weight="bold")
            ax.text(x+.078,.40,"C",ha="center",va="center",color=RCOL["C"],weight="bold")
            ax.text(x-.033,.76,"A",ha="right",va="center",color=RCOL["A"],weight="bold")
    ax.text(.5,.98,"Expression profile",ha="center",va="top",fontsize=8,color="#555")

def bracket(ax,x1,x2,y,h,text):
    ax.plot([x1,x1,x2,x2],[y,y+h,y+h,y],c="#333",lw=.55,clip_on=False)
    ax.text((x1+x2)/2,y+h+.015,text,ha="center",va="bottom",fontsize=5.2)

def panel_c(fig,gs,d):
    tests=[]; plotrows=[]; axes={}
    roles=["P","C","A"]; positions=np.arange(3)
    for ri,mode in enumerate(["TD","PD","TRD"]):
        for ci,cl in enumerate(CLASSES):
            ax=fig.add_subplot(gs[ri,ci]); axes[(ri,ci)]=ax
            sub=d[(d.duplication_type==mode)&(d.two_method_agreement==cl)]
            vals={"P":np.log2(pd.to_numeric(sub.parent_mean_TPM)+1).to_numpy(),
                  "C":np.log2(pd.to_numeric(sub.child_mean_TPM)+1).to_numpy(),
                  "A":np.log2(pd.to_numeric(sub.ancestor_mean_TPM)+1).to_numpy()}
            bp=ax.boxplot([vals[r] for r in roles],positions=positions,widths=.54,
                          patch_artist=True,showfliers=False,
                          medianprops={"color":"white","lw":.8},
                          whiskerprops={"lw":.55},capprops={"lw":.55},boxprops={"lw":.55})
            for patch,r in zip(bp["boxes"],roles):
                patch.set_facecolor(RCOL[r]); patch.set_alpha(.92)
            for r in roles:
                plotrows += [{"duplication_type":mode,"trajectory":SHORT[cl],"role":r,
                              "log2_TPM_plus_1":v} for v in vals[r]]
            for r1,r2 in [("P","C"),("P","A"),("C","A")]:
                try: stat,pv=wilcoxon(vals[r1],vals[r2],zero_method="wilcox")
                except ValueError: stat,pv=np.nan,1.0
                tests.append({"duplication_type":mode,"trajectory":SHORT[cl],
                              "comparison":f"{r1} vs {r2}","n":len(sub),
                              "statistic":stat,"p_raw":pv})
            ax.set_xticks(positions,roles)
            if ri==0: ax.set_title(SHORT[cl],fontsize=9,weight="bold",pad=5)
            if ci==0:
                ax.set_ylabel(r"$\log_2$(TPM + 1)")
                ax.text(-.42,.5,mode,transform=ax.transAxes,
                        ha="right",va="center",fontsize=9,weight="bold")
                ax.text(.04,.88,f"n={len(sub)}",transform=ax.transAxes,
                        ha="left",va="top",fontsize=6.5,color="#555555")
            else: ax.set_yticklabels([])
            clean(ax)
    q=bh([r["p_raw"] for r in tests])
    for r,z in zip(tests,q): r["p_BH"]=z; r["significance"]=star(z)
    top=max(ax.get_ylim()[1] for ax in axes.values())
    for ax in axes.values(): ax.set_ylim(0,top+1.18)
    for ri in range(3):
        for ci in range(4):
            ax=axes[(ri,ci)]; block=tests[(ri*4+ci)*3:(ri*4+ci+1)*3]; base=top+.08
            bracket(ax,0,1,base,.07,block[0]["significance"])
            bracket(ax,1,2,base,.07,block[2]["significance"])
            bracket(ax,0,2,base+.54,.07,block[1]["significance"])
    handles=[Line2D([0],[0],marker="s",color="none",markerfacecolor=RCOL[r],
                    markeredgecolor="none",markersize=6,label=r) for r in roles]
    axes[(0,1)].legend(handles=handles,frameon=False,ncol=3,loc="upper center",
                       bbox_to_anchor=(.5,1.35),columnspacing=1.1)
    return pd.DataFrame(tests),pd.DataFrame(plotrows)
def main():
    raw=pd.read_csv(SRC,sep="\t")
    d=raw[(raw.analysis_status=="classified") & raw.two_method_agreement.isin(CLASSES)].copy()
    fig=plt.figure(figsize=(7.2,9.4))
    outer=fig.add_gridspec(2,1,height_ratios=[.92,3.05],
                          left=.16,right=.985,top=.972,bottom=.05,hspace=.25)
    top=outer[0].subgridspec(1,2,width_ratios=[.82,1.38],wspace=.25)
    axa=fig.add_subplot(top[0,0]); tab,pct=panel_a(axa,d)
    axb=fig.add_subplot(top[0,1]); panel_b(axb)
    cgs=outer[1].subgridspec(3,4,wspace=.10,hspace=.22)
    tests,long=panel_c(fig,cgs,d)
    fig.text(.018,.973,"A",fontsize=13,weight="bold",va="top")
    fig.text(.44,.973,"B",fontsize=13,weight="bold",va="top")
    fig.text(.018,outer[1].get_position(fig).y1+.018,"C",fontsize=13,weight="bold",va="top")
    base=OUT/"Figure2_hapA_two_method_consensus"
    fig.savefig(base.with_suffix(".svg")); fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"),dpi=300)
    fig.savefig(base.with_suffix(".tiff"),dpi=600,pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)
    tab.to_csv(OUT/"Figure2A_two_method_counts.tsv",sep="\t")
    pct.to_csv(OUT/"Figure2A_two_method_percentages.tsv",sep="\t")
    tests.to_csv(OUT/"Figure2C_two_method_all_pairwise_wilcoxon_BH.tsv",sep="\t",index=False)
    long.to_csv(OUT/"Figure2C_two_method_plot_data.tsv",sep="\t",index=False)
    summary={"two_method_consensus_n":len(d),
             "by_duplication_type":d.duplication_type.value_counts().to_dict(),
             "by_trajectory":d.two_method_agreement.value_counts().to_dict(),
             "tests":len(tests)}
    (OUT/"Figure2_two_method_analysis_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    im=Image.open(base.with_suffix(".tiff"))
    assert im.info["dpi"][0]>=599
    assert "<text" in base.with_suffix(".svg").read_text(encoding="utf-8")
    print(json.dumps(summary,indent=2))

if __name__=="__main__": main()















