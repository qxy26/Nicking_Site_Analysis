"""Generate terminal nucleotide logos and frequency bar charts."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from circular_fragments import plt,TextPath,PathPatch,Affine2D,BASES,COLORS

def plot_figures(output):
    out=Path(output)
    for end,offset in [('first6',0),('last6',6)]:
        labels=np.arange(1,7) if offset==0 else np.arange(-6,0)
        for kind,ylabel in [('frequency','Frequency'),('information_bits','bits')]:
            df=pd.read_csv(out/f'end6_{kind}.csv').iloc[offset:offset+6].copy(); df['Position']=labels; df.to_csv(out/f'{end}_{kind}.csv',index=False)
            mat=df[list(BASES)].to_numpy(); maximum=mat.sum(1).max(); limit=maximum*1.10 if maximum>0 else 1e-3
            fig,ax=plt.subplots(figsize=(6.5,3.7))
            for i,row in enumerate(mat):
                bottom=0
                for j in np.argsort(row):
                    h=row[j]
                    if h<=0: continue
                    path=TextPath((0,0),BASES[j],size=1); b=path.get_extents(); transform=Affine2D().translate(-b.x0,-b.y0).scale(.8/b.width,h/b.height).translate(i+.1,bottom)
                    ax.add_patch(PathPatch(path,transform=transform+ax.transData,color=COLORS[j])); bottom+=h
            ax.set(xlim=(0,6),ylim=(0,limit),ylabel=ylabel,xlabel='Position',title='First 6 bases' if offset==0 else 'Last 6 bases'); ax.set_xticks(np.arange(6)+.5,labels); ax.spines[['top','right']].set_visible(False); fig.tight_layout()
            for ext in ['png','svg']: fig.savefig(out/f'{end}_{kind}_logo.{ext}',dpi=240)
            plt.close(fig)
        df=pd.read_csv(out/'end6_counts.csv').iloc[offset:offset+6].copy(); df['Position']=labels; df.to_csv(out/f'{end}_counts.csv',index=False)
        freq=pd.read_csv(out/f'{end}_frequency.csv')[list(BASES)].to_numpy(); fig,ax=plt.subplots(figsize=(6.5,3.7)); x=np.arange(6)
        for j,b in enumerate(BASES): ax.bar(x+(j-1.5)*.2,freq[:,j],.2,color=COLORS[j],label=b)
        ax.set_xticks(x,labels); ax.set(xlabel='Position',ylabel='Frequency',ylim=(0,freq.max()*1.15 if freq.max()>0 else 1)); ax.legend(frameon=False,ncol=4); ax.spines[['top','right']].set_visible(False); fig.tight_layout()
        for ext in ['png','svg']: fig.savefig(out/f'{end}_frequency_bars.{ext}',dpi=240)
        plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True); a=p.parse_args(); plot_figures(a.output)
