"""Circular full-length ungapped fragment mapping and endpoint analysis."""
import argparse, json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D

BASES='ACGT'
COLORS=['#38a169','#3182ce','#d69e2e','#e53e3e']

def find_hits(query, reference, threshold=.95):
    """Return every circular window in forward-reference coordinates, 0-based."""
    n=len(reference); m=len(query)
    if m>n: return []
    import math
    target=reference+reference
    hits=[]
    allowed=math.floor(m*(1-threshold)+1e-10)
    for strand,s in [('+',query),('-',str(Seq(query).reverse_complement()))]:
        candidates=set()
        # With <=k mismatches at least one of k+1 disjoint seeds is exact.
        pieces=min(m,allowed+1)
        for part in range(pieces):
            left=part*m//pieces; right=(part+1)*m//pieces
            seed=s[left:right]; pos=target.find(seed)
            while pos>=0:
                start=pos-left
                if 0<=start<n: candidates.add(start)
                pos=target.find(seed,pos+1)
        if allowed>=m: candidates=set(range(n))
        for start in sorted(candidates):
            count=sum(a==b for a,b in zip(s,target[start:start+m]))
            if count/m+1e-12>=threshold: hits.append((strand,start,count/m))
    return hits

def plot_ends(sequences,out):
    counts=np.zeros((12,4),int)
    for s in sequences:
        for i,b in enumerate(s[:6]+s[-6:]): counts[i,BASES.index(b)]+=1
    total=len(sequences); freq=counts/total if total else counts.astype(float)
    entropy=-(freq*np.log2(np.where(freq>0,freq,1))).sum(1)
    bits=freq*(2-entropy)[:,None] if total else freq.copy()
    labels=[f"5p_{i}" for i in range(1,7)]+[f"3p_{i}" for i in range(-6,0)]
    for name,matrix in [('counts',counts),('frequency',freq),('information_bits',bits)]:
        df=pd.DataFrame(matrix,columns=list(BASES)); df.insert(0,'Position',labels); df.to_csv(out/f'end6_{name}.csv',index=False)
    for name,matrix,limit,ylabel in [('frequency_logo',freq,1,'Frequency'),('information_logo',bits,2,'bits')]:
        fig,ax=plt.subplots(figsize=(12,4))
        for i,row in enumerate(matrix):
            bottom=0
            for j in np.argsort(row):
                height=row[j]
                if height<=0: continue
                p=TextPath((0,0),BASES[j],size=1); b=p.get_extents()
                transform=Affine2D().translate(-b.x0,-b.y0).scale(.85/b.width,height/b.height).translate(i+.075,bottom)
                ax.add_patch(PathPatch(p,transform=transform+ax.transData,color=COLORS[j])); bottom+=height
        ax.set(xlim=(0,12),ylim=(0,limit),ylabel=ylabel,xlabel='Fragment ends (original 5′ to 3′ orientation)'); ax.set_xticks(np.arange(12)+.5,labels,rotation=45); ax.axvline(6,color='gray',ls='--'); fig.tight_layout()
        for ext in ['png','svg']: fig.savefig(out/f'end6_{name}.{ext}',dpi=200)
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,4)); x=np.arange(12)
    for j,b in enumerate(BASES): ax.bar(x+(j-1.5)*.2,freq[:,j],.2,label=b,color=COLORS[j])
    ax.set_xticks(x,labels,rotation=45); ax.set(ylim=(0,1),ylabel='Frequency'); ax.axvline(5.5,color='gray',ls='--'); ax.legend(); fig.tight_layout()
    for ext in ['png','svg']: fig.savefig(out/f'end6_frequency_bars.{ext}',dpi=200)
    plt.close(fig)

def analyze(fragments,reference,output,threshold=.95):
    if not 0<=threshold<=1: raise ValueError('Identity must be between 0 and 1')
    refs=list(SeqIO.parse(reference,'fasta'))
    if len(refs)!=1: raise ValueError('Reference FASTA must contain exactly one plasmid')
    ref=str(refs[0].seq).upper(); n=len(ref)
    if not n or set(ref)-set(BASES): raise ValueError('Reference must be nonempty ACGT DNA')
    records=list(SeqIO.parse(fragments,'fasta'))
    if not records: raise ValueError('No fragments found')
    if len({r.id for r in records})!=len(records): raise ValueError('Duplicate fragment IDs: make IDs unique first')
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    rows=[]; matched=[]; first=[]; last=[]; rejected=[]; short=[]
    for index,r in enumerate(records,1):
        q=str(r.seq).upper(); invalid=not q or bool(set(q)-set(BASES)); hits=[] if invalid else find_hits(q,ref,threshold)
        if hits:
            matched.append(SeqRecord(Seq(q),id=r.id,description=f'mapping_count={len(hits)}'))
            if len(q)>=6:
                first.append(SeqRecord(Seq(q[:6]),id=r.id,description='first6 original_orientation'))
                last.append(SeqRecord(Seq(q[-6:]),id=r.id,description='last6 original_orientation'))
            else: short.append(r.id)
        else: rejected.append({'Fragment_ID':r.id,'Length':len(q),'Reason':'non_ACGT_or_empty' if invalid else ('longer_than_plasmid' if len(q)>n else 'below_identity_threshold')})
        for strand,start,identity in hits:
            end=(start+len(q)-1)%n
            five=start if strand=='+' else end
            three=end if strand=='+' else start
            # Boundaries are indexed by the forward-reference base immediately before the cut.
            boundary5=start if strand=='+' else (end+1)%n
            boundary3=(end+1)%n if strand=='+' else start
            rows.append([r.id,start+1,end+1,five+1,three+1,strand,len(q),identity,start+len(q)>n,boundary5,boundary3])
        if index%1000==0 or index==len(records): print(f'{index}/{len(records)} ({index/len(records):.1%})',flush=True)
    columns=['Fragment_ID','Start_bp','End_bp','Fragment_5prime_bp','Fragment_3prime_bp','Strand','Length','Identity','Wraps_origin','Cut_5prime_after_bp','Cut_3prime_after_bp']
    pd.DataFrame(rows,columns=columns).to_csv(out/'fragment_matches.csv',index=False)
    pd.DataFrame(rejected,columns=['Fragment_ID','Length','Reason']).to_csv(out/'unmatched_fragments.csv',index=False)
    SeqIO.write(matched,out/'matched_fragments_both_strands.fasta','fasta')
    SeqIO.write(first,out/'first6.fasta','fasta'); SeqIO.write(last,out/'last6.fasta','fasta'); combined=[]
    for label,group in [('first6',first),('last6',last)]:
        combined.extend(SeqRecord(r.seq,id=r.id+'|'+label,description='') for r in group)
    SeqIO.write(combined,out/'ends6_combined.fasta','fasta')
    plot_ends([str(r.seq) for r in matched if len(r)>=6],out)
    summary={'reference_length':n,'identity_threshold':threshold,'input_fragments':len(records),'matched_unique_fragments':len(matched),'alignment_rows':len(rows),'unmatched_fragments':len(rejected),'end6_fragments':len(first),'shorter_than_6_excluded_from_end_analysis':short,'mapping_mode':'full-query ungapped circular exhaustive; all qualifying sites','end_statistics_weight':'one per matched input fragment, original 5prime-to-3prime orientation','cut_interpretation':'fragment boundary candidates, not experimentally proven cleavage sites'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    from plot_figures import plot_figures
    plot_figures(out)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--fragments',required=True); p.add_argument('--reference',required=True); p.add_argument('--output',default='circular_results'); p.add_argument('--identity',type=float,default=.95); args=p.parse_args()
    print(json.dumps(analyze(args.fragments,args.reference,args.output,args.identity),indent=2))
