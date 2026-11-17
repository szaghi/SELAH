import colorsys
def lum(h):
    r,g,b=[int(h[i:i+2],16)/255 for i in (1,3,5)]
    f=lambda c: c/12.92 if c<=0.03928 else ((c+0.055)/1.055)**2.4
    return 0.2126*f(r)+0.7152*f(g)+0.0722*f(b)
def cr(a,b):
    la,lb=sorted([lum(a),lum(b)],reverse=True); return (la+0.05)/(lb+0.05)
def shade(h,k):
    r,g,b=[int(h[i:i+2],16)/255 for i in (1,3,5)]
    hh,l,s=colorsys.rgb_to_hls(r,g,b); l=max(0,min(1,l*k))
    r,g,b=colorsys.hls_to_rgb(hh,l,s); return '#%02x%02x%02x'%tuple(round(c*255) for c in (r,g,b))
def fit(h,bg,target,dark):
    k=1.0
    while cr(h if k==1 else shade(h,k),bg)<target:
        k = k*0.97 if not dark else k*1.03
    return h if k==1 else shade(h,k)
WASHI="#f3ede1"; NIGHT="#110f0d"
BASE=[("akatsuki","Akatsuki","暁","dawn","#cf531c","#ff8f4f"),
      ("yamabuki","Yamabuki","山吹","noon gold","#b57c00","#ffc23a"),
      ("yuyake","Yūyake","夕焼","sunset","#c0283c","#ff6f82"),
      ("tasogare","Tasogare","黄昏","twilight","#7a46b8","#c4a0ff"),
      ("ruri","Ruri","瑠璃","night lapis","#2456b8","#79a9ff"),
      ("tokiwa","Tokiwa","常盤","evergreen, the everlasting","#1d7a5c","#46d6a4")]
PAL={}
for key,name,kan,gloss,a,an in BASE:
    acc=fit(a,WASHI,3.6,False); link=fit(acc,WASHI,5.2,False)
    accn=fit(an,NIGHT,7.0,True); linkn=fit(shade(accn,1.08),NIGHT,9.0,True)
    PAL[key]=dict(key=key,name=name,kanji=kan,gloss=gloss,acc=acc,link=link,accn=accn,linkn=linkn,
       c=dict(acc=cr(acc,WASHI),link=cr(link,WASHI),accn=cr(accn,NIGHT),linkn=cr(linkn,NIGHT)))
if __name__=="__main__":
    for p in PAL.values(): print(p)
