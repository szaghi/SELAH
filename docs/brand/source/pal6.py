from pal import fit, cr, shade, WASHI, NIGHT
BASE=[("tsuyukusa","Tsuyukusa","露草","dayflower blue, the colour of morning dew","#1f6fc0","#6fb8ff"),
      ("kakitsubata","Kakitsubata","杜若","iris by the pond, blue-violet","#4b4fb0","#a3a7ff"),
      ("fujinezumi","Fujinezumi","藤鼠","wisteria grey, a quiet violet","#6a6596","#bdb7f0"),
      ("uguisu","Uguisu","鶯","bush-warbler olive","#6b7a1e","#c6d95c"),
      ("sango","Sango","珊瑚","coral","#d0573f","#ff8f78"),
      ("kuchiba","Kuchiba","朽葉","fallen-leaf amber","#9a5b13","#f2a64c"),
      ("neonhisui","Neon Hisui","翡翠","neon jade","#00a383","#2bffc8"),
      ("neonbotan","Neon Botan","牡丹","neon peony","#b5177d","#ff5ad6")]
PAL6={}
for key,name,kan,gloss,a,an in BASE:
    acc=fit(a,WASHI,3.6,False); link=fit(acc,WASHI,5.2,False)
    accn=fit(an,NIGHT,7.0,True); linkn=fit(shade(accn,1.08),NIGHT,9.0,True)
    PAL6[key]=dict(key=key,name=name,kanji=kan,gloss=gloss,acc=acc,link=link,accn=accn,linkn=linkn,
       c=dict(acc=cr(acc,WASHI),link=cr(link,WASHI),accn=cr(accn,NIGHT),linkn=cr(linkn,NIGHT)))
if __name__=="__main__":
    for p in PAL6.values(): print(p['name'],p['acc'],p['link'],p['accn'],p['linkn'],{k:round(v,1) for k,v in p['c'].items()})
