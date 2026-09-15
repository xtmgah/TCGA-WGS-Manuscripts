"""Read-only IWA forensic inventory; compiled schemas are kept alongside evidence."""
from pathlib import Path
import sys, zipfile, snappy, json, csv, hashlib, argparse
from collections import Counter
cli=argparse.ArgumentParser(description='Read-only IWA inventory; manual figure assignment remains in keynote_panel_map.tsv.')
cli.add_argument('keynote',type=Path)
cli.add_argument('--schema-dir',type=Path,required=True)
cli.add_argument('--output-dir',type=Path,required=True)
args=cli.parse_args()
ROOT=args.output_dir.resolve()
ROOT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(args.schema_dir.resolve()))
import TSPArchiveMessages_pb2 as TSP
import TSDArchives_pb2 as TSD
import TSWPArchives_pb2 as TSWP
import KNArchives_pb2 as KN
from google.protobuf.json_format import MessageToDict

def varint(b, p):
    value=shift=0
    while True:
        v=b[p];p+=1;value|=(v&127)<<shift
        if not v&128:return value,p
        shift+=7

def parse_iwa(b):
    p=0;chunks=[]
    while p<len(b):
        kind=b[p];n=int.from_bytes(b[p+1:p+4],'little');p+=4
        chunk=b[p:p+n];p+=n
        chunks.append(snappy.uncompress(chunk) if kind==0 else chunk)
    raw=b''.join(chunks);p=0;objects=[]
    while p<len(raw):
        n,p=varint(raw,p);a=TSP.ArchiveInfo();a.ParseFromString(raw[p:p+n]);p+=n
        for mi in a.message_infos:
            body=raw[p:p+mi.length];p+=mi.length
            objects.append({'id':a.identifier,'type':mi.type,'body':body,
                            'refs':list(mi.object_references),'data_refs':list(mi.data_references)})
    assert p==len(raw)
    return objects

def decoded(row,klass):
    m=klass();m.ParseFromString(row['body']);return m

def fields(b):
    p=0;out=[]
    while p<len(b):
        key,p=varint(b,p);f,w=key>>3,key&7
        if w==0:v,p=varint(b,p)
        elif w==2:
            n,p=varint(b,p);v=b[p:p+n];p+=n
        elif w in (1,5):
            n=8 if w==1 else 4;v=b[p:p+n];p+=n
        else:raise ValueError(w)
        out.append((f,w,v))
    return out

def write_tsv(name,rows):
    if not rows:return
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    src=args.keynote.resolve()
    z=zipfile.ZipFile(src)
    iwas={n:parse_iwa(z.read(n)) for n in z.namelist() if n.endswith('.iwa')}
    objects={o['id']:dict(o,iwa=n) for n,oo in iwas.items() for o in oo}
    metadata=decoded(iwas['Index/Metadata.iwa'][0],TSP.PackageMetadata)
    data={d.identifier:MessageToDict(d,preserving_proto_field_name=True) for d in metadata.datas}
    write_tsv('metadata_media.tsv',[{'data_id':k,'original_plot_filename':v.get('preferred_file_name',''),
                                    'package_file_name':v.get('file_name',''),
                                    'sha256':hashlib.sha256(z.read('Data/'+v['file_name'])).hexdigest() if 'Data/'+v.get('file_name','') in z.namelist() else ''}
                                   for k,v in data.items()])
    nodes={k:decoded(o,KN.SlideNodeArchive) for k,o in objects.items() if o['type']==4}
    doc=decoded(next(o for o in iwas['Index/Document.iwa'] if o['type']==1),KN.DocumentArchive)
    show=decoded(objects[doc.show.identifier],KN.ShowArchive)
    if show.slideTree.HasField('rootSlideNode'):
        roots=[nodes[show.slideTree.rootSlideNode.identifier]]
    else:
        # Current Keynote encodes ordered SlideNode references in SlideTree field 2.
        tree=next(v for f,w,v in fields(objects[doc.show.identifier]['body']) if f==3)
        ordered_ids=[next(vv for ff,ww,vv in fields(v) if ff==1) for f,w,v in fields(tree) if f==2]
        roots=[nodes[k] for k in ordered_ids]
    slides=[]
    def traverse(n):
        if n.HasField('slide'):slides.append((n.slide.identifier,n.isHidden))
        for c in n.children:traverse(nodes[c.identifier])
    for n in roots:traverse(n)
    print('Slides:',slides,'canvas',show.size.width,show.size.height)
    inventories=[];media=[];texts=[];slideinfo=[]
    for num,(sid,hidden) in enumerate(slides,1):
        sr=objects[sid];slide=decoded(sr,KN.SlideArchive)
        slideinfo.append({'keynote_slide':num,'slide_id':sid,'iwa_file':sr['iwa'],'is_hidden':hidden})
        seen=set()
        def visit(oid,parent='',parent_pos=(0,0)):
            if oid in seen:return
            seen.add(oid);o=objects.get(oid)
            if o is None:return
            k=o['type'];m=None;drawable=None;children=[];text=''
            if k==3005:
                m=decoded(o,TSD.ImageArchive);drawable=m.super
            elif k==3008:
                m=decoded(o,TSD.GroupArchive);drawable=m.super;children=[c.identifier for c in m.children]
            elif k in (2011,7):
                m=decoded(o,TSWP.ShapeInfoArchive) if k==2011 else decoded(o,KN.PlaceholderArchive).super
                drawable=m.super.super
                st=objects.get(m.containedStorage.identifier)
                if st:text=''.join(decoded(st,TSWP.StorageArchive).text)
            elif k==3004:
                m=decoded(o,TSD.ShapeArchive);drawable=m.super
            row={'keynote_slide':num,'slide_id':sid,'object_id':oid,'object_type':k,'parent_object_id':parent,
                 'iwa_file':o['iwa'],'x':None,'y':None,'width':None,'height':None,'angle':None,'text':text}
            if drawable:
                g=drawable.geometry
                row.update(x=g.position.x+parent_pos[0],y=g.position.y+parent_pos[1],width=g.size.width,height=g.size.height,angle=g.angle)
            inventories.append(row)
            if text:texts.append(row)
            if k==3005:
                for did in o['data_refs']:
                    dd=data.get(did,{})
                    name=dd.get('preferred_file_name','')
                    media.append(dict(row,data_id=did,original_plot_filename=name,package_file_name=dd.get('file_name',''),is_preview='-small.' in name))
            for c in children:visit(c,oid,(row['x'] or 0,row['y'] or 0))
        for d in slide.drawables:visit(d.identifier)
    write_tsv('slide_order.tsv',slideinfo)
    write_tsv('drawable_objects.tsv',inventories)
    write_tsv('media_object_inventory.tsv',media)
    write_tsv('text_objects.tsv',texts)
    # Keep all raw objects/references, including non-visible styles and notes.
    with (ROOT/'all_iwa_objects.jsonl').open('w') as f:
        for o in objects.values():f.write(json.dumps({k:v for k,v in o.items() if k!='body'})+'\n')
    for s in slideinfo:
        n=s['keynote_slide'];print(n,[x['text'].replace('\n',' ') for x in texts if x['keynote_slide']==n])
        for m in media:
            if m['keynote_slide']==n and not m['is_preview']:print(' ',m['object_id'],m['x'],m['y'],m['width'],m['height'],m['original_plot_filename'])

if __name__=='__main__':main()
