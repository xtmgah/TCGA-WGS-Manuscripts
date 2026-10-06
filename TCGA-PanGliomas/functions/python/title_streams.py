"""Translate selected PDF text runs without reserializing drawing commands."""
from io import BytesIO
import copy
import pypdf
from pypdf.generic import TextStringObject,ByteStringObject,NameObject
from pypdf.generic._data_structures import read_non_whitespace,read_until_regex,read_object
from pypdf._cmap import get_encoding


def parse(data):
 stream=BytesIO(data);items=[];args=[];start=None
 while True:
  peek=read_non_whitespace(stream)
  if not peek:break
  stream.seek(-1,1);pos=stream.tell()
  if peek==b'%':
   while peek not in (b'\r',b'\n',b''):peek=stream.read(1)
   continue
  if start is None:start=pos
  if peek.isalpha() or peek in [b"'",b'"']:
   op=read_until_regex(stream,NameObject.delimiter_pattern);assert op!=b'BI','Unexpected inline raster'
   items.append((args,op,start,stream.tell()));args=[];start=None
  else:args.append(read_object(stream,None,None))
 return items


def shift(doc,path,targets,fitz):
 reader=pypdf.PdfReader(path);found=set();patches={};cache={};seen=set();page=doc[0];to_page=page.transformation_matrix
 # PDF numbers do not support scientific notation, including tiny translations.
 def matrix_bytes(m):return (' '.join(format(v,'.12f').rstrip('0').rstrip('.') or '0' for v in m)+' Tm').encode()
 def decode(args,op,font):
  key=font.indirect_reference.idnum
  if key not in cache:cache[key]=get_encoding(font)
  enc,cmap=cache[key];text=''
  for val in (args[0] if op==b'TJ' else [args[-1]]):
   if not isinstance(val,(TextStringObject,ByteStringObject)):continue
   raw=val.original_bytes if isinstance(val,TextStringObject) else bytes(val)
   rawtext=raw.decode(enc,errors='replace') if isinstance(enc,str) else ''.join(enc.get(c,chr(c)) for c in raw)
   text+=''.join(cmap.get(c,c) for c in rawtext)
  return text
 def walk(stream,res,state,stack=None):
  ref=stream.indirect_reference;key=ref.idnum;res=stream.get('/Resources',res);data=stream.get_data();ops=parse(data);state=copy.copy(state);stack=[] if stack is None else stack;tm=fitz.Matrix(1,1);tlm=fitz.Matrix(1,1);positioned=True;pending=None
  def commit(begin,end,before,after):
   patch=(begin,end,before,after)
   if (key,begin) in seen:assert patch in patches[key]
   else:patches.setdefault(key,[]).append(patch);seen.add((key,begin))
  for j,(args,op,start,end) in enumerate(ops):
   if op==b'q':stack.append(copy.copy(state))
   elif op==b'Q':state=stack.pop()
   elif op==b'cm':state['ctm']=fitz.Matrix(*map(float,args))*state['ctm']
   elif op==b'Do':
    sub=res['/XObject'][args[0]].get_object()
    if sub.get('/Subtype')=='/Form':
     child=copy.copy(state);child['ctm']=fitz.Matrix(*map(float,sub.get('/Matrix',[1,0,0,1,0,0])))*state['ctm'];walk(sub,res,child)
   elif op==b'BT':tm=fitz.Matrix(1,1);tlm=fitz.Matrix(1,1);positioned=True
   elif op==b'Tf':state['font']=res['/Font'][args[0]].get_object();state['size']=float(args[1])
   elif op==b'Tm':tm=fitz.Matrix(*map(float,args));tlm=fitz.Matrix(tm);positioned=True
   elif op in [b'Td',b'TD']:
    tx,ty=map(float,args);tlm=fitz.Matrix(1,0,0,1,tx,ty)*tlm;tm=fitz.Matrix(tlm);positioned=True
    if op==b'TD':state['leading']=-ty
   elif op==b'TL':state['leading']=float(args[0])
   elif op==b'T*':tlm=fitz.Matrix(1,0,0,1,0,-state.get('leading',0))*tlm;tm=fitz.Matrix(tlm);positioned=True
   elif op in [b'Tj',b'TJ',b"'",b'"']:
    text=decode(args,op,state['font']);origin=fitz.Point(tm.e,tm.f)*state['ctm']*to_page
    following=[o for _,o,_,_ in ops[j+1:]]
    reset=next((o for o in following if o in [b'Tm',b'Td',b'TD',b'T*',b'BT',b'ET',b'Tj',b'TJ',b"'",b'"']),None)
    if pending is not None:
     matches=[(k,t) for k,t in enumerate(targets) if ''.join(text.split())==''.join(t['text'].split()) and abs(t['origin_pt'][1]-pending['y'])<.1 and abs(t['dx_pt']-pending['dx'])<.01]
     assert len(matches)==1,('Unexpected inline continuation',text,pending)
     k,target=matches[0];found.add(k);target['source_xref']=key;target['source_operation_offset']=start;target['inline_with_previous_title_run']=True
     if reset not in [b'Tj',b'TJ',b"'",b'"']:
      commit(pending['start'],end,pending['before'],pending['after']);pending=None
     positioned=False;continue
    for k,target in enumerate(targets):
     if ''.join(text.split())!=''.join(target['text'].split()):continue
     if abs(origin.x-target['origin_pt'][0])>.1 or abs(origin.y-target['origin_pt'][1])>.1:continue
     assert positioned and op in [b'Tj',b'TJ'],('Unsupported continued text run',target)
     # The original line matrix is restored after the show operation. Every
     # following text show must first position itself, so no glyph advance is lost.
     inv=~state['ctm'];v=fitz.Point(target['dx_pt'],-target.get('dy_pt',0))*inv-fitz.Point(0,0)*inv;new=fitz.Matrix(tm);new.e+=v.x;new.f+=v.y
     before=b'\n'+matrix_bytes(new)+b'\n';after=b'\n'+matrix_bytes(tm)+b'\n'
     if reset in [b'Tj',b'TJ',b"'",b'"']:pending={'start':start,'before':before,'after':after,'dx':target['dx_pt'],'y':target['origin_pt'][1]}
     else:commit(start,end,before,after)
     found.add(k);target['source_xref']=key;target['source_operation_offset']=start;target['source_tm']=list(tm);target['shifted_tm']=list(new)
    positioned=False
  assert pending is None,'Unfinished inline title'
  return state
 state={'ctm':fitz.Matrix(1,1),'font':None,'size':0};stack=[]
 for content in reader.pages[0]['/Contents'] if isinstance(reader.pages[0]['/Contents'],list) else [reader.pages[0]['/Contents']]:
  state=walk(content.get_object(),reader.pages[0]['/Resources'],state,stack)
 assert len(found)==len(targets),('Unmatched title runs',[t for k,t in enumerate(targets) if k not in found])
 for key,items in patches.items():
  data=doc.xref_stream(key)
  for start,end,before,after in sorted(items,reverse=True):data=data[:start]+before+data[start:end]+after+data[end:]
  doc.update_stream(key,data)
 return targets
