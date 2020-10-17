import json
import pandas
from os.path import join as pjoin
from subprocess import call
import os
from tqdm import tqdm
import tempfile
import gzip
from datetime import datetime
from numpy import mean,median
import re


root = "/home/moritz/data/SRAprov/"
stot = lambda x : datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f')
data = pjoin(root, "data")

log_files = []
for v in tqdm(os.walk(root)):
    for vv in v[2]:
        if vv.endswith(".log.gz"):
            log_files += [pjoin(root,v[0], vv)]

logs = dict()
for l in tqdm(log_files):
    with gzip.open(l) as handle:
        logs[l.split("/")[-1][:-7]] =  [l.decode().strip() for l in handle if l.decode().startswith("===")]

times = {k : ( " ".join(l[0].split()[-3:-1]), " ".join(l[0-2].split()[-3:-1])) for k,l in logs.items()}
times = {kk : (stot(k[1])-stot(k[0])).seconds/3600 for kk, k in times.items()}


checkms = dict()
for l in tqdm(log_files):
    with gzip.open(l.replace('.log', '.bins.json')) as handle:
        checkms.update(json.load(handle))

with open(pjoin(data, "dbs", "all_checkms.tsv"), "w") as handle:
    handle.writelines(['Bin Id\tCompleteness\tContamination\n'] + [k + "\t" + str(v['Completeness']) + "\t" + str(v['Contamination']) +"\n" for k,v in checkms.items()] )

bin_count = {c.split('.')[0] : 0 for c in checkms.keys()}
for c in checkms:
    bin_count[c.split(".")[0]] += 1

with open(pjoin(data, "dbs", "sra_data.json")) as handle:
    sra_md = json.load(handle)

sra_md = pandas.DataFrame.from_dict({v['SRA_ID'] : v for k,v in sra_md.items()}, orient = "index")

del sra_md['TARGETED_LOCI']
del sra_md['LIBRARY_CONSTRUCTION_PROTOCOL']
del sra_md['study_id']
del sra_md['LIBRARY_SOURCE']
del sra_md['LIBRARY_STRATEGY']
del sra_md['Platform']
del sra_md['SRA_ID']
del sra_md['Statistics']
attrs = sra_md.sample_attributes.to_dict()

del sra_md['sample_attributes']
del sra_md['sample_id']
del sra_md['sample_name']

tt = [attrs[i].get('lat_lon', "NA") for i in sra_md.index]
sra_md['coord'] = ["NA" if t in ['missing', 'not collected', "Missing", "N/A", 'not applicable','Not Applicable', 'Not applicable', "NULL", "Unknown"] else t  for t in tt]

def tolat_lon(lat, lon):
    if lat >= 0 :
        lat = str(lat) + " N "
    else :
        lat = str(-lat) + " S "
    if lon >= 0 :
        lon = str(lon) + " E"
    else :
        lon = str(-lon) + " W"
    return lat + lon



for k,v in attrs.items():
    if 'Latitude Start' in v and 'Latitude End' not in v:
        sra_md['coord'][k] = tolat_lon(float(attrs[k]['Latitude Start']), float(attrs[k]['Longitude Start']))

for k,v in attrs.items():
    if 'Latitude Start' in v and 'Latitude End' in v:
        lat = (float(attrs[k]['Latitude Start']) + float(attrs[k]['Latitude End']))/2
        lon = (float(attrs[k]['Longitude Start']) + float(attrs[k]['Longitude End']))/2
        sra_md['coord'][k] = tolat_lon(lat, lon)

clean_char = set("0123456789,.-")
to_dec = lambda v : v[0] + v[1]/60 + v[2]/3600

for k,v in attrs.items():
    if 'geographic location (latitude)' in v and 'geographic location (longitude)' in v:
        lat = v['geographic location (latitude)'].replace(" DD", "")
        lon = v['geographic location (longitude)'].replace(" DD", "")
        if all([c in clean_char for c in lat + lon]):
            lat = float(lat.replace(",", "."))
            lon = float(lon.replace(",", "."))
            sra_md['coord'][k] = tolat_lon(lat, lon)
        else :
            lat_tri = [float(l) for l in re.split("[^0-9.,]",lat) if l != ""]
            lon_tri = [float(l) for l in re.split("[^0-9.,]",lon) if l != ""]
            if len(lat_tri) == 3 and len(lon_tri) == 3:
                if ("N" in lat or "S" in lat) and ("E" in lon or "W" in lon):
                    lat_str = str(to_dec(lat_tri)) + (" N " if "N" in lat else " S ")
                    lon_str = str(to_dec(lon_tri)) + (" E " if "E" in lat else " W ")
                    sra_md['coord'][k] = lat_str + lon_str
            elif k.startswith("ERR133318"):
                lon = lon + " E"
                sra_md['coord'][k] = lat + " " + lon
            elif k == 'ERR1299101':
                lon = lon.replace("'", "")
                sra_md['coord'][k] = tolat_lon(float(lat), float(lon))



for id, lat_lon in sra_md.coord[sra_md.coord.apply(lambda x : len(x.split()) != 4 and x != 'NA')].items():
    lat_lon_tri = [float(l) for l in re.split("[^0-9.]",lat_lon) if l != ""]
    print("none decimal e.g. 50032.74 50008.79 coords")
    if len(lat_lon_tri) ==6:
        lat_tri = lat_lon_tri[0:3]
        lon_tri = lat_lon_tri[3:]
    else :
        lat_tri = lat_lon_tri[0:2] + [0]
        lon_tri = lat_lon_tri[2:] + [0]
    if ("N" in lat_lon or "S" in lat_lon) and ("E" in lat_lon or "W" in lat_lon):
        lat_str = str(to_dec(lat_tri)) + (" N " if "N" in lat_lon else " S ")
        lon_str = str(to_dec(lon_tri)) + (" E " if "E" in lat_lon else " W ")
        sra_md['coord'][id] = lat_str + lon_str
    else :
        sra_md['coord'][id] = tolat_lon(to_dec(lat_tri), to_dec(lon_tri))



lat_long2lat = lambda x: (float(x.split()[0]) * (-1 if "S" in x else 1) ) if x != "NA" else x
lat_long2lon = lambda x: (float(x.split()[2]) * (-1 if "W" in x else 1) ) if x != "NA" else x

sra_md['lat'] = sra_md.coord.apply(lat_long2lat)
sra_md['lon'] = sra_md.coord.apply(lat_long2lon)

id2tax = {v[1] : v[0]  for k,v in sra_md[['taxon', 'taxon_id']].iterrows() if v[0] != ""}
sra_md['taxon'] = [id2tax[v] for v in sra_md.taxon_id]

bad_study = sra_md.loc['ERR4193663']['study']

lats = sra_md.loc[sra_md.study == bad_study, 'lat']
lons = sra_md.loc[sra_md.study == bad_study, 'lon']
sra_md.loc[sra_md.study == bad_study, 'lat'] = lons
sra_md.loc[sra_md.study == bad_study, 'lon'] = lats
sra_md.loc[sra_md.study == bad_study, 'coord'] =  [tolat_lon(float(b),float(a)) for a, b in zip(lats, lons)]


sra_md.to_csv(pjoin(data, "dbs", "sra_table.csv"), index_label= "SRA_ID")

"mOTUlize.py  -o data/dbs/mOTUs.json -k data/dbs/all_checkms.tsv -c 24 --similarities simis3.tsv  -F all_bins/* --MC 70"

with open(pjoin(data, "dbs", "mOTUs.json") ) as handle:
    full_mOTUs = json.load(handle)

get_goods = lambda g : [gg['name'] for gg in g['genomes'] if gg['checkm_complet'] > 70]
