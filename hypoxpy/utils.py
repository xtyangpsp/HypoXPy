#This file contains utility functions for the package.
import pandas as pd
import numpy as np
import os,json,shutil,glob,time
from importlib import resources as impresources
from hypoxpy import templates
from obspy import UTCDateTime
import warnings
from tqdm import tqdm
from datetime import datetime
#####
def pick_vars(*names, scope):
    """
    Pick and pass variable names to dictionary items.
    """
    return {name: scope[name] for name in names}

def basename_list():
    """
    Output a list of available template basenames for HypoXPy.
    """
    return ['hypoinv','hypodd']
#
def cat_files(pattern, fout):
    
    with open(fout, 'w') as out:
        for fname in sorted(glob.glob(pattern)):
            with open(fname) as f:
                shutil.copyfileobj(f, out)
#
def tic():
    return time.perf_counter()

def toc(t0, label):
    dt = time.perf_counter() - t0
    print(f"[TIMER] {label}: {dt:.2f} s")

def get_filelist(dir=None,extension=None,pattern=None,sort=True):
    """
    Get list of files with absolute path, by specifying the format extension.
    Modified from SeisGo.utils.get_filelist()

    ===========PARAMETERS=============
    dir: directory containing the files.
    extension: file extension (the ending format tag), for example "h5" for asdf file.
    pattern: pattern to use in searching. Wildcards are NOT considered here.
    sort: (optional) to sort the list, default is True.

    ============RETURN=============
    flist: the list of file names with paths.
    """
    if dir is None:
        dir="."
    if extension is None:
        flist=[os.path.join(dir,f) for f in os.listdir(dir)]
    else:
        flist=[os.path.join(dir,f) for f in os.listdir(dir) if f[-len(extension):].lower()==extension.lower()]
    if pattern is not None:
        flist2=[]
        for f in flist:
            if f.find(pattern)>=0: flist2.append(f)
        flist=flist2
    if sort:
        return  sorted(flist)
    else:
        return  flist
##    
def stainfo_json2csv(infile,outfile=None,informat=None):
    """
    Convert seismic station information from json database format to CSV format (with the option of return a Pandas.DataFrame object).
    
    ====== PARAMETERS =======
    infile: input file name of the json database. [required]
    outfile: output file name. Default is None (return the dataframe object). If specified (not None), the dataframe will be saved as a csv file.
    informat: input json format, either from EQTransformer (eqt) or GAMMA.
    ====RETURN====
    outdata: [optional] return the object when outfile is not specified (None)

    """
    if informat is None:
        informat='EQT'
        warnings.warn('default input station format (eqt or gamma) not set. will use: '+informat)
    # allocate arrays/list for output columns
    net_all=[]
    sta_all=[]
    chan_all=[]
    lat_all=[]
    lon_all=[]
    ele_all=[]
    if informat.lower() == 'eqt':
        indata=pd.read_json(infile,orient='index')
        indata['station']=indata.index
        indata=indata.reset_index(drop=True)
        for i in range(len(indata)):
            # [lat[i], lon[i],ele[i]]=coords_all[i]
            for j in range(len(indata.channels[i])):
                net_all.append(indata.network[i])
                sta_all.append(indata.station[i])
                chan_all.append(indata.channels[i][j])
                lat_all.append(indata.coords[i][0])
                lon_all.append(indata.coords[i][1])
                ele_all.append(indata.coords[i][2])
    elif informat.lower() == 'gamma':
        with open(infile, "r") as file:
            indata = json.load(file)
        #
        for station,details in indata.items():
            net,sta,loc,ch=station.split('.')
            # print(net,sta,loc,ch)
            component=details['component']
            for i in range(len(component)):
                chan=ch+component[i]
                net_all.append(net)
                sta_all.append(sta)
                chan_all.append(chan)
                lon_all.append(details['longitude'])
                lat_all.append(details['latitude'])
                ele_all.append(details['elevation(m)'])
    ##
    outdata=pd.DataFrame(list(zip(net_all,sta_all,chan_all,lat_all,lon_all,ele_all)),columns=['network','station','channel',
                                                                             'latitude','longitude','elevation'])
    if outfile is not None: #save to file.
        fhead=os.path.split(outfile)[0]
        if len(fhead)>0:
            if not os.path.isdir(fhead):os.makedirs(fhead)
        outdata.to_csv(outfile, index=False)
    #return the pandas.DataFrame object.
    return outdata
#
def reformat_stainfo_hypoinverse(infile,outfile, informat='csv',channel=True,channel_default='HHZ',
                     force_channel_type=None,lat_code='N', lon_code='W',rename_component_dict=None,
                     ignore_component=False):
    """
    Reformat station information file for hypoinver run. Input CSV needs to have at least the following columns:
    network,station,channel,latitude,longitude,elevation.

    Additional columns will be ignored.

    ====== PARAMETERS =======
    infile: input file name (including path).
    outfile: output file name (including path).
    informat: input file format, choose from 'csv', 'json', 'json-eqt','json-gamma'. Default is 'csv'.
    channel: is channel column available. Default is True. Otherwise, channel_default is used.
    channel_default: when channel is not available (channel=False), this value will be used. Default='HHZ'
    force_channel_type: treat all channel types (i.e., EH?, BH?, or HH?) as the same. Default=None (use true channel).
    rename_component_dict: dictionary used to rename component label, e.g., HH1 to HHN. Default is {'1':'N','2':'E'}
    lat_code: code for latitudes. Default is 'N'.
    lon_code: code for longitudes. Default is 'W'. lat_code and lon_code are used to format the coordinates.
    ignore_component: if True, only save the channel type information, such as BH instead of BHZ. Default False.
    """
    if force_channel_type is not None:
        if len(force_channel_type) != 2:
            raise ValueError('force_channel_type has to be two characters. Wrong value: '+force_channel_type)
    if informat.lower() == 'csv':
        if infile[-4:].lower() =='json':
            raise ValueError(infile+' may be a json file. change informat argument to json.')
        indata=pd.read_csv(infile)
    elif informat.lower() == 'json' or informat.lower() == 'json-eqt':
        if infile[-3:].lower() =='csv':
            raise ValueError(infile+' may be a csv file. change informat argument to csv.')
        indata=stainfo_json2csv(infile)
    elif informat.lower() == 'json-gamma':
        if infile[-3:].lower() =='csv':
            raise ValueError(infile+' may be a csv file. change informat argument to csv.')
        indata=stainfo_json2csv(infile,informat='gamma')
    else:
        raise ValueError('input file format of %s not recoganized. Use "csv" or "json".'%(informat))
    #
    if rename_component_dict is not None: 
        rename_keys=list(rename_component_dict.keys())

    fhead=os.path.split(outfile)[0]
    if len(fhead)>0:
        if not os.path.isdir(fhead):os.makedirs(fhead)
    fout = open(outfile,'w')
    for i in range(len(indata)):
        net,sta,lat,lon,ele=[indata['network'][i],indata['station'][i],indata['latitude'][i],indata['longitude'][i],indata['elevation'][i]]
        lat, lon, ele = abs(lat), abs(lon), int(ele)
        lat_deg = int(lat)
        lat_min = 60*(lat-int(lat))
        lon_deg = int(lon)
        lon_min = 60*(lon-int(lon))
        lat = '{:2} {:7.4f}{}'.format(lat_deg, lat_min, lat_code)
        lon = '{:3} {:7.4f}{}'.format(lon_deg, lon_min, lon_code)
        if not channel:
            chan=channel_default
        else:
            chan=indata['channel'][i]
        if force_channel_type is not None:
            chan=force_channel_type+chan[-1]
        if rename_component_dict is not None: #force to rename component label.
            if chan[-1] in rename_keys:
                chan=chan[0:2]+str(rename_component_dict[chan[-1]])
        if ignore_component:
            chan=chan[0:2] #drop the component information.
        # hypoinverse format 2
        fout.write("{:<5} {}  {}  {}{}{:4}\n".format(sta, net, chan,lat, lon, ele))
    fout.close()

#
def reformat_stainfo_hypodd(infile,outfile,informat="csv",combine_net_sta=True,elev_unit="m"):
    """
    Reformat station information file for HypoDD station.dat.

    Input must contain at least:
        network, station, latitude, longitude, elevation

    Output format (HypoDD):
        STA LAT LON ELEV(km)

    ===== PARAMETERS =====
    infile : str
        Input station file
    outfile : str
        Output HypoDD station.dat
    informat : str
        'csv', 'json', or 'json-gamma'
    combine_net_sta : bool
        If True, station name = NET.STA (strongly recommended)
    elev_unit : str
        'm' or 'km'. Default assumes meters.
    """

    # ---- read input ----
    if informat.lower() == "csv":
        indata = pd.read_csv(infile)
    elif informat.lower() in ("json", "json-eqt"):
        indata = stainfo_json2csv(infile)
    elif informat.lower() == "json-gamma":
        indata = stainfo_json2csv(infile, informat="gamma")
    else:
        raise ValueError(f"Unrecognized informat: {informat}")

    # ---- output directory ----
    fhead = os.path.split(outfile)[0]
    if fhead and not os.path.isdir(fhead):
        os.makedirs(fhead)

    # ---- write station.dat ----
    donelist=[] #to avoid duplicate stations.
    with open(outfile, "w") as fout:
        for _, row in indata.iterrows():

            net = str(row["network"]).strip()
            sta = str(row["station"]).strip()

            if combine_net_sta:
                sta_name = f"{net}.{sta}"
            else:
                sta_name = sta
            if sta_name in donelist:
                continue
            donelist.append(sta_name)
            lat = float(row["latitude"])
            lon = float(row["longitude"])
            ele = float(row["elevation"])

            # elevation to km
            if elev_unit.lower() == "m":
                ele /= 1000.0

            fout.write(
                f"{sta_name:<7} {lat:9.4f} {lon:10.4f} {ele:7.3f}\n"
            )
#
def reformat_stainfo_legacy(fin,fout):
    """ 
    Reformat input station file for HypoDD. This is a legacy function moidified from hypo-interface-py. I kept it here for backward compatibility.
    We use net.sta, lat, lon format, instead of sta, lat, lon, to keep network info.
    Parameters
    ----------
    fin : str
        Input station file path.
    fout : str
        Output station file path.   
    """
    foutid=open(fout,'w')
    
    done_list = []
    f=open(fin); 
    lines=f.readlines(); 
    f.close()
    for line in lines:
        codes = line.split(',')
        netsta = codes[0].strip()
        if netsta in done_list: continue
        lat, lon = [float(code) for code in codes[1:3]]
        foutid.write('{} {} {}\n'.format(netsta, lat, lon))
        done_list.append(netsta)
    foutid.close()
def get_template_list(basename,pattern='',fullpath=False):
    """ 
    Get list of available template files for a specified basename.
    ==========PARAMETERS=============
    basename: the template basename, e.g., 'hypoinv' or 'hypodd

    ==========RETURN=============
    templatetail: list of available template file tails for the specified basename.
    """
    if basename not in basename_list():
        raise ValueError('basename %s NOT recognized. Available basenames are: %s'
                         %(basename,str(basename_list())))
    
    templatedir=impresources.files(templates)
    templatelist=get_filelist(templatedir,pattern='%s_template_%s'%(basename,pattern),sort=True)
    if fullpath:
        return templatelist
    else:
        templatetail=[]
        for tf in templatelist:
            ftail=os.path.split(tf)[1]
            templatetail.append(ftail)
        #
        return templatetail
#
def load_template(template_name=None):
    """
    Load template file lines into a list.
    ==========PARAMETERS=============
    template_name: template file name. If None, an error will be raised.
    ==========RETURN=============
    lines: list of lines from the template file.
    """
    if template_name is None:
        raise ValueError('template_name NOT specified. run get_template_list() to get a list of available templates.')
    else:
        if os.path.isfile(template_name):
            inp_file=template_name
        else:
            inp_file = (impresources.files(templates) / template_name)
        f=open(inp_file)
        lines=f.readlines()
        f.close()
        
        return lines
# format output
def write_csv(fout, line, evid, lat_code, lon_code, mag_dict=None):
    """
    Write hypoinverse output lines into CSV format with earthquake source parameters.

    ==========PARAMETERS==============
    fout: output CSV file name.
    line: hypoinverse output earthquake parameter line.
    evid: evid for record. 
    lat_code and lon_code: the character labeling the latitude and longitude values. E.g., lat_code='N', lon_code='W' for 
        geographical locations in the northern hemisphere with west longitude.
    mag_dict: a dictionary containing the magnitude of each event id in the formm of {'id',magnitude}. Default None.
    """
    #
    grd_ele = 0 #the original code from Hypo_Interface_Py added the correction for grid elevation. Force to 0 now. 
    #Not clear why this was added. Will work on this. Be aware of this point. Noted by Xiaotao for HypoInvPy
    
    # mag = 0.0 #force magnitude to 0 for now. Will add magnitude information later.
    codes = line.split()
    date, hrmn, sec = codes[0:3]
    dtime = date + hrmn + sec.zfill(5)
    lat_deg = float(line[20:22])
    lat_min = float(line[23:28])
    lat = lat_deg + lat_min/60 if lat_code=='N' else -lat_deg - lat_min/60
    lon_deg = float(line[29:32])
    lon_min = float(line[33:38])
    lon = lon_deg + lon_min/60 if lon_code=='E' else -lon_deg - lon_min/60
    dep = float(line[38:45])
    mag = float(line[45:52])
    if mag_dict is not None:
        mag = mag_dict[str(evid)]
    
    fout.write('{},{:.4f},{:.4f},{:.1f},{:.2f},{}\n'.format(dtime, lat, lon, dep+grd_ele, mag, evid))
# read fpha with evid
def load_hypo_phasedata(fin,separator=','):
    """
    Load phase data from a phase file with event id (evid) at the event header. 
    Each evid element consists of multiple lines (event header + picks).
    ==========PARAMETERS=============
    fin: input phase file name.
    separator: separator used in the phase file. Default is comma (,).
    ==========RETURN=============
    pha_dict: a dictionary containing phase data for each event id (evid).
    """
    pha_dict = {}
    evid = None

    with open(fin) as f:
        for line in f:
            codes = line.split(separator)
            if len(codes[0]) >= 14:
                evid = codes[-1].strip()
                pha_dict[evid] = []
            elif evid is not None:
                pha_dict[evid].append(line)
    return pha_dict
#
def resolve_event_depth_column(df,label=None):
    """
    Resolve depth column and unit automatically.

    Returns
    -------
    label : str
        Column name to use
    scale : float
        Scale factor to convert to km
    """
    # 1. User explicitly specifies
    if label is not None:
        if label not in df.columns:
            raise KeyError(f"Depth column '{label}' not found in dataframe.")
        return label, 1.0

    # 2. Column name based detection
    if "depth_km" in df.columns:
        return "depth_km", 1.0

    if "depth(m)" in df.columns:
        return "depth(m)", 1e-3

    if "depth_m" in df.columns:
        return "depth_m", 1e-3

    if "depth" in df.columns:
        # 3. Infer from values
        vals = df["depth"].dropna().values
        if len(vals) == 0:
            raise ValueError("Depth column exists but contains no valid values.")

        median = np.nanmedian(vals)

        if median > 1000:  # almost certainly meters
            warnings.warn(
                "Depth column 'depth' appears to be in meters. Converting to km."
            )
            return "depth", 1e-3
        elif median < 100:  # likely km
            warnings.warn(
                "Depth column 'depth' assumed to be in kilometers."
            )
            return "depth", 1.0
        else:
            raise ValueError(
                "Ambiguous depth units in 'depth' column. "
                "Please specify depth column explicitly."
            )

    # 4. Nothing usable
    raise KeyError(
        "No recognizable depth column found. "
        "Expected one of: depth_km, depth_m, depth(m), depth."
    )

# functions for GaMMA to Hypoinverse phase file conversion
def format_pick_line_gamma2hypoinverse(pick, default_component=None): #gamma to hypoinverse
    """
    Formats a single pick line for Hypoinverse phase file.
    Includes network, station, channel, component info.
    ==========PARAMETERS=============
    pick: a single pick record (Pandas Series).
    default_component: default component to use if not specified in picks.
    ==========RETURN=============
    formatted pick line in Hypoinverse format.
    """
    try:
        net, sta, cha, comp = pick["id"].split(".")
    except ValueError:
        warnings.warn(f"Invalid pick id: {pick['id']}")
        return None

    if not comp and default_component:
        comp = default_component

    phase = str(pick["type"]).upper()
    if phase not in ("P", "S"):
        return None

    network_code, station_code, comp_code, channel_code = pick['id'].split('.')
    phase = pick['type']
    phase_weight = min(max(int((1-pick['prob']) / (1 - 0.3) * 4) - 1, 0), 3)
    picktime = datetime.strptime(pick["timestamp"], "%Y-%m-%dT%H:%M:%S.%f")
    pickmin = picktime.strftime("%Y%m%d%H%M")
    picksec = picktime.strftime("%S%f")[:-4]
    phase_amplitude=pick['phase_amplitude']

    if default_component is not None and len(comp_code)<1:
        comp_code = default_component

    templine = f"{station_code:<5}{network_code:<2}  {channel_code:<2}{comp_code:<1}"
    
    if phase.upper() == 'P':
        pick_line = f"{templine:<13} P {phase_weight:<1d}{pickmin} {picksec}                    {phase_amplitude:<7}"
    elif phase.upper() == "S":
        pick_line = f"{templine:<13}   4{pickmin} {'':<12}{picksec} S {phase_weight:<1d}    {phase_amplitude:7f}"
    else:
        warnings.warn("Phase Type Error: "+phase.upper())
        return None
    
    return pick_line

#
def format_event_line_hypoinverse(ev, evid, evdep_label, evdep_scale):
    """
    Format a HYPOINVERSE event header line (fixed-width), HYPOINVERSE ARCHIVE-2000 FILE, NO SHADOWS.
    Magnitude is set to 0.0 to avoid negative magnitude issues, which caurse HYPOINVERSE to fail due to length restriction.
    """
    fixmag = True  # force magnitude to 0.0 to avoid negative magnitude issues.
    if fixmag:
        mag = 0.0
    else:
        mag = ev.get("magnitude", 0.0)
    # get event parameters 

    ot = ev["time"]

    lat = ev["latitude"]
    lon = ev["longitude"]
    depth = ev[evdep_label] * evdep_scale
    
    # ---- Latitude ----
    south = "S" if lat < 0 else " "
    lat = abs(lat)
    lat_degree = int(lat)
    lat_minute = (lat - lat_degree) * 60 * 100

    # ---- Longitude ----
    west = "W" if lon < 0 else " "
    lon = abs(lon)
    lon_degree = int(lon)
    lon_minute = (lon - lon_degree) * 60 * 100

    event_time = datetime.strptime(ot, "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y%m%d%H%M%S%f")[:-4]

    event_line = f"{event_time}{abs(lat_degree):2d}{south}{abs(lat_minute):4.0f}{abs(lon_degree):3d}{west}{abs(lon_minute):4.0f}{depth:5.0f}" 
    event_final = event_line + f"{mag:3.1f}{' ':97}{evid:10}"

    return event_final
####
def convert_gamma2hypoinverse(eventfile,pickfile,outfile="phase_output.phs",qc=False,separator=",",default_component=None,
    verbose=False,evid_label="event_id",mapping_evid=True,evid_label_mapped=None,save_cleaned_data=False,cleaned_eventfile=None,
    cleaned_pickfile=None):
    """
    Convert GAMMA event + pick CSV files to HYPOINVERSE phase format.
    ==========PARAMETERS=============
    eventfile: input event CSV file from GaMMA.
    pickfile: input pick CSV file from GaMMA.
    outfile: output Hypoinverse phase file name.
    qc: if True, only events with both P and S picks are written. Default is False.
    separator: separator used in the input CSV files. Default is comma (,).
    default_component: default component to use if not specified in picks.
    verbose: if True, print progress messages. Default is False.
    evid_label: column name for event ID in the input event and pick files. Default is "event_id".
    mapping_evid: if True, remap event IDs to integers starting from 1 for HypoDD compatibility. Default is True. 
                If False, original event IDs are used (may cause issues if IDs are not integers or not starting from 1).
    evid_label_mapped: evid label after remapping. This is only effective when mapping_evid is True.
    save_cleaned_data: if True, save cleaned event and pick data after removing invalid entries. Default is False.
                Set to True automatically if mapping_evid is True.
    cleaned_eventfile: output file name for cleaned event data. Default is None (auto-generate).
    cleaned_pickfile: output file name for cleaned pick data. Default is None (auto-generate).
    ==========RETURN=============
    None
    """
    events = pd.read_csv(eventfile, sep=separator)
    picks = pd.read_csv(pickfile, sep=separator)

    evdep_label, evdep_scale = resolve_event_depth_column(events)

    if "event_idx" in picks.columns:
        picks = picks[picks["event_idx"] >= 0]
    if "event_index" in events.columns:
        events = events[events["event_index"] >= 0]

    picks = picks.loc[:, ~picks.columns.duplicated()]

    # ---- remap event IDs ----
    if mapping_evid:
        save_cleaned_data = True
        if evid_label_mapped is None:
            evid_label_mapped = evid_label + "_mapped"

        evid_map = {
            evid: i + 1 for i, evid in enumerate(events[evid_label].unique())
        }
        events[evid_label_mapped] = events[evid_label].map(evid_map)
        picks[evid_label_mapped] = picks[evid_label].map(evid_map)
    else:
        evid_label_mapped = evid_label

    if save_cleaned_data:
        cleaned_eventfile = cleaned_eventfile or eventfile.replace(".csv", "_cleaned.csv")
        cleaned_pickfile = cleaned_pickfile or pickfile.replace(".csv", "_cleaned.csv")
        events.to_csv(cleaned_eventfile, index=False)
        picks.to_csv(cleaned_pickfile, index=False)

    picks_by_event = picks.groupby(evid_label_mapped)

    with open(outfile, "w") as f:
        for _, ev in tqdm(events.iterrows(), total=len(events),desc='Converting GaMMA to HYPOINVERSE phase file'):
            evid = ev[evid_label_mapped]
            if evid not in picks_by_event.groups:
                continue

            header = format_event_line_hypoinverse(
                ev, evid, evdep_label, evdep_scale
            )

            phase_lines = []
            has_p = has_s = False

            for _, pk in picks_by_event.get_group(evid).iterrows():
                line = format_pick_line_gamma2hypoinverse(pk, default_component)
                if line is None:
                    continue
                phase_lines.append(line + "\n")
                if pk["type"].upper() == "P":
                    has_p = True
                elif pk["type"].upper() == "S":
                    has_s = True

            if qc and not (has_p and has_s):
                continue

            f.write(header+"\n")
            f.writelines(phase_lines)
            f.write("\n")

    if verbose:
        print(f"[DONE] HYPOINVERSE phase file written → {outfile}")
######### end of functions for GaMMA to Hypoinverse phase file conversion

######### functions for GaMMA to HypoDD phase file conversion
# -------------------------------
# Format one HypoDD pick line
# -------------------------------
def format_pick_line_gamma2hypodd(pick,origin_time,combine_net_sta=True):
    """
    Convert a single GAMMA pick to HypoDD phase format.
    Station name = NET + STA. Pick line is based on the formatting in convert_legacyphase2hypodd().
    ==========PARAMETERS=============
    pick: a single pick record (Pandas Series) from GAMMA output.
    origin_time: origin time of the event (datetime object).
    combine_net_sta: if True, station name = NET.STA. Default is True [recommended].
    ==========RETURN=============
    formatted pick line in HypoDD format.
    """

    # Parse pick ID
    try:
        net, sta, _, _ = pick["id"].split(".")
    except ValueError:
        warnings.warn(f"Invalid pick id format: {pick['id']}")
        return None
    if combine_net_sta:
        station = f"{net}.{sta}"
    else:
        station = sta
    phase = pick["type"].upper()

    if phase not in ("P", "S"):
        return None

    # Pick time
    t = datetime.strptime(pick["timestamp"], "%Y-%m-%dT%H:%M:%S.%f")

    picktime = UTCDateTime(t)
    phase_dt = picktime - UTCDateTime(origin_time)

    weight = 1.0  # default weight

    pickstr="{:<7}{:s}{:6.3f}  {:2f}   {:s}\n".format(
        station, ' ' * 6,
        phase_dt,
        weight,
        phase
    )

    return pickstr
# ---------------------------------------
# Main GAMMA → HypoDD conversion function
# ---------------------------------------
def convert_gamma2hypodd(eventfile,pickfile,outfile="hypodd.phase",separator=",",
                         combine_net_sta=True,qc=False,verbose=False,evid_label="event_id",
                         mapping_evid=True,evid_label_mapped=None,save_cleaned_data=False, 
                         cleaned_eventfile=None, cleaned_pickfile=None):
    """
    Convert GAMMA event + pick CSVs to HypoDD phase file.
    ==========PARAMETERS=============
    eventfile: input event CSV file from GAMMA.
    pickfile: input pick CSV file from GAMMA.
    outfile: output HypoDD phase file name.
    separator: separator used in the input CSV files. Default is comma (,).
    qc: if True, only events with both P and S picks are written. Default is False.
    combine_net_sta: if True, station name = NET.STA. Default is True [recommended].
    verbose: if True, print progress messages. Default is False.
    evid_label: column name for event ID in the input event and pick files. Default is "event_id".
    mapping_evid: if True, remap event IDs to integers starting from 1 for HypoDD compatibility. Default is True. 
                If False, original event IDs are used (may cause issues if IDs are not integers or not starting from 1).
    evid_label_mapped: evid label after remapping. This is only effective when mapping_evid is True.
    save_cleaned_data: if True, save cleaned event and pick data after removing invalid entries. Default is False.
                Set to True automatically if mapping_evid is True.
    cleaned_eventfile: output file name for cleaned event data. Default is None (auto-generate).
    cleaned_pickfile: output file name for cleaned pick data. Default is None (auto-generate).
    ==========RETURN=============
    None
    """

    # Read input
    events = pd.read_csv(eventfile, sep=separator)
    picks = pd.read_csv(pickfile, sep=separator)
    evdep_label,evdep_scale=resolve_event_depth_column(events)
    #remove rows with event_idx negative. These are invalid picks/events.
    if "event_idx" in picks.columns:
        picks = picks[picks["event_idx"] >= 0]

    if "event_index" in events.columns:
        events = events[events["event_index"] >= 0]

    # Clean columns
    picks = picks.loc[:, ~picks.columns.duplicated()]
    if mapping_evid:
        save_cleaned_data = True
        # remap event IDs to integers to be compatible with HypoDD
        unique_evids = events[evid_label].unique()
        evid_map = {evid: idx + 1 for idx, evid in enumerate(unique_evids)} # start from 1 for HypoDD compatibility

        # add a new column for mapped event IDs
        if evid_label_mapped is None:
            evid_label_mapped = evid_label + "_mapped"
        #apply mapping
        events[evid_label_mapped] = events[evid_label].map(evid_map)
        picks[evid_label_mapped] = picks[evid_label].map(evid_map)
    else:
        evid_label_mapped = evid_label

    if save_cleaned_data:
        # Save cleaned event and pick data
        if cleaned_eventfile is None:
            cleaned_eventfile = os.path.splitext(eventfile)[0] + "_cleaned.csv"
        if cleaned_pickfile is None:
            cleaned_pickfile = os.path.splitext(pickfile)[0] + "_cleaned.csv"
        events.to_csv(cleaned_eventfile, index=False)
        picks.to_csv(cleaned_pickfile, index=False)
        print(f"Cleaned event data saved to {cleaned_eventfile}")
        print(f"Cleaned pick data saved to {cleaned_pickfile}")
    # Group picks by event
    picks_by_event = picks.groupby(evid_label_mapped)

    with open(outfile, "w") as f:
        #loop over events
        for event in tqdm(events.to_dict(orient="records"), desc="Converting GAMMA to HypoDD phase file"):

            evid = event[evid_label_mapped]
            if evid not in picks_by_event.groups:
                continue

            # Event header line
            ot = datetime.strptime(
                event["time"], "%Y-%m-%dT%H:%M:%S.%f"
            )
            lat = event["latitude"]
            lon = event["longitude"]
            dep = event[evdep_label]*evdep_scale  # convert m to km 
            mag = event.get("magnitude", 0.0)
            f.write(
                "# {:4d} {:02d} {:02d}  {:02d} {:02d} {:06.3f}  "
                "{:7.4f} {:9.4f}  {:6.2f} {:4.2f} 0.00  0.00  0.00 {:9d}\n".format(
                    ot.year,
                    ot.month,
                    ot.day,
                    ot.hour,
                    ot.minute,
                    ot.second + ot.microsecond / 1e6,
                    lat,
                    lon,
                    dep,
                    mag,
                    evid,
                )
            )

            # Pick lines
            lines = []
            has_p = has_s = False

            for _, pick in picks_by_event.get_group(evid).iterrows():
                line = format_pick_line_gamma2hypodd(pick, ot,combine_net_sta=combine_net_sta)
                if line:
                    lines.append(line)
                    if pick["type"].upper() == "P":
                        has_p = True
                    elif pick["type"].upper() == "S":
                        has_s = True

            # Optional QC: require both P and S
            if not qc or (has_p and has_s):
                f.writelines(lines)
                if verbose:
                    print(f"Event {evid}: {len(lines)} picks written")
            elif verbose:
                print(f"Event {evid} skipped (QC failed)")
##### end of functions for GaMMA to HypoDD phase file conversion
def convert_hypoinverse2hypodd(phase_file_in,phase_file_out, verbose=False):
    """
    Reformat input hypoinverse phase file for HypoDD run.

    Parameters
    ----------
    phase_file_in : str
        Input phase file path.
    phase_file_out : str
        Output phase file path.
    verbose : bool
        If True, print progress messages. Default is False.    

    Returns
    -------
    None
    """
    #this is a untested function. Use with caution.
    # print warning message here
    if verbose:
        print("[WARNING] convert_hypoinverse2hypodd() is untested. Use with caution.")
    fout_dir = os.path.split(phase_file_out)[0]
    os.makedirs(fout_dir, exist_ok=True)


    evid_list = []

    # --------------------------------------------------
    # Read original phase file
    # --------------------------------------------------
    with open(phase_file_in) as f:
        lines = f.readlines()

    with open(phase_file_out, 'w') as fout:
        for line in lines:
            codes = line.strip().split(',')

            # ------------------------------------------
            # Event header line
            # ------------------------------------------
            if len(codes[0]) >= 14:
                ot = UTCDateTime(codes[0])
                lat, lon, dep, mag = [float(code) for code in codes[1:5]]
                evid = int(codes[-1])

                write_event = True
                evid_list.append(evid)

                # Format time
                date = '{:4} {:2} {:2}'.format(ot.year, ot.month, ot.day)
                time = '{:2} {:2} {:5.2f}'.format(
                    ot.hour, ot.minute, ot.second + ot.microsecond / 1e6
                )

                # Format location
                loc = '{:7.4f} {:9.4f}  {:6.2f} {:4.2f}'.format(
                    lat, lon, dep, mag
                )

                fout.write(
                    '# {} {}  {}  0.00  0.00  0.00  {:>9}\n'.format(
                        date, time, loc, evid
                    )
                )

            # ------------------------------------------
            # Station pick lines
            # ------------------------------------------
            else:
                if not write_event:
                    continue

                netsta = codes[0].strip()
                wp, ws = 1.0, 1.0

                # P pick
                if codes[1] != '-1':
                    tp = UTCDateTime(codes[1])
                    ttp = tp - ot
                    fout.write(
                        '{:<7}{}{:6.3f}  {:6.3f}   P\n'.format(
                            netsta, ' ' * 6, ttp, wp
                        )
                    )

                # S pick
                if codes[2] != '-1':
                    ts = UTCDateTime(codes[2])
                    tts = ts - ot
                    fout.write(
                        '{:<7}{}{:6.3f}  {:6.3f}   S\n'.format(
                            netsta, ' ' * 6, tts, ws
                        )
                    )

    # --------------------------------------------------
    print(f"[INFO] Wrote single phase file: {phase_file_out}")
    print(f"[INFO] Number of events: {len(evid_list)}")
#
def convert_legacyphase2hypodd(phase_file_in,phase_file_out, dep_corr=5,
                       time_range=None,lat_range=None,lon_range=None):
    """
    Reformat input phase file for HypoDD run. This is a legacy function modified from hypo-interface-py.
    I kept it here for backward compatibility.

    Applies:
      - time window filtering
      - lat/lon window filtering
      - depth correction

    Parameters
    ----------
    phase_file_in : str
        Input phase file path.
    phase_file_out : str
        Output phase file path.
    dep_corr : float
        Depth correction to avoid air quakes. Default: 5 km.
    time_range : str
        Time range filter in 'YYYY-MM-DDTHH:MM:SS-YYYY-MM-DDTHH:MM:SS' format. Default: None (no filter).
    lat_range : tuple
        Latitude range filter as (lat_min, lat_max). Default: None (no filter).
    lon_range : tuple
        Longitude range filter as (lon_min, lon_max). Default: None (no filter).    

    Returns
    -------
    phase_file_out : str
        Output phase file path.
    evid_list : np.ndarray
        Array of event IDs included in the output phase file.
    """
    fout_dir = os.path.split(phase_file_out)[0]
    os.makedirs(fout_dir, exist_ok=True)

    #subset filters
    if time_range is not None:
        ot_min, ot_max = [UTCDateTime(date) for date in time_range.split('-')]
        subset_time = True
    if lat_range is not None:
        lat_min, lat_max = lat_range
    else:
        lat_min, lat_max = -90.0, 90.0
    if lon_range is not None:
        lon_min, lon_max = lon_range
    else:
        lon_min, lon_max = -180.0, 180.0

    evid_list = []

    # --------------------------------------------------
    # Read original phase file
    # --------------------------------------------------
    with open(phase_file_in) as f:
        lines = f.readlines()

    with open(phase_file_out, 'w') as fout:
        for line in lines:
            codes = line.strip().split(',')

            # ------------------------------------------
            # Event header line
            # ------------------------------------------
            if len(codes[0]) >= 14:
                ot = UTCDateTime(codes[0])
                lat, lon, dep, mag = [float(code) for code in codes[1:5]]
                dep += dep_corr
                evid = int(codes[-1])

                # Filters
                if subset_time:
                    if not (ot_min < ot < ot_max):
                        write_event = False
                        continue
                if not (lat_min <= lat <= lat_max and lon_min <= lon <= lon_max):
                    write_event = False
                    continue

                write_event = True
                evid_list.append(evid)

                # Format time
                date = '{:4} {:2} {:2}'.format(ot.year, ot.month, ot.day)
                time = '{:2} {:2} {:5.2f}'.format(
                    ot.hour, ot.minute, ot.second + ot.microsecond / 1e6
                )

                # Format location
                loc = '{:7.4f} {:9.4f}  {:6.2f} {:4.2f}'.format(
                    lat, lon, dep, mag
                )

                fout.write(
                    '# {} {}  {}  0.00  0.00  0.00  {:>9}\n'.format(
                        date, time, loc, evid
                    )
                )

            # ------------------------------------------
            # Station pick lines
            # ------------------------------------------
            else:
                if not write_event:
                    continue

                netsta = codes[0].strip()
                wp, ws = 1.0, 1.0

                # P pick
                if codes[1] != '-1':
                    tp = UTCDateTime(codes[1])
                    ttp = tp - ot
                    fout.write(
                        '{:<7}{}{:6.3f}  {:6.3f}   P\n'.format(
                            netsta, ' ' * 6, ttp, wp
                        )
                    )

                # S pick
                if codes[2] != '-1':
                    ts = UTCDateTime(codes[2])
                    tts = ts - ot
                    fout.write(
                        '{:<7}{}{:6.3f}  {:6.3f}   S\n'.format(
                            netsta, ' ' * 6, tts, ws
                        )
                    )

    # --------------------------------------------------
    print(f"[INFO] Wrote single phase file: {phase_file_out}")
    print(f"[INFO] Number of events: {len(evid_list)}")

    return phase_file_out, np.array(evid_list)