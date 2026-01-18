#This file contains utility functions for the package.
import pandas as pd
import os,json,shutil,glob,time
from importlib import resources as impresources
from hypoxpy import templates
import warnings
from tqdm import tqdm
from datetime import datetime
#####

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
    with open(outfile, "w") as fout:
        for _, row in indata.iterrows():

            net = str(row["network"]).strip()
            sta = str(row["station"]).strip()

            if combine_net_sta:
                sta_name = f"{net}.{sta}"
            else:
                sta_name = sta

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
def get_template_list(basename):
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
    templatelist=get_filelist(templatedir,pattern='%s_template_'%(basename))
    templatetail=[]
    for tf in templatelist:
        ftail=os.path.split(tf)[1]
        templatetail.append(ftail)
    #
    return(templatetail)
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
# functions for GaMMA to Hypoinverse phase file conversion
def format_event_line_gamma2hypoinverse(event): #gamma to hypoinverse
    """
    Formats a single event line for Hypoinverse phase file.
    Keeps event ID and magnitude.
    """
    event_time = datetime.strptime(event["time"], "%Y-%m-%dT%H:%M:%S.%f")
    lat_deg = int(event["latitude"])
    lon_deg = int(event["longitude"])
    lat_min = (abs(event["latitude"]) - abs(lat_deg)) * 60 * 100
    lon_min = (abs(event["longitude"]) - abs(lon_deg)) * 60 * 100
    south = "S" if event["latitude"] < 0 else " "
    east = "E" if event["longitude"] >= 0 else " "
    depth = event["depth(m)"] / 1000
    mag = event["magnitude"]

    line = (f"{event_time.strftime('%Y%m%d%H%M%S%f')[:-4]}"
            f"{abs(lat_deg):2d}{south}{abs(lat_min):4.0f}"
            f"{abs(lon_deg):3d}{east}{abs(lon_min):4.0f}"
            f"{depth:5.0f}{mag:3.1f}{' ':97}{event['event_index']:10}\n")
    return line

def format_pick_line_gamma2hypoinverse(pick, default_component=None): #gamma to hypoinverse
    """
    Formats a single pick line for Hypoinverse phase file.
    Includes network, station, channel, component info.
    """
    try:
        network_code, station_code, comp_code, channel_code = pick['id'].split('.')
    except ValueError:
        warnings.warn(f"Pick ID format invalid: {pick['id']}")
        return None

    if default_component and not comp_code:
        comp_code = default_component

    phase = pick['type'].upper()
    weight = min(max(int((1 - pick['prob']) / (1 - 0.3) * 4) - 1, 0), 3)
    pick_time = datetime.strptime(pick["timestamp"], "%Y-%m-%dT%H:%M:%S.%f")
    pick_min = pick_time.strftime("%Y%m%d%H%M")
    pick_sec = pick_time.strftime("%S%f")[:-4]
    amp = pick.get('phase_amplitude', 0)

    templine = f"{station_code:<5}{network_code:<2}  {channel_code:<2}{comp_code:<1}"

    if phase == 'P':
        return f"{templine:<13} P {weight:<1d}{pick_min} {pick_sec}                    {amp:<7}"
    elif phase == 'S':
        return f"{templine:<13}   4{pick_min} {'':<12}{pick_sec} S {weight:<1d}    {amp:7f}"
    else:
        warnings.warn(f"Unknown phase type: {phase}")
        return None

def convert_gamma2hypoinverse(eventfile, pickfile, outfile='phase_output.phs',
                  qc=False, separator=',', default_component=None, verbose=False):
    """
    Converts GaMMA event and pick files to Hypoinverse phase file format.
    """
    # Read files
    events = pd.read_csv(eventfile, sep=separator)
    picks = pd.read_csv(pickfile, sep=separator)
    
    # Remove duplicate columns
    picks = picks.loc[:, ~picks.columns.duplicated()]
    # Consistent column name for event index
    picks.rename(columns={"event_idx": "event_index"}, inplace=True)
    # Group picks by event index for faster access
    picks_eventwise = picks.groupby("event_index").groups

    with open(outfile, 'w') as f:
        for _, event in tqdm(events.iterrows(), total=len(events)):
            lines = [format_event_line_gamma2hypoinverse(event)]
            has_p = has_s = False

            # Get picks for this event
            for idx in picks_eventwise.get(event["event_index"], []):
                pick_line = format_pick_line_gamma2hypoinverse(picks.iloc[idx], default_component)
                if pick_line:
                    lines.append(pick_line + "\n")
                    phase_type = picks.iloc[idx]['type'].upper()
                    if phase_type == 'P': has_p = True
                    elif phase_type == 'S': has_s = True

            # Write lines if QC passes
            if not qc or (has_p and has_s):
                f.writelines(lines)
                f.write("\n")
                if verbose:
                    print(f"Event {event['event_index']} saved with {len(lines)-1} picks.")
            elif verbose:
                print(f"Event {event['event_index']} skipped (QC failed).")
######### end of functions for GaMMA to Hypoinverse phase file conversion


# -------------------------------
# Format one HypoDD pick line
# -------------------------------
def format_pick_line_gamma2hypodd(pick):
    """
    Convert a single GAMMA pick to HypoDD phase format.
    Station name = NET + STA
    """

    # Parse pick ID
    try:
        net, sta, _, _ = pick["id"].split(".")
    except ValueError:
        warnings.warn(f"Invalid pick id format: {pick['id']}")
        return None

    station = f"{net}.{sta}"
    phase = pick["type"].upper()

    if phase not in ("P", "S"):
        return None

    # Pick time
    t = datetime.strptime(pick["timestamp"], "%Y-%m-%dT%H:%M:%S.%f")

    # HypoDD time fields
    yyyy = t.year
    mm = t.month
    dd = t.day
    hh = t.hour
    mi = t.minute
    sec = t.second + t.microsecond / 1e6

    # Weight from GAMMA probability
    # prob ~ [0.3, 1.0] → weight [0,4]
    prob = pick.get("prob", 1.0)
    weight = min(max(int((1 - prob) / (1 - 0.3) * 4), 0), 4)

    evid = int(pick["event_index"])

    return (
        f"{station:<7s} {phase} "
        f"{yyyy:4d} {mm:02d} {dd:02d} "
        f"{hh:02d} {mi:02d} {sec:06.3f} "
        f"{weight:d} {evid:d}"
    )


# ---------------------------------------
# Main GAMMA → HypoDD conversion function
# ---------------------------------------
def convert_gamma2hypodd(eventfile,
                            pickfile,
                            outfile="hypodd.phase",
                            separator=",",
                            qc=False,
                            verbose=False):
    """
    Convert GAMMA event + pick CSVs to HypoDD phase file.
    """

    # Read input
    events = pd.read_csv(eventfile, sep=separator)
    picks = pd.read_csv(pickfile, sep=separator)

    # Clean columns
    picks = picks.loc[:, ~picks.columns.duplicated()]
    picks.rename(columns={"event_idx": "event_index"}, inplace=True)

    # Group picks by event
    picks_by_event = picks.groupby("event_index")

    with open(outfile, "w") as f:
        for _, event in tqdm(events.iterrows(), total=len(events)):

            evid = event["event_index"]
            if evid not in picks_by_event.groups:
                continue

            lines = []
            has_p = has_s = False

            for _, pick in picks_by_event.get_group(evid).iterrows():
                line = format_pick_line_gamma2hypodd(pick)
                if line:
                    lines.append(line + "\n")
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
