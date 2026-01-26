#This module contains core functions for running HypoInverse interface.
#Import needed packages first.
import pandas as pd
import os,glob,shutil
import numpy as np
from hypoxpy import utils
import subprocess
from pathlib import Path
from obspy import UTCDateTime
from tqdm import tqdm
#
#
def generate_parfile(config,pardir='input',outdir='output',template=None,magline=None):
    """
    Generate parameter files for running hypoinverse, with given HypoInvConfig object.

    ========PARAMETERS=========
    config: A HypoInvConfig object containing all controling parameters.
    pardir: output directory to store the parameter files. Default: input.
    outdir: hypoinverse relocation output directory. Default: output.
    template: provide template or use the built-in template file for parameters.
    magline: the line in the parameter file to compute magnitude. Defalt is None.

    ========RETURNS===========
    filelist: list of the parameter files, including path.
    """
    if not os.path.isdir(pardir):
        os.makedirs(pardir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    filelist=[]
    for ztr in config.ztrlist:
        # set control file
        fhyp = os.path.join(pardir,'%s-%s.hypoinv.par'%(config.namebase, ztr))
        filelist.append(fhyp)

        #save parameters by modifying the template parameters.
        fout=open(fhyp,'w')
        if template is None:
            lines=utils.load_template(config.template_parfile)
        else:
            lines=utils.load_template(template)
        for line in lines:
            # loc params
            if line[0:3]=='ZTR': line = "ZTR %s F \n"%ztr
            if line[0:3]=='RMS': line = "RMS %s \n"%config.rms_weight
            if line[0:3]=='DI1': line = "DI1 %s \n"%config.dist_initial
            if line[0:3]=='DIS': line = "DIS %s \n"%config.dist_weight
            if line[0:3]=='WET': line = "WET %s \n"%config.weight_code
            # i/o paths
            if line[0:3]=='STA': line = "STA '%s' \n"%config.station_file
            if line[0:3]=='PHS': line = "PHS '%s' \n"%config.phase_file
            if line[0:5]=='CRE 1': line = "CRE 1 '%s' %s T \n"%(config.pmodel, config.ref_ele)
            if line[0:5]=='CRE 2': line = "CRE 2 '%s' %s T \n"%(config.smodel, config.ref_ele)
            if line[0:3]=='POS': line = "POS %s \n"%(config.poisson)
            if line[0:3]=='SUM': line = "SUM '%s/%s-hypoinv_%s.sum' \n"%(outdir,config.namebase, ztr)
            if line[0:3]=='MIN': line = "MIN %d \n"%(config.min_nsta)
            if line[0:3]=='PRT': 
                line = "PRT '%s/%s-%s.ptr' \n"%(outdir,config.namebase, ztr) if config.get_prt else ''
            if line[0:3]=='ARC': 
                line = "ARC '%s/%s-%s.arc' \n"%(outdir,config.namebase, ztr) if config.get_arc else ''
            #if line[0:3]=='H71': line = "H71 1 1 3" #use hypoinverse summary output format (first integer)
            if line[0:3]=='STO': 
                continue
            else:
                fout.write(line)
        #get magnitude
        if magline is not None:
            # line="MAG 1 T 1 1\n"
            line = "MFL '%s/%s-%s.mag' \n"%(outdir,config.namebase, ztr)
            fout.write(line)
            if magline[-1] != '\n':
                magline += '\n'
            fout.write(magline)
        # Always stop HypoInverse
        fout.write("STO\n")
        fout.close()
    #
    return filelist

def summary_to_catalog(df, lat_code, lon_code, mag_dict=None, evid_label="event_id"):
    """
    Convert hypoinverse summary dataframe to a catalog dataframe.
    """
    rows = []

    for _, row in df.iterrows():
        line = row["line"]
        evid = row["evid"]

        # origin time → UTCDateTime ISO format
        date, hrmn, sec = line.split()[0:3]
        sec = sec.zfill(5)
        ot = UTCDateTime(f"{date}{hrmn}{sec}")
        time_str = ot.strftime("%Y-%m-%dT%H:%M:%S.%f")
        # latitude
        lat_deg = float(line[20:22])
        lat_min = float(line[23:28])
        lat = lat_deg + lat_min / 60 if lat_code == "N" else -lat_deg - lat_min / 60

        # longitude
        lon_deg = float(line[29:32])
        lon_min = float(line[33:38])
        lon = lon_deg + lon_min / 60 if lon_code == "E" else -lon_deg - lon_min / 60

        # depth & magnitude
        dep = float(line[38:45])
        mag = float(line[45:52])

        if mag_dict is not None and str(evid) in mag_dict:
            mag = mag_dict[str(evid)]

        rows.append(
            {
                "time": time_str,
                "latitude": lat,
                "longitude": lon,
                "depth_km": dep,
                "magnitude": mag,
                evid_label: evid,
            }
        )

    return pd.DataFrame(rows)

#####
#####
class HypoInvConfig(object):
    """
    Container class to store key configuration parameters for running HypoInverse.
    ======== PARAMETERS TO INITIATE =========
    binpath: path to the hypoinverse binary. Default None (use system PATH).
    indir: input directory for hypoinverse files. Default 'input'.
    outdir: output directory for hypoinverse files. Default 'output'.
    phase_file: input phase file in hypoinverse format.
    station_file: input station file in hypoinverse format.
    pmodel: P-wave velocity model file in CRE format.
    smodel: S-wave velocity model file in CRE format.
    poisson: Poisson's ratio to compute S-wave velocity model from P-wave velocity model.
    namebase: tag for the relocation run. Default 'hyp'.
    hypoinv_bin: path to the hypoinverse binary. Default 'hyp1.40'.
    get_prt: whether to output the .prt file. Default False.
    get_arc: whether to output the .arc file. Default False.
    lat_code: latitude code for hypoinverse format. Default 'N'.
    lon_code: longitude code for hypoinverse format. Default 'W'.
    ref_ele: reference elevation for the velocity model. Default 0.0.
    grd_ele: ground elevation for the stations. Default 0.0.
    ztrlist: list of initial depths for the relocation run. Default np.arange(0
    ,20,1).
    rms_weight: RMS weighting parameters. Default '4 0.3 1 3'.
    dist_initial: initial distance weighting parameters. Default '1 50 1 2'.
    dist_weight: distance weighting parameters. Default '4 20 1 3'.
    weight_code: weight code parameters. Default '1 0.6 0.3 0.2'.
    min_nsta: minimum number of stations for a valid event. Default 4.
    =============================
    """
    def __init__(self,binpath=None,indir='input',outdir='output',phase_file=None,station_file=None,pmodel=None,
                 smodel=None,poisson=1.73,namebase='hyp',hypoinv_bin='hyp1.40',get_prt=False,get_arc=False,
               lat_code='N',lon_code='W',ref_ele=0.0,grd_ele=0.0,ztrlist = np.arange(0,20,1),
               rms_weight='4 0.3 1 3',dist_initial = '1 50 1 2',dist_weight = '4 20 1 3',
               weight_code='1 0.6 0.3 0.2',min_nsta=4,template_parfile=None):
        self.type="HypoInvConfig object"
        if binpath is None:
            HYPBIN = shutil.which("hyp1.40")
            if HYPBIN is None:
                raise RuntimeError("hyp1.40 not found in PATH")
            binpath = HYPBIN # default path to hypoinverse binary, assuming it is in the system PATH
        self.binpath = os.path.join(binpath,'hyp1.40') #hard-coded hypoinverse program. 
        self.indir = indir
        self.outdir = outdir
        self.namebase = namebase
        self.hypoinv_bin=hypoinv_bin
        # i/o paths
        self.station_file = station_file
        self.phase_file = phase_file
        self.get_prt = get_prt
        self.get_arc = get_arc
        # geo ref
        self.lat_code = lat_code
        self.lon_code = lon_code
        self.ref_ele = ref_ele # ref ele for CRE mod (max sta ele)
        self.grd_ele = grd_ele # typical station elevation
        # loc params
        self.ztrlist = ztrlist #initial depth for the relocation run.
        self.p_weight = 0 # weight code index
        self.s_weight = 1
        self.min_nsta=min_nsta
        self.rms_weight = rms_weight
        self.dist_initial = dist_initial
        self.dist_weight = dist_weight
        self.weight_code = weight_code
        if template_parfile is not None:
            self.template_parfile = template_parfile
        else:
            self.template_parfile = utils.get_template_list('hypoinv')[1]
            print(f"[INFO] Using built-in template parameter file: {self.template_parfile}")
        self.pmodel = pmodel #'input/velo_p_eg.cre'
        self.smodel = smodel #[None, 'input/velo_s_eg.cre'][1]
        self.poisson = poisson #1.73 # provide smod or pos

        #
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if not os.path.exists(self.indir):
            os.makedirs(self.indir)
    def __str__(self):
        lines = [
            f"{self.type}",
            "-" * len(self.type),

            # Executable / paths
            f"HYPOINVERSE binary : {self.binpath}",
            f"Input directory   : {self.indir}",
            f"Output directory  : {self.outdir}",
            f"Name base         : {self.namebase}",

            # Input files
            f"Phase file        : {self.phase_file}",
            f"Station file     : {self.station_file}",
            f"Template parfile  : {self.template_parfile}",
            f"P velocity model : {self.pmodel}",
            f"S velocity model : {self.smodel}",
            f"Poisson ratio    : {self.poisson}",

            # Geographic reference
            f"Latitude code    : {self.lat_code}",
            f"Longitude code   : {self.lon_code}",
            f"Reference elev.  : {self.ref_ele}",
            f"Grid elev.       : {self.grd_ele}",

            # Location / weighting parameters
            f"Initial depths   : {list(self.ztrlist)}",
            f"Min # stations   : {self.min_nsta}",
            f"RMS weight       : {self.rms_weight}",
            f"Initial dist.    : {self.dist_initial}",
            f"Dist. weight     : {self.dist_weight}",
            f"Weight code      : {self.weight_code}",

            # Output controls
            f"Write PRT file   : {self.get_prt}",
            f"Write ARC file   : {self.get_arc}",
        ]

        return "\n".join(lines)
    def __repr__(self):
        return (
            f"HypoInvConfig("
            f"phase_file={self.phase_file}, "
            f"station_file={self.station_file}, "
            f"template_parfile={self.template_parfile}, "
            f"pmodel={self.pmodel}, "
            f"smodel={self.smodel}, "
            f"outdir={self.outdir})"
        )

    #-------------------------------------------------
    # core function to run hypoinverse
    #-------------------------------------------------
    def run(self,parfilelist):
        """
        Run hypoinverse for a list of parameter files.
        ======== PARAMETERS ==========
        parfilelist: list of parameter files for hypoinverse.
        cleanup: whether to remove intermediate files after running hypoinverse. Default True. 
        ========
        """
        for fhyp in tqdm(parfilelist,"Running %s by depth grids"%(self.binpath)):
            log_file = Path(fhyp).with_suffix(".log")
            with open(log_file, "w") as log:
                p = subprocess.Popen(
                    [self.binpath],
                    stdin=subprocess.PIPE,
                    stdout=log,
                    stderr=subprocess.STDOUT,  # capture errors too
                    encoding="utf-8"
                )

                s = f"@{fhyp}\n"
                p.communicate(s)

    ###
    def finalize(self,mag_dict=None,evid_label="event_id",cleanup=False,out_good=None,out_bad=None):
                #-------------------------------------------------
        cat_good,cat_bad = self._merge_summary(mag_dict=mag_dict,evid_label=evid_label,
                            out_good=out_good,out_bad=out_bad)
        # remove intermidiate files
        if cleanup:
            files_to_remove = []

            files_to_remove.extend(glob.glob('fort.*'))
            files_to_remove.extend(glob.glob(f'{self.indir}/{self.namebase}-*.hypoinv.par'))
            files_to_remove.extend(glob.glob(f'{self.indir}/{self.namebase}-*.hypoinv.log'))
            files_to_remove.extend(glob.glob(f'{self.outdir}/{self.namebase}-hypoinv_*.sum'))

            for fname in files_to_remove:
                try:
                    os.remove(fname)
                except FileNotFoundError:
                    pass
                except Exception as e:
                    print(f"[WARN] Could not remove {fname}: {e}")
        #
        return cat_good,cat_bad
    ####
    def _merge_summary(self,mag_dict=None,evid_label="event_id",
                       out_good=None,out_bad=None):
        """
        Extract earthquake parameters based on final quality after merging all summary files.
        Uses pandas DataFrame internally.

        PARAMETERS
        ----------
        mag_dict : dict or str, optional
            Magnitude dictionary {event_id: magnitude} or CSV catalog path.
        evid_label : str
            Event ID column name if mag_dict is a CSV.
        """

        filelist = glob.glob(f"{self.outdir}/{self.namebase}-hypoinv_*.sum")

        if out_good is None:
            out_good = f"{self.outdir}/{self.namebase}_good.csv"
        if out_bad is None:
            out_bad  = f"{self.outdir}/{self.namebase}_bad.csv"

        # --------------------------------------------------
        # Load magnitude dictionary if provided
        # --------------------------------------------------
        mag_dict_use = None
        if mag_dict is not None:
            if isinstance(mag_dict, str):
                events = pd.read_csv(mag_dict)
                mag_dict_use = dict(
                    zip(events[evid_label].astype(str), events["magnitude"])
                )
            elif isinstance(mag_dict, dict):
                mag_dict_use = mag_dict
            else:
                raise ValueError("mag_dict must be a dict or CSV filename.")
        # --------------------------------------------------
        # Read and concatenate all summary files
        # --------------------------------------------------
        records = []

        for fsum in filelist:
            with open(fsum) as f:
                for line in f:
                    evid = line.split()[-1]

                    codes = line.split()
                    is_loc = 0 if ("-" in codes or "#" in codes) else 1

                    records.append(
                        {
                            "evid": evid,
                            "line": line,
                            "is_loc": is_loc,
                            "qua": line[80:81],
                            "azm": float(line[56:59]),
                            "npha": 1.0 / float(line[52:55]),
                            "rms": float(line[64:69]),
                        }
                    )

        df = pd.DataFrame(records)

        # --------------------------------------------------
        # Sort using Hypoinverse priority logic
        # --------------------------------------------------
        df = df.sort_values(
            by=["evid", "qua", "azm", "npha", "rms"],
            ascending=[True, True, True, True, True],
        )

        # --------------------------------------------------
        # Pick best location per event
        # --------------------------------------------------
        good_rows = []
        bad_rows = []

        for evid, g in df.groupby("evid"):
            g_loc = g[g["is_loc"] == 1]

            if len(g_loc) > 0:
                good_rows.append(g_loc.iloc[0])
            else:
                bad_rows.append(g.iloc[0])

        df_good = pd.DataFrame(good_rows)
        df_bad  = pd.DataFrame(bad_rows)

        cat_good = summary_to_catalog(df_good,self.lat_code,self.lon_code,mag_dict=mag_dict_use,evid_label=evid_label)
        cat_bad  = summary_to_catalog(df_bad,self.lat_code,self.lon_code,mag_dict=mag_dict_use,evid_label=evid_label)

        # --------------------------------------------------
        # Save with headers
        # --------------------------------------------------
        cat_good.to_csv(out_good, index=False)
        cat_bad.to_csv(out_bad, index=False)

        print(
            f"Earthquakes are saved in:\n"
            f"  {out_good} (good)\n"
            f"  {out_bad} (bad)"
        )

        return cat_good,cat_bad




        

