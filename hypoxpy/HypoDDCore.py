#this package contains functions to run hypodd earthquake relocation codes.
#import needed packages.
""" 
    Download hypoDD at https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Or use the copy included in this python interface package.
"""
from pathlib import Path
import os,shutil, glob, subprocess, warnings
import numpy as np
import pandas as pd
from obspy import UTCDateTime
from hypoxpy import utils
warnings.filterwarnings("ignore")

#
#
def run_ph2dt(config,verbose=True):
    """
    Call ph2dt once for the entire dataset (no grids).

    config: HypoDDConfig object
    """

    if verbose: print("running ph2dt...")

    # --------------------------------------------------
    # Prepare input phase file
    # --------------------------------------------------
    # Expect a single reformatted phase file:
    #   input/{namebase}.pha
    if not os.path.exists(config.phase_file):
        raise FileNotFoundError(f"Phase file not found: {config.phase_file}")

    # costumize ph2dt input file
    ph2dt_inp = '%s/ph2dt_%s.inp'%(config.indir, config.namebase)
    fout = open(ph2dt_inp,'w')
    f=open(config.ph2dt_inp_template); lines=f.readlines(); f.close()
    for line in lines:
        #skip lines starting with "*"
        if line[0] == "*":
            pass
        else:
            if 'input/phase.dat' in line: line = '%s \n'%(config.phase_file)
            if 'input/station.dat' in line: line = '%s \n'%(config.station_file)
        fout.write(line)
    fout.close()
    # --------------------------------------------------
    # Run ph2dt
    # --------------------------------------------------
    out_file = f'{config.outdir}/{config.namebase}.ph2dt'
    os.makedirs(config.outdir, exist_ok=True)

    with open(out_file, 'w') as f:
        subprocess.run(
            [f'{config.bin_ph2dt}', ph2dt_inp],
            stdout=f,
            stderr=subprocess.STDOUT,
            check=False
        )

    # --------------------------------------------------
    # Collect outputs
    # --------------------------------------------------
    # ph2dt produces:
    #   event.sel
    #   dt.ct
    if not os.path.exists('event.sel') or not os.path.exists('dt.ct'):
        raise RuntimeError("ph2dt failed: missing output files")

    shutil.move('event.sel', f'{config.indir}/event.dat')
    shutil.move('dt.ct', f'{config.indir}/dt.ct')

    # --------------------------------------------------
    # move log file to the config.indir directory
    # --------------------------------------------------
    if os.path.exists('ph2dt.log'):
        shutil.move('ph2dt.log', f'{config.outdir}/{config.namebase}_ph2dt.log')

    if verbose: print(f"[INFO] ph2dt completed: {out_file}")

    return out_file


#
def update_phasefile_with_reloc(infile,catalog,outfile):
    """
    Update phase file with relocated event info from catalog.
    Parameters
    ----------
    infile : str
        Input phase file path.
    catalog : str
        Catalog file path with relocated event info.
    outfile : str
        Output phase file path.
    """
    #load input phase data
    pha_dict = utils.load_hypo_phasedata(infile)

    # read catalog
    with open(catalog) as f:
        evlines = f.readlines()
    
    # collect output catalog and phase file
    with open(outfile, 'w') as out_pha:
        for evline in evlines:
            evid = evline.strip().split(',')[-1]
            out_pha.write(evline)

            for pha_line in pha_dict[evid]:
                out_pha.write(pha_line)
        #
#
class HypoDDConfig(object):
    """
    Configuration class for HypoDD.
    =============================
    Parameters
    ----------
    binpath : str
        Path to hypoDD binaries. Default: system PATH)
    indir : str
        Input directory. Default: 'input'
    outdir : str
        Output directory. Default: 'output'
    namebase : str
        Base name for input/output files.
    station_file : str
        Path to station file (after reformatting for HypoDD).
    phase_file : str
        Path to phase file (after reformatting for HypoDD).
    dep_corr : float
        Depth correction to avoid air quakes. Default: 5 km.
    hypodd_inp_template : str
        Template file for hypoDD input.inp
    ph2dt_inp_template : str
        Template file for ph2dt input.inp
    =============================   
    Operations
    ----------
    run()
        Run HypoDD once for the entire dataset.
    collect_outputs(pha_dict,cleanup=True)
        Collect and clean up HypoDD outputs.
    help()
        Print help message.
    =============================
    """
    def __init__(self,binpath=None,indir='input',outdir='output',namebase=None,station_file=None, phase_file=None,
               dep_corr = 5,hypodd_inp_template=None,ph2dt_inp_template=None):
        self.type="HypoDDConfig object"
        # i/o paths
        # phase_file: needs to be the file after reformatted to be used by ph2dt and hypoDD.
        if binpath is None:
            HYPODD = shutil.which("hypoDD")
            if HYPODD is None:
                raise RuntimeError("hypoDD not found in PATH")
            self.bin_hypodd = HYPODD # default path to hypoDD binaries, assuming it is in the system PATH

            PH2DT = shutil.which("ph2dt")
            if PH2DT is None:
                raise RuntimeError("ph2dt not found in PATH")
            self.bin_ph2dt = PH2DT
        else:
            self.bin_hypodd = os.path.join(binpath,'hypoDD')
            self.bin_ph2dt = os.path.join(binpath,'ph2dt')
        self.binpath = binpath
        self.indir = indir
        self.outdir = outdir
        self.namebase = namebase
        self.station_file = station_file
        if phase_file is not None:
            self.phase_file = phase_file
        else:
            self.phase_file = f'{self.indir}/{self.namebase}.pha' #this is the phase file after reformatting to be used by ph2dt and hypoDD.
        self.hypodd_inp_template = hypodd_inp_template
        self.ph2dt_inp_template = ph2dt_inp_template

        if self.hypodd_inp_template is None:
            self.hypodd_inp_template = f'{self.indir}/template_hypodd_par.inp'
            template_name = utils.get_template_list('hypodd',pattern='hypodd_par.inp',fullpath=True)[0]
            shutil.copyfile(template_name,self.hypodd_inp_template)
        if self.ph2dt_inp_template is None:
            self.ph2dt_inp_template = f'{self.indir}/template_ph2dt_par.inp'
            template_name = utils.get_template_list('hypodd',pattern='ph2dt_par.inp',fullpath=True)[0]
            shutil.copyfile(template_name,self.ph2dt_inp_template)
            
        # run ph2dt & hypoDD 
        self.dep_corr = dep_corr # avoid air quake, modify velo_mod accordingly

        # create output directory if not exists yet
        if not os.path.exists(self.indir):
            os.makedirs(self.indir)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
    #
    def __str__(self):
        lines = [
            f"{self.type}",
            "-" * len(self.type),

            # Executables
            f"HypoDD binary     : {self.bin_hypodd}",
            f"ph2dt binary      : {self.bin_ph2dt}",
            f"Binary base path  : {self.binpath}",

            # I/O
            f"Input directory   : {self.indir}",
            f"Output directory  : {self.outdir}",
            f"Name base         : {self.namebase}",

            # Input files
            f"Station file      : {self.station_file}",
            f"Phase file        : {self.phase_file}",

            # Templates
            f"HypoDD template   : {self.hypodd_inp_template}",
            f"ph2dt template    : {self.ph2dt_inp_template}",

            # Relocation settings
            f"Depth correction  : {self.dep_corr}",
        ]

        return "\n".join(lines)

    #
    def __repr__(self):
        return (
            f"HypoDDConfig("
            f"namebase={self.namebase}, "
            f"phase_file={self.phase_file}, "
            f"station_file={self.station_file}, "
            f"outdir={self.outdir})"
        )

    """
    Core runner for HypoDD (single-run version)
    """
    def help(self):
        print("HypoDDConfig object:")
        print("  Attributes:")
        print("    binpath: path to hypoDD binaries")
        print("    indir: input directory")
        print("    outdir: output directory")
        print("    namebase: base name for input/output files")
        print("    station_file: path to station file (after reformatting for HypoDD)")
        print("    phase_file: path to phase file (after reformatting for HypoDD)")
        print("    dep_corr: depth correction to avoid air quakes")
        print("    hypodd_inp_template: template file for hypoDD input.inp")
        print("    ph2dt_inp_template: template file for ph2dt input.inp")
        print("  Methods:")
        print("    run(cleanup): Run HypoDD once for the entire dataset.")
        print("    help(): Print this help message.")

    # --------------------------------------------------
    # Public API
    # --------------------------------------------------
    def run(self,verbose=True):
        """
        Run HypoDD once for the entire dataset. 

        Workflow:
            1. run ph2dt
            2. write hypoDD input files
            3. run hypoDD core
 
        Parameters
        ----------
        cleanup : bool
            Whether to clean up intermediate files. Default: True.

        """

        # -------------------------------
        # run ph2dt
        # -------------------------------
        run_ph2dt(self,verbose=verbose)

        self._write_hypodd_input_files()

        if verbose:
            print('running hypodd...')
            self._run_hypodd_core()
        
        print('finished hypoDD flow. Call HypoDDConfig.finalize() to collect the outputs.')
    #
    def finalize(self,cleanup=False,out_catalog_file=None):
        cat_out=self._collect_outputs(out_catalog_file=out_catalog_file)

        #clean up
        if cleanup:
            reloc_grids = glob.glob(f'{self.outdir}/hypoDD_{self.namebase}.reloc.*')
            for f in reloc_grids:
                if os.path.exists(f):
                    os.remove(f)

            files_to_remove = [
                'event.dat',
                'station.sel',
                'hypoDD.log',
                f'{self.indir}/event.dat',
                f'{self.indir}/dt.ct',
                f'{self.outdir}/{self.namebase}.ph2dt',
                f'{self.outdir}/hypoDD_{self.namebase}.loc',
                f'{self.outdir}/hypoDD_{self.namebase}.res',
                f'{self.outdir}/hypoDD_{self.namebase}.sta',
                f'{self.outdir}/hypoDD_{self.namebase}.src',
            ]
            for f in files_to_remove:
                if os.path.exists(f):
                    os.remove(f)

            print("[INFO] Cleaned up intermediate files.")
        #
        return cat_out

    # --------------------------------------------------
    # Internal helpers
    # --------------------------------------------------
    def _write_hypodd_input_files(self):
        """
        Write all HypoDD input files for a single run.
        Assumes ph2dt has already been executed.
        """
        # Input file names
        os.makedirs(self.indir, exist_ok=True)
        fout = open(f'{self.indir}/hypoDD_{self.namebase}.inp','w')
        f=open(self.hypodd_inp_template); lines=f.readlines(); f.close()
        for line in lines:
            #skip lines starting with *
            if line[0] == "*":
                pass
            else:
                if 'dt.ct' in line: line = f'{self.indir}/dt.ct \n'
                if 'dt.cc' in line: line = f'{self.indir}/dt.cc \n'  #
                if 'event.dat' in line: line = f'{self.indir}/event.dat \n'
                if 'station.dat' in line: line = f'{self.station_file} \n'

                if 'hypoDD.reloc' in line: line = f'{self.outdir}/hypoDD_{self.namebase}.reloc \n'
                if 'hypoDD.loc' in line: line = f'{self.outdir}/hypoDD_{self.namebase}.loc \n'
                if 'hypoDD.sta' in line: line = f'{self.outdir}/hypoDD_{self.namebase}.sta \n'
                if 'hypoDD.res' in line: line = f'{self.outdir}/hypoDD_{self.namebase}.res \n'
                if 'hypoDD.src' in line: line = f'{self.outdir}/hypoDD_{self.namebase}.src \n'
            
            fout.write(line)
        fout.close()

    def _run_hypodd_core(self):
        """Execute HypoDD binary"""
        exe = self.bin_hypodd
        inp = Path(self.indir) / f'hypoDD_{self.namebase}.inp'

        if not os.path.exists(exe):
            raise FileNotFoundError(f"hypoDD binary not found: {exe}")

        cmd = [str(exe), str(inp)]

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if proc.returncode != 0:
            raise RuntimeError(
                "HypoDD failed:\n" + proc.stderr
            )

        # Write stdout log
        logf = os.path.join(self.outdir, f'{self.namebase}_hypoDD.log')
        with open(logf, 'w') as f:
            f.write(str(proc.stdout))

    def _collect_outputs(self, out_catalog_file=None):
        """
        Collect, clean, and save HypoDD output catalog.

        Parameters
        ----------
        out_catalog_file : str
            Path to save final catalog CSV. Default: {self.outdir}/{self.namebase}.ctlg

        Returns
        -------
        pandas.DataFrame
            DataFrame containing final relocated events.
        """
        if out_catalog_file is None:
            out_catalog_file = f'{self.outdir}/{self.namebase}.ctlg'

        freloc = f'{self.outdir}/hypoDD_{self.namebase}.reloc'
        if not os.path.exists(freloc):
            warnings.warn(f"HypoDD reloc file not found: {freloc}")
            return pd.DataFrame(
                columns=["origin_time", "latitude", "longitude", "depth_km", "magnitude", "event_id"]
            )

        records = []
        with open(freloc) as f:
            lines = f.readlines()

        for line in lines:
            codes = line.split()
            if len(codes) < 17:
                continue  # skip incomplete lines

            try:
                evid = codes[0]

                # parse location
                lat = float(codes[1])
                lon = float(codes[2])
                dep = round(float(codes[3]) - self.dep_corr, 2)

                # magnitude
                mag = float(codes[16])

                # origin time
                year, mon, day, hour, mnt, sec = codes[10:16]
                sec = '59.999' if sec == '60.000' else sec
                ot = UTCDateTime(f'{year}{mon:0>2}{day:0>2}{hour:0>2}{mnt:0>2}{sec:0>6}')
                time_str = ot.strftime("%Y-%m-%dT%H:%M:%S.%f")
                records.append({
                    "time": time_str,
                    "latitude": lat,
                    "longitude": lon,
                    "depth_km": dep,
                    "magnitude": mag,
                    "event_id": evid
                })

            except Exception:
                continue

        df = pd.DataFrame.from_records(records)

        # save CSV with header
        df.to_csv(out_catalog_file, index=False)

        print(f"[INFO] Wrote output catalog: {out_catalog_file} with {len(df)} events.")

        return df


