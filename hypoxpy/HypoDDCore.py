#import needed packages.
""" 
    Download hypoDD at https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Or use the copy included in this python interface package.
"""
from pathlib import Path
import os,shutil, glob, subprocess, warnings
import numpy as np
from obspy import UTCDateTime
from hypoxpy import utils
warnings.filterwarnings("ignore")

def reformat_stationfile(fin,fout):
    """ 
    Reformat input station file for HypoDD.
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
#
#
def run_ph2dt(config):
    """
    Call ph2dt once for the entire dataset (no grids).

    config: HypoDDConfig object
    """

    print("running ph2dt...")

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
        if 'input/phase.dat' in line: line = '%s \n'%(config.phase_file)
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

    print(f"[INFO] ph2dt completed: {out_file}")

    return out_file

# 
def reformat_phasefile(config, phase_file_in,phase_file_out=None, out_dir=None,
                       time_range=None,lat_range=None,lon_range=None):
    """
    Reformat input phase file for HypoDD into ONE single phase file.
    Output:
      {fout_dir}/{config.namebase}.pha

    Still applies:
      - time window filtering
      - lat/lon window filtering
      - depth correction

    Parameters
    ----------
    config : HypoDDConfig
    fout_dir : str
        Output directory
    """
    if phase_file_out is not None:
        outfile = phase_file_out
        #update config.phase_file accordingly
        config.phase_file = phase_file_out
    else:
        outfile = config.phase_file 
    if out_dir is not None:
        fout_dir = out_dir
    else:
        fout_dir = config.indir

    os.makedirs(fout_dir, exist_ok=True)

    dep_corr = config.dep_corr

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

    with open(outfile, 'w') as fout:
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
    print(f"[INFO] Wrote single phase file: {outfile}")
    print(f"[INFO] Number of events: {len(evid_list)}")

    return outfile, np.array(evid_list)
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
               dep_corr = None,hypodd_inp_template=None,ph2dt_inp_template=None):

        # i/o paths
        # phase_file: needs to be the file after reformatted to be used by ph2dt and hypoDD.
        if binpath is None:
            self.bin_hypodd = 'hypoDD' # default path to hypoDD binaries, assuming it is in the system PATH
            self.bin_ph2dt = 'ph2dt'
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

        if self.hypodd_inp_template is None or self.ph2dt_inp_template is None:
            if self.hypodd_inp_template is None:
                self.hypodd_inp_template = f'{self.indir}/template_hypodd_par.inp'
                template_name = utils.get_template_list('hypodd',pattern='hypodd_par.inp',fullpath=True)[0]
                shutil.copyfile(template_name,self.hypodd_inp_template)
            if self.ph2dt_inp_template is None:
                self.ph2dt_inp_template = f'{self.indir}/template_ph2dt_par.inp'
                template_name = utils.get_template_list('hypodd',pattern='ph2dt_par.inp',fullpath=True)[0]
                shutil.copyfile(template_name,self.ph2dt_inp_template)
            
        # run ph2dt & hypoDD 
        self.dep_corr = dep_corr or 5 # avoid air quake, modify velo_mod accordingly

        # create output directory if not exists yet
        if not os.path.exists(self.indir):
            os.makedirs(self.indir)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
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
    def run(self, cleanup=True):
        """
        Run HypoDD once for the entire dataset. 

        Workflow:
            1. run ph2dt
            2. write hypoDD input files
            3. run hypoDD core
            4. collect outputs
 
        Parameters
        ----------
        cleanup : bool
            Whether to clean up intermediate files. Default: True.
        =============
        Returns
        -------
        outfile_catalog : str
            Path to output catalog file.
        """

        # -------------------------------
        # run ph2dt
        # -------------------------------
        run_ph2dt(self)

        self._write_input_files()
        self._run_hypodd_core()
        outcatalog = self._collect_outputs(cleanup=cleanup)

        return outcatalog
        #return outcatalog

    # --------------------------------------------------
    # Internal helpers
    # --------------------------------------------------
    def _write_input_files(self):
        """
        Write all HypoDD input files for a single run.
        Assumes ph2dt has already been executed.
        """
        # Input file names
        os.makedirs(self.indir, exist_ok=True)
        fout = open(f'{self.indir}/hypoDD_{self.namebase}.inp','w')
        f=open(self.hypodd_inp_template); lines=f.readlines(); f.close()
        for line in lines:
            if 'dt.ct' in line: line = f'{self.indir}/dt.ct \n'
            if 'event.dat' in line: line = f'{self.indir}/event.dat \n'
            if 'hypoDD.reloc' in line: line = f'{self.outdir}/hypoDD_{self.namebase}.reloc \n'
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

    def _collect_outputs(self,cleanup=True):
        """
        Collect and clean up HypoDD outputs.

        Parameters
        ----------
        pha_dict : dict
            Dictionary of phase lines keyed by event ID.
        cleanup : bool
            Whether to clean up intermediate files. Default: True.
        Returns
        -------
        outfile_catalog : str
            Path to output catalog file.
        """
        outfile_catalog = f'{self.outdir}/{self.namebase}.ctlg'
        # collect output catalog and phase file
        with open(outfile_catalog, 'w') as out_ctlg:

            freloc = f'{self.outdir}/hypoDD_{self.namebase}.reloc'
            if not os.path.exists(freloc):
                return

            with open(freloc) as f:
                lines = f.readlines()

            for line in lines:
                codes = line.split()
                evid = codes[0]

                # location
                lat, lon, dep = codes[1:4]
                try:
                    dep = round(float(dep) - self.dep_corr, 2)
                    mag = float(codes[16])
                except Exception:
                    continue

                # origin time
                year, mon, day, hour, mnt, sec = codes[10:16]
                sec = '59.999' if sec == '60.000' else sec

                ot = UTCDateTime(
                    f'{year}{mon:0>2}{day:0>2}{hour:0>2}{mnt:0>2}{sec:0>6}'
                )

                out_ctlg.write(f'{ot},{lat},{lon},{dep},{mag},{evid}\n')
        #
        print(f"[INFO] Wrote output catalog: {outfile_catalog} with {len(lines)} events.")  

        #
        """Remove intermediate files generated by HypoDD."""    
        if cleanup:
            # cleanup intermediate files
            reloc_grids = glob.glob(f'{self.outdir}/hypoDD_{self.namebase}.reloc.*')
            for f in reloc_grids:
                if os.path.exists(f):
                    os.remove(f)
            print("[INFO] Cleaned up intermediate files.")

        return outfile_catalog

