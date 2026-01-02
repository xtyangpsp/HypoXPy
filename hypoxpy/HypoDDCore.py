#import needed packages.
""" 
    Download hypoDD at https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Or use the copy included in this python interface package.
"""
from pathlib import Path
import os,shutil, glob
import numpy as np
from obspy import UTCDateTime
from pathlib import Path
import subprocess
import warnings
warnings.filterwarnings("ignore")

def cat_files(pattern, fout):
    
    with open(fout, 'w') as out:
        for fname in sorted(glob.glob(pattern)):
            with open(fname) as f:
                shutil.copyfileobj(f, out)

# read fpha with evid
def load_phasedata(fin):
    pha_dict = {}
    evid = None

    with open(fin) as f:
        for line in f:
            codes = line.split(',')
            if len(codes[0]) >= 14:
                evid = codes[-1].strip()
                pha_dict[evid] = []
            elif evid is not None:
                pha_dict[evid].append(line)

    return pha_dict

#
def run_ph2dt(config):
    """
    Call ph2dt once for the entire dataset (no grids).

    config: HypoDDConfig object
    """

    print("run ph2dt (single run)")

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
            [f'{config.binpath}/ph2dt', ph2dt_inp],
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

def reformat_stationfile(fin,fout):
    """ 
    Reformat input station file for HypoDD
    """
    foutid=open(fout,'w')
    
    done_list = []
    f=open(fin); 
    lines=f.readlines(); 
    f.close()
    for line in lines:
        codes = line.split(',')
        net, sta = codes[0].split('.')
        if sta in done_list: continue
        lat, lon = [float(code) for code in codes[1:3]]
        foutid.write('{} {} {}\n'.format(sta, lat, lon))
        done_list.append(sta)
    foutid.close()
#

# grid params
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

    outfile = os.path.join(fout_dir, f"{config.namebase}.pha")

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

                sta = codes[0].split('.')[1]
                wp, ws = 1.0, 1.0

                # P pick
                if codes[1] != '-1':
                    tp = UTCDateTime(codes[1])
                    ttp = tp - ot
                    fout.write(
                        '{:<5}{}{:6.3f}  {:6.3f}   P\n'.format(
                            sta, ' ' * 6, ttp, wp
                        )
                    )

                # S pick
                if codes[2] != '-1':
                    ts = UTCDateTime(codes[2])
                    tts = ts - ot
                    fout.write(
                        '{:<5}{}{:6.3f}  {:6.3f}   S\n'.format(
                            sta, ' ' * 6, tts, ws
                        )
                    )

    # --------------------------------------------------
    print(f"[INFO] Wrote single phase file: {outfile}")
    print(f"[INFO] Number of events: {len(evid_list)}")

    return outfile, np.array(evid_list)
#
#
class HypoDDConfig(object):
  def __init__(self,binpath=None,indir='input',outdir='output',namebase=None,station_file=None, phase_file=None,
               dep_corr = None,hypodd_inp_template='hypoDD.inp',ph2dt_inp_template='ph2dt.inp'):

    # i/o paths
    # phase_file: needs to be the file after reformatted to be used by ph2dt and hypoDD.
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
    # run ph2dt & hypoDD 
    self.dep_corr = dep_corr or 5 # avoid air quake, modify velo_mod accordingly

    # create output directory if not exists yet
    if not os.path.exists(indir):
        os.makedirs(indir)
    if not os.path.exists(outdir):
      os.makedirs(outdir)
#
class HypoDDCore(object):
    """
    Core runner for HypoDD (single-run version)

    Responsibilities:
      - prepare HypoDD input files
      - run hypoDD once for the entire dataset
      - no grids, no parallelism, no merging
    """

    def __init__(self, config, evid_list, pha_dict):
        self.config = config
        self.evid_list = evid_list
        self.pha_dict = pha_dict

        self.binpath = Path(self.config.binpath)
        self.outdir = Path(self.config.outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # Public API
    # --------------------------------------------------
    def run(self):
        """Run HypoDD once for the entire dataset"""
        self._write_input_files()
        self._run_hypodd()
        self._collect_outputs()

    # --------------------------------------------------
    # Internal helpers
    # --------------------------------------------------
    def _write_input_files(self):
        """
        Write all HypoDD input files for a single run.
        Assumes ph2dt has already been executed.
        """
        # Input file names
        os.makedirs('input', exist_ok=True)
        fout = open('input/hypoDD_%s.inp'%(self.config.namebase),'w')
        f=open(self.config.hypodd_inp_template); lines=f.readlines(); f.close()
        for line in lines:
            if 'dt.ct' in line: line = 'input/dt.ct \n'
            if 'event.dat' in line: line = 'input/event.dat \n'
            if 'hypoDD.reloc' in line: line = 'output/hypoDD_%s.reloc \n'%(self.config.namebase)
            fout.write(line)
        fout.close()

    def _run_hypodd(self):
        """Execute HypoDD binary"""
        exe = self.binpath / 'hypoDD'
        inp = Path('input') / f'hypoDD_{self.config.namebase}.inp'

        if not exe.exists():
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
        logf = self.outdir / f'{self.config.namebase}_hypoDD.log'
        with open(logf, 'w') as f:
            f.write(str(proc.stdout))

    def _collect_outputs(self):
        """
        Collect and clean up HypoDD outputs.
        Produces:
          - {outdir}/{namebase}.ctlg
          - {outdir}/{namebase}.pha
          - {outdir}/{namebase}_full.pha
        """
        # clean up outputs
        with open(f'{self.config.outdir}/{self.config.namebase}.ctlg', 'w') as out_ctlg, \
            open(f'{self.config.outdir}/{self.config.namebase}.pha', 'w') as out_pha, \
            open(f'{self.config.outdir}/{self.config.namebase}_full.pha', 'w') as out_pha_full:

            freloc = f'{self.config.outdir}/hypoDD_{self.config.namebase}.reloc'
            if not os.path.exists(freloc):
                return

            with open(freloc) as f:
                lines = f.readlines()

            for line in lines:
                codes = line.split()
                evid = codes[0]

                if int(evid) not in self.evid_list:
                    continue

                pha_lines = self.pha_dict[evid]
                # location
                lat, lon, dep = codes[1:4]
                try:
                    dep = round(float(dep) - self.config.dep_corr, 2)
                    mag = float(codes[16])
                except Exception:
                    continue

                # origin time
                year, mon, day, hour, mnt, sec = codes[10:16]
                sec = '59.999' if sec == '60.000' else sec

                ot = UTCDateTime(
                    f'{year}{mon:0>2}{day:0>2}{hour:0>2}{mnt:0>2}{sec:0>6}'
                )

                out_ctlg.write(f'{ot},{lat},{lon},{dep},{mag}\n')
                out_pha.write(f'{ot},{lat},{lon},{dep},{mag}\n')
                out_pha_full.write(f'{ot},{lat},{lon},{dep},{mag},{evid}\n')

                for pha_line in pha_lines:
                    out_pha.write(pha_line)
                    out_pha_full.write(pha_line)

