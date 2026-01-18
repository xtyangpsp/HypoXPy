#import needed packages.
""" 
    Download hypoDD at https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Or use the copy included in this python interface package.
"""
from pathlib import Path
import os,shutil, glob, subprocess, warnings
import numpy as np
from obspy import UTCDateTime
warnings.filterwarnings("ignore")

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
    run(pha_dict)
        Run HypoDD once for the entire dataset.
    
    """
    def __init__(self,binpath=None,indir='input',outdir='output',namebase=None,station_file=None, phase_file=None,
               dep_corr = None,hypodd_inp_template='hypoDD.inp',ph2dt_inp_template='ph2dt.inp'):

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
        # run ph2dt & hypoDD 
        self.dep_corr = dep_corr or 5 # avoid air quake, modify velo_mod accordingly

        # create output directory if not exists yet
        if not os.path.exists(indir):
            os.makedirs(indir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        #
        #

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
        print("    run(pha_dict): Run HypoDD once for the entire dataset.")
        print("    help(): Print this help message.")

    # --------------------------------------------------
    # Public API
    # --------------------------------------------------
    def run(self,pha_dict,cleanup=True):
        """
        Run HypoDD once for the entire dataset. 

        Workflow:
            1. run ph2dt
            2. write hypoDD input files
            3. run hypoDD core
            4. collect outputs

        Parameters
        ----------
        pha_dict : dict
            Dictionary of phase data, keyed by event ID.
        cleanup : bool
            Whether to clean up intermediate files. Default: True.
        Returns
        -------
        outfile_catalog : str
            Path to output catalog file.
        outfile_phase : str
            Path to output phase file.
        outfile_phase_full : str
            Path to full output phase file.
        """

        # -------------------------------
        # run ph2dt
        # -------------------------------
        run_ph2dt(self)

        self._write_input_files()
        self._run_hypodd_core()
        self._collect_outputs(pha_dict)

        # cleanup intermediate files
        if cleanup:
            reloc_grids = glob.glob(f'{self.outdir}/hypoDD_{self.namebase}.reloc.*')
            for f in reloc_grids:
                if os.path.exists(f):
                    os.unlink(f)

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

    def _collect_outputs(self,pha_dict):
        """
        Collect and clean up HypoDD outputs.
        Produces:
          - {outdir}/{namebase}.ctlg
          - {outdir}/{namebase}.pha
          - {outdir}/{namebase}_full.pha
        """
        evid_list = list(pha_dict.keys())
        outfile_catalog = f'{self.outdir}/{self.namebase}.ctlg'
        outfile_phase = f'{self.outdir}/{self.namebase}.pha'
        outfile_phase_full = f'{self.outdir}/{self.namebase}_full.pha'
        # clean up outputs
        with open(outfile_catalog, 'w') as out_ctlg, \
            open(outfile_phase, 'w') as out_pha, \
            open(outfile_phase_full, 'w') as out_pha_full:

            freloc = f'{self.outdir}/hypoDD_{self.namebase}.reloc'
            if not os.path.exists(freloc):
                return

            with open(freloc) as f:
                lines = f.readlines()

            for line in lines:
                codes = line.split()
                evid = codes[0]

                if int(evid) not in evid_list:
                    continue

                pha_lines = pha_dict[evid]
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

                out_ctlg.write(f'{ot},{lat},{lon},{dep},{mag}\n')
                out_pha.write(f'{ot},{lat},{lon},{dep},{mag}\n')
                out_pha_full.write(f'{ot},{lat},{lon},{dep},{mag},{evid}\n')

                for pha_line in pha_lines:
                    out_pha.write(pha_line)
                    out_pha_full.write(pha_line)
        #
        print(f"[INFO] Wrote output catalog: {outfile_catalog}")
        print(f"[INFO] Wrote output phase file: {outfile_phase}")
        print(f"[INFO] Wrote full output phase file: {outfile_phase_full}")

        return outfile_catalog, outfile_phase, outfile_phase_full
