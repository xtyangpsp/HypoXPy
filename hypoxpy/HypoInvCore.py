#This module contains core functions for running HypoInvPy interface.
#Import needed packages first.
import pandas as pd
import os,glob
import numpy as np
from hypoxpy import utils
import subprocess
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
        fhyp = os.path.join(pardir,'%s-%s.hyp'%(config.namebase, ztr))
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
            if line[0:3]=='SUM': line = "SUM '%s/%s-%s.sum' \n"%(outdir,config.namebase, ztr)
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
                 smodel=None,poisson=1.73,
               namebase='hyp',hypoinv_bin='hyp1.40',get_prt=False,get_arc=False,
               lat_code='N',lon_code='W',ref_ele=0.0,grd_ele=0.0,
               ztrlist = np.arange(0,20,1),rms_weight='4 0.3 1 3',dist_initial = '1 50 1 2',
               dist_weight = '4 20 1 3',weight_code='1 0.6 0.3 0.2',min_nsta=4):
        if binpath is None:
            binpath = 'hyp1.40' # default path to hypoinverse binary, assuming it is in the system PATH
        self.binpath = binpath
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
        self.template_parfile = utils.get_template_list('hypoinv')[1]
        self.pmodel = pmodel #'input/velo_p_eg.cre'
        self.smodel = smodel #[None, 'input/velo_s_eg.cre'][1]
        self.poisson = poisson #1.73 # provide smod or pos

        #
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if not os.path.exists(self.indir):
            os.makedirs(self.indir)

    #-------------------------------------------------
    # core function to run hypoinverse
    #-------------------------------------------------
    def run(self,parfilelist, cleanup=True, merge_summary=False):
        """
        Run hypoinverse for a list of parameter files.
        ======== PARAMETERS ==========
        parfilelist: list of parameter files for hypoinverse.
        cleanup: whether to remove intermediate files after running hypoinverse. Default True. 
                If merge_summary is True, the summary files will be removed after merging is done.

        merge_summary: whether to merge summary files after running hypoinverse. Default False.
                User can also call merge_summary() function separately.
        ========
        """
        
        for fhyp in parfilelist:
            # 2. run hypoinverse
            p = subprocess.Popen([self.binpath], stdin=subprocess.PIPE,encoding='utf-8')
            s = "@{}".format(fhyp) + '\n'
            p.communicate(s)
        #
        if merge_summary:
            self.merge_summary(cleanup=cleanup)
        #-------------------------------------------------
        # remove intermidiate files
        if cleanup:
            for fname in glob.glob('fort.*'): os.unlink(fname)
            for fname in glob.glob('input/%s-*.hyp'%self.namebase): os.unlink(fname)

    ####
    def merge_summary(self,filelist=None,mag_dict=None,cleanup=True):
        """
        Extract earthquake parameters based on final quality after merging all summary files.

        ======PARAMETERS======
        filelist: list of summary files to merge. Default None (all summary files in output dir).
        mag_dict: magnitude dictionary in the form of {'id',mag} for all events. Default None (use hypoinverse output).
                mag_dict could also be specified as a catalog file in csv format. E.g., the catalog from GAMMA.
        """
        if filelist is None:
            filelist=glob.glob('%s/%s-*.sum'%(self.outdir,self.namebase))
        file_good='%s/%s_good.csv'%(self.outdir,self.namebase)
        file_bad='%s/%s_bad.csv'%(self.outdir,self.namebase)

        if mag_dict is not None:
            if isinstance(mag_dict,str):
                events=pd.read_csv(mag_dict)
                mag_dict_use=dict()
                for i in range(len(events.time)):
                    event = events.iloc[i]
                    # print(event['event_index'].astype(str))
                    mag_dict_use[event['event_index'].astype(str)] = event['magnitude']
                #
            elif isinstance(mag_dict,dict):
                mag_dict_use = mag_dict
            else:
                raise ValueError('mag_dict is wrong in type. CSV catalog or dictionary.')
        else:
            mag_dict_use = mag_dict
        fout_bad = open(file_bad,'w')
        fout_good = open(file_good,'w')
        
        # read sum files
        sum_dict = {}
        for fsum in filelist:
            f=open(fsum); sum_lines=f.readlines(); f.close()
            for sum_line in sum_lines:
                evid = sum_line.split()[-1]
                if evid not in sum_dict: sum_dict[evid] = [sum_line]
                else: sum_dict[evid].append(sum_line)
        
        # merge sum lines
        for evid, sum_lines in sum_dict.items():
            sum_list = []
            dtype = [('line','O'),('is_loc','O'),('qua','O'),('azm','O'),('npha','O'),('rms','O')]
            for sum_line in sum_lines:
                codes = sum_line.split()
                is_loc = 1 # whether loc reliable
                if '-' in codes or '#' in codes: is_loc = 0
                qua = sum_line[80:81]
                npha = 1 / float(sum_line[52:55])
                azm  = float(sum_line[56:59])
                rms  = float(sum_line[64:69])
                sum_list.append((sum_line, is_loc, qua, azm, npha, rms))
            sum_list = np.array(sum_list, dtype=dtype)
            sum_list = np.sort(sum_list, order=['qua','azm','npha','rms'])
            sum_list_loc = sum_list[sum_list['is_loc']==1]
            num_loc = len(sum_list_loc)
            # if no reliable loc
            if num_loc==0:
                sum_list_loc = sum_list
                utils.write_csv(fout_bad, sum_list_loc[0]['line'], evid,self.lat_code,self.lon_code,mag_dict=mag_dict_use)
            else:
                utils.write_csv(fout_good, sum_list_loc[0]['line'], evid,self.lat_code,self.lon_code,mag_dict=mag_dict_use)

        fout_bad.close()
        fout_good.close()

        print('Earthquakes are saved in: '+file_good+' and '+file_bad+' for good and bad sources.')
        
        # remove summary files.
        if cleanup:
            for fname in glob.glob(self.outdir+'/'+self.namebase+'*.sum'): os.unlink(fname)



        

