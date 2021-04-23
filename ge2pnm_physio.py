
import os
import numpy
import tkinter as tk
from tkinter import filedialog as fd



def combine_physio(resp_file=None, ppg_file=None, scantrig_file=None, ppgtype=None, output_file=None, overwrite=0):

    #Compile various physiology files into a single file, as required by FSL
    
    if resp_file is None:
        raise RuntimeError('No resp_file passed!')

    if ppg_file is None:
        raise RuntimeError('No card_file passed!')

    if scantrig_file is None:
        raise RuntimeError('No scantrig_file passed!')

    if os.path.exists(output_file):
        if overwrite:
            print('Output file exists and overwrite set.')
            print('Deleting: {}'.format(output_file))
            os.remove(output_file)
        else:
            raise RuntimeError('Output file exists and overwrite not set.')


    #Read in the three files
    with open(ppg_file, 'r') as ppg_open:
        ppg_contents = []
        for line in ppg_open:
            ppg_contents.append(line.strip())

    with open(scantrig_file, 'r') as scantrig_open:
        scantrig_contents = []
        for line in scantrig_open:
            scantrig_contents.append(line.strip())

    with open(resp_file, 'r') as resp_open:
        resp_contents = []
        for line in resp_open:
            resp_contents.append(line.strip())

    #Make sure the number of entries in each file matches the others
    #If the number of entries in the cardiac trigger file is within 4 of the entries in the respiratory file,
    #make the trigger file match the respiratory file.
    if (len(ppg_contents) != len(resp_contents)) or (len(scantrig_contents) != len(resp_contents)):
        print('Number of entries in the physiology files do not match!')
        print('PPG entires: {}'.format(len(ppg_contents)))
        print('Scantrig entires: {}'.format(len(scantrig_contents)))
        print('Respiratory entires: {}'.format(len(resp_contents)))
        if abs((len(ppg_contents) - len(resp_contents)) < 5):
            print('The number of entries in the ppg file is within 4 of the number of entries in the respiratory file.')
            if ppgtype=='Triggers':
                print('Matching the length of the ppg trigger file to the resp file.')
                if len(ppg_contents) < len(resp_contents):
                    while len(ppg_contents) < len(resp_contents):
                        ppg_contents.append('0')
                elif len(ppg_contents) > len(resp_contents):
                    while len(ppg_contents) > len(resp_contents):
                        ppg_contents.pop()
            else:
                print('Trimming the longer one to match the shorter one.')
                if len(ppg_contents) < len(resp_contents):
                    while len(ppg_contents) < len(resp_contents):
                        resp_contents.pop()
                elif len(ppg_contents) > len(resp_contents):
                    while len(ppg_contents) > len(resp_contents):
                        ppg_contents.pop()
        else:
            print('EXITTING!')
            raise OSError('Number of entries in the physiology files do not match!')
    if len(ppg_contents) != len(resp_contents):
        print('Length of cardiac trigger time-course still does not match that of the respiratory time-course after attempted correction!')
        print('EXITTING!')
        raise OSError('Number of entries in the physiology files do not match!')

    #Write the combined physiology file
    with open(output_file, 'w') as f:
        for element in range(len(resp_contents)):
            f.write('{} {} {}\n'.format(ppg_contents[element], resp_contents[element], scantrig_contents[element]))

    print('Combined physiology file written: {}'.format(output_file))

    return output_file



def create_scantrig(physio_data_file=None, sampling_rate=None, TR=None, overwrite=0):
    """Create a scanner trigger file that contains 1s at the beginnings of TRs
        and 0s elsewhere.
    """

    print('-----Starting: create_scantrig-----')

    if physio_data_file is None:
        raise RuntimeError('No file passed!')
    else:
        print('Setting input file to: '+str(physio_data_file))
        physio_data_file = str(physio_data_file)

    if os.path.exists(physio_data_file) != 1:
        raise RuntimeError('Input file cannot be found: {}'.format(physio_data_file))

    if sampling_rate is None:
        raise RuntimeError('No sampling rate passed!')
    else:
        sampling_rate = float(sampling_rate)

    if TR is None:
        raise RuntimeError('No TR passed!')
    else:
        TR = float(TR)

    with open(physio_data_file, 'r') as f:
        contents = []
        for line in f:
            contents.append(line)

    num_entries = len(contents)
    scantrig_time = numpy.zeros(num_entries)

    #Make every entry at the beginning of a TR 1
    out_file_path = os.path.split(physio_data_file)[0]
    out_file_name = os.path.join(out_file_path, 'Scanner_Trig.txt')

    if os.path.exists(out_file_name):
        if overwrite:
            print('Overwrite set; deleting existing: {}'.format(out_file_name))
        else:
            raise RuntimeError('Scanner trig file exists and overwrite not set!')

    with open(out_file_name, 'w') as out_file_open:

        tr = float(TR)
        count_max = (sampling_rate * tr) - 1
        count = 1
        out_file_open.write('1\n')
        master_count = 1
        for element in range(1,len(scantrig_time)):
            if count > count_max:
                if (num_entries - master_count) > (sampling_rate / 2.0):   #Don't write another 1 if we're very close to the end of the file
                  out_file_open.write('1\n')
                  count = 1
                else:
                  out_file_open.write('0\n')
                  count = 1
            else:
                out_file_open.write('0\n')
                count = count + 1
            master_count = master_count + 1

    print('Scanner trigger file written.')
    print('Returning file name: {}'.format(out_file_name))

    print('-----Finishing: create_scantrig-----')

    return out_file_name




def cut_thirty_seconds(samp_per_sec=None, input_file=None, overwrite=0):
    """Remove the first 30 seconds from a physiology file.

    The GE physiology files start from 30 seconds before the beginning of
    the fMRI scan and go until the end of the scan.

    ARGUMENTS:
        (input_file):  '/PATH/FILE'; file to cut
        (samp_per_sec):  INT; sampling frequency of the input physiology file
    """

    print('-----Starting: cut_thirty_seconds-----')

    #The first 30 seconds of each respiratory and cardiac file need to be removed because they correspond
    #to time before the start of the scan acquisition.

    if input_file is None:
        raise RuntimeError('No input file passed')
    else:
        print('Setting input file to: '+str(input_file))
        input_file = str(input_file)

    if not os.path.exists(input_file):
        raise RuntimeError('Input file cannot be found!')

    read_okay = os.access(input_file, os.R_OK)
    if not read_okay:
        raise RuntimeError('Input file cannot be read from!')

    if samp_per_sec is None:
        raise RuntimeError('No sampline frequency (samp_per_sec) passed!')

    #Read input file
    input_list = []
    with open(input_file, 'r') as input_file_open:
        for line in input_file_open:
            input_list.append(line)

    #Calculate number of elements to remove
    cut_num = 30 * int(samp_per_sec)
    print('Number of time points being removed: '+str(cut_num))

    #Write left-over elements to a new file
    new_file_list = input_list[cut_num:]
    new_file_name = input_file+'_cut'

    if os.path.exists(new_file_name):
        if overwrite:
            print('Overwrite set; deleting existing: {}'.format(new_file_name))
        else:
            raise RuntimeError('Cut file exists and overwrite not set!')

    print('Writing new file: '+str(new_file_name))
    with open(new_file_name, 'w') as new_file_open:
        for element in new_file_list:
            new_file_open.write(element)

    print('-----Finishing: cut_thirty_seconds-----')

    return new_file_name


def cut_presteady_physio(cut_num=None, tr=None, samp_per_sec=None, input_file=None, overwrite=0):
    """Remove the pre-steady-state time from a physiology file.  The amount of
    time removed will depend on how many TRs were removed from the fMRI data.


    ARGUMENTS:
        (input_file):  '/PATH/FILE'; file to cut
        (samp_per_sec):  INT; sampling frequency of the input physiology file
        (cut_num):  INT; the number of volumes removed from the fMRI data.  The safest thing
                         to do is just leave it at None; it will then read in the class
                         variable value.
        (tr): FLOAT; the TR, in seconds, of the fMRI data.  The safest thing to do is leave it set to None;
                     as with cut_num, it will then read in the class variable value.
    """

    print('-----Starting: cut_presteady_physio-----')

    if input_file is None:
        raise RuntimeError('No input file passed!')
    else:
        print('Setting input file to: '+str(input_file))
        input_file = str(input_file)

    if os.path.exists(input_file) != 1:
        raise RuntimeError('Input file cannot be found: {}'.format(input_file))

    read_okay = os.access(input_file, os.R_OK)
    if read_okay != 1:
        raise RuntimeError('Input file cannot be read from: {}'.format(input_file))

    if samp_per_sec is None:
        raise RuntimeError('No sampling frequency passed!')

    if cut_num is None:
        raise RuntimeError('No number of cut time-points passed!')

    if tr is None:
        raise RuntimeError('No TR passed!')
    else:
        tr = float(tr)

    print('Removing pre-steady-state time from physiology file: {}'.format(input_file))
    print('Samples per second: {}'.format(samp_per_sec))
    print('Number of volumes removed from fMRI data: {}'.format(cut_num))
    print('fMRI data TR (sec): {}'.format(tr))

    #Read input file
    input_list = []
    with open(input_file, 'r') as input_file_open:
        for line in input_file_open:
            input_list.append(line)

    #Calculate number of elements to remove
    remove_num = int(cut_num * tr * int(samp_per_sec))
    print('Number of time points being removed: '+str(remove_num))

    #Write left-over elements to a new file
    new_file_list = input_list[remove_num:]
    new_file_name = input_file+'_poststeady'

    if os.path.exists(new_file_name):
        if overwrite:
            print('Overwrite set; deleting existing: {}'.format(new_file_name))
        else:
            raise RuntimeError('Trimmed file exists and overwrite not set!')

    print('Writing new file: '+str(new_file_name))
    with open(new_file_name, 'w') as new_file_open:
        for element in new_file_list:
            new_file_open.write(element)

    return new_file_name



def convert_trig(input_file=None, example_data_file=None, overwrite=0):
    """Change the PPG trigger file from a list of trigger recording points
       to a time-course of '1's and '0's.

       The subsampled PPG data file must be supplied to set the length of the new trigger
       file.
       The input_file must be the 'subsampled' PPG trigger file.
    """

    print('-----Starting: convert_trig-----')

    #The PPG trigger file originally comes as a list of the recording time points
    #at which a trigger occured.  In order to work with FSL, this needs to be
    #changed into a time-course file where time points without a trigger are set
    #to 0 and points with a trigger are set to 1.

    #Check inputs
    if input_file is None:
        raise RuntimeError('No input file passed!')
    else:
        print('Setting input file to: '+str(input_file))
        input_file = str(input_file)

    if example_data_file is None:
        raise RuntimeError('No PPG data file passed!')
    else:
        print('Setting PPG data file to: '+str(example_data_file))
        input_file = str(input_file)

    #Read PPG trigger file
    with open(input_file, 'r') as f:
        trig_list = []
        for line in f:
            trig_list.append(int(line.strip()))

    #Read PPG data file
    with open(example_data_file, 'r') as f:
        data_list = []
        for line in f:
            data_list.append(line.strip())

    #Create an array of zeroes with length equal to the PPG data file
    print('Number of time-points in subsampled PPG data file: {}'.format(len(data_list)))
    trig_timecourse = numpy.zeros(len(data_list))

    #Set time points with a PPG trigger to 1
    for element in trig_list:
        if (element-1) < len(trig_timecourse):
            trig_timecourse[element-1] = 1
        else:
            print('Trigger file entry is outside the length of the example data file!')
            print('trig_list element: {}'.format(element))
            print('len(trig_timecourse): {}'.format(len(trig_timecourse)))
            raise RuntimeError('Trigger entry outside length of example data file!')

    #Write trigger time-course file
    trigger_file_name = str(input_file)+'_timecourse'

    if os.path.exists(trigger_file_name):
        if overwrite:
            print('Overwrite set; deleting existing: {}'.format(trigger_file_name))
        else:
            raise RuntimeError('Trigger timecourse file exists and overwrite not set!')

    with open(trigger_file_name, 'w') as trigger_file:
        for element in trig_timecourse:
            trigger_file.write(str(element)+'\n')

    print('Trigger time-course file written: {}'.format(trigger_file_name))
    
    print('-----Finishing: convert_trig-----')

    return trigger_file_name



def resample_trig(input_file=None, subsample_factor=None, overwrite=0):
    """Subsamples the passed physio. trigger file by a passed factor.

    The default sampling frequencies of the GE respiratory and cardiac
    physiology files are 25 and 100 Hz, respectively.

    ARGUMENTS:
        (input_file):  '/PATH/FILE'; file to subsample
        (subsample_factor):  INT; factor by which to subsample the input file

    RETURN:
        (subsamp_file_name): STR; output file name containing subsampled trigger times
    """
    print('-----Starting: resample_ppgtrig-----')

    #Check inputs
    if subsample_factor is None:
        print('No subsample_factor passed!')
        print('DEFAULTING TO 4!')
        subsample_factor = 4
    else:
        subsample_factor = int(subsample_factor)
        print('Subsample factor passed: '+str(subsample_factor))

    if input_file is None:
        print('No input file passed!')
        print('RETURNING 0!')
        return 0
    else:
        input_file = str(input_file)
        print('Input file passed: '+str(input_file))

    if os.path.exists(input_file) != 1:
        print('Input file cannot be found!')
        print('input_file: '+str(input_file))
        print('RETURNING 0!')
        return 0

    read_okay = os.access(input_file, os.R_OK)
    if read_okay != 1:
        print('Input file cannot be read from!')
        print('input_file: '+str(input_file))
        print('RETURNING 0!')
        return 0

    #Read in input PPG data
    ppgtrig_list = []

    with open(input_file, 'rt') as fp:
        for line in fp:
            ppgtrig_list.append(line)

    #Subsample ppg_list
    ppg_subsamp_list = []
    for element in ppgtrig_list:
        ppg_subsamp_list.append(int(round(float(element) / float(subsample_factor))))

    #Write subsampled PPG data file
    subsamp_file_name = str(input_file)+'_subsamp'+str(subsample_factor)

    if os.path.exists(subsamp_file_name):
        if overwrite:
            print('Overwrite set; deleting existing: {}'.format(subsamp_file_name))
        else:
            raise RuntimeError('Subsampled trigger file exists and overwrite not set!')

    with open(subsamp_file_name, 'w') as ppg_subsamp_file:
        for element in ppg_subsamp_list:
            ppg_subsamp_file.write(str(element)+'\n')

    print('Subsampled file written: '+str(subsamp_file_name))

    print('-----Finishing: resample_ppgtrig-----')
    
    return subsamp_file_name



def resample_data(input_file=None, subsample_factor=None, overwrite=0):
    """Subsamples the passed file by a passed factor.

    The default sampling frequencies of the GE respiratory and cardiac
    physiology files are 25 and 100, respectively.

    ARGUMENTS:
        (input_file):  '/PATH/FILE'; file to subsample
        (subsample_factor):  INT; factor by which to subsample the input file

    RETURN:
        (subsamp_file_name): STR; output file written containing subsampled data
    """
    print('-----Starting: resample_data-----')

    #Check inputs
    if subsample_factor is None:
        raise RuntimeError('No subsample_factor passed!')
    else:
        subsample_factor = int(subsample_factor)
        print('Subsample factor passed: '+str(subsample_factor))

    if input_file is None:
        raise RuntimeError('No input file passed!')
    else:
        input_file = str(input_file)
        print('Input file passed: '+str(input_file))

    if os.path.exists(input_file) != 1:
        raise RuntimeError('Input file cannot be found: {}'.format(input_file))

    read_okay = os.access(input_file, os.R_OK)
    if read_okay != 1:
        raise RuntimeError('Input file cannot be read from: {}'.format(input_file))

    #Read in input PPG data
    ppg_list = []

    with open(input_file, 'rt') as fp:
        for line in fp:
            ppg_list.append(line)

    #Subsample ppg_list
    count = subsample_factor-1
    ppg_subsamp_list = []
    for element in ppg_list:
        if count == (subsample_factor-1):
            ppg_subsamp_list.append(element)
            count = 0
        else:
            count = count+1

    #Write subsampled PPG data file
    subsamp_file_name = str(input_file)+'_subsamp'+str(subsample_factor)

    if os.path.exists(subsamp_file_name):
        if overwrite:
            print('Overwrite set; deleting existing: {}'.format(subsamp_file_name))
        else:
            raise RuntimeError('Subsampled data file exists and overwrite not set!')

    with open(subsamp_file_name, 'w') as ppg_subsamp_file:
        for element in ppg_subsamp_list:
            ppg_subsamp_file.write(element)

    print('Subsampled file written: '+str(subsamp_file_name))

    print('-----Finishing: resample_ppg-----')

    return subsamp_file_name



def convert_ge_to_pnm(ppg_input, resp_input, ppgsamprate, respsamprate, ppgtype, tr, trim_flag, cut_trs, delete_flag, output_file, overwrite=0):

    try:
        print('----convet_ge_to_pnm----')
        print('ppg_input: {}'.format(ppg_input))
        print('resp_input: {}'.format(resp_input))
        print('ppgsamprate: {}'.format(ppgsamprate))
        print('respsamprate: {}'.format(respsamprate))
        print('ppgtype: {}'.format(ppgtype))
        print('tr: {}'.format(tr))
        print('trim_flag: {}'.format(trim_flag))
        print('cut_trs: {}'.format(cut_trs))
        print('delete_flag: {}'.format(delete_flag))
        print('output_file: {}'.format(output_file))
        print('overwrite: {}'.format(overwrite))
        print('------------------------')
        ppg_input = str(ppg_input)
        resp_input = str(resp_input)
        ppgsamprate = float(ppgsamprate)
        respsamprate = float(respsamprate)
        ppgtype = str(ppgtype)
        tr = float(tr)
        trim_flag = int(trim_flag)
        cut_trs = int(cut_trs)
        delete_flag = int(delete_flag)
        output_file = str(output_file)
        overwrite = int(overwrite)

    except Exception as ex:
        print('Something is wrong with inputs to covert_ge_to_pnm!')
        print(ex)
        raise RuntimeError

    intermediate_files = []
    ##TODO: check inputs
    if not os.path.exists(ppg_input):
        raise RuntimeError('ppg_input file not found: {}'.format(ppg_input))
    if not os.path.exists(resp_input):
        raise RuntimeError('resp_input file not found: {}'.format(resp_input))

    working_ppg = ppg_input
    working_resp = resp_input

    #Resample the higher-Hz file to match the lower-Hz file
    if ppgsamprate > respsamprate:
        try:
            samprate = respsamprate
            subsamp_factor = int(ppgsamprate/respsamprate)
            #The PPG file might be triggers
            if ppgtype == 'Recorded Data':
                working_ppg = resample_data(input_file=working_ppg, subsample_factor=subsamp_factor, overwrite=overwrite)
            elif ppgtype == 'Triggers':
                working_ppg = resample_trig(input_file=working_ppg, subsample_factor=subsamp_factor, overwrite=overwrite)
            intermediate_files.append(working_ppg)
        except Exception as ex:
            print('Error resampling PPG file!')
            print(ex)
            raise RuntimeError
    elif respsamprate > ppgsamprate:
        try:
            samprate = ppgsamprate
            subsamp_factor = int(respsamprate/ppgsamprate)
            #PNM cannot accept a resp. data file of triggers, so it will always contain recorded data.
            working_resp = resample_data(input_file=working_resp, subsample_factor=subsamp_factor, overwrite=overwrite)
            intermediate_files.append(working_resp)
        except Exception as ex:
            print('Error resampling Resp. file!')
            print(ex)
            raise RuntimeError
    else:
        print('Sampling factors are equal, no resampling needed...')


    #Convert the cardiac trigger file to a list of 1s and 0s instead of a list
    #of time points.
    try:
        if ppgtype == 'Recorded Data':
            print('Cardiac input is recorded data...')
        elif ppgtype == 'Triggers':
            print('Cardiac input is trigger times...')
            print('Converting to a list of 1s and 0s...')
            working_ppg = convert_trig(input_file=working_ppg, example_data_file=working_resp, overwrite=overwrite)
    except Exception as ex:
        print('Error converting cardiac trigger file!')
        print(ex)
        raise RuntimeError

    intermediate_files.append(working_ppg)


    #Remove the first 30 seconds of each file, if required
    if trim_flag:
        try:
            print('Removing first 30 seconds from PPG file...')
            working_ppg = cut_thirty_seconds(samp_per_sec=samprate, input_file=working_ppg, overwrite=overwrite)
            working_resp = cut_thirty_seconds(samp_per_sec=samprate, input_file=working_resp, overwrite=overwrite)
        except Exception as ex:
            print('Error removing first 30 seconds of a file!')
            print(ex)
            raise RuntimeError

        intermediate_files.append(working_ppg)
        intermediate_files.append(working_resp)

    #Remove any additional early time points
    if cut_trs > 0:
        try:
            print('Removing initial TRs: {}'.format(cut_trs))
            working_ppg = cut_presteady_physio(cut_num=cut_trs, tr=tr, samp_per_sec=samprate, input_file=working_ppg, overwrite=overwrite)
            working_resp = cut_presteady_physio(cut_num=cut_trs, tr=tr, samp_per_sec=samprate, input_file=working_resp, overwrite=overwrite)
        except Exception as ex:
            print('Error removing initial TRs!')
            print(ex)
            raise RuntimeError

        intermediate_files.append(working_ppg)
        intermediate_files.append(working_resp)

    #Create scanner trigger file
    try:
        scan_trig_file = create_scantrig(physio_data_file=working_resp, sampling_rate=samprate, TR=tr, overwrite=overwrite)
        intermediate_files.append(scan_trig_file)
    except Exception as ex:
        print('Error creating scanner trigger file!')
        print(ex)
        raise RuntimeError

    #Merge data and write output
    try:
        combined_file = combine_physio(resp_file=working_resp, ppg_file=working_ppg, scantrig_file=scan_trig_file, ppgtype=ppgtype, output_file=output_file, overwrite=overwrite)
    except Exception as ex:
        print('Error creating combined output file!')
        print(ex)
        raise RuntimeError

    #Delete intermediate files if requested
    if delete_flag:
        print('Deleting intermediate files...')
        for element in intermediate_files:
            if os.path.exists(element):
                print('Deleting: {}'.format(element))
                os.remove(element)


class physio_gui():
    def __init__(self):
        self.window = tk.Tk()
        self.window.title('Write PNM Input')

        self.frame_files = tk.Frame()
        #Create frame of input file information
        self.ppg_file_lbl = tk.Label(self.frame_files, text='PPG File:')
        self.ppgfile_var = tk.StringVar(self.window)
        self.ppg_file_entry = tk.Entry(self.frame_files, textvariable=self.ppgfile_var)
        self.ppgdatatype_var = tk.StringVar(self.window)
        self.ppg_file_btn = tk.Button(self.frame_files, text='Open', command=self.gui_open_ppg)
        self.ppgdatatype_var.set('Recorded Data')
        ppgtypechoices = ['Recorded Data', 'Triggers']
        self.ppg_datatype_drop = tk.OptionMenu(self.frame_files, self.ppgdatatype_var, *ppgtypechoices)
        self.resp_file_lbl = tk.Label(self.frame_files, text='Resp. File:')
        self.respfile_var = tk.StringVar(self.window)
        self.resp_file_entry = tk.Entry(self.frame_files, textvariable=self.respfile_var)
        self.resp_file_btn = tk.Button(self.frame_files, text='Open', command=self.gui_open_resp)
        self.out_dir_lbl = tk.Label(self.frame_files, text='Output Dir.:')
        self.outdir_var = tk.StringVar(self.window)
        self.out_dir_entry = tk.Entry(self.frame_files, textvariable=self.outdir_var)
        self.out_dir_btn = tk.Button(self.frame_files, text='Select', command=self.gui_select_outdir)
        self.out_file_lbl = tk.Label(self.frame_files, text='Output Filename:')
        self.outfile_var = tk.StringVar(self.window)
        self.out_file_entry = tk.Entry(self.frame_files, textvariable=self.outfile_var)
        self.ppg_file_lbl.grid(row=0, column=0, sticky='w')
        self.ppg_file_entry.grid(row=0, column=1, sticky='w')
        self.ppg_file_btn.grid(row=0, column=2, sticky='w')
        self.ppg_datatype_drop.grid(row=0, column=3, sticky='w')
        self.resp_file_lbl.grid(row=1, column=0, sticky='w')
        self.resp_file_entry.grid(row=1, column=1, sticky='w')
        self.resp_file_btn.grid(row=1, column=2, sticky='w')
        self.out_dir_lbl.grid(row=2, column=0, sticky='w')
        self.out_dir_entry.grid(row=2, column=1, sticky='w')
        self.out_dir_btn.grid(row=2, column=2, sticky='w')
        self.out_file_lbl.grid(row=2, column=3, sticky='w')
        self.out_file_entry.grid(row=2, column=4, sticky='w')

        self.frame_samprates = tk.Frame()
        #Create frame for sampling rates
        self.ppg_samprate_lbl = tk.Label(self.frame_samprates, text='PPG Sampling Rate (Hz):')
        self.ppgsamprate_var = tk.StringVar(self.window)
        self.ppg_samprate_entry = tk.Entry(self.frame_samprates, textvariable=self.ppgsamprate_var)
        self.resp_samprate_lbl = tk.Label(self.frame_samprates, text='Resp. Sampling Rate (Hz):')
        self.respsamprate_var = tk.StringVar(self.window)
        self.resp_samprate_entry = tk.Entry(self.frame_samprates, textvariable=self.respsamprate_var)
        self.ppg_samprate_lbl.grid(row=0, column=0, sticky='w')
        self.ppg_samprate_entry.grid(row=0, column=1, sticky='w')
        self.resp_samprate_lbl.grid(row=1, column=0, sticky='w')
        self.resp_samprate_entry.grid(row=1, column=1, sticky='w')

        self.frame_options = tk.Frame()
        #Create frame for some additional options/inputs
        self.tr_lbl = tk.Label(self.frame_options, text='TR (sec):')
        self.tr_var = tk.StringVar()
        self.tr_entry = tk.Entry(self.frame_options, textvariable=self.tr_var)
        self.cut_var = tk.IntVar()
        self.cut_var.set(1)
        self.cut_cbtn = tk.Checkbutton(self.frame_options, text='Remove Pre-Scan 30 sec', variable=self.cut_var, onvalue=1, offvalue=0)
        self.delete_var = tk.IntVar()
        self.delete_var.set(1)
        self.delete_cbtn = tk.Checkbutton(self.frame_options, text='Delete Intermediate Files', variable=self.delete_var, onvalue=1, offvalue=0)
        self.dummytr_lbl = tk.Label(self.frame_options, text='Remove Initial TRs:')
        self.dummytr_var = tk.StringVar(self.window)
        self.dummytr_var.set('0')
        self.dummytr_entry = tk.Entry(self.frame_options, textvariable=self.dummytr_var)
        self.overwrite_var = tk.IntVar()
        self.overwrite_var.set(0)
        self.overwrite_cbtn = tk.Checkbutton(self.frame_options, text='Overwrite Existing Files', variable=self.overwrite_var, onvalue=1, offvalue=0)
        self.tr_lbl.grid(row=0, column=0, sticky='w')
        self.tr_entry.grid(row=0, column=1, sticky='w')
        self.cut_cbtn.grid(row=0, column=3, sticky='w')
        self.dummytr_lbl.grid(row=1, column=0, sticky='w')
        self.dummytr_entry.grid(row=1, column=1, sticky='w')
        self.delete_cbtn.grid(row=2, column=0, sticky='w')
        self.overwrite_cbtn.grid(row=2, column=1, sticky='w')

        self.frame_responses = tk.Frame()
        #Create frame for response buttons
        self.go_btn = tk.Button(self.frame_responses, text='Go', command=self.gui_run_btn)
        self.quit_btn = tk.Button(self.frame_responses, text='Quit', command=self.gui_quit_btn)
        self.go_btn.grid(row=0, column=0, sticky='w')
        self.quit_btn.grid(row=0, column=1, sticky='w')

        #Put all the frames into the GUI window
        self.frame_files.pack(fill=tk.X, expand=True)
        self.frame_samprates.pack(fill=tk.X, expand=True)
        self.frame_options.pack(fill=tk.X, expand=True)
        self.frame_responses.pack(fill=tk.X, expand=True)


    def run_loop(self):
        self.window.mainloop()


    def gui_open_ppg(self):
        print('Opening file selection dialogue...')
        #Open file-selection dialogue to get a file path and name
        ppg_filename = fd.askopenfilename(parent=self.window)
        self.ppgfile_var.set(ppg_filename)


    def gui_open_resp(self):
        print('Opening file selection dialogue...')
        #Open file-selection dialogue to get a file path and name
        resp_filename = fd.askopenfilename(parent=self.window)
        self.respfile_var.set(resp_filename)

    def gui_select_outdir(self):
        print('Opening directory selection dialogue...')
        output_dir = fd.askdirectory(parent=self.window)
        self.outdir_var.set(output_dir)


    def gui_run_btn(self):
        intermediate_files = []
        delete_intermediates = 1
        #Pull variables from the GUI and check them
        try:
            ppg_input = str(self.ppgfile_var.get())
            resp_input = str(self.respfile_var.get())
            ppgsamprate = float(self.ppgsamprate_var.get())
            respsamprate = float(self.respsamprate_var.get())
            trim_flag = int(self.cut_var.get())
            cut_trs = int(self.dummytr_var.get())
            tr = float(self.tr_var.get())
            ppgtype = str(self.ppgdatatype_var.get())
            delete_flag = int(self.delete_var.get())
            outdir = str(self.outdir_var.get())
            outfile = str(self.outfile_var.get())
            output_file = os.path.join(outdir, outfile)
            overwrite = int(self.overwrite_var.get())

            if not os.path.exists(ppg_input):
                raise RuntimeError('PPG file not found: {}'.format(ppg_input))

            if not os.path.exists(resp_input):
                raise RuntimeError('Resp. file not found: {}'.format(resp_input))

            if trim_flag not in [1,0]:
                raise RuntimeError('Something went wrong with the trim flag: {}'.format(trim_flag))

            if not os.path.exists(outdir):
                raise RuntimeError('Output dir not found: {}'.format(outdir))

            working_ppg = ppg_input
            working_resp = resp_input

        except Exception as ex:
            print('Something went wrong checking variables')
            print(ex)
            self.window.destroy()
            raise RuntimeError

        try:
            #Run the conversion
            convert_ge_to_pnm(ppg_input, resp_input, ppgsamprate, respsamprate, ppgtype, tr, trim_flag, cut_trs, delete_flag, output_file, overwrite=overwrite)
        except Exception as ex:
            print('Error occurred during conversion: {}'.format(ex))
            self.window.destroy()
            raise RuntimeError

        print('Conversion successful!')



    def gui_quit_btn(self):

        #Close the GUI
        self.window.destroy()


if __name__ == '__main__':
    obj = physio_gui()
    obj.run_loop()
