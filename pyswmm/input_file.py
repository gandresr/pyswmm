import fileinput, re
import os

def set_sim_time(swmm_file, start_date, start_time, last_date, last_time):
    for line in fileinput.input(swmm_file, inplace=True):
        if 'REPORT_START_DATE' in line:
            print('REPORT_START_DATE ' + start_date)
            continue
        if 'REPORT_START_TIME' in line:
            print('REPORT_START_TIME ' + start_time)
            continue
        if 'START_DATE' in line:
            print('START_DATE ' + start_date)
            continue
        if 'START_TIME' in line:
            print('START_TIME ' + start_time)
            continue
        if 'END_DATE' in line:
            print('END_DATE ' + last_date)
            continue
        if 'END_TIME' in line:
            print('END_TIME ' + last_time)
            continue
        print(line.strip('\n'))


def set_rainfile(swmm_file, rain_gage, rain_file, station, units='MM'):
    in_rain_gages = False
    rain_file = rain_file.replace('/','\\')
    fl = get_first_last(rain_file)
    sd = fl[0].decode('utf-8').strip().split()
    ld = fl[1].decode('utf-8').strip().split()
    start_date = sd[2] + '/' + sd[3] + '/' + sd[1]
    start_time = ':'.join(sd[4:6]) + ':00'
    last_date =  ld[2] + '/' + ld[3] + '/' + ld[1]
    last_time = ':'.join(ld[4:6]) + ':00'

    set_sim_time(swmm_file, start_date, start_time, last_date, last_time)

    for line in fileinput.input(swmm_file, inplace=True):
        if in_rain_gages:
            l = re.findall(r'(?:[^\s,"]|"(?:\\.|[^"])*")+', line)
            if str(rain_gage) == l[0]:
                l[5] = '"' + rain_file + '"'
                l[6] = station
                l[7] = units
                print('\t'.join(l))
                in_rain_gages = False
                continue
        if '[RAINGAGES]' in line:
            in_rain_gages = True
        print(line.strip('\n'))

def get_first_last(ffile):
    with open(ffile, "rb") as f:
        first = f.readline()        # Read the first line.
        f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
        while f.read(1) != b"\n":   # Until EOL is found...
            f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
        last = f.readline()         # Read last line.
    return (first, last)

# swmm_file = 'C:/Users/ga.riano949/Documents/GitHub/mpc/Paper ACC 2018/models/10 states/model_pumps.inp'
# rain_gage = 1
# rain_file = 'C:/Users/ga.riano949/Documents/GitHub/mpc/Paper ACC 2018/NOAA data/precipitation_data/COOP_021574.dat'
# set_rainfile(swmm_file, rain_gage, rain_file, 'COOP_021574')
