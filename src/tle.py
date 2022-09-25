from math import pi, sqrt
import string

class TLE:

        
    def _checksum(self, line):
        # stolen from https://space.stackexchange.com/questions/5358/what-does-check-sum-tle-mean
        L = line.strip()
        cksum = 0
        for i in range(68):
            c = L[i]
            if c == ' ' or c == '.' or c == '+' or c in string.ascii_letters:
                continue
            elif c == '-':
                cksum = cksum + 1
            else:
                cksum = cksum + int(c)

        cksum %= 10

        return cksum
    
    def __init__(self, name: str, raw_line_1: str, raw_line_2: str, perform_check=False, body_mu = 3.986004418 * (10 ** 14)):
        self.name = name.strip()
        # split_line_1 = raw_line_1.strip().split(' ')
        # split_line_2 = raw_line_2.strip().split(' ')
        # print(name)
        # print(raw_line_1)
        # print(raw_line_2)
        # print()
        # split_line_1 = list(filter(lambda x: x!='', split_line_1))
        # split_line_2 = list(filter(lambda x: x!='', split_line_2))
        # print(split_line_1)
        # print(split_line_2)

        self.body_mu = body_mu # default to Earth's mu

        #first line
        self.catalogue_num = raw_line_1[2:7]
        self.classification = raw_line_1[7:8]
        self.international_designator_year = raw_line_1[9:11]
        self.international_designator_launch_num = raw_line_1[11:14]
        self.international_designator_piece = raw_line_1[14:17].strip()
        self.epoch_year = raw_line_1[18:20]
        self.epoch = raw_line_1[20:32]
        self.d_mean_motion =  raw_line_1[33:43].strip()
        self.d2_mean_motion =  raw_line_1[44:52].strip()
        self.perturbation_term =  raw_line_1[53:61].strip()
        self.ephemeris_type = raw_line_1[62]
        self.element_set_number =  raw_line_1[64:68].strip()
        self.check_sum_line_1 =  int(raw_line_1[68])

        #second line
        if self.catalogue_num != raw_line_2[2:7].strip():
            raise Exception("Catalogue Number Mismatched Between Two Lines: {} vs {}".format(self.catalogue_num, raw_line_2[2:7].strip()))
            
        self.inclination = float(raw_line_2[8:16].strip())
        self.raan = float(raw_line_2[17:25].strip())
        self.eccentricity = float(raw_line_2[26:33].strip())
        self.arg_pe = float(raw_line_2[34:42].strip())
        self.mean_anomaly = float(raw_line_2[43:51].strip())
        self.mean_motion = float(raw_line_2[52:63].strip())
        self.revolutions = int(raw_line_2[63:68].strip())
        self.check_sum_line_2 = int(raw_line_2[68].strip())

        #check sum
        if perform_check:
            if self._checksum(raw_line_1) != self.check_sum_line_1:
                raise Exception("line 1 checksum doesn't match!!")
            if self._checksum(raw_line_2) != self.check_sum_line_2:
                raise Exception("line 2 checksum doesn't match!!")

        #pre computation:
        self.n_rad_per_sec = self.mean_motion * (2 * pi / 86400)
        self.period = 2 * self.n_rad_per_sec * pi
        # self.semi_major_axis = ((self.body_mu * self.period ** 2) / (4 * pi ** 2)) ** (1. / 3) # we got this, but ioncorrect
        # print(self.name)
        # print(self.body_mu)
        # print(self.semi_major_axis)
        
        self.semi_major_axis = (self.body_mu ** (1./ 3)) / (((2 * self.mean_motion * pi)) / 86400  ) ** (2./3) # correct
        # print(self.mean_motion)
        # print(self.semi_major_axis)
        # print()

        
    def get_inclination_rad(self):
        return self.inclination * pi / 180

    def get_semi_major_axis(self):
        return self.semi_major_axis





