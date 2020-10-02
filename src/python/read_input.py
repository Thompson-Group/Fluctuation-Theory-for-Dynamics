class user_input:

    def __init__(self,inputfile):
        ip = open(inputfile)

        self.molec=[]
        count = 0
        molshift=0
        for i, line in enumerate(ip):
            if i == 6:
                self.prog=str(line.strip())
            if i == 10:
                # 11th line
                self.num_molecs=int(line)
                molshift=int(line)
            if i >= 12 and count < molshift:
                count += 1
                self.molec.append(str(line.strip()))
            if i == 12+molshift+1:
                self.start_config=int(line.strip())
            if i == 12+molshift+3:
                self.end_config=int(line.strip())
            if i == 12+molshift+5:
                self.sep_config=int(line.strip())
            if i == 12+molshift+7:
                self.timestep=float(line.strip())
            if i == 12+molshift+9:
                self.num_times=int(line.strip())
            if i == 12+molshift+11:
                self.nblocks=int(line.strip())
            if i == 12+molshift+13:
                self.segsplit=int(line.strip())
            if i == 12+molshift+15:
                self.num_files=int(line.strip())
            if i == 12+molshift+17:
                self.nve_length=int(line.strip())
            if i == 12+molshift+19:
                self.num_rpj=int(line.strip())
            if i == 12+molshift+21:
                self.cab=str(line.strip())
            if i == 12+molshift+23:
                self.constraint=float(line.strip())
        ip.close()
