'''
Reads data from Agilent MSO6054A
Tested on Windows 10
'''

import visa
import time
import argparse
# start of Read_Agilent_01


__description__ = 'Read_Agilent_01'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--output_file', type=str, required=False, default='temp.txt', help='ascii output file')
PARSER.add_argument('-n', '--num_events',type=int, required=False, default=10, help='Number of events to analyze')

def main(**kwargs):

        file  = open(kwargs['output_file'], "w")

        rm = visa.ResourceManager()
        MSO6054A = rm.open_resource('USB0::0x0957::0x1744::MY47180105::0::INSTR')
        #selezionare il SETUP desiderato:
        MSO6054A.write(':RECall:SETup:STARt "%s"' % ('/usb0/setup/s10'))
        #MSO6054A.write(':WAVeform:POINts:MODE %s' % ('RAW'))
        MSO6054A.write(':WAVeform:POINts %s,%d' % ('RAW', 2000))
        MSO6054A.write(':WAVeform:FORMat %s' % ('ASCii'))
        
        num_events = kwargs['num_events']
        
        for i in range(0,num_events):
                ##################
                # START SEQUENCE #
                ##################
                MSO6054A.write(':RUN')
                MSO6054A.write(':STOP')
                rm = visa.ResourceManager()
                temp_values = MSO6054A.query_ascii_values(':WAVeform:XORigin?')
                XORigin = temp_values[0]

                temp_values = MSO6054A.query_ascii_values(':WAVeform:XINCrement?')
                XINCrement = temp_values[0]

                temp_values = MSO6054A.query_ascii_values(':WAVeform:YORigin?')
                YORigin = temp_values[0]

                temp_values = MSO6054A.query_ascii_values(':WAVeform:POINts?')
                points = int(temp_values[0])

                temp_values = MSO6054A.query_binary_values(':WAVeform:DATA?','s',False)
                binary_block_data = temp_values[0].decode().split(',')
                
                
                
                if(i%10==0):
                        print('Reading Event '+str(i))
                
                ################
                # END SEQUENCE #
                ################
                
                file.write('NewWaveform\n')
                file.write(str(points))
                file.write('\nTime\tVolt\n')
                x=0.
                y=0.
                for i in range (0, points):
                        x = XORigin + i * XINCrement
                        y = binary_block_data[i]
                        file.write(str(x) + '\t' + str(y) + '\n')


        MSO6054A.close()
        rm.close()
        file.close() 

if __name__ == '__main__':
        args = PARSER.parse_args()
        main(**args.__dict__)

