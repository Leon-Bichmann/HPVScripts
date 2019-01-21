from Bio.SeqUtils import six_frame_translations

seq=seq.upper()
sft=six_frame_translations(seq)
sft_split=sft.split('\n')

frame1=''.join([sft_split[r] for r in range(5,475,10)]).replace(' ','').replace('*','')
frame2=''.join([sft_split[r] for r in range(6,475,10)]).replace(' ','').replace('*','')
frame3=''.join([sft_split[r] for r in range(7,475,10)]).replace(' ','').replace('*','')
frame4=''.join([sft_split[r] for r in range(10,475,10)]).replace(' ','').replace('*','')
frame5=''.join([sft_split[r] for r in range(11,475,10)]).replace(' ','').replace('*','')
frame6=''.join([sft_split[r] for r in range(12,475,10)]).replace(' ','').replace('*','')