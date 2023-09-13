set -uex

#
# Input files are obtained from Revigo
#

#sexual_up:
python ~/bin/CirGO.py -inputFile sexual_up_BP_revigo.csv -outputFile sexual_up_BP_cirgo.svg -legend 'Name and Proportion of Biological Process (inner ring)'

python ~/bin/CirGO.py -inputFile sexual_up_MF_revigo.csv -outputFile sexual_up_MF_cirgo.svg -legend 'Name and Proportion of Molecular Function (inner ring)'

python ~/bin/CirGO.py -inputFile sexual_up_CC_revigo.csv -outputFile sexual_up_CC_cirgo.svg -legend 'Name and Proportion of Cellular Component (inner ring)'

#sexual_down:
python ~/bin/CirGO.py -inputFile sexual_down_BP_revigo.csv -outputFile sexual_down_BP_cirgo.svg -legend 'Name and Proportion of Biological Process (inner ring)'

python ~/bin/CirGO.py -inputFile sexual_down_MF_revigo.csv -outputFile sexual_down_MF_cirgo.svg -legend 'Name and Proportion of Molecular Function(inner ring)'

#asexual_down:
python ~/bin/CirGO.py -inputFile asexual_down_BP_revigo.csv -outputFile asexual_down_BP_cirgo.svg -legend 'Name and Proportion of Biological Process (inner ring)'

python ~/bin/CirGO.py -inputFile asexual_down_MF_revigo.csv -outputFile asexual_down_MF_cirgo.svg -legend 'Name and Proportion of Molecular Function (inner ring)'

python ~/bin/CirGO.py -inputFile asexual_down_CC_revigo.csv -outputFile asexual_down_CC_cirgo.svg -legend 'Name and Proportion of Cellular Component (inner ring)'
