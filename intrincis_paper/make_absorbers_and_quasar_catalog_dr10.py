import numpy

"""
This file is a substitute to make_absorbers_and_quasar_catalog_dr9 file.
It will use Britt's absorber catalog for CIV absorber and new spAll for Dr10.
"""


def load_britt_file(absorber_file):
    names = ['qsonum','sdssj','zS10','zHW10','ra_degdec_deg','plate','fiber','mjd','imag','flux_20cm','snr_20cm',
        'FIRST_targ','FIRST_sep','BALflag1','BALflag2','BALYA','sysnum','zsys','grade','beta','lamlo','lamhi',
        'wlim_1216_1241','wlim_1241_1400','wlim_1400_1550','wlim_1550_1909','wlim_1909_2799',
        'wlim_2799_3969','wlim_3969_8200','wlim_8200_9200','NVa','NVa_err','OI','OI_err','SiIIa',
        'SiIIa_err','CII','CII_err','SiIVa','SiIVa_err','SiIVb','SiIVb_err','SiIIb','SiIIb_err',
        'CIVa','CIVa_err','CIVb','CIVb_err','CI','CI_err','AlII','AlII_err','AlIIIa','AlIIIa_err','AlIIIb','AlIIIb_err',
        'ZnII','ZnII_err','CrII','CrII_err','FeIIa','FeIIa_err','FeIIb','FeIIb_err','FeIIc','FeIIc_err','FeIId','FeIId_err',
        'MnII','MnII_err','FeIIe','FeIIe_err','MgIIa','MgIIa_err','MgIIb','MgIIb_err','MgI','MgI_err',
        'TiII','TiII_err','CaIIa','CaIIa_err','CaIIb','CaIIb_err','NaI','NaI_err']
    #61 column, T-> text(3), M-R-> Int, G-I -> Int, B->Text(25), A-> Int
    formats = [numpy.int, 'S25',]
    formats.extend(4 * [numpy.float64,])
    formats.extend(3 * [numpy.int,])
    formats.extend(3 * [numpy.float64,])
    formats.extend(7 * [numpy.int,])
    formats.extend(['S4',])
    formats.extend(66 * [numpy.float64,])

    dtype={'names':names, 'formats': formats}

    return numpy.loadtxt(absorber_file, skiprows=1, dtype=dtype)

def main():
    test_load_britt_file()

def test_load_britt_file():
    filename='/home/vikas/intrinsic_code/intrinsic/Britts_022113_shorten.dat'
    data = load_britt_file(filename)
    print data


if __name__=='__main__':
    main()
