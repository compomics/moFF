import pytest

@pytest.fixture
def data_mbr():
    import pandas
    return  pandas.read_csv('absence_peak_data/mbr_output/B002413_Ap_22cm_Yeast_171215184201_match.txt',sep="\t")

@pytest.fixture
def data_mbr_2():
    import pandas
    d1=pandas.read_csv('absence_peak_data/mbr_output/B002413_Ap_22cm_Yeast_171215184201_match.txt', sep="\t")
    d2=pandas.read_csv('absence_peak_data/mbr_output/B002417_Ap_22cm_iRT_PRC-Hans_equimolar_100fmol_match.txt', sep="\t")
    d3=pandas.read_csv('absence_peak_data/mbr_output/B002419_Ap_22cm_iRT_PRC-Hans_equimolar_100fmol_inYeast_match.txt',sep="\t")
    d4=pandas.read_csv('absence_peak_data/mbr_output/B002421_Ap_22cm_iRT_PRC-Hans_equimolar_100fmol_match.txt',sep="\t")
    return ( d1,d2,d3,d4)



def test_mbr_iRT(data_mbr):
    iRT_match = data_mbr[data_mbr['prot']=='Biognosys']
    assert (iRT_match['peptide'].unique().shape[0] ==  11 ), 'not right unique iRT peptide matched'
    assert (iRT_match['rt'].mean()/ 60 == 36.38888249800042), 'not right mean iRT rt matched peptide'
    assert (iRT_match.shape[0]== 15) ,'not right size '


def test_general_mbr(data_mbr_2):
    assert (data_mbr_2[0].shape[0] ==  8356 ), 'not right unique iRT peptide matched'
    assert (data_mbr_2[1].shape[0] ==  8427 ), 'not right unique iRT peptide matched'
    assert (data_mbr_2[2].shape[0] ==  8264 ), 'not right unique iRT peptide matched'
    assert (data_mbr_2[3].shape[0] ==  8424 ), 'not right unique iRT peptide matched'
