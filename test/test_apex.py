import pytest


@pytest.fixture
def data_apex():
    import pandas
    return  pandas.read_csv('test/B002417_Ap_22cm_iRT_PRC-Hans_equimolar_100fmol_moff_result.txt',sep="\t")


def test_apex(data_apex):
    assert all([a == b for a, b in zip(data_apex.columns, ['peptide', 'mod_peptide', 'prot', 'Type', 'Raw file', 'Experiment',
       'mz', 'charge', 'another m/z', 'mass', 'rt', 'PEP', 'Reverse',
       'intensity', 'rt_peak', 'lwhm', 'rwhm', '5p_noise', '10p_noise', 'SNR',
       'log_L_R', 'log_int'])]),"wrong  column names"
    assert not(data_apex.iloc[:, 13].isnull().all()), 'missing value on intensity '
    assert not(data_apex.iloc[:, 14].isnull().all()), 'missing value on rt_peak'
    assert not(data_apex.iloc[:, 15].isnull().all()), 'missing value on  lwhm'
    assert not(data_apex.iloc[:, 16].isnull().all()), 'missing value on rwhm'
    assert not(data_apex.iloc[:, 17].isnull().all()), 'missing value on 5p_noise'
    assert not(data_apex.iloc[:, 18].isnull().all()), 'missing value on 10p_noise'
    assert not(data_apex.iloc[:, 19].isnull().all()), 'missing value on SNR'
    assert not(data_apex.iloc[:, 20].isnull().all()), 'missing value on log_L_R'
    assert not(data_apex.iloc[:, 21].isnull().all()), 'missing value on log_int'
    assert data_apex.shape[0]== 106 , "wrong data size"
    assert data_apex[data_apex.log_L_R == -1].shape[0] ==10," worng number of record with log_L_R "
    assert data_apex.log_int.mean() == 19.892951681679005, "wrong log_int mean "