from  conda/miniconda3

RUN apt-get update
RUN conda install python=3.6
RUN conda install numpy pandas scikit-learn
RUN conda install -c conda-forge pynumpress
RUN conda install simplejson
RUN pip install  pymzML pyteomics brain-isotopic-distribution
RUN apt-get -y install git
RUN echo "deb http://download.mono-project.com/repo/debian wheezy-apache24-compat main" | tee -a /etc/apt/sources.list.d/mono-xamarin.list
RUN apt-get update
RUN apt-get install -y mono-complete
RUN git clone  -b master  --single-branch https://github.com/compomics/moff /moFF
WORKDIR /moFF

