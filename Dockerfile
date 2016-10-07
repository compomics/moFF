FROM continuumio/anaconda

RUN apt-get update
RUN conda install numpy pandas scikit-learn 
RUN pip install argparse pymzML
RUN apt-get install git
RUN git clone https://github.com/compomics/moff /moff
WORKDIR /moff
