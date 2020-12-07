#FROM python:3

FROM continuumio/miniconda3
RUN conda create -n env python=3.8
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
RUN conda install -c anaconda pip

RUN mkdir src
WORKDIR src/
COPY . .

RUN pip install -r requirements.txt
RUN conda install jupyter
RUN python3 transmissionTools.py
WORKDIR /src/sample_data
RUN conda install -c bioconda parsnp -y
RUN conda install -c bioconda harvesttools -y
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install parsnp=1.5 -y
# RUN wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz
# RUN tar -xvf parsnp-Linux64-v1.2.tar.gz
# Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "notebook", "--port=8899", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

# RUN pip3 install -r requirements.txt
# RUN pip3 install jupyter
# RUN python3 transmissionTools.py
# WORKDIR /src/sample_data
# RUN wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz
# RUN tar -xvf parsnp-Linux64-v1.2.tar.gz
# # Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
# ENV TINI_VERSION v0.6.0
# ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
# RUN chmod +x /usr/bin/tini
# ENTRYPOINT ["/usr/bin/tini", "--"]
# CMD ["jupyter", "notebook", "--port=8899", "--no-browser", "--ip=0.0.0.0", "--allow-root"]