# set base image (host OS)
FROM python:3.9

# set the working directory in the container
WORKDIR /igem_tau_2021

# copy the dependencies file to the working directory
COPY requirements.txt .

# install dependencies
RUN pip install -r requirements.txt

# install meme
COPY meme-5.4.1.tar.gz .
RUN tar -xzf meme-5.4.1.tar.gz && cd meme-5.4.1 && \
     ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && cd scripts && cpan XML::Parser::Expat && cd .. && make install && cd ..

ENV PATH /root/meme/bin:/root/meme/libexec/meme-5.4.1:$PATH

# Install communique
COPY GUI ./GUI
COPY modules ./modules
COPY promoters_not_for_user ./promoters_not_for_user
COPY static ./static
COPY templates ./templates
COPY flaskgui.py .
COPY main.py .

ENV PYTHONPATH /igem_tau_2021
