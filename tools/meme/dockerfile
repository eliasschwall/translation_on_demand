# useing the tidyverse rstudio image as our base
FROM rocker/tidyverse:4.4

# for the meme suite analysis we need to add some depenencies
RUN cpan XML::Parser XML::Simple HTML::Template HTML::TreeBuilder JSON Sys::Info
RUN sudo apt-get update
RUN sudo apt-get install apt-utils
RUN sudo apt-get install -y imagemagick ghostscript libxslt1-dev libglpk40

# now we install meme itself
RUN mkdir meme_installation && \
    cd meme_installation && \
    wget https://meme-suite.org/meme/meme-software/5.5.7/meme-5.5.7.tar.gz && \
    tar zxf meme-5.5.7.tar.gz && \
    rm meme-5.5.7.tar.gz && \
    cd meme-5.5.7 && \
    ./configure --prefix=$(pwd)/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && \
    make test && \
    make install
ENV PATH=meme_installation/meme-5.5.7/meme/bin:meme_installation/meme-5.5.7/meme/libexec/meme-5.5.7:$PATH

# next we install the memes R package needed for the analysis
RUN R -e "BiocManager::install('memes')"
# and we check if R the memes package finds the meme installation
RUN R -e "memes::check_meme_install('/meme_installation/meme-5.5.7/meme/bin/')"
# lastly we install other R packages we will need for the analysis 
RUN R -e "BiocManager::install(c('DESeq2', 'biomaRt', 'Biostrings', 'clusterProfiler', 'org.Mm.eg.db', 'org.Hs.eg.db', 'ComplexHeatmap'))"


