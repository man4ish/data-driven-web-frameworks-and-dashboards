FROM perl:5.32

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libmariadb-dev \
    mariadb-client \
    pkg-config \
    libssl-dev \
    libexpat1-dev \
    zlib1g-dev \
    perl \
    cpanminus \
    default-jre \
    git \
    wget \
    curl \
    python3 \
    python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

    #ENV MYSQL_LIB=/usr/lib/x86_64-linux-gnu
    #ENV MYSQL_INCLUDE=/usr/include/mysql
    
    #RUN cpanm --verbose Bio::Perl \
    #    && cpanm --verbose DBI JSON::XS \
    #    && cpanm --force --verbose DBD::mysql

    WORKDIR /app

    COPY . /app
    
    #ENV PERL5LIB=/app/bin/ensembl/modules:/app/bin/ensembl-variation/modules:/app/bin/ensembl-funcgen/modules:/app/bin/ensembl-external/modules:/app/bin/ensembl-vep-release-114/modules
    
    EXPOSE 8000
    
    # Install Python dependencies if requirements.txt is present
    RUN if [ -f requirements.txt ]; then pip3 install --no-cache-dir -r requirements.txt; fi

    CMD ["python3", "app.py"]    