FROM centos:centos7
RUN yum -y clean all
RUN yum -y update

RUN yum -y install epel-release
RUN yum -y install gcc gcc-c++ python-devel
RUN yum -y install python-matplotlib numpy scipy

ADD calculate_phist.py /calculate_phist.py
ADD cross_multinomial.py /cross_multinomial.py

ENTRYPOINT ["python", "/cross_multinomial.py"]