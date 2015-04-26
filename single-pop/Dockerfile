FROM centos:centos7
RUN yum -y clean all
RUN yum -y update

RUN yum -y install epel-release
RUN yum -y install gcc gcc-c++ python-devel
RUN yum -y install python-matplotlib numpy scipy

ADD calculate_phist.py /calculate_phist.py
ADD single_pop.py /single_pop.py

ENTRYPOINT ["python", "/single_pop.py"]