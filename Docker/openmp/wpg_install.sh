wget --no-check-certificate -nc -O /opt/wpg.zip https://github.com/samoylv/WPG/archive/feature/openmp.zip
unzip -o /opt/wpg.zip -d /opt 
cd /opt/WPG-feature-openmp/ && make all
# rm -rf /opt/WPG-develop/build