wget --no-check-certificate -nc -O /opt/wpg.zip https://github.com/samoylv/WPG/archive/develop.zip
unzip -o /opt/wpg.zip -d /opt 
cd /opt/WPG-develop/ && make all
# rm -rf /opt/WPG-develop/build