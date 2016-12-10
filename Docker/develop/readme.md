To run jupyter notebook

1) Run ```docker run -it --rm -u `id -u`:`id -g` -p 8088:8888 -v ~/tmp/data:/data buzmakov/wpg3```
where ~/tmp/data - your local folder to save data

2) Open browser on http://localhost:8088

WPG tutorials can be found here /opt/WPG-develop/samples/Tutorials

To run your python script put it in ~/tmp/data (or other directory mount as /data in starting docker command) and run ```docker run -it --rm  -p 6611:8888 -v ~/tmp/data:/data buzmakov/wpg3 python /data/my_script``` . WPG is already in python path.
