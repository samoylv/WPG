To run jupyter notebook

1) Run ```docker run -it --rm -p 8088:8888 -v ~/tmp/data:/data buzmakov/wpg3:openmp```
where ~/tmp/data - your local folder to save data. This will take some time do download wpg-docker contatiner for the fitst time (~ 2 GB)

2) Open browser on http://localhost:8088

WPG tutorials can be found here /opt/WPG-feature-openmp/samples/Tutorials

To run your python script put it in ~/tmp/data (or other directory mount as /data in starting docker command) and run ```docker run -it --rm   -u `id -u`:`id -g` -v ~/tmp/data:/data buzmakov/wpg3:openmp python /data/my_script``` . WPG is already in python path.
