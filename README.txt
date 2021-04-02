To connect to the GMU omics server:
ssh kjastrze@omics.gmu.edu 
input password

To upload an app update
1.  Zip the folder
zip -r APP.zip APP
copy the zipped folder to the swerver
scp APP.zip kjastrze@omics.gmu.edu:/home/kjastrze


ONCE ON THE server

unzip the newly uploaded app
unzip app.zip

cd APP

to install the required libraries:
pip install --upgrade --user pip
pip install -U -r requirements.txt
to run the app if it goes down:

# kill the current process
ps -ef|grep python

# get rid of the previous output file
rm nohup.out

nohup gunicorn ngampp:server -b :8080 &