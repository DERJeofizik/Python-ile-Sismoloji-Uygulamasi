#%matplotlib inline
import os, glob
# Veri Okunmasi
from obspy import read
st = read('2014-03-10T05_17_30.008400Z.PB.B045.EH1.SAC',format = 'SAC')
st += read('2014-03-10T05_17_30.008400Z.PB.B045.EH2.SAC',format = 'SAC')
st += read('2014-03-10T05_17_30.008400Z.PB.B045.EHZ.SAC',format = 'SAC')

# Veri Icerigi
print(st[0].stats)
print(st[1].stats)
print(st[2].stats)

# Deprem Dalgasinin Cizdirilmesi
st.plot()

# Filtreleme
st_filt = st.copy()
st_filt.detrend()
st_filt[0].filter("bandpass", freqmin=0.1, freqmax=10)
st_filt.plot()

# Kesme
st_kesik0 = st.copy()
st_kesik0[0].trim(st[0].stats.starttime + 0,st[0].stats.starttime + 60)
st_kesik0.plot()
st_kesik1 = st.copy()
st_kesik1[0].trim(st[0].stats.starttime + 60,st[0].stats.starttime + 120)
st_kesik1.plot()

# Birlestirme
st_birlestik = st_kesik0
st_birlestik += st_kesik1
st_birlestik.merge(method=1,fill_value = 0)
st_birlestik.plot()

''' Ivme = cm/s^2
Hiz = cm/s
Yer Degistirme = cm '''
# Integral Alma
st_hiz = st.copy()
st_hiz.filter("bandpass", freqmin=0.1, freqmax=10)
st_hiz.integrate(method='cumtrapz')
st_hiz.plot()
#Turev Alma
st_ivme = st_hiz.copy()
st_ivme.differentiate(method='gradient')
st_ivme.plot()

# Deprem Ilk Varis Belirlenmesi
from obspy.signal.trigger import plot_trigger, classic_sta_lta
df = st[0].stats.sampling_rate
cft = classic_sta_lta(st[0].data, int(2 * df), int(6 * df))
plot_trigger(st[0], cft, 2.95, 0.2,show=True)

# Spectrogram
st[0].spectrogram(log=True, title=str(st[0].stats.station) + ' ' + str(st[0].stats.channel) +  ' Spectrogram Grafigi ',show=True)

# Fourier Transformu
import scipy.fftpack
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
st_alcakgecis = st.copy()
st_alcakgecis[0].filter('lowpass',freq=3, corners=2, zerophase=True)
dt = st[0].stats.delta
npts = st[0].stats.npts
yf = scipy.fftpack.fft(st_alcakgecis[0])
xf = np.linspace(0.0, 1.0/(2.0*dt), npts/2)
fig, ax = plt.subplots()
ax.plot(xf, 2.0/npts * np.abs(yf[:npts//2]))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frekanslar')
plt.ylabel('Buyukluk')
fname = str(st[0].stats.starttime.year) + str(st[0].stats.starttime.month) + str(st[0].stats.starttime.day) + '_fft.jpg'
savefig(fname, dpi=50, bbox_inches='tight')
plt.show()

# Kuzey-Dogu Istasyonlarini Radyal ve Transverse'e Dondurme
from obspy.signal.rotate import rotate_ne_rt
# Baslangic Bitis Farkinin Ortadan Kaldirilmasi
t_bas_fark = st[0].stats.starttime - st[1].stats.starttime
t_son_fark = st[0].stats.endtime - st[1].stats.endtime
if t_bas_fark < 0:
  st[0].trim(st[0].stats.starttime,st[1].stats.endtime)
elif t_bas_fark > 0:
  st[1].trim(st[0].stats.starttime,st[1].stats.endtime)
if t_son_fark < 0:
  st[1].trim(st[1].stats.starttime,st[0].stats.endtime)
elif t_son_fark > 0:
  st[0].trim(st[0].stats.starttime,st[1].stats.endtime)
trace_kuzey = st[0].copy()
trace_dogu = st[1].copy()
rt = rotate_ne_rt(trace_kuzey.data,trace_dogu.data, st[0].stats.sac.baz)
#Cizdirme
plt.plot(rt[1])
plt.plot(st[1])
plt.show()

# Beach Ball Cizimi
from obspy.imaging.beachball import beachball
# Gosterim: Strike,Dip,Rake
mt = [228, 79, -2]
beachball(mt, size=200, linewidth=2, facecolor='r')

# Harita Cizimi
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
# Harita Altilgi
m = Basemap(projection='cyl', llcrnrlon = -130, llcrnrlat = 37, urcrnrlon = -120, urcrnrlat = 45,
            lat_0=42, lon_0=-120, resolution='i')
# Enlem Boylam Cizgileri
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-90., 90., 3), labels=[1, 1, 0, 0], fmt="%.2f",dashes=[2, 2])
m.drawmeridians(np.arange(-180., 180., 3.), labels=[0, 0, 1, 1], fmt="%.2f",dashes=[2, 2])
m.drawmapboundary()
# Istasyon
lons = [st[0].stats.sac.stlo]
lats = [st[0].stats.sac.stla]
names = [str(st[0].stats.station)]
x, y = m(lons, lats)
m.scatter(x, y, 200, color="r", marker="v", edgecolor="k", zorder=3)
for i in range(len(names)):
    plt.text(x[i], y[i], names[i], va="top", family="monospace", weight="bold")
m.plot(x, y, 'vk', markersize=5)
# Beach Ball
x, y = m(st[0].stats.sac.evlo, st[0].stats.sac.evla)
focmecs = [228, 79, -2]
ax = plt.gca()
b = beach(focmecs, xy=(x, y), width=0.5, linewidth=1, alpha=0.85)
b.set_zorder(10)
ax.add_collection(b)
plt.show()

# Alet Cevabini silme
cwd = os.getcwd()
folders = glob.glob(cwd + '/' + '2014-03-10_MW6.9_Off_Coast_Of_Northern_California.743736')
# Yeni SAC Dosyasi Kaydecek Fonksiyon
def sacwriter(time, stanw, stanm, stach, tr):
    folder = os.getcwd() + '/rrsac/'
    tr.write(folder + time + '.' + stanw + '.' +  stanm + '.' + stach + '.SAC', format="SAC")
    st2 = read(folder + time + '.' + stanw + '.' + stanm + '.' + stach + '.SAC')
    st2[0].write(folder + time + '.' + stanw + '.' + stanm + '.' + stach + '.SAC', format="SAC")
    return
for folder in folders:
  os.chdir(folder)
  # Alet Cevabi Dosyalarini Listeleme
  resps = glob.glob('RESP.*')
  # SAC Dosyalarini Listeleme
  sacs = glob.glob('*.SAC')
  for resp in resps:
    # SAC Dosyalari ile Alet Cevap Dosyalarini Eslestirme
    for sac in sacs:
      st = read(sac)
      STNW = st[0].stats.network
      STNM = st[0].stats.station
      STCH = st[0].stats.channel
      time = str(st[0].stats.starttime)
      STLA = st[0].stats.sac.stla
      STLO = st[0].stats.sac.stlo
      EVLA = st[0].stats.sac.evla
      EVLO = st[0].stats.sac.evlo
      # UNITS: ACC - Ivme, VEL - Hiz, DISP - Yer Degistirme
      seedresp = {'filename': resp, 'units': 'ACC'}
      if resp == 'RESP.' + STNW + '.' + STNM + '..' + STCH:
	print(resp,sac)
	pre_filt = [0.1, 0.5, 40, 50]
	st.simulate(paz_remove=None, pre_filt=pre_filt, seedresp=seedresp)
	# Yeni Dosyalari Kaydetme
	sacwriter(time, STNW, STNM, STCH, st[0])

# Toplu Deprem Indirme
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
client = Client()
def makedir(dirname):
  try:
    os.makedirs(dirname)
  except:
    pass
  
starttime = UTCDateTime("2018-08-13")
endtime = UTCDateTime("2018-09-13")
cat = client.get_events(starttime=starttime, endtime=endtime,minlatitude = 36, maxlatitude = 42, minlongitude = 26, maxlongitude = 45, minmagnitude=5,maxdepth=500)
cat.plot('local')

for i, event in enumerate(cat):
  # Yil, ay, gun, saat, dakika, saniye
  t = str(cat[i].origins[0].time)
  time = t[0:19]
  eventname = cat[i].event_descriptions[0].text + time
  # Bosluklari nokta ile dolduruyoruz
  eventname = str(eventname)
  a = eventname.replace(" ", "")
  b = a.replace('.','')
  event_dir = '/event'
  makedir(event_dir)
  cwd = os.getcwd()
  completeName = cwd + '/event/'
  cat[i].write(completeName + b + '.xml', "QUAKEML")
  
from obspy import read_events
from obspy.clients.fdsn.mass_downloader import (CircularDomain, MassDownloader,
                                                Restrictions, RectangularDomain)
def download(eqname, t0, lat, lon, min_length=360):
    domain = RectangularDomain(36, 42, 26, 45)

    restrictions = Restrictions(starttime=t0,endtime=t0+min_length,chunklength_in_sec=60*7,network="*", station="*", location="", channel="BH*",
        reject_channels_with_gaps=True,minimum_length=0.0,minimum_interstation_distance_in_m=100.0)

    waveform_dir = "{}/waveforms".format(eqname)
    stationxml_dir = "{}/stations".format(eqname)
    makedir(waveform_dir)
    makedir(stationxml_dir)

    mdl = MassDownloader(providers=["http://eida.koeri.boun.edu.tr:8080"])
    mdl.download(domain, restrictions,
                 mseed_storage=waveform_dir,
                 stationxml_storage=stationxml_dir)

path = (os.getcwd()+'/event/')
xmls = os.listdir(path)

for filename in xmls:
  os.chdir(path)
  eqname, _ = os.path.split(filename)[1].split(".")
  cat = read_events(filename)
  event = cat[0]
  origin = event.preferred_origin()
  download(eqname, origin.time, origin.latitude, origin.longitude)
