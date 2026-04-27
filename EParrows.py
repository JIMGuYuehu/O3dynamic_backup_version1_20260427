#! env python

import xarray as xr
import numpy as np
import aostools_functions as ac
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.style.use(['seaborn-whitegrid','seaborn-poster','plotstyle.mplstyle'])
# plt.style.use('fivethirtyeight')

a0 = 6376.0e3
g  = 9.81


era_dir = 'data/'

def GetAxSize(fig,ax,dpi=False):
	"""get width and height of a given axis.
	   output is in inches if dpi=False, in dpi if dpi=True
	"""
	bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	width, height = bbox.width, bbox.height
	if dpi:
		width *= fig.dpi
		height *= fig.dpi
	return width, height

def set_size(w,h, ax=None):
	""" w, h: width, height in inches """
	if not ax: ax=plt.gca()
	l = ax.figure.subplotpars.left
	r = ax.figure.subplotpars.right
	t = ax.figure.subplotpars.top
	b = ax.figure.subplotpars.bottom
	figw = float(w)/(r-l)
	figh = float(h)/(t-b)
	ax.figure.set_size_inches(figw, figh)

def CorrectAspect(x,y):
	# x-coordinate in kilometers = a*phi
	delta_x = 6.4e3*np.deg2rad(x[-1]-x[0])
	# delta_y = y[-1] - y[0]
	delta_y = max(y) - min(y)
	# figure aspect ratio
	width,height = GetAxSize(plt.gcf(),plt.gca())
	return delta_y/height,delta_x/width

def Palmer(F_phi,F_p,pres='level',lat='lat',xlim=None):
	# Palmer JAS (1981)
	# F_phi = a*cosphi*(-u'v')
	# F_p   = a*cosphi*(f*v'theta'/theta_p)
	# In Palmer (1981), F_z is used, so that we have to
	#  convert F_p to F_z here
	H = 6.4e3
	p0 = 1000
	p = F_p[pres]
	z = -H*np.log(p/p0)
	F_z = -F_p * H/p
	a = a0
	Omega = 2*np.pi/86400.
	cosphi = np.cos(np.deg2rad(F_phi[lat]))
	# Fhat_phi = 2*pi*a*cosphi*exp(-z/H)*F_phi
	# Fhat_p   = 2*pi*a*cosphi*exp(-z/H)*F_z
	factr = 2*np.pi*a*cosphi*np.exp(-z/H)
	# Palmer also adjusts for aspect ratio with a factor c
	if xlim is None:
		xlim = [F_phi[lat][0],F_phi[lat][-1]]
	dx,dy = CorrectAspect(xlim,z*1e-3) # z in km
	c = dx/dy
	print('Palmer c = {0}'.format(c.values))
	return (factr*c*F_phi).transpose(),(factr*F_z).transpose()

def Baldwin(F_phi,F_p,pres='level',lat='latitude',xlim=[-90,90]):
	# Baldwin et al JAS (1985), Figure 6E
	# F_phi = a*cosphi*(-u'v')
	# F_p   = a*cosphi*(f*v'theta'/theta_p)
	# log(p) coordinates. Note here that
	#  they say they use ln(p) as vertical coordinate,
	#  but in fact they use log10(p)
	#  They do not say what value they use for H either,
	#   so we have to guess here that it's the same as Palmer (1981)
	H = 6.4e3
	p0 = 1000.
	z = -H*np.log(F_p[pres]/p0)
	# They use theta_z instead of theta_p, so we need to convert
	F_z = -F_p*H/p0
	# Then, they multiply by exp(z/H) = p0/0
	# Fhat_phi = F_phi*exp(z/H) = F_phi*p0/p
	# Fhat_p   = F_z*exp(z/H) = F_z*p0/p
	# in addition,
	factr = p0/F_p[pres]
	# also correct for aspect ratio
	dx,dy = CorrectAspect(xlim,z*1.e-3)
	return factr*dx*F_phi,factr*dy*F_z

def Taguchi(F_phi,F_p,pres='level',lat='latitude'):
	# F is assumed in geostrophic scaling form from
	#  Andrews et al (1987). We assume this uses F_z instead of F_p
	# F_phi = rho*a*cosphi*(-u'v' + ...)
	# F_z   = rho*a*cosphi*(-u'w(z)' + ...)
	#  as rho is a constant and multiplies both components,
	#  we can safely set it to unity.
	# From Taguchi and Hartman (2006) we do not get any information
	#  on whether they take care of aspect ratio etc. Probably not.
	#  They only mention a division by sqrt(1000/pressure).
	# They *almost certainly* do not divide by sqrt() but multiply by sqrt().
	fact = np.sqrt(1000/F_p[pres])
	# Benefit of the doubt: assume they do this correctly
	# dx,dy = CorrectAspect(F_phi[lat],z*1.e-3)
	dx,dy = CorrectAspect(F_phi[lat],np.log10(F_p[pres]))
	return fact*dx*F_phi,fact*dy*F_p

def ESRL(F_phi,F_p,pres='level',lat='latitude',multfact=1):
	# From esrl.noaa.gov website
	# F_phi = a*cosphi*(-u'v')
	# F_p   = a*cosphi*(f*v'theta'/theta_p)
	# aspect ratio scale factors are simply
	cosphi = np.cos(np.deg2rad(F_phi[lat]))
	Fhat_phi = cosphi*F_phi/a0
	Fhat_p   = cosphi*F_p
	dx = 1/np.pi
	dy = 1./1000.
	if multfact == 1:
		factr = np.sqrt(1000./F_p[pres])
	else:
		factr = xr.where(F_p[pres] > 100,1,multfact)
	return dx*factr*Fhat_phi,-dy*factr*Fhat_p


def Reading(F_phi,F_p,pres='level',lat='latitude'):
	# from met.reading.ac.uk
	return ESRL(F_phi,F_p,pres,lat)


def Andrews(F_phi,F_p,pres='level',lat='latitude',xlim=[-90,90]):
	# Andrews et al JAS (1983)
	#  F_phi = a*cosphi*(-u'v' + ...)
	#  F_p   = a*cosphi*(-u'w' + ...)
	p0 =1013.25
	# note that they don't use H, or set H=1. That doesn't work.
	#  so we use the same H as all others here.
	H = 6.4e3
	z = -H*np.log(F_p[pres]/p0)
	# "To represent the vector F in a rectangular
	#  plot with phi and z as coordinates, Fphi and Fz
	#  are separately stretched in an obvious way."
	# They do not say what this obvious way is.
	# We give them the benefit of the doubt and do this correctly
	dx,dy = CorrectAspect([-90,90],z*1e-3)
	# What they do say is the following scaling:
	Fhat_phi = F_phi*F_p[pres]/p0
	Fhat_p   = -H*F_p/p0
	# note that this is equivalent to only rescaling F_p by -H/p
	# Fhat_phi = F_phi
	# Fhat_p   = -H*F_p/F_p[pres]
	return (dx*Fhat_phi).transpose(),(dy*Fhat_p).transpose(),z

def Fz(F_phi,F_p,pres='level',lat='latitude'):
	# using F_z instead of F_p makes things slightly more complicated
	# input F_phi,F_p are the same as for all other funcitons, i.e. in pressure
	#  coordinates
	H = 6.4e3
	p0 = 1000
	z = -H*np.log(F_p[pres]/p0)
	rho = F_p[pres]*100/g/H # rho_0 = p[Pa]/(gH)
	cosphi = np.cos(np.deg2rad(F_phi[lat]))
	# The correct scaling here is this:
	# Fhat_phi = 2*np.pi*a0*cosphi*rho*F_phi
	# That ultimately, that will yield momentum change.
	#  But usually, we are interested in acceleration, i.e.
	#  speed per mass. So we don't multiply by density
	Fhat_phi = 2*np.pi*a0*cosphi*F_phi
	F_z = -F_p*H/F_p[pres]
	# Same scaling argument as for Fhat_phi: no density
	Fhat_z   = 2*np.pi*a0**2*cosphi*F_z
	width,height = GetAxSize(plt.gcf(),plt.gca())
	dx = width/(np.pi/180*(F_phi[lat][0]-F_phi[lat][-1]))
	dy = height/(max(z)-min(z))
	return (dx*Fhat_phi).transpose(),(dy*Fhat_z).transpose()


def ReadERA(dates,lskip=2,pskip=2,divs=False):
	if isinstance(dates,str):
		import datetime as dtt
		formt = '%Y-%m-%d'
		prev_day = dtt.datetime.strftime(dtt.datetime.strptime(dates,formt) - dtt.timedelta(days=1),formt)
		next_day = dtt.datetime.strftime(dtt.datetime.strptime(dates,formt) + dtt.timedelta(days=1),formt)
		sel_dates = slice(prev_day,next_day)
	else:
		sel_dates = dates
	era_u_v_t = xr.open_mfdataset([era_dir+'ERA5_daily.1979.u.nc',era_dir+'ERA5_daily.1979.v.nc',era_dir+'ERA5_daily.1979.t.nc']).sel(time=sel_dates,latitude=slice(None,None,lskip),level=slice(None,None,pskip))
	ep1,ep2,div1,div2 = ac.ComputeEPfluxDiv(era_u_v_t.latitude.values,era_u_v_t.level.values,era_u_v_t.u.values,era_u_v_t.v.values,era_u_v_t.t.values)
	f1 = xr.DataArray(ep1,coords=[era_u_v_t.time,era_u_v_t.level,era_u_v_t.latitude],name='f1')
	f2 = xr.DataArray(ep2,coords=[era_u_v_t.time,era_u_v_t.level,era_u_v_t.latitude],name='f2')
	cosphi = np.cos(np.deg2rad(era_u_v_t.latitude))
	F_phi = a0*cosphi*f1
	F_p = a0*cosphi*f2
	if isinstance(dates,str):
		out = [f1.sel(time=dates),f2.sel(time=dates),F_phi.sel(time=dates),F_p.sel(time=dates)]
	else:
		out = [f1.mean(dim='time'),f2.mean(dim='time'),F_phi.mean(dim='time'),F_p.mean(dim='time')]
	if divs:
		div1 = xr.DataArray(div1,coords=f1.coords,name='div1')
		div2 = xr.DataArray(div2,coords=f1.coords,name='div2')
		if isinstance(dates,str):
			div1 = div1.sel(time=dates)
			div2 = div2.sel(time=dates)
		else:
			div1 = div1.mean(dim='time')
			div2 = div2.mean(dim='time')
		out.append(div1)
		out.append(div2)
	return out

if __name__ == '__main__':
	quiv_args = {'angles':'uv','scale_units':'inches','linewidths':1,'edgecolors':'k','color':'none'}
	for scale in ['log']:#['linear','log']:
		Fx = {}
		Fy = {}
		# Baldwin et al 1985, Figure 6:
		method = 'Baldwin'
		f1,f2,F_phi,F_p,div1,div2 = ReadERA(slice('1979-01-05','1979-03-03'),pskip=1,divs=True)
		lat = F_phi.sel(latitude=slice(80,0)).latitude
		pres = np.array([50,70,100,150,200,250,300,400,500,700,850])
		fig,ax = plt.subplots()
		set_size(5.2,5.2,ax)
		# Reference method
		Fx[method],Fy[method] = Baldwin(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres),xlim=[0,80])
		ax.quiver(lat,pres,Fx[method],Fy[method],scale=1.5e10,**quiv_args)
		# Jucker (2020) method
		_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
			,fig,ax,xlim=[0,80],ylim=[1000,50],yscale=scale)
		ax.get_children()[1].set_facecolor('r')
		# if scale == 'linear':
		#     ax.invert_yaxis()
		plt.title('Baldwin et al 1985')
		ax.set_xlabel('latitude')
		ax.set_ylabel('pressure [hPa]')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/Baldwin85.{0}.pdf'.format(scale),bbox_inches='tight')

		# # Correct method, but with F_z
		# method = 'Fz'
		# # pres = F_phi.level
		# pres = F_phi.sel(level=np.array([850,500,300,200,150,100,70,50,30,20,10,7,5,2,1])[::-1]).level
		# fig,ax = plt.subplots()
		# Fx[method],Fy[method] = Fz(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres),lat='latitude')
		# H = 6.4e3
		# p0 = 1000
		# z = -H*np.log(pres/p0)
		# # make sure the z- and p-axes have the same extent
		# ax.set_ylim(z[-1],z.interp({'level':1},kwargs={'fill_value':"extrapolate"}))
		# if scale == 'log':
		# 	axp = ax.twinx()
		# 	# (div1+div2).plot.contourf(ax=axp,levels=np.linspace(),robust=True)
		# 	# _,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
		# 	# 	,fig,axp,xlim=[0,80],ylim=[1000,1],yscale=scale,newax=False)
		# else:
		# 	axp = ax
		# (div1+div2).plot.contourf(ax=axp,levels=np.linspace(-15,15,16),robust=True)
		# _,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
		# 	,fig,axp,xlim=[0,80],ylim=[850,1],yscale=scale,newax=False)
		# if scale == 'log':
		# 	ax.quiver(lat,z,Fx[method],Fy[method],angles='uv')
		# 	ax.set_ylabel('height [m]')
		# else:
		# 	ax.quiver(lat,pres,Fx[method],Fy[method],angles='uv')
		# ax.set_xlabel('latitude')
		# plt.title('Fz')
		# fig.savefig('figures/Fz.{0}.pdf'.format(scale),bbox_inches='tight')

		# Andrews et al (1983)
		# they analyse model data, which I don't have.
		#  I therefore use the same as the Baldwin et al (1985) analysis.
		method = 'Andrews'
		f1,f2,F_phi,F_p = ReadERA(slice('1979-01-05','1979-03-03'),pskip=1)
		pres = np.array([850,500,300,200,150,100,70,50,30,20,10,7,5,2])[::-1]
		# pres = F_p.sel(level=slice(1,850)).level
		lat = F_phi.latitude
		# pres = np.array([50,70,100,150,200,250,300,400,500,700,850])
		fig,ax = plt.subplots()
		set_size(8.8,5.3,ax)
		# Reference method
		Fx[method],Fy[method],z = Andrews(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres))
		if scale == 'log':
			ax.quiver(lat,z,Fx[method],Fy[method],scale=3e9,**quiv_args)
			# Jucker (2020) method
			_,_,axp = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
				,fig,ax,xlim=[-90,90],ylim=[1000,1],yscale=scale,newax=True,scale=1e16)
			axp.get_children()[0].set_facecolor('r')
			axp.yaxis.set_major_formatter(ticker.ScalarFormatter())
			ax.set_ylim(0,6.4e3*np.log(1000))
			ax.set_ylabel('height [m]')
		else:
			ax.quiver(lat,pres,Fx[method],Fy[method],scale=3e9,**quiv_args)
			# Jucker (2020) method
			_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres),fig,ax,xlim=[-90,90],ylim=[1000,1],yscale=scale,newax=False,scale=1e16)
			ax.get_children()[1].set_facecolor('r')
		plt.title('Andrews et al 1983')
		ax.set_xlabel('latitude')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/Andrews83.{0}.pdf'.format(scale),bbox_inches='tight')

		# Taguchi and Hartmann (2006)
		# they analyse model data which I don't have access to, so we will use the
		# SSW of 1979-02-22 which happened during the Baldwin et al (1985) analysis
		#  There is too little information in the paper to do this properly.
		# method = 'Taguchi'
		# f1,f2,F_phi,F_p = ReadERA(slice('1979-02-10','1979-02-22'))
		# lat = F_phi.sel(latitude=slice(85,0)).latitude
		# pres = F_p.sel(level=slice(1,975)).level
		# fig,ax = plt.subplots()
		# set_size(5.2,5.2,ax)
		# Fx[method],Fy[method] = Taguchi(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres))
		# _,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
		#     ,fig,ax,xlim=[0,85],ylim=[1000,1],yscale='log')
		# ax.quiver(lat,pres,Fx[method],Fy[method],angles='uv')
		# plt.title('Taguchi and Hartmann 2006')
		# ax.set_xlabel('latitude')
		# ax.set_ylabel('pressure [hPa]')
		# fig.savefig('figures/Taguchi06.pdf',bbox_inches='tight')

		# We do the same thing with ESRL, as they seem to use about the same method (maybe?)
		method = 'ESRL'
		f1,f2,F_phi,F_p = ReadERA(slice('1979-01-05','1979-03-03'),pskip=1)
		lat = F_phi.sel(latitude=slice(85,-85)).latitude
		# pres = F_p.sel(level=slice(10,975)).level
		# NCEP standard levels
		pres = np.array([925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10])[::-1]
		# first, don't do special scaling
		fig,ax = plt.subplots()
		set_size(8.7,5.1)
		# Reference method
		Fx[method],Fy[method] = ESRL(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres),multfact=0.99)
		ax.quiver(lat,pres,Fx[method],Fy[method],scale=2e2,**quiv_args)
		# Jucker (2020) method
		_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
			,fig,ax,xlim=[-85,85],ylim=[1000,10],yscale=scale,scale=2e16)
		ax.get_children()[1].set_facecolor('r')
		plt.title('PSL')
		ax.set_xlabel('latitude')
		ax.set_ylabel('pressure [hPa]')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/PSL_raw.{0}.pdf'.format(scale),bbox_inches='tight')
		# then, use sqrt(1000/p) scaling
		fig,ax = plt.subplots()
		set_size(8.7,5.1)
		# Reference method
		Fx[method],Fy[method] = ESRL(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres))
		ax.quiver(lat,pres,Fx[method],Fy[method],scale=2e2,**quiv_args)
		# Jucker (2020) method
		_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
			,fig,ax,xlim=[-85,85],ylim=[1000,10],yscale=scale,scale=2e16)
		ax.get_children()[1].set_facecolor('r')
		plt.title(u'PSL ($\sqrt{1000/p}$)')
		ax.set_xlabel('latitude')
		ax.set_ylabel('pressure [hPa]')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/PSL_sqrt.{0}.pdf'.format(scale),bbox_inches='tight')
		# then, use arbitrary scaling above  100hPa
		fig,ax = plt.subplots()
		set_size(8.7,5.1)
		# Reference method
		Fx[method],Fy[method] = ESRL(F_phi.sel(latitude=lat,level=pres),F_p.sel(latitude=lat,level=pres),multfact=5)
		ax.quiver(lat,pres,Fx[method],Fy[method],scale=2e2,**quiv_args)
		# Jucker (2020) method
		_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres).level,f1.sel(latitude=lat,level=pres),f2.sel(latitude=lat,level=pres)
			,fig,ax,xlim=[-85,85],ylim=[1000,10],yscale=scale,scale=2e16)
		ax.get_children()[1].set_facecolor('r')
		plt.title('PSL (x5 above 100hPa)')
		ax.set_xlabel('latitude')
		ax.set_ylabel('pressure [hPa]')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/PSL_mult.{0}.pdf'.format(scale),bbox_inches='tight')


		# Palmer 1981
		method = 'Palmer'
		f1,f2,F_phi,F_p = ReadERA('1979-02-21',pskip=1)
		lat = F_phi.sel(latitude=slice(85,20)).latitude
		# here we want to use the same aspect ratio as in the paper
		# Figure 7 (top)
		# Palmer has a special equation to define the pressure levels
		m = np.arange(10)+1
		pres = 10**((10-m[1:])/3.)[::-1]
		# pres = F_p.sel(level=slice(1,975)).level
		fig,ax = plt.subplots()
		set_size(6.2,4.2,ax)
		# Reference method
		Fx[method],Fy[method] = Palmer(F_phi.sel(latitude=lat,level=pres,method='nearest'),F_p.sel(latitude=lat,level=pres,method='nearest'),lat='latitude',xlim=[85,20])
		ax.quiver(lat,pres,Fx[method],Fy[method],scale=3e13,**quiv_args)
		# we are plotting from +85 to +20 instead of 20-85 (for some reason)
		_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres,method='nearest').level,f1.sel(latitude=lat,level=pres,method='nearest'),f2.sel(latitude=lat,level=pres,method='nearest')
			,fig,ax,xlim=[85,20],ylim=[1000,1],yscale=scale)
		ax.get_children()[1].set_facecolor('r')
		plt.title('Palmer 1981')
		ax.set_xlabel('latitude')
		ax.set_ylabel('pressure [hPa]')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/Palmer81a.{0}.pdf'.format(scale),bbox_inches='tight')
		# Figure 7 (bottom)
		#pres = F_p.sel(level=slice(1,50)).level
		pres = pres[:6]
		fig,ax = plt.subplots()
		set_size(6.2,4.2,ax)
		# Reference method
		Fx[method],Fy[method] = Palmer(F_phi.sel(latitude=lat,level=pres,method='nearest'),F_p.sel(latitude=lat,level=pres,method='nearest'),lat='latitude',xlim=[85,20])
		ax.quiver(lat,pres,Fx[method],Fy[method],scale=1e13,**quiv_args)
		# we are plotting from +85 to +20 instead of 20-85 (for some reason)
		_,_ = ac.PlotEPfluxArrows(f1.sel(latitude=lat).latitude,f1.sel(level=pres,method='nearest').level,f1.sel(latitude=lat,level=pres,method='nearest'),f2.sel(latitude=lat,level=pres,method='nearest')
			,fig,ax,xlim=[85,20],ylim=[50,1],yscale=scale)
		ax.get_children()[1].set_facecolor('r')
		plt.title('Palmer 1981')
		ax.set_xlabel('latitude')
		ax.set_ylabel('pressure [hPa]')
		ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
		fig.savefig('figures/Palmer81b.{0}.pdf'.format(scale),bbox_inches='tight')
