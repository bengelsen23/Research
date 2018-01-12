;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"

;***********************************************


; use this script to generate the Weisman analytic sounding and
; hodograph, and then dump the profiles into an input_sounding file
; that can be read by CM1
;
; adapted from code from CM1 (in base.F), as modified by G. Marion


 begin

  case = "defaultBR"


; some settings
  fact = 1       	; Scaling factor for low-level winds (for increasing 0-2 SRH) --- base value 1
  umax1_val = 7	        ; Hodograph radius (vmax) (default = 7 m/s)
  umax2_val = 31 	; Maximum u-wind (default = 31 m/s)
  udep1=2000.0   	; Level where v=umax for quarter-circle hodograph
  udep2=6000.0   	; Level of maximum winds - above here, constant at umax2
  thexp = 1.25     	; Exponent for potential temperature calculation (lapse rate) --- base value 1.25 for WK82 
  vertgrid = 500 	; Vertical grid spacing (originally 500) 
  theta_tr = 343.       ; Tropopause potential temperature - base value 343 for WK82

  T_tr     = 213        ; Tropopause temperature - base value 213 for WK82 (for now, get from sounding - need to work on this)

; some constants
  rpi = atan(1.)*4.
  radcnvrt = 180./rpi
  Rd = 287.04  ; J/kg/K
  g = 9.81   ; m/s
  Rv = 461.5
  cp = 1005.7
  cpd = 1005.7
  cpv = 1870.0
  cpl = 4190.0
  cpi = 2118.636
  cv  = cp-Rd
  p00 = 100000
  T00 = 273.15
  m2km = 0.001
  km2m = 1000.
  L = 2.501*10^6
  epsilon = Rd/Rv
  xlv = L
  xls = 2836017.0
 
  kappa = (cp-cv)/cp
  kappa_d = Rd/cp
  rp00 = 1/p00
  reps = Rv/Rd
  eps = epsilon
  rddcp = kappa_d
  cpdrd = cp/Rd
  cpdg = cp/g


; height coordinates; assuming 1/2 grid heights
  z      = ispan(250, 20000, vertgrid) ; Change vertical grid spacing here (20 km as model top recommended)
  zflt   = int2flt(z)
  th     = new(dimsizes(z),float)
  thv    = new(dimsizes(z),float)
  p      = new(dimsizes(z),float)
  pmb    = new(dimsizes(z),float)
  qv     = new(dimsizes(z),float)
  rh     = new(dimsizes(z),float)
  pi     = new(dimsizes(z),float)
  T      = new(dimsizes(z),float)
  Tc     = new(dimsizes(z),float)
  u      = new(dimsizes(z),float)
  v      = new(dimsizes(z),float)
  qvgkg  = new(dimsizes(z),float)
  es     = new(dimsizes(z),float)
  el     = new(dimsizes(z),float)
  Td     = new(dimsizes(z),float)
  Tdc    = new(dimsizes(z),float)
  V      = new(dimsizes(z),float)

  z_tr     = 12000.
  qv_0     = 0.014
  p_0      = 100000.
  dimsizez = dimsizes(z)-1
  th_0     = 300.
  tr_lev   = (z_tr-(vertgrid/2))/vertgrid

; parameters
  th_sfc = th_0
  qv_pbl = qv_0
  z_tr = z_tr
  t_tr = T_tr
  th_tr = theta_tr
  p_sfc = p_0
  pmb_sfc = p_sfc/100

  pi_sfc  = (p_sfc/p00)^(Rd/cp)
  qv_sfc  = qv_0 
  thv_sfc = th_sfc*(1+qv_sfc*reps)/(1+qv_sfc)

  qv(:) = 0
  rh(:) = 0
  pi(:) = 0
  ; calculate
  ; WK82 pp 506
  do k = 0,dimsizez
    if (z(k).lt.z_tr)
      th(k) = th_sfc+(th_tr-th_sfc)*((z(k)/z_tr)^(thexp))
      rh(k) = 1.0-0.75*((z(k)/z_tr)^1.25)
    else
      th(k) = th_tr*exp((g/(cp*t_tr))*(z(k)-z_tr))
      rh(k) = 0.25
    end if
  end do 

  do n = 0,19
    do k = 0,dimsizez
      thv(k) = th(k)*(1+reps*qv(k))/(1+qv(k))
    end do
    pi(0) = pi_sfc-g*z(0)/(cp*0.5*(thv_sfc+thv(0)))
    do k = 1,dimsizez
      pi(k) = pi(k-1)-g*(z(k)-z(k-1))/(cp*0.5*(thv(k)+thv(k-1)))
    end do
    do k = 0,dimsizez
      p(k)   = p00*(pi(k))^(cp/Rd)
      pmb(k) = p(k)/100
      T(k)   = th(k)*(p00/p(k))^(-1*kappa_d)
      Tc(k)  = T(k)-273.15
      es(k)  = 611.2*exp(17.67*Tc(k)/(Tc(k)+243.5))
    end do
    do k = 0,dimsizez
      qv(k)    = rh(k)*epsilon*es(k)/(p(k)-es(k))
      qvgkg(k) = qv(k)*1000
      el(k)    = log((qv(k)/epsilon)*p(k)/(100/(1+(qv(k)/epsilon))))
      Td(k)    = T00+(243.5*el(k)-440.8)/(19.48-el(k))
      Tdc(k)   = Td(k)-273.15
    end do
    do k = 0,dimsizez
      if(qv(k).gt.qv_pbl)
        qv(k)    = qv_pbl
        qvgkg(k) = qv(k)*1000
        el(k)    = log((qv(k)/epsilon)*p(k)/(100/(1+(qv(k)/epsilon))))
        Td(k)    = T00+(243.5*el(k)-440.8)/(19.48-el(k))
        Tdc(k)   = Td(k)-273.15
      end if
    end do
  end do



; wind profile - default is quarter circle
   udep1=2000.0
   udep2=6000.0
   umax1=umax1_val
   umax2=umax2_val
     
   do k = 0,dimsizez
     if (z(k).le.udep1) then
       ANGLE = 90.0*(z(k)/udep1)*(rpi/180.0)
       u(k)     = umax1-(umax1*fact)*cos(ANGLE)
       v(k)     = umax1*sin(ANGLE)
     else if (z(k).gt.udep1 .and. z(k).le.udep2) then
       u(k)     = umax1+(z(k)-udep1)*(umax2-umax1)/(udep2-udep1)
       v(k)     = umax1
     else
       u(k)     = umax2
       v(k)     = umax1
     end if
     end if
   end do

; wind profile - full-circle
; udep1=2000.0
; udep2=6000.0
; umax1=umax1_val
; umax2=umax2_val
   
; do k = 0,dimsizez
;   if (z(k).le.udep2) then
;     ANGLE    = 360.0*(z(k)/udep2)*((rpi)/180.0)
;     u(k)     = umax1-((umax1*fact)*cos(ANGLE))
;     v(k)     = umax1*sin(ANGLE)
;   ;else if (z(k).gt.udep1 .and. z(k).le.udep2) then
    ;  u(k)     = umax1+(z(k)-udep1)*(umax2-umax1)/(udep2-udep1)
    ;  v(k)     = umax1
;   else
;     u(k)     = u(0)
;     v(k)     = v(0)
;   end if
;   ;end if
; end do

;*************************
; create plot
;*************************
;  wks  = gsn_open_wks ("ncgm", "skewt"+case+"")
   wks  = gsn_open_wks ("pdf", "skewt"+case+"")


; --- Create background skew-T and plot sounding + winds----------------

  dataOpts           = True    ; options describing data and plotting
  dataOpts@WspdWdir  = False     ; wind speed and dir [else: u,v]
  dataOpts@hemisphere    = "NH"
  dataOpts@colTemperature = "red"
  dataOpts@colDewPt       = "darkgreen"  

  skewtOpts                 = True
  skewtOpts@DrawFahrenheit  = False
  skewtOpts@DrawStandardAtm = False

  skewt_bkgd = skewT_BackGround (wks, skewtOpts)
  skewt_data = skewT_PlotData   (wks, skewt_bkgd, pmb,Tc,Tdc,z,u,v,dataOpts)
  draw (skewt_bkgd)
  draw (skewt_data)
  frame(wks)

; now, dump back out the edited file in the input_sounding format
  slist = [/pmb_sfc, th_sfc, qv_pbl/]
  dlist = [/zflt, th, qvgkg, u, v/]
  fileout1 = "input_sounding_"+case+""

  write_table(fileout1,"w",slist,"%12.4f%12.4f%12.4f")
  write_table(fileout1,"a",dlist,"%12.4f%12.4f%12.4f%12.4f%12.4f")

  windlist = [/u,v/]
  fileout2 = "winds_"+case+""
  write_table(fileout2,"w",windlist,"%12.4f%12.4f")

 end
