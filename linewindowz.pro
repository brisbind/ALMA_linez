;-------------------------------------------------------------
;+
; NAME:
;       HOR
; PURPOSE:
;       Plot a horizontal line on a graph at specified y value.
; CATEGORY:
; CALLING SEQUENCE:
;       hor, y
; INPUTS:
;       y = Y value of horizontal line.  Scalar or array.    in
; KEYWORD PARAMETERS:
;       Keywords:
;         /DEVICE means work in device coordinates.
;         /NORMALIZED means work in normalized coordinates.
;           Default is data coordinates.
;         LINESTYLE=s.  Linestyle (def=!p.linestyle).
;         COLOR=c.      Line color (def=!p.color).
;         THICKNESS=t   Line thickness (def=!p.thick).
;         FILL=clr        Optional color to fill between line pairs.
;           Fills between lines 0 and 1, 2 and 3, and so on.
;         POINTER=pt      Draw arrowhead pointers at left and right
;           instead of lines.  Arrowhead dimensions may be given as
;           fraction of screen or plot window size, the value of
;           pt is height, or [height, width].  For /pointer the
;           default used is [.03,.03].
;         /LEFT  used with POINTER to plot left pointers only.
;         /RIGHT  used with POINTER to plot right pointers only.
; OUTPUTS:
; COMMON BLOCKS:
; NOTES:
;       Notes: see ver.
; MODIFICATION HISTORY:
;       R. Sterner, 2 Aug, 1989.
;       R. Sterner, 21 May, 1992 --- fixed for log X axes.
;       R. Sterner,  3 Nov, 1992 --- Added /device.
;       R. Sterner, 20 Jun, 1993 --- Added /normalized.
;       R. Sterner,  1994 Feb  2 --- Added THICK.
;       R. Sterner, 1994 Jun 3 --- Added FILL.
;       R. Sterner, 1994 Jun 16 --- Added POINTER, /TOP, /BOTTOM.
;
; Copyright (C) 1989, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	pro hor, y, help=hlp, device=device, linestyle=ls, color=clr, $
	  thickness=thk, normalized=norm, fill=fill, pointer=pt, $
	  left=left, right=right
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Plot a horizontal line on a graph at specified y value.'
	  print,' hor, y'
	  print,'   y = Y value of horizontal line.  Scalar or array.    in'
	  print,' Keywords:'
          print,'   /DEVICE means work in device coordinates.'
          print,'   /NORMALIZED means work in normalized coordinates.'
	  print,'     Default is data coordinates.'
	  print,'   LINESTYLE=s.  Linestyle (def=!p.linestyle).'
	  print,'   COLOR=c.      Line color (def=!p.color).'
	  print,'   THICKNESS=t   Line thickness (def=!p.thick).'
          print,'   FILL=clr        Optional color to fill between line pairs.'
          print,'     Fills between lines 0 and 1, 2 and 3, and so on.'
          print,'   POINTER=pt      Draw arrowhead pointers at left and right'
          print,'     instead of lines.  Arrowhead dimensions may be given as'
          print,'     fraction of screen or plot window size, the value of'
          print,'     pt is height, or [height, width].  For /pointer the'
          print,'     default used is [.03,.03].'
	  print,'   /LEFT  used with POINTER to plot left pointers only.'
	  print,'   /RIGHT  used with POINTER to plot right pointers only.'
	  print,' Notes: see ver.'
	  return
	end
 
	yy = y
	n = n_elements(yy)
	if n_elements(ls) eq 0 then ls = !p.linestyle
	if n_elements(clr) eq 0 then clr = !p.color
	if n_elements(thk) eq 0 then thk = !p.thick
        ;------  Handle pointers  -----------
        if n_elements(pt) eq 0 then pt = 0
        pflag = 0
        if pt(0) gt 0 then begin
          if pt(0) eq 1 then pt=.03
          ht = pt(0)
          wd = pt(n_elements(pt)-1)
          if n_elements(pt) eq 1 then wd = ht/2.
          pflag = 1
        endif
	lflag=0
	rflag=0
	if keyword_set(left) then lflag=1
	if keyword_set(right) then rflag=1
	if (lflag+rflag) eq 0 then begin
	  lflag=1
	  rflag=1
	endif
 
        ;--------  Device  ------------
        if keyword_set(device) then begin
          xx = [0,!d.x_size-1]
          for i = 0, n-1 do begin
            ;--------  Filled line pairs  -----------
            if n_elements(fill) ne 0 then begin
              if (i mod 2) eq 0 then begin
                x1 = xx(0) & x2 = xx(1)
                y1 = yy(i) & y2 = yy((i+1)<(n-1))
                polyfill, /dev, [x1,x2,x2,x1],[y1,y1,y2,y2],color=fill
              endif
            ;---------  Single lines  ---------------
            endif else if pflag eq 0 then begin
              plots, /device,xx,[0,0]+yy(i),linestyle=ls,color=clr,thick=thk
            ;---------  Pointers  -------------------
            endif else begin
              dx = round((!d.x_size-1)*ht)
              x1 = [0,dx,0]
              x2 = !d.x_size-1 - [0,dx,0]
              dy = round((!d.y_size-1)*wd/2.)
              dy = [-dy,0,dy]+yy(i)
              if lflag then polyfill,/dev,x1,dy,col=clr
              if rflag then polyfill,/dev,x2,dy,col=clr
            endelse
          endfor
	;---------  Normalized  ----------
        endif else if keyword_set(norm) then begin
          xx = [0,1]
          for i = 0, n-1 do begin
            ;--------  Filled line pairs  -----------
            if n_elements(fill) ne 0 then begin
              if (i mod 2) eq 0 then begin
                x1 = xx(0) & x2 = xx(1)
                y1 = yy(i) & y2 = yy((i+1)<(n-1))
                polyfill, /norm, [x1,x2,x2,x1],[y1,y1,y2,y2],color=fill
              endif
            ;---------  Single lines  ---------------
            endif else if pflag eq 0 then begin
              plots, /norm,xx,[0,0]+yy(i),linestyle=ls,color=clr,thick=thk
            ;---------  Pointers  -------------------
            endif else begin
              dx = ht
              x1 = [0,dx,0]
              x2 = [1,1-dx,1]
              dy = wd/2.
              dy = [-dy,0,dy]+yy(i)
              if lflag then polyfill,/norm,x1,dy,col=clr
              if rflag then polyfill,/norm,x2,dy,col=clr
            endelse
          endfor
	;----------  Data  -------------
        endif else begin
	  xx = [!x.range, !x.crange]
	  for i = 0, n-1 do begin
	    ;------  Linear X axis  ------------
	    if !x.type eq 0 then begin
              ;--------  Filled line pairs  -----------
              if n_elements(fill) ne 0 then begin
                if (i mod 2) eq 0 then begin
                  x1 = min(xx)  &  x2 = max(xx)
                  y1 = yy(i) & y2 = yy((i+1)<(n-1))
                  polyfill, [x1,x2,x2,x1],[y1,y1,y2,y2],color=fill,noclip=0
                endif
              ;---------  Single lines  ---------------
              endif else if pflag eq 0 then begin
	        oplot,[min(xx),max(xx)],[1.,1.]*yy(i),linestyle=ls,$
		  color=clr, thick=thk
              ;---------  Pointers  -------------------
              endif else begin
                dx = (!x.crange(1)-!x.crange(0))*ht
                x1 = [0,dx,0]+!x.crange(0)
                x2 = [0,-dx,0]+!x.crange(1)
                dy = (!y.crange(1)-!y.crange(0))*wd/2.
                if !y.type eq 0 then dy=[-dy,0,dy]+yy(i) else $
                  dy=10^([-dy,0,dy]+alog10(yy(i)))
                if lflag then polyfill,x1,dy,col=clr
                if rflag then polyfill,x2,dy,col=clr
              endelse
	    ;------  Log X axis  ------------
	    endif else begin
              ;--------  Filled line pairs  -----------
              if n_elements(fill) ne 0 then begin
                if (i mod 2) eq 0 then begin
                  y1 = yy(i) & y2 = yy((i+1)<(n-1))
                  x1 = min(xx)  &  x2 = max(xx)
                  polyfill, 10^[x1,x2,x2,x1],[y1,y1,y2,y2],color=fill,noclip=0
                endif
              ;---------  Single lines  ---------------
              endif else if pflag eq 0 then begin
	        oplot,10^[min(xx),max(xx)],[1.,1.]*yy(i),linestyle=ls,$
		  color=clr, thick=thk
              ;---------  Pointers  -------------------
              endif else begin
                dy = (!y.crange(1)-!y.crange(0))*wd/2.
                if !y.type eq 0 then dy=[-dy,0,dy]+yy(i) else $
                  dy=10^([-dy,0,dy]+alog10(yy(i)))
                dx = (!x.crange(1)-!x.crange(0))*ht
                x1 = 10^([0,dx,0]+!x.crange(0))
                x2 = 10^([0,-dx,0]+!x.crange(1))
                if lflag then polyfill,x1,dy,col=clr
                if rflag then polyfill,x2,dy,col=clr
              endelse
	    endelse  ; !x.type.
	  endfor
	endelse
 
	return
	end











;-------------------------------------------------------------
;+
; NAME:
;       VER
; PURPOSE:
;       Plot a vertical line on a graph at specified x value.
; CATEGORY:
; CALLING SEQUENCE:
;       ver, x
; INPUTS:
;       x = X value of vertical line. Scalar or array.    in
; KEYWORD PARAMETERS:
;       Keywords:
;         /DEVICE means work in device coordinates.
;         /NORMALIZED means work in normalized coordinates.
;           Default is data coordinates.
;         LINESTYLE=s.    Linestyle (def=!p.linestyle).
;         COLOR=c.        Line Color (def=!p.color).
;         THICKNESS=thk   Line thickness (def=!p.thick).
;         FILL=clr        Optional color to fill between line pairs.
;           Fills between lines 0 and 1, 2 and 3, and so on.
;         POINTER=pt      Draw arrowhead pointers at top and bottom
;           instead of lines.  Arrowhead dimensions may be given as
;           fraction of screen or plot window size, the value of
;           pt is height, or [height, width].  For /pointer the
;           default used is [.03,.03].
;         /BOTTOM  used with POINTER to plot bottom pointers only.
;         /TOP  used with POINTER to plot top pointers only.
; OUTPUTS:
; COMMON BLOCKS:
; NOTES:
;       Note: see hor.
; MODIFICATION HISTORY:
;       R. Sterner, 2 Aug, 1989.
;       R. Sterner, 21 May, 1992 --- fixed for log Y axes.
;       R. Sterner,  3 Nov, 1992 --- Added /device.
;       R. Sterner, 27 Jan, 1993 --- dropped reference to array.
;       R. Sterner 20 Jun, 1993 --- added /norm.
;       R. Sterner 1994 Feb 2 --- Add THICK.
;       R. Sterner, 1994 Jun 3 --- Added FILL.
;       R. Sterner, 1994 Jun 16 --- Added POINTER.
;
; Copyright (C) 1989, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	pro ver, x, help=hlp, device=device, linestyle=ls, color=clr, $
	  normalized=norm, thickness=thk, fill=fill, pointer=pt, $
	  top=top, bottom=bot
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Plot a vertical line on a graph at specified x value.'
	  print,' ver, x'
	  print,'   x = X value of vertical line. Scalar or array.    in'
	  print,' Keywords:'
	  print,'   /DEVICE means work in device coordinates.'
	  print,'   /NORMALIZED means work in normalized coordinates.'
	  print,'     Default is data coordinates.'
	  print,'   LINESTYLE=s.    Linestyle (def=!p.linestyle).'
	  print,'   COLOR=c.        Line Color (def=!p.color).'
	  print,'   THICKNESS=thk   Line thickness (def=!p.thick).'
	  print,'   FILL=clr        Optional color to fill between line pairs.'
	  print,'     Fills between lines 0 and 1, 2 and 3, and so on.'
	  print,'   POINTER=pt      Draw arrowhead pointers at top and bottom'
	  print,'     instead of lines.  Arrowhead dimensions may be given as'
	  print,'     fraction of screen or plot window size, the value of'
	  print,'     pt is height, or [height, width].  For /pointer the'
	  print,'     default used is [.03,.03].'
          print,'   /BOTTOM  used with POINTER to plot bottom pointers only.'
          print,'   /TOP  used with POINTER to plot top pointers only.'
	  print,' Note: see hor.'
	  return
	end
 
	xx = x
	n = n_elements(xx)
	if n_elements(ls) eq 0 then ls = !p.linestyle
	if n_elements(clr) eq 0 then clr = !p.color
	if n_elements(thk) eq 0 then thk = !p.thick
	;------  Handle pointers  -----------
	if n_elements(pt) eq 0 then pt = 0
	pflag = 0
	if pt(0) gt 0 then begin
	  if pt(0) eq 1 then pt=.03
	  ht = pt(0)
	  wd = pt(n_elements(pt)-1)
	  if n_elements(pt) eq 1 then wd = ht/2.
	  pflag = 1
	endif
        bflag=0
        tflag=0
        if keyword_set(bot) then bflag=1
        if keyword_set(top) then tflag=1
        if (bflag+tflag) eq 0 then begin
          bflag=1
          tflag=1
        endif
 
	;--------  Device  ------------
	if keyword_set(device) then begin
	  yy = [0,!d.y_size-1]
	  for i = 0, n-1 do begin
	    ;--------  Filled line pairs  -----------
	    if n_elements(fill) ne 0 then begin
	      if (i mod 2) eq 0 then begin
	        x1 = xx(i) & x2 = xx((i+1)<(n-1))
		y1 = yy(0) & y2 = yy(1)
		polyfill, /dev, [x1,x2,x2,x1],[y1,y1,y2,y2],color=fill
	      endif
	    ;---------  Single lines  ---------------
	    endif else if pflag eq 0 then begin
 	      plots,/device,[0,0]+xx(i),yy,linestyle=ls,color=clr,thick=thk
	    ;---------  Pointers  -------------------
	    endif else begin
	      dx = round((!d.x_size-1)*wd/2.)
	      dx = [-dx,0,dx]+xx(i)
	      dy = round((!d.y_size-1)*ht)
	      y1 = [0,dy,0]
	      y2 = !d.y_size-1 - [0,dy,0]
	      if bflag then polyfill,/dev,dx,y1,col=clr
	      if tflag then polyfill,/dev,dx,y2,col=clr
	    endelse
	  endfor
	;---------  Normalized  ----------
	endif else if keyword_set(norm) then begin
	  yy = [0,1]
	  for i = 0, n-1 do begin
	    ;--------  Filled line pairs  -----------
	    if n_elements(fill) ne 0 then begin
	      if (i mod 2) eq 0 then begin
	        x1 = xx(i) & x2 = xx((i+1)<(n-1))
		y1 = yy(0) & y2 = yy(1)
		polyfill, /norm, [x1,x2,x2,x1],[y1,y1,y2,y2],color=fill
	      endif
	    ;---------  Single lines  ---------------
	    endif else if pflag eq 0 then begin
 	      plots,/norm,[0,0]+xx(i),yy,linestyle=ls,color=clr,thick=thk
	    ;---------  Pointers  -------------------
	    endif else begin
	      dx = wd/2.
	      dx = [-dx,0,dx]+xx(i)
	      dy = ht
	      y1 = [0,dy,0]
	      y2 = [1,1-dy,1]
	      if bflag then polyfill,/norm,dx,y1,col=clr
	      if tflag then polyfill,/norm,dx,y2,col=clr
	    endelse
	  endfor
	;----------  Data  -------------
	endif else begin
	  yy = [!y.range, !y.crange]
	  for i = 0, n-1 do begin
	    ;--------  Linear Y axis  ----------
	    if !y.type eq 0 then begin
	      ;--------  Filled line pairs  -----------
	      if n_elements(fill) ne 0 then begin
	        if (i mod 2) eq 0 then begin
	          x1 = xx(i) & x2 = xx((i+1)<(n-1))
	  	  y1 = min(yy)  &  y2 = max(yy)
  		  polyfill, [x1,x2,x2,x1],[y1,y1,y2,y2],color=fill,noclip=0
	        endif
	      ;---------  Single lines  ---------------
	      endif else if pflag eq 0 then begin
 	        oplot,[1.,1.]*xx(i),[min(yy),max(yy)],linestyle=ls,$
		  color=clr, thick=thk
	      ;---------  Pointers  -------------------
	      endif else begin
	        dx = (!x.crange(1)-!x.crange(0))*wd/2.
		if !x.type eq 0 then dx=[-dx,0,dx]+xx(i) else $
		  dx=10^([-dx,0,dx]+alog10(xx(i)))
	        dy = (!y.crange(1)-!y.crange(0))*ht
		y1 = [0,dy,0]+!y.crange(0)
		y2 = [0,-dy,0]+!y.crange(1)
	        if bflag then polyfill,dx,y1,col=clr
	        if tflag then polyfill,dx,y2,col=clr
	      endelse
	    ;--------  Log Y axis  ----------
	    endif else begin
	      ;--------  Filled line pairs  -----------
	      if n_elements(fill) ne 0 then begin
	        if (i mod 2) eq 0 then begin
	          x1 = xx(i) & x2 = xx((i+1)<(n-1))
	  	  y1 = min(yy)  &  y2 = max(yy)
  		  polyfill, [x1,x2,x2,x1],10^[y1,y1,y2,y2],color=fill,noclip=0
	        endif
	      ;---------  Single lines  ---------------
	      endif else if pflag eq 0 then begin
 	        oplot,[1.,1.]*xx(i),10^[min(yy),max(yy)],linestyle=ls,$
		  color=clr, thick=thk
	      ;---------  Pointers  -------------------
	      endif else begin
	        dx = (!x.crange(1)-!x.crange(0))*wd/2.
		if !x.type eq 0 then dx=[-dx,0,dx]+xx(i) else $
		  dx=10^([-dx,0,dx]+alog10(xx(i)))
	        dy = (!y.crange(1)-!y.crange(0))*ht
		y1 = 10^([0,dy,0]+!y.crange(0))
		y2 = 10^([0,-dy,0]+!y.crange(1))
	        if bflag then polyfill,dx,y1,col=clr
	        if tflag then polyfill,dx,y2,col=clr
	      endelse
	    endelse  ; !y.type.
	  endfor
	endelse
 
	return
	end





;----------------------------------------------------------------------------
function arrayreplicate, array, dimensions, PUTITINFRONT=putitinfront

;a function to make a massive matrix with repeating entries of a given array

;INPUT
;array: the array you want to replicate
;dimensions: an array containing dimensions of desired output matrix 
;(but not including the array's dimension itself)
;KEYWORD
;putitinfront: the default is to make the final matrix dimension be the same 
;as the array.  If putitinfront is activated then the first matrix dimension 
;will be that of the array
;OUTPUT
;a matrix with desired dimensions
;EXAMPLE
;array has m entries, dimensions=[n,p,q,...,z]
;output is a matrix with dimensionality [n,p,q,...,z,m]
;if PUTITINFRONT is active then output has dimensionality  [m,n,p,q,...,z]

if keyword_set(PUTITINFRONT) then result = rebin(array, $
	[n_elements(array),dimensions],/sample) else result = transpose(rebin(array, $
	[n_elements(array),reverse(dimensions)]))
;stop
return, result
end
;----------------------------------------------------------------------------









pro linewindowz, linerest, transitionz,UNITS=units, LINENAME=linename, ZRANGE=zrange, $
SETXRANGE=setxrange,SPECIALWINDOWS=specialwindows,OUTFILE=outfile, PWVFILE=pwvfile $
                 , SHOW=show

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;       LINEWINDOWZ
;
; PURPOSE:
;
;       Solve for and plot the redshift range for which a set of rest
;       wavelengths all fall within a specified set of observing
;       window boundaries.
;
; EXPLANATION:
;
;       The program takes a set of rest wavelengths and a set of
;       observing window boundaries - most commonly this will be the
;       frequencies covered by the receiver bands. It searches
;       redshift space (over a range that can be specified by the
;       user) and reports the range of redshifts for which all of the
;       transitions can be observed. That is, if you start with a set
;       of transitions then this tells you the redshift range over
;       which ALMA can observe all of them.
;
;       This is useful for planning high redshift line observations
;       with ground based sub-mm/mm observatories like ALMA and want
;       to target a specific set of lines.  The program computes the
;       redshift ranges over which the full set of lines is accessible
;       and makes a .PS plot showing a visual representation of the
;       observed frame lines as a function of redshift with observing
;       windows and optionally atmospheric transmission overlaid.
;
; CALLING SEQUENCE:
;
;       LINEWINDOWZ, linerest, [transitionz, UNITS= , LINENAME= , ZRANGE= , 
;           SETXRANGE= , SPECIALWINDOWS= , OUTFILE= , PWVFILE= 
;
; INPUTS:
;
;       LINEREST - An array of rest frame wavelengths or frequencies of interest. By 
;                default this is interpreted as microns. (See UNITS.)
;       
; OPTIONAL OUTPUTS:
;
;       TRANSITIONZ - An array giving the redshifts at which the set of lines 
;                   transitions between being observable and not (or vice versa) 
;
; OPTIONAL INPUTS:
;
;       UNITS - string identifying the units of linerest. This also affects the 
;             interpretation of SETXRANGE, and SPECIALWINDOWS if supplied. By default 
;             these are interpreted as microns. Supported settings are 'GHz','microns', 
;             and equivalently, 'um'
;       LINENAME - a string array containing labels to identify lines on the output 
;                plot. By default they will be labelled by their supplied wavelength/
;                frequency.
;       
;       ZRANGE - A two element array listing the minimum and maximum redshift ranges of 
;              interest. By default this is set to [2,7]
;       
;       SETXRANGE - A two element array setting the x-range displayed on the output 
;                 plot, by default in microns (See UNITS)
;       
;       SPECIALWINDOWS - A two by N array where N is the number of observing windows of 
;                      interest. By default this is set to ALMA bands 3,4,6,7,8, and 9. 
;                      (Note that bands 6 and 7 abut each other and are considered as a
;                      single band for the purposes of this code.) 
;                      This keyword may be useful for avoiding specific ranges of poor 
;                      atmospheric transmission or if interested in other observing 
;                      facilities. If this keyword is supplied it is interpreted to be 
;                      in microns, unless otherwise specified by UNITS input.
;       
;       OUTFILE - string giving file name for output plot. Defaults to 
;               'linetransplot.ps'
;       
;       PWVFILE - string giving path and file name for an ascii data of atmospheric 
;               transmission. Data must be a two column file containing frequency (in 
;               GHz,) and transmissivity. Useful data sets can be downloaded from the 
;               ALMA Atmosphere Model calculator available on the ALMA science portal: 
;               http://almascience.eso.org/about-alma/weather/atmosphere-model
;               It will be read in using READCOL procedure. If no file is supplied the 
;               output plot will not show atmospheric transmission.
;
; EXAMPLES:
;       Consider redshifts over which the [CII] 157.7 micron line, [OI] 63, and the 
;       12CO (16-15) line are all observable
;
;       IDL> LINEWINDOWZ, [63.18372,157.74093,162.81163], transitionz
;       
;       On output, transitionz will be a two element array, [5.59011, 6.88174] 
;       indicating that all three lines are accessible within that redshift range.
;       The output plot 'linetransplot.ps' will be created in the current directory. 
;       The plot shows Redshift on the y-axis and observed frequency and wavelength on 
;       the x-axes. The lines of interest are plotted as black curved lines, identified 
;       by their wavelengths. Dashed vertical lines indicate boundaries of ALMA 
;       observing windows. Dashed horizontal lines indicate redshift range over which 
;       the lines fall within observing windows.
;       
;       Consider redshifts over which just [CII] 157.7 and 12CO (16-15) are available, 
;       but specify observing windows to avoid band 9, and also consider a larger 
;       redshift range of 1<z<8. To specify the observing windows in GHz, switch to GHz 
;       for the lines as well.
;       
;       IDL> LINEWINDOWZ, [1900.537,1841.345], transitionz, $
;            LINENAME=['[CII]','CO(16-15)'], ZRANGE=[1,8], SPECIALWINDOWS=[[84,116], $
;            [125,163],[211,373],[385,500]], UNITS='GHz', /SHOW
;
;       On output transitionz is a four element array indicating the lines are visible 
;       between redshifts 2.801 and 3.782 as well as 4.095 and 7.727. The plot has been 
;       scaled to the appropriate redshift, and the line labels are given as specified 
;       by LINENAME.
;
; PROCEDURES CALLED
;       This procedure makes use of procedures in the Goddard IDL Astrolib, 
;       and also VER and HOR, created by Johns Hopkins University/Applied Physics 
;       Laboratory and currently available in the University of Washington IDL library:
;       http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library09.html?VER
;       http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library09.html?HOR
;
; REVISION HISTORY:
;       Written         D. Brisbin                 November, 2013
;
; This software may be used, copied, or redistributed as long as it is not
; sold. This routine is provided as is without any express or implied 
; warranties and its results do not reflect the opinions or recommendations
; of NRAO.
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;tweakable parameters
nz=10000; how many test locations do we want to evaluate the lines at (z resolution will be ~ Delta(zrange)/zbins)



;defaults
if keyword_set(OUTFILE) eq 0 then outfile="linetransplot.ps"
if keyword_set(ZRANGE) eq 0 then zrange=[2,7]*1e0 else zrange=zrange*1e0
if keyword_set(SETXRANGE) eq 0 then myxrange=[50,1000.] else myxrange=SETXRANGE*1e0
if keyword_set(SPECIALWINDOWS) eq 0 then obswindows=[[84,116],[125,163],[211,373],[385,500],[602,720]] else obswindows=specialwindows
if keyword_set(LINENAME) eq 0 then linename=strtrim(string(linerest),2)



;;;constants
c=2.99792458d14	;um/s
dims=size(obswindows,/dimensions)
if N_elements(dims) eq 1 then nbands=1 else nbands=dims[1]
;nbands=(size(obswindows,/dimensions))[1]; Note that in this code, by default ALMA bands 6 and 7 are counted as one band since they abut each other



;The default units are microns. If the user inputs GHz convert it here
;(This block could be coded more simply, but I'm trying to make it flexible to build in many more unit types in the future)
;some conversion notes: this program was built with sub-mm astronomers in mind so by default the input lines are assumed to be in microns
;However, since atmospheric transmission plots usually increase in frequency left to right, by the time we get to plotting, we've switched
;to GHz, so by default the code expects myxrange and obswindows to be given in GHz and linerest to be given in microns.
;This block makes everything copacetic as long as the user is consistent (ie linerest, SETXRANGE, and SPECIALWINDOWS -if given- are all in the same units)

unitpossibilities=['microns','um','GHz']

if keyword_set(UNITS) then begin
	unittype=where(UNITS eq unitpossibilities,count)
	if count eq 0 then begin
		print, "acceptable unit options: ",unitpossibilities
		print, "(Default is microns)"
		stop
	endif
	if unittype eq 0 OR unittype eq 1 then begin
		if keyword_set(SETXRANGE) then myxrange=c/myxrange*1e-9
		if keyword_set(SPECIALWINDOWS) then obswindows=c/obswindows*1e-9
	endif
	if unittype eq 2 then begin
		linerest=c/(linerest*1e9)
	endif
endif else begin
	if keyword_set(SETXRANGE) then myxrange=c/myxrange*1e-9
	if keyword_set(SPECIALWINDOWS) then obswindows=c/obswindows*1e-9
endelse
;make sure, if the user has input SPECIALWINDOWS or SETXRANGE, that the order is correct 
;(lower number followed by higher number in frequency units)
if keyword_set(SETXRANGE) then myxrange=myxrange[sort(myxrange)]
if keyword_set(SPECIALWINDOWS) then begin
	reversebounds=where(obswindows[1,*]-obswindows[0,*] lt 0,reversecount)
	if reversecount ne 0 then begin
		for i=0, reversecount-1 do begin
			obswindows[*,reversebounds[i]]=reverse(obswindows[*,reversebounds[i]])
		endfor
	endif
endif
lineorder=sort(linerest)
linerest=linerest[lineorder]
linename=linename[lineorder]



;create a grid of z values to test
zstep=(zrange[1]-zrange[0])/(nz-1.0)
redshifts=findgen(nz)*zstep+zrange[0]



;shift the lines to redshifts in the grid and convert to using frequency units (GHz)
nline=n_elements(linerest)
gridline=arrayreplicate(linerest,nz)
gridz=arrayreplicate(redshifts,nline,/putitinfront)
gridline=gridline*(1.0+gridz)
;convert to frequencies
gridnu=(c/gridline)/(1e9)



;check to see if each line is visible
;a bit of logic follows to first make a grid of dimension (number of z bins x number of lines x number of windows) 
;which equals 1 where a line at a particular redshift falls in a particular band (and 0 if not)
evaluationgrid=intarr(nz,nline,nbands)
for i=0, nbands-1 do begin
	evaluationgrid[*,*,i]=gridnu ge obswindows[0,i] AND gridnu le obswindows[1,i]
endfor
;we now collapse the evaluationgrid along the window axis, resulting in a grid (number of z bins x number of lines) which equals 1 if the line is 
;visible in ANY band
if nbands gt 1 then evaluationgrid=total(evaluationgrid,3)



;this block adds up the evaluationgrid along the second axis to produce one integer for each zbin. If that integer equals the number of lines input,
;it means that all lines are visible at that redshift in some band (but not necessarily all in the same band.)
if nline gt 1 then n_linesvisible=total(evaluationgrid,2) else n_linesvisible=evaluationgrid; (if there's only 1 line to worry about, 
		;no need to sum evaluation grid.)
goodz=n_linesvisible eq nline; an array of length nz, equal to 1 if all lines are observable, 0 if not



;This block pulls out contiguous chunks of the array goodz to determine at what redshifts things transition from observable to unobservable (or vice-versa)
;foo=goodz[0:nz-2] eq goodz[1:nz-1]; will equal 0 wherever there's a change from observable to unobservable and vice versa
;changes=where(foo eq 0,n_transitionz)
;the z values that we want to note are then roughly:
;transitionz=redshifts[changes]
;actually, given the bin-by-bin way we compared things above, the best estimates for transition z's is midway between two bins:
;if n_transitionz ne 0 then transitionz=(redshifts[changes]+redshifts[changes+1])/2.0

reg = label_region(goodz)
n_transitionz = max(reg)*2
transitionz = fltarr(n_transitionz)
for i = 1, max(reg) do begin
   ind = where(reg eq i)
   transitionz[2*(i-1)] = min(redshifts[ind])
   transitionz[2*(i-1)+1] = max(redshifts[ind])
endfor

;Make a PS file

set_plot, 'ps'
device, file = outfile
;device, /inches,/landscape;, xsize = 7.0, ysize = 10., xoffset = 0.5, yoffset = 0.5
device, /inches,/landscape, xsize = 10.0, ysize = 7.5, xoffset = 0.4, yoffset = 11.
!P.Multi = [0,1,1]
device, /color
TVLCT, [0,255,0,0,0,255,255], [0,0,127,0,255,255,0], [0,0,0,255,255,0,255]

;plotting parameters
th=3
lth=7
ct=4
cs=1.3

nu=myxrange
trans=[-1.0,-1.0]
ystyle=5
;read in and plot atmospheric transmission model
;;NOTE: this transmission model is currently used only for aesthetics.
if keyword_set(PWVFILE) then begin 
	readcol,pwvfile, format="f,f", nu, trans,/silent
	ystyle=9
endif

plot, nu, trans, xrange=myxrange,xstyle=9,yrange=[0,1.1],ystyle=ystyle,xtitle="Freq. (GHz)",ytitle="Transmission",/nodata, $
	xthick=th, ythick=th, charthick=ct, charsize=cs

; ... shaded background
for i = 0, n_elements(obswindows)/2-1 do begin $
   xlo = obswindows[2*i]
   xhi = obswindows[2*i+1]
   ylo = 0
   yhi = 1.1
   polyfill, [xlo, xlo, xhi, xhi, xlo], [ylo, yhi, yhi, ylo, ylo], color=200, /clip
endfor

oplot, nu, trans,color=3

;plot the obswindows

; ... lines
ver, obswindows, linestyle=2, color=1, thick=lth

;set up the wavelength axis
umxrange=c/(!X.CRANGE*1e9)

;non_linear_axis,
xtickn=[300,350,400,500,600,800,1000,1500,2500]
xtickv=c/xtickn/1e9
xticks=n_elements(xtickn)
xtickn=strtrim(string(xtickn),2)
AXIS, XAXIS=1, xticks=xticks, xtickv=xtickv, xtickn=xtickn, XSTYLE = 1, xTITLE = 'Wavelength (microns)', ymargin=[5,5], $
		xthick=th, ythick=th, charthick=ct, charsize=cs

;set up the redshift axis
zyrange=!Y.CRANGE*(zrange[1]-zrange[0])/((!y.crange)[1]-(!y.crange)[0])+zrange[0]
AXIS, yAXIS=1, yRANGE = zyrange, ySTYLE = 1, yTITLE = 'Redshift', xmargin=[5,5],/save, $
		xthick=th, ythick=th, charthick=ct, charsize=cs;,ticklen=1,ygridstyle=1

if keyword_set(PWVFILE) eq 0 then AXIS, yAXIS=0, yRANGE = zyrange, ySTYLE = 1, xmargin=[5,5],/save, $
		xthick=th, ythick=th, charthick=ct, charsize=cs;,ticklen=1,ygridstyle=1

;plot redshifted line tracks
for i=0, nline-1 do begin

        for j=0, n_transitionz/2-1 do begin
           ind = where(gridz[*,i] ge transitionz[2*j] and $
                       gridz[*,i] le transitionz[2*j+1], ct)
           if ct gt 0 then $
              oplot, gridnu[ind,i], gridz[ind,i],color=100, thick=th*5.
        endfor

	oplot, gridnu[*,i],gridz[*,i], thick=th

endfor

;create line lables
;a kludgey loop to stagger the line labels:
for i=0, nline-1 do begin
	alternate=i mod 2
	index=round(nz*([0.97,0.94])[alternate])
	xyouts, [gridnu[index,i]],[redshifts[index]],linename[i], charsize=0.75, charthick=2.5
endfor

;plot and label redshift demarcations where lines are observable
if n_transitionz ne 0 then begin
	hor, transitionz, linestyle=2,thick=4;, color=4
	xyouts, myxrange[1]+(myxrange[1]-myxrange[0])/20.0,transitionz,string(format='(f7.4)',transitionz),charsize=.75, charthick=2.5;,string(format='f6.4',transitionz)

	print, "Valid redshift ranges to observe all lines:"
        print, "(Over the range ", min(zrange), " to ", max(zrange), ")"
	for i=0, n_transitionz/2-1 do begin
		print, transitionz[2*i], ' to ', transitionz[2*i+1]
	endfor
endif else begin
	if goodz[0] eq 1 then validity="valid." else validity="invalid."
	print, "No transition z's found."
	print, "Entire z range tested is "+validity
endelse


device, /close
set_plot, 'x', _extra = ex

if keyword_set(show) then $
   spawn, 'gv '+outfile+' &'

;stop



end
