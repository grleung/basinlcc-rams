!##############################################################################
Subroutine thermo ()

use mem_grid
use mem_basic
use mem_micro
use mem_radiate
use mem_scratch
use micphys, only:level
use node_mod, only:mzp,mxp,myp

implicit none

if (level .le. 1) then

   CALL drythrm (mzp,mxp,myp,1,mxp,1,myp  &
      ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1)   &
      ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv(1,1,1),level)

elseif (level .eq. 2) then

   CALL satadjst (mzp,mxp,myp,1,mxp,1,myp  &
      ,basic_g(ngrid)%pp(1,1,1)  ,scratch%scr1(1)             &
      ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1) &
      ,basic_g(ngrid)%pi0(1,1,1)                              &
      ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv(1,1,1)    &
      ,micro_g(ngrid)%rcp(1,1,1) ,scratch%scr2(1))

elseif (level .eq. 3) then

   CALL wetthrm3 (mzp,mxp,myp,1,mxp,1,myp  &
     ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1)  &
     ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv   (1,1,1)  &
     ,micro_g(ngrid)%rcp(1,1,1) ,micro_g(ngrid)%rrp  (1,1,1)  &
     ,micro_g(ngrid)%rpp(1,1,1) ,micro_g(ngrid)%rsp  (1,1,1)  &
     ,micro_g(ngrid)%rap(1,1,1) ,micro_g(ngrid)%rgp  (1,1,1)  &
     ,micro_g(ngrid)%rhp(1,1,1) ,micro_g(ngrid)%q6   (1,1,1)  &
     ,micro_g(ngrid)%q7(1,1,1)  ,micro_g(ngrid)%rdp  (1,1,1)  &
     ,basic_g(ngrid)%pi0(1,1,1) ,basic_g(ngrid)%pp   (1,1,1)  &
     ! GRL 2024-03-22 added variables for RCEMIP
     ,basic_g(ngrid)%pres(1,1,1),basic_g(ngrid)%tmpt (1,1,1)  &
     ,basic_g(ngrid)%rhl(1,1,1) ,basic_g(ngrid)%rhi(1,1,1) &
     ,basic_g(ngrid)%rsatvl(1,1,1),basic_g(ngrid)%rsatvi(1,1,1)  &
     ,basic_g(ngrid)%thte(1,1,1),basic_g(ngrid)%tcon (1,1,1)  &
     ,basic_g(ngrid)%cfrac(1,1,1),basic_g(ngrid)%clrflag(1,1) &
     ,radiate_g(ngrid)%fthrd(1,1,1),radiate_g(ngrid)%fthrdlw (1,1,1)  &
     ,radiate_g(ngrid)%fthrdsw(1,1,1),radiate_g(ngrid)%clrhr (1,1,1)  &
     ,radiate_g(ngrid)%clrhrlw(1,1,1),radiate_g(ngrid)%clrhrsw (1,1,1)  &
     ,radiate_g(ngrid)%swup(1,1,1),radiate_g(ngrid)%swdn (1,1,1)  &
     ,radiate_g(ngrid)%lwup(1,1,1),radiate_g(ngrid)%lwdn (1,1,1)  &
     ,radiate_g(ngrid)%clrswup(1,1,1),radiate_g(ngrid)%clrswdn (1,1,1)  &
     ,radiate_g(ngrid)%clrlwup(1,1,1),radiate_g(ngrid)%clrlwdn (1,1,1)  &
     )

elseif (level .eq. 4) then
      CALL wetthrm3_bin (mzp,mxp,myp,1,mxp,1,myp  &
     ,basic_g(ngrid)%pi0(1,1,1) ,basic_g(ngrid)%pp   (1,1,1)  &
     ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1)  &
     ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv   (1,1,1)  &
     ,micro_g(ngrid)%ffcd(1,1,1,1) ,micro_g(ngrid)%ffic  (1,1,1,1)  &
     ,micro_g(ngrid)%ffip(1,1,1,1) ,micro_g(ngrid)%ffid  (1,1,1,1)  &
     ,micro_g(ngrid)%ffsn(1,1,1,1) ,micro_g(ngrid)%ffgl  (1,1,1,1)  &
     ,micro_g(ngrid)%ffhl(1,1,1,1))
else

   stop 'Thermo option not supported...LEVEL out of bounds'

endif

return
END SUBROUTINE thermo

!##############################################################################
Subroutine drythrm (m1,m2,m3,ia,iz,ja,jz,thil,theta,rt,rv,level)

! This routine calculates theta and rv for the case where no condensate is
! allowed.

implicit none

integer m1,m2,m3,ia,iz,ja,jz,i,j,k,level
real thil(m1,m2,m3),theta(m1,m2,m3),rt(m1,m2,m3),rv(m1,m2,m3)

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         theta(k,i,j) = thil(k,i,j)
      enddo
      if (level .eq. 1) then
         do k = 1,m1
            rv(k,i,j) = rt(k,i,j)
         enddo
      endif
   enddo
enddo

return
END SUBROUTINE drythrm

!##############################################################################
Subroutine satadjst (m1,m2,m3,ia,iz,ja,jz  &
   ,pp,p,thil,theta,pi0,rtp,rv,rcp,rvls)

! This routine diagnoses theta, rv, and rcp using a saturation adjustment
! for the case when water is in the liquid phase only

use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real :: pp(m1,m2,m3),p(m1,m2,m3),thil(m1,m2,m3),theta(m1,m2,m3)  &
   ,pi0(m1,m2,m3),rtp(m1,m2,m3),rv(m1,m2,m3)  &
   ,rcp(m1,m2,m3),rvls(m1,m2,m3)
real :: t(m1,m2,m3)
real, external :: rslf
integer :: i,j,k,iterate
real :: picpi,til,tt

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         picpi = (pi0(k,i,j) + pp(k,i,j)) * cpi
         p(k,i,j) = p00 * picpi ** 3.498
         til = thil(k,i,j) * picpi
         t(k,i,j) = til

         do iterate = 1,20
            rvls(k,i,j) = rslf(p(k,i,j),t(k,i,j))
            rcp(k,i,j) = max(rtp(k,i,j) - rvls(k,i,j),0.)
            tt = 0.7 * t(k,i,j) + 0.3 * til  &
               * (1. + alvl * rcp(k,i,j)  &
               / (cp * max(t(k,i,j),253.)))
            if (abs(tt - t(k,i,j)) .le. 0.001) go to 1
            t(k,i,j) = tt
         enddo
1             continue
         rv(k,i,j) = rtp(k,i,j) - rcp(k,i,j)
         theta(k,i,j) = t(k,i,j) / picpi
      enddo
   enddo
enddo

return
END SUBROUTINE satadjst

!##############################################################################
Subroutine wetthrm3 (m1,m2,m3,ia,iz,ja,jz                   &
   ,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rap,rgp,rhp,q6,q7,rdp  &
   ,pi0,pp,pres,tmpt,rhl,rhi,rsatvl,rsatvi,thte,tcon,cfrac,clrflag       &
   ,fthrd,fthrdlw,fthrdsw,clrhr,clrhrlw,clrhrsw             &
   ,swup,swdn,lwup,lwdn,clrswup,clrswdn,clrlwup,clrlwdn     &
   )

! This routine calculates theta and rv for "level 3 microphysics"
! given prognosed theta_il, cloud, rain, pristine ice, snow, aggregates,
! graupel, hail, q6, and q7.

use rconstants
use micphys, only:tair,til,qhydm,rliq,rice,jnmb

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: tcoal,fracliq,tairstr,es,tcon_cond,sat_cond,pi,tmpl,rsatv,tlcl
real, dimension(m1) :: picpi
real, dimension(m1,m2,m3) :: pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rap &
                            ,rgp,rhp,q6,q7,rdp,pres,tmpt,rhl,rhi,rsatvl,rsatvi,thte,tcon,cfrac &
                            ,fthrd,fthrdlw,fthrdsw,clrhr,clrhrlw,clrhrsw,clrswdn,clrswup,clrlwdn,clrlwup,swdn,swup,lwdn,lwup
real, dimension(m2,m3)  :: clrflag
real, external :: rsif,rslf

tcon_cond = 1.e-5 !condition for cloud fraction to be 1

do j = ja,jz
   do i = ia,iz

         do k = 1,m1
            picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
            tair(k) = theta(k,i,j) * picpi(k)
            til(k) = thp(k,i,j) * picpi(k)
            rliq(k) = 0.
         rice(k) = 0.
      enddo

      if (jnmb(1) .ge. 1) then
         do k = 1,m1
            rliq(k) = rliq(k) + rcp(k,i,j)
         enddo
      endif

      if (jnmb(2) .ge. 1) then
         do k = 1,m1
            rliq(k) = rliq(k) + rrp(k,i,j)
         enddo
      endif

      if (jnmb(3) .ge. 1) then
         do k = 1,m1
            rice(k) = rice(k) + rpp(k,i,j)
         enddo
      endif

      if (jnmb(4) .ge. 1) then
         do k = 1,m1
            rice(k) = rice(k) + rsp(k,i,j)
         enddo
      endif

      if (jnmb(5) .ge. 1) then
         do k = 1,m1
            rice(k) = rice(k) + rap(k,i,j)
         enddo
      endif

      if (jnmb(6) .ge. 1) then
         do k = 1,m1
            CALL qtc (q6(k,i,j),tcoal,fracliq)
            rliq(k) = rliq(k) + rgp(k,i,j) * fracliq
            rice(k) = rice(k) + rgp(k,i,j) * (1. - fracliq)
         enddo
      endif

      if (jnmb(7) .ge. 1) then
         do k = 1,m1
            CALL qtc (q7(k,i,j),tcoal,fracliq)
            rliq(k) = rliq(k) + rhp(k,i,j) * fracliq
            rice(k) = rice(k) + rhp(k,i,j) * (1. - fracliq)
         enddo
      endif

      if (jnmb(8) .ge. 1) then
         do k = 1,m1
            rliq(k) = rliq(k) + rdp(k,i,j)
         enddo
      endif

      do k = 1,m1
         qhydm(k) = alvl * rliq(k) + alvi * rice(k)
         rv(k,i,j) = rtp(k,i,j) - rliq(k) - rice(k)
      enddo
      
      do k = 1,m1
         if (tair(k) .gt. 253.) then
            tairstr = 0.5 * (til(k)  &
               + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
         else
            tairstr = til(k) * (1. + qhydm(k) * cp253i)
         endif
         theta(k,i,j) = tairstr / picpi(k)
      enddo

      ! GRL 2024-03-22 adding RCEMIP output

      ! clearsky flag starts as 1 (yes clearsky), then gets set to 0 if any of points within column are cloudy
      clrflag(i,j) = 1.

      do k = 1,m1
         pi = pi0(k,i,j) + pp(k,i,j)
         tmpt(k,i,j) = theta(k,i,j)*pi*cpi
         pres(k,i,j) = ( pi/cp )**cpor*p00 

         ! calculate saturation and relative humidity with respect to liquid and ice
         rsatvl(k,i,j) = rslf(pres(k,i,j),tmpt(k,i,j))
         rsatvi(k,i,j) = rsif(pres(k,i,j),tmpt(k,i,j))

         rhl(k,i,j) = 100. * rv(k,i,j)/rsatvl(k,i,j)
         rhi(k,i,j) = 100. * rv(k,i,j)/rsatvi(k,i,j)

         if (tmpt(k,i,j)>=273.15) then
            rsatv = rsatvl(k,i,j)
         else
            rsatv = rsatvi(k,i,j)
         endif 

         ! define alterantive cloud condition as 1% of saturation mixing ratio relative to liquid or ice (use lower threshold, which is ice, above 0C)
         sat_cond = 0.01*rsatv

         ! using eq 22 and 43 of Bolton 1980, where eq 43 is converted to take rv in kg/kg rather than g/kg
         if (rhl(k,i,j)>=1.e-28) then 
            tlcl = (1/((1/(tmpt(k,i,j)-55.))+(log(rhl(k,i,j)/100.)/2840.)))+55.
            thte(k,i,j)=tmpt(k,i,j) * ((p00/pres(k,i,j))**(0.2854*(1-(.28*rv(k,i,j))))) * exp(((3.376/tlcl)-.00254) * 1.e3 * rv(k,i,j) * (1+(0.81*rv(k,i,j))))
         else
            ! if rh is very low then the above equation converges to this, so just using this to avoid small number errors
            thte(k,i,j)=tmpt(k,i,j) * ((p00/pres(k,i,j))**0.2854) 
         endif 

         tcon(k,i,j) = rtp(k,i,j) - rv(k,i,j)

         if ((tcon(k,i,j) .ge. tcon_cond) .or. (tcon(k,i,j) .ge. sat_cond)) then
            cfrac(k,i,j) = 1.
            clrflag(i,j) = 0.
         else
            cfrac(k,i,j) = 0.
         endif
      enddo
      

      ! GRL 2024-03-25 add terms for clear-sky radiative heating
      ! check that column is clrsky (cfrac=0 throughout column), then copy heating rate terms

      if (clrflag(i,j) .eq. 1) then
         do k = 1,m1
            clrhr(k,i,j) = fthrd(k,i,j)
            clrhrlw(k,i,j) = fthrdlw(k,i,j)
            clrhrsw(k,i,j) = fthrdsw(k,i,j)
            clrswdn(k,i,j) = swdn(k,i,j)
            clrswup(k,i,j) = swup(k,i,j)
            clrlwdn(k,i,j) = lwdn(k,i,j)
            clrlwup(k,i,j) = lwup(k,i,j)
         enddo
      else
         do k = 1,m1
            clrhr(k,i,j) = 0.
            clrhrlw(k,i,j) = 0.
            clrhrsw(k,i,j) = 0.
            clrswdn(k,i,j) = 0.
            clrswup(k,i,j) = 0.
            clrlwdn(k,i,j) = 0.
            clrlwup(k,i,j) = 0.
         enddo
      endif



   enddo
enddo

return
END SUBROUTINE wetthrm3

!##############################################################################
Subroutine wetthrm3_bin (m1,m2,m3,ia,iz,ja,jz  &
   ,pi0,pp,thp,theta,rtp,rv,ffcd,ffic,ffip,ffid,ffsn,ffgl,ffhl)

! This routine calculates theta and rv for "level 4 microphysics"

use rconstants
use micro_prm, only:nkr,iceprocs
use micphys, only:ipris,igraup,ihail,tair,til,qhydm

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,ngrid,izero,i,j,k
real, dimension(m1,m2,m3) :: pi0,pp,thp,theta,rtp,rv
real, dimension(m1,m2,m3,nkr) :: ffcd,ffic,ffip,ffid,ffsn,ffgl,ffhl
real, dimension(m1) :: picpi
real, allocatable, dimension(:,:,:) :: rliq,rice
real, allocatable, dimension(:,:,:,:) :: ice_bins
real :: tcoal,fracliq,tairstr

allocate(ice_bins(m1,m2,m3,nkr))
allocate(rliq(m1,m2,m3))
allocate(rice(m1,m2,m3))

izero=0
CALL sum_bins (ffcd,rliq,m1,m2,m3,1,nkr,izero)

if(iceprocs.eq.1) then
   ice_bins = ffsn
   if (ipris == 1 .or. ipris >= 4) ice_bins = ice_bins + ffic
   if (ipris == 2 .or. ipris >= 4) ice_bins = ice_bins + ffip
   if (ipris >= 3) ice_bins = ice_bins + ffid
   if (igraup > 0) ice_bins = ice_bins + ffgl
   if (ihail > 0) ice_bins = ice_bins + ffhl
   izero=0
   CALL sum_bins (ice_bins,rice,m1,m2,m3,1,nkr,izero)
else
   rice=0.
endif

do j = ja,jz
   do i = ia,iz

      do k = 1,m1
         picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
         tair(k) = theta(k,i,j) * picpi(k)
         til(k) = thp(k,i,j) * picpi(k)

         qhydm(k) = alvl * rliq(k,i,j) + alvi * rice(k,i,j)
         rv(k,i,j) = rtp(k,i,j) - rliq(k,i,j) - rice(k,i,j)
         if (tair(k) .gt. 253.) then
            tairstr = 0.5 * (til(k)  &
               + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
         else
            tairstr = til(k) * (1. + qhydm(k) * cp253i)
         endif
         theta(k,i,j) = tairstr / picpi(k)
      enddo

   enddo
enddo

deallocate(ice_bins,rliq,rice)

return
END SUBROUTINE wetthrm3_bin
