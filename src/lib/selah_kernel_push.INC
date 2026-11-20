   ! Loops of the particle push of SELAH (selah_pushers), one per pusher, for the layout of the macros PP(k) (component k
   ! of particle i), EE(k) and BB(k) (the fields at particle i). Included once per layout. Every loop is one particle per
   ! iteration, scalars only (no array temporaries on the device); the caller's arrays, device resident in a device
   ! build. Notation: qh = q dt/(2m), t = qh B (/gamma), a cross product written out component by component.
   select case(pusher)
   case(PUSHER_BORIS) ! non-relativistic: components iv..iv+2 hold v
      !$acc parallel loop independent DEVICEVAR(p, e, b) private(qh, m1, m2, m3, t1, t2, t3, s1, s2, s3, w1, w2, w3, tt)
      !$omp OMPLOOP DEVICEPTR(p, e, b) private(qh, m1, m2, m3, t1, t2, t3, s1, s2, s3, w1, w2, w3, tt)
      do i=1, n
         qh = 0.5_R_P * dt * PP(iq) / PP(im)
         m1 = PP(iv) + qh*EE(1) ; m2 = PP(iv+1) + qh*EE(2) ; m3 = PP(iv+2) + qh*EE(3)
         t1 = qh*BB(1) ; t2 = qh*BB(2) ; t3 = qh*BB(3)
         tt = 2._R_P / (1._R_P + t1*t1 + t2*t2 + t3*t3)
         s1 = tt*t1 ; s2 = tt*t2 ; s3 = tt*t3
         w1 = m1 + (m2*t3 - m3*t2) ; w2 = m2 + (m3*t1 - m1*t3) ; w3 = m3 + (m1*t2 - m2*t1)
         m1 = m1 + (w2*s3 - w3*s2) ; m2 = m2 + (w3*s1 - w1*s3) ; m3 = m3 + (w1*s2 - w2*s1)
         PP(iv) = m1 + qh*EE(1) ; PP(iv+1) = m2 + qh*EE(2) ; PP(iv+2) = m3 + qh*EE(3)
         PP(ix)   = PP(ix)   + dt * PP(iv)
         PP(ix+1) = PP(ix+1) + dt * PP(iv+1)
         PP(ix+2) = PP(ix+2) + dt * PP(iv+2)
      enddo
   case(PUSHER_BORIS_RELATIVISTIC) ! components iv..iv+2 hold u = gamma v
      !$acc parallel loop independent DEVICEVAR(p, e, b) &
      !$acc private(qh, m1, m2, m3, t1, t2, t3, s1, s2, s3, w1, w2, w3, tt, g)
      !$omp OMPLOOP DEVICEPTR(p, e, b) private(qh, m1, m2, m3, t1, t2, t3, s1, s2, s3, w1, w2, w3, tt, g)
      do i=1, n
         qh = 0.5_R_P * dt * PP(iq) / PP(im)
         m1 = PP(iv) + qh*EE(1) ; m2 = PP(iv+1) + qh*EE(2) ; m3 = PP(iv+2) + qh*EE(3)
         g = sqrt(1._R_P + (m1*m1 + m2*m2 + m3*m3) * c2i)
         t1 = qh*BB(1)/g ; t2 = qh*BB(2)/g ; t3 = qh*BB(3)/g
         tt = 2._R_P / (1._R_P + t1*t1 + t2*t2 + t3*t3)
         s1 = tt*t1 ; s2 = tt*t2 ; s3 = tt*t3
         w1 = m1 + (m2*t3 - m3*t2) ; w2 = m2 + (m3*t1 - m1*t3) ; w3 = m3 + (m1*t2 - m2*t1)
         m1 = m1 + (w2*s3 - w3*s2) ; m2 = m2 + (w3*s1 - w1*s3) ; m3 = m3 + (w1*s2 - w2*s1)
         PP(iv) = m1 + qh*EE(1) ; PP(iv+1) = m2 + qh*EE(2) ; PP(iv+2) = m3 + qh*EE(3)
         g = sqrt(1._R_P + (PP(iv)**2 + PP(iv+1)**2 + PP(iv+2)**2) * c2i)
         PP(ix)   = PP(ix)   + dt * PP(iv)   / g
         PP(ix+1) = PP(ix+1) + dt * PP(iv+1) / g
         PP(ix+2) = PP(ix+2) + dt * PP(iv+2) / g
      enddo
   case(PUSHER_VAY) ! components iv..iv+2 hold u = gamma v
      !$acc parallel loop independent DEVICEVAR(p, e, b) &
      !$acc private(qh, m1, m2, m3, t1, t2, t3, w1, w2, w3, tt, g, us, sg)
      !$omp OMPLOOP DEVICEPTR(p, e, b) private(qh, m1, m2, m3, t1, t2, t3, w1, w2, w3, tt, g, us, sg)
      do i=1, n
         qh = 0.5_R_P * dt * PP(iq) / PP(im)
         g = sqrt(1._R_P + (PP(iv)**2 + PP(iv+1)**2 + PP(iv+2)**2) * c2i)
         w1 = PP(iv) / g ; w2 = PP(iv+1) / g ; w3 = PP(iv+2) / g            ! v^{n-1/2}
         m1 = PP(iv)   + 2._R_P*qh*EE(1) + qh*(w2*BB(3) - w3*BB(2))          ! u'
         m2 = PP(iv+1) + 2._R_P*qh*EE(2) + qh*(w3*BB(1) - w1*BB(3))
         m3 = PP(iv+2) + 2._R_P*qh*EE(3) + qh*(w1*BB(2) - w2*BB(1))
         t1 = qh*BB(1) ; t2 = qh*BB(2) ; t3 = qh*BB(3)                      ! tau
         tt = t1*t1 + t2*t2 + t3*t3
         us = (m1*t1 + m2*t2 + m3*t3) * ci                               ! u* = u'.tau/c
         sg = 1._R_P + (m1*m1 + m2*m2 + m3*m3) * c2i - tt                ! gamma(u')^2 - |tau|^2
         g = sqrt(0.5_R_P * (sg + sqrt(sg*sg + 4._R_P * (tt + us*us))))  ! gamma^{n+1/2}
         t1 = t1/g ; t2 = t2/g ; t3 = t3/g
         tt = 1._R_P / (1._R_P + t1*t1 + t2*t2 + t3*t3)
         us = m1*t1 + m2*t2 + m3*t3
         PP(iv)   = tt * (m1 + us*t1 + (m2*t3 - m3*t2))
         PP(iv+1) = tt * (m2 + us*t2 + (m3*t1 - m1*t3))
         PP(iv+2) = tt * (m3 + us*t3 + (m1*t2 - m2*t1))
         g = sqrt(1._R_P + (PP(iv)**2 + PP(iv+1)**2 + PP(iv+2)**2) * c2i)
         PP(ix)   = PP(ix)   + dt * PP(iv)   / g
         PP(ix+1) = PP(ix+1) + dt * PP(iv+1) / g
         PP(ix+2) = PP(ix+2) + dt * PP(iv+2) / g
      enddo
   case(PUSHER_HIGUERA_CARY) ! components iv..iv+2 hold u = gamma v
      !$acc parallel loop independent DEVICEVAR(p, e, b) &
      !$acc private(qh, m1, m2, m3, t1, t2, t3, w1, w2, w3, tt, g, us, sg)
      !$omp OMPLOOP DEVICEPTR(p, e, b) private(qh, m1, m2, m3, t1, t2, t3, w1, w2, w3, tt, g, us, sg)
      do i=1, n
         qh = 0.5_R_P * dt * PP(iq) / PP(im)
         m1 = PP(iv) + qh*EE(1) ; m2 = PP(iv+1) + qh*EE(2) ; m3 = PP(iv+2) + qh*EE(3) ! u-
         t1 = qh*BB(1) ; t2 = qh*BB(2) ; t3 = qh*BB(3)                             ! tau
         tt = t1*t1 + t2*t2 + t3*t3
         us = (m1*t1 + m2*t2 + m3*t3) * ci                                      ! u* = u-.tau/c
         sg = 1._R_P + (m1*m1 + m2*m2 + m3*m3) * c2i - tt                       ! gamma(u-)^2 - |tau|^2
         g = sqrt(0.5_R_P * (sg + sqrt(sg*sg + 4._R_P * (tt + us*us))))         ! gamma(u bar)
         t1 = t1/g ; t2 = t2/g ; t3 = t3/g
         tt = 1._R_P / (1._R_P + t1*t1 + t2*t2 + t3*t3)
         us = m1*t1 + m2*t2 + m3*t3
         w1 = tt * (m1 + us*t1 + (m2*t3 - m3*t2))                               ! u bar
         w2 = tt * (m2 + us*t2 + (m3*t1 - m1*t3))
         w3 = tt * (m3 + us*t3 + (m1*t2 - m2*t1))
         PP(iv)   = w1 + (w2*t3 - w3*t2) + qh*EE(1)                               ! u bar + u bar x t + qh E
         PP(iv+1) = w2 + (w3*t1 - w1*t3) + qh*EE(2)
         PP(iv+2) = w3 + (w1*t2 - w2*t1) + qh*EE(3)
         g = sqrt(1._R_P + (PP(iv)**2 + PP(iv+1)**2 + PP(iv+2)**2) * c2i)
         PP(ix)   = PP(ix)   + dt * PP(iv)   / g
         PP(ix+1) = PP(ix+1) + dt * PP(iv+1) / g
         PP(ix+2) = PP(ix+2) + dt * PP(iv+2) / g
      enddo
   endselect
