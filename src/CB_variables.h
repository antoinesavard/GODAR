!=======================================================================
!   Common block Variablesn dynamic variables
!=======================================================================

    double precision ::           &   
                x       (n),      & ! x positions
                y       (n),      & ! y positions
                r       (n),      & ! radii
                h       (n),      & ! thicknesses
                u       (n),      & ! u velocity in x axis
                v       (n),      & ! v velocity in y axis
                f       (n,2),    & ! perp/para forces collision frame
                d       (n,n),    & ! perp damping force
                s       (n,n),    & ! storage array for increments force
                tfx     (n),      & ! total force in x
                tfy     (n),      & ! total force in y
                sx      (n,n),    & !
                sy      (n,n),    & !
                fsx     (n,n),    & !
                fsy     (n,n),    & !
                teta    (n),      & !
                ome     (n),      & !
                m       (n),      & !
                fw      (n,n),    & !
                fwx     (n,n),    & !
                fwy     (n,n),    & !
                w       (n,n),    & !
                wx      (n,n),    & !
                wy      (n,n),    & !
                dx      (n,n),    & !
                dy      (n,n),    & !
                dafx    (n),      & !
                dafy    (n),      & !
                wafx    (n),      & !
                wafy    (n),      & !
                l       (n,n),    & !
                lx      (n,n),    & !
                ly      (n,n),    & !
                mass    (n)         ! mass

    common/variables/        &
                x       ,    &
                y       ,    &
                r       ,    &
                h       ,    &
                u       ,    &
                v       ,    &
                f       ,    &
                d       ,    &
                s       ,    &
                tfx     ,    &
                tfy     ,    &
                sx      ,    &
                sy      ,    &
                fsx     ,    &
                fsy     ,    &
                teta    ,    &
                ome     ,    &
                m       ,    &
                fw      ,    &
                fwx     ,    &
                fwy     ,    &
                w       ,    &
                wx      ,    &
                wy      ,    &
                dx      ,    &
                dy      ,    &
                dafx    ,    &
                dafy    ,    &
                wafx    ,    &
                wafy    ,    &
                l       ,    &
                lx      ,    &
                ly      ,    &
                mass            