!=======================================================================
!   Common block Variablesn dynamic variables
!=======================================================================

    double precision ::           &   
                xy      (n,4),    &
                vu      (n,2),    &
                f       (n,n),    &
                d       (n,n),    &
                s       (n,n),    &
                tfx     (n),      &
                tfy     (n),      &
                sx      (n,n),    &
                sy      (n,n),    &
                fsx     (n,n),    &
                fsy     (n,n),    &
                teta    (n),      &
                ome     (n),      &
                m       (n),      &
                fw      (n,n),    &
                fwx     (n,n),    &
                fwy     (n,n),    &
                w       (n,n),    &
                wx      (n,n),    &
                wy      (n,n),    &
                dx      (n,n),    &
                dy      (n,n),    &
                dafx    (n),      &
                dafy    (n),      &
                wafx    (n),      &
                wafy    (n),      &
                l       (n,n),    &
                lx      (n,n),    &
                ly      (n,n)   

    common/variables/        &
                xy      ,    &
                vu      ,    &
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
                ly          