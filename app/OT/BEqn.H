{
    std::array<std::shared_ptr<Field<vector> >, 4> dB; 
    for( int i=0; i<4; i++ )
    {   
        reuseTmp<vector> dBt( mesh );
        dB[i]= dBt();
    }

    reuseTmp<vector> Boldt( mesh );
    std::shared_ptr<Field<vector> > Bold( Boldt() );


    *Bold = B;

    for( int rk=0; rk<4; rk++ )
    {
        dB[rk] =
        time.dt() *
        (
           ex::curl( U^B ) //decompose?
           +ex::laplacian( DB, B ) 
        );
	
	Bstar = Bold;
	for( int j=0; j<=rk; j++ )
        {   
	      Bstar = Bstar+settings::RK4a()[rk][j]*dB[j];
        } 

    
	//pB poisson loop
//        Bstar.correctBoundaryConditions();
//        pBEqn.rhs(-ex::div(Bstar) / time.dt() / settings::RK4c()[rk] );
//        pBEqn.solve();
//        pB.correctBoundaryConditions();
//        B = Bstar + settings::RK4c()[rk]*time.dt() * ex::grad(pB);
//        B.correctBoundaryConditions();        

    }

Bstar.correctBoundaryConditions();
pBEqn.rhs(-ex::div(Bstar)/time.dt());
pBEqn.solve();
pB.correctBoundaryConditions();
B = Bstar + ex::grad(pB)*time.dt();
B.correctBoundaryConditions();


//    if( time.writelog() )
//    {   
       auto divB = ex::div(B);
       auto maxdivB = tools::maximum(divB);

        if( parallelCom::master() )
        {   
            std::cout<<"Magnetic field divergence: " << maxdivB<<std::endl;
        }   
//    }

}
