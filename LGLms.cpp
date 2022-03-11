//
// Created by ohmy on 2022-03-08.
//

#include "LGLms.h"

void LGLms::generate_constraints() {
    int nx = X.size1();
    double h = T/(N-1);
    J_ = 0;

    W = MX::vertcat({W, X(all,0)});
    for (int k = 0; k < N; ++k) {
        // Create collocation states
        auto o = ocp.variable(nx, n);


        auto xo = MX::horzcat({ X(all,k), o});
        W = MX::horzcat({W, o});
        MX u = MX::repmat(U(all, k), 1, n+1);
        // Defect constraints
        auto x_dot_ = f({{o}, {U(all, k)}});
        auto x_dot = MX::vertcat(x_dot_);
        auto F = 0.5*h  * x_dot ;
        auto g_d = mtimes(D(Slice(1,n+1), all), transpose(xo)) - transpose(F);

        // Shooting gap
        auto g_s = X(all,k+1) - o(all,n-1);

        g.push_back(g_s);
        g.push_back(MX::vec(g_d));



        // Cost
        auto l_k_ = J({{xo}, {U(all,k)}});
        auto l_k = MX::vertcat(l_k_);
        J_ += 0.5* h * mtimes(l_k , w); // quadrature

    }
    std::cout << W.size1() << std::endl;
    std::cout << W.size2() << std::endl;
    std::cout << X.size1() << std::endl;
    std::cout << X.size2() << std::endl;
}


//void LGLms::generate_constraints() {
//    int nx = X.size1();
//    double h = T/((N-2)/3-1);
//    J_ = 0;
//
//    int u_k = 0;
//    for (int k = 0; k < N-1; k+=n) {
//
//        // Create collocation states
//        auto o = X(all,Slice(k+1, k+n+1));
//        auto xo = MX::horzcat({ X(all,Slice(k, k+n+1))});
//        MX u = MX::repmat(U(all, k), 1, n+1);
//        // Defect constraints
//        auto x_dot_ = f({{o}, {U(all, u_k)}});
//        auto x_dot = MX::vertcat(x_dot_);
//        std::cout << x_dot.size1() << std::endl;
//        std::cout << x_dot.size2() << std::endl;
//        auto F = 0.5*h  * x_dot ;
//        auto g_d = mtimes(D(Slice(1,n+1), all), transpose(xo)) - transpose(F);
//        std::cout << g_d.size1() << std::endl;
//        std::cout << g_d.size2() << std::endl;
//        auto h = MX::vec(g_d);
//        std::cout << h.size1() << std::endl;
//        std::cout << h.size2() << std::endl;
//        // Shooting gap
//        auto g_s = X(all,k+1) - o(all,n-1);
//
//        //g.push_back(g_s);
//        g.push_back(MX::vec(g_d));
//
//
//
//        // Cost
//        auto l_k_ = J({{xo}, {U(all,u_k)}});
//        auto l_k = MX::vertcat(l_k_);
//        std::cout << l_k.size1() << std::endl;
//        std::cout << l_k.size2() << std::endl;
//        std::cout << w.size1() << std::endl;
//        std::cout << w.size2() << std::endl;
//        J_ += MX::sum2(l_k) ;//0.5* h * mtimes(l_k , w); // quadrature
//
//        u_k++;
//    }
//
//}

MX LGLms::integrated_cost(MX t0, MX tf, int N) {
    return J_;
}

void LGLms::create_LGL_params(int degree){
    if(degree == 1){
        this->tau = DM( {-1, 1 });
        this->D = DM ({
            {-0.500000000000000, 0.500000000000000, 0},
            {-0.500000000000000, 0.500000000000000, 0}   });
        this->D = this->D(Slice(0,2));
        this->w = DM({1, 1});

    }
    if(degree == 2){
        this->tau = DM( {-1, 0, 1 });
        this->D = DM ({
            {-1.50000000000000,	2,	-0.500000000000000},
            {-0.500000000000000,	0,	0.500000000000000},
            {0.500000000000000,	-2,	1.50000000000000} });
        this->w = DM({0.333333333333333, 1.333333333333333, 0.333333333333333});

    }
    if(degree == 3){
        this->tau = DM( {-1, -0.447213595499958, 0.447213595499958, 1 });
        this->D = DM ({
                        {-3,	4.04508497187474,	-1.54508497187474,	0.500000000000000},
                        {-0.809016994374947,	0, 1.11803398874990, -0.309016994374947},
                        {0.309016994374947,	-1.11803398874990,	0,	0.809016994374947},
                        {-0.500000000000000,	1.54508497187474,	-4.04508497187474,	3}  });
        this->w = DM({0.166666666666667, 0.833333333333333, 0.833333333333333, 0.166666666666667});
    }

    if( degree == 4 ){
        this->tau = DM({-1, -0.654653670707977, 0, 0.654653670707977, 1});
        this->D = DM( { {-5, 6.75650248872424, -2.66666666666667, 1.41016417794243, -0.500000000000000},
                        {-1.24099025303098,	0, 1.74574312188794 , -0.763762615825973, 0.259009746969017},
                        {0.375000000000000,	-1.33658457769545,	0,	1.33658457769545,	-0.375000000000000},
                        {-0.259009746969017,	0.763762615825973,	-1.74574312188794,	0,	1.24099025303098},
                        {0.500000000000000,	-1.41016417794243,	2.66666666666667,	-6.75650248872424,	5} });
        this->w = DM({0.100000000000000, 0.544444444444444, 0.711111111111111, 0.544444444444444, 0.100000000000000});
    }

    if( degree == 5 ){
        this->tau = DM({-1, -0.765055323929465, -0.285231516480645, 0.285231516480645, 0.765055323929465, 1});
        this->D = DM({ {-7.50000000000000,	10.1414159363197,	-4.03618727030535,	2.24468464817617,	-1.34991331419049,	0.500000000000000},
                       {-1.78636494833910,	2.22044604925031e-16,	2.52342677742946,	-1.15282815853593,	0.653547507429800,	-0.237781177984231},
                       {0.484951047853569,	-1.72125695283023,	0,	1.75296196636787,	-0.786356672223241,	0.269700610832039},
                       {-0.269700610832039,	0.786356672223241,	-1.75296196636787,	0,	1.72125695283023,	-0.484951047853569},
                       {0.237781177984231,	-0.653547507429800,	1.15282815853593,	-2.52342677742946,	0,	1.78636494833910},
                       {-0.500000000000000,	1.34991331419049,	-2.24468464817617,	4.03618727030535,	-10.1414159363197,	7.50000000000000} });
        this->w = DM({0.0666666666666667, 0.378474956297847, 0.554858377035487, 0.554858377035487, 0.378474956297847, 0.0666666666666667});
    }
}




//    auto a = casadi::collocation_points(3, "legendre");
//    DM C, B, D;
//    casadi::collocation_coeff(a, C,D,B);
//    std::cout << DM({a}) << std::endl;
//    std::cout << std::endl;
//    std::cout << C << std::endl;
//    std::cout << std::endl;
//    std::cout << D << std::endl;
//    std::cout << std::endl;
//    std::cout << B << std::endl;
//
//    for (int k = 0; k < N; ++k) {
//        // Create collocation states
//        auto o = ocp.variable(nx, n);
//
//
//        auto xo = MX::horzcat({  o});
//        MX u = MX::repmat(U(all, k), 1, n+1);
//        // Defect constraints
//        auto x_dot_ = f({{xo}, {U(all, k)}});
//        auto x_dot = MX::vertcat(x_dot_);
//        std::cout << x_dot.size1() << std::endl;
//        std::cout << x_dot.size2() << std::endl;
//        auto F = 0.5*h  * x_dot ;
//        auto g_d = mtimes(C, transpose(xo)) - transpose(F);
//        std::cout << g_d.size1() << std::endl;
//        std::cout << g_d.size2() << std::endl;
//        auto h = MX::vec(g_d);
//        std::cout << h.size1() << std::endl;
//        std::cout << h.size2() << std::endl;
//        // Shooting gap
//        auto g_s = X(all,k+1) - o(all,n-1);
//
//        g.push_back(g_s);
//        g.push_back(MX::vec(g_d));
//
//
//
//        // Cost
//        auto l_k_ = J({{xo}, {U(all,k)}});
//        auto l_k = MX::vertcat(l_k_);
//        std::cout << l_k.size1() << std::endl;
//        std::cout << l_k.size2() << std::endl;
//        std::cout << w.size1() << std::endl;
//        std::cout << w.size2() << std::endl;
//        J_ += MX::sum2(l_k) ;//0.5* h * mtimes(l_k , w); // quadrature
//    }

/*for (int k = 0; k < N; ++k) {
    // Create collocation states
    auto o = ocp.variable(nx, n);
    // Shooting gap
    auto g_s = X(all,k) - o(all,0);

    auto xo = MX::horzcat({ o, X(all,k+1)});
    //MX u = MX::repmat(U(all, k), 1, n+1);
    // Defect constraints
    auto x_dot_ = f({{xo}, {U(all, k)}});
    auto x_dot = MX::vertcat(x_dot_);
    std::cout << x_dot.size1() << std::endl;
    std::cout << x_dot.size2() << std::endl;
    auto F = 0.5*h  * x_dot ;
    auto g_d = mtimes(D, transpose(xo)) - transpose(F);

    //g.push_back(g_s);
    g.push_back(MX::vec(g_d));

    // Cost
    auto l_k_ = J({{xo}, {U(all,k)}});
    auto l_k = MX::vertcat(l_k_);
    std::cout << l_k.size1() << std::endl;
    std::cout << l_k.size2() << std::endl;
    std::cout << w.size1() << std::endl;
    std::cout << w.size2() << std::endl;
    J_ += 0.5* h * mtimes(l_k , w); // quadrature
}*/