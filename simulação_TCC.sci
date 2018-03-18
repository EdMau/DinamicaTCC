


function z = gV(t,V,K,J,T,M,g,gama,D) //EDO usada para o cálculo da velocidade do foguete no referencial do solo
    z = (T/M)*cos(K)-g*sin(gama)-(D/M);
endfunction

function z = gK(t,V,K,J) //EDO usada para o cálculo do ângulo de ataque pt.1
    z = J;
endfunction

function z = gJ(t,V,K,J,T,I,L,VW,gama,Xcg,d1) //EDO usada para o cálculo do ângulo de ataque pt.2
    z = (T/I)*Xcg*sin(k)-(L/I)*d1-(VW/I)*d1*cos(gama);
endfunction

function z = lift(t,densidade,Al,Cl,V) //força lateral (lift)
    z = 0.5*densidade*Al*Cl*(V^2);
endfunction

function z = ar_frontal(t,A0,K) //área forntal efetiva do foguete (área frontal com relação ao ângulo de ataque)
    z = A0*cos(k);
endfunction

function z = drag(t,densidade,Af,Cd,V) //força de arrasto (drag)
    z = 0.5*densidade*Af*Cd*(V^2);
endfunction

function z = V_ar(t,V,VW,K) //velocidade do foguete no referencial do ar (alinhada com o eixo do foguete)
    z = V*cos(K)+sqrt((V*cos(K))^2 + VW^2);
endfunction

function z = ang_gama(t,K,V_air,VW) //ângulo entre o vetor velocidade com referencial no solo e o horizonte
    z = asin((sin(K)*V_air)/VW);
endfunction

function z = cg(t,M0_mot,M_s,Xcg0,M_mot)  //centro de gravidade do foguete
    z = ((M0_mot+M_s)*Xcg0)/(M_mot+M_s);
endfunction

function z = massa_motor(t,M0_mot,taxa_mot) //massa do motor em um certo instante de tempo
    z = M0_mot-(taxa_mot*t);
endfunction

function z = massa_foguete(t,M_s,M0_mot,M_mot) //massa do foguete em um certo instante de tempo
    z = M_s+M0_mot-M_mot;
endfunction

function z = empuxo(t,V_eje,taxa,Pe_noz,Patm,A_noz) //empuxo do foguete
    z = V_eje*taxa+(Pe_noz-Patm)*A_noz;
endfunction

function  z = dist_1(t,Xcg,Xcp) //distância entre o centro de gravidade e o centro de pressão do foguete
    z = Xcg-Xcp;
endfunction

function z = altitude(t,V,gama) //altitude do foguete (y)
    z = t*V*sin(gama);
endfunction

function z = distancia(t,V,gama) //distância x percorrida pelo foguete
    z = t*V*cos(gama);
endfunction

function z = densid(t,densid_0,g,y,Temp) //densidade do ar para uma altitude - constante dos gases para pressão em kPa
    z = densid_0^(-y*g*1000/8.314*Temp);
endfunction

function z = pressao_atm(t,densidade,g,y) //pressão atmosférica para determinada altitude 
    z =  101325+densidade*g*y;
endfunction

function [t,V,K,J,X,Y] = rk4DinamicaFoguete(a,b,h,Af0,Al0,Cl,Cd,Inercia,M0_motor,taxa_motor,M_seco,Xcg0,Xcp,Vvento,Temp,Vexaus,Pbocal,Abocal)
    t = a:h:b;
    n = length(t);
    ForcaLift = lift(1,1.2928,Al0,Cl,0);
    ForcaDrag = drag(1,1.2928,Af0,Cd,0);
    MassaMotor = massa_motor(0,M0_motor,taxa_motor);
    MassaFoguete = massa_foguete(1,M_seco,M0_motor,MassaMotor);
    CentroGravidade = cg(1,M0_motor,M_seco,Xcg0,MassaMotor);
    Comprimento1 = dist_1(1,CentroGravidade,Xcp);
    Velocidade_ref_vent = V_ar(1,0,0,0);
    Gama = %pi/2;
    AnguloAtaque = 0;
    Altitude = 0;
    Distancia = 0;
    DensidadeAr = densid(1,1.2928,9.81,Altitude,Temp);
    PressaoAtmosferica = pressao_atm(1,1.2928,9.81,Altitude);
    Empuxo = empuxo(1,Vexaus,taxa_motor,Pbocal,PressaoAtmosferica,Abocal);
    V(1)=0;
    K(1)=0;
    J(1)=0;
    
    for i = 1:n-1
        // gV(t,V,K,J,T,M,g,gama,D)
        // gK(t,V,K,J)
        // gJ(t,V,K,J,T,I,L,VW,gama,Xcg,d1) 
        
        k1V = gV(t(i),V(i),K(i),J(i),Empuxo,MassaFoguete,9.81,Gama,ForcaDrag);
        k1K = gK(t(i),V(i),K(i),J(i));
        k1J = gJ(t(i),V(i),K(i),J(i),Empuxo,Inercia,ForcaLift,Vvento,Gama,CentroGravidade,Comprimento1);
        
        k2V = gV(t(i)+h/2,V(i)+k1V*h/2,K(i)+k1K*h/2,J(i)+k1J*h/2,Empuxo,MassaFoguete,9.81,Gama,ForcaDrag);
        k2K = gK(t(i)+h/2,V(i)+k1V*h/2,K(i)+k1K*h/2,J(i)+k1J*h/2);
        k2J = gJ(t(i)+h/2,V(i)+k1V*h/2,K(i)+k1K*h/2,J(i)+k1J*h/2,Empuxo,Inercia,ForcaLift,Vvento,Gama,CentroGravidade,Comprimento1);
        
        k3V = gV(t(i)+h/2,V(i)+k2V*h/2,K(i)+k2K*h/2,J(i)+k2J*h/2,Empuxo,MassaFoguete,9.81,Gama,ForcaDrag);
        k3K = gK(t(i)+h/2,V(i)+k2V*h/2,K(i)+k2K*h/2,J(i)+k2J*h/2);
        k3J = gJ(t(i)+h/2,V(i)+k2V*h/2,K(i)+k2K*h/2,J(i)+k2J*h/2,Empuxo,Inercia,ForcaLift,Vvento,Gama,CentroGravidade,Comprimento1);
        
        k4V = gV(t(i)+h,V(i)+k3V*h,K(i)+k3K*h,J(i)+k3J*h,Empuxo,MassaFoguete,9.81,Gama,ForcaDrag);
        k4K = gK(t(i)+h,V(i)+k3V*h,K(i)+k3K*h,J(i)+k3J*h);
        k4J = gJ(t(i)+h,V(i)+k3V*h,K(i)+k3K*h,J(i)+k3J*h,Empuxo,Inercia,ForcaLift,Vvento,Gama,CentroGravidade,Comprimento1);
        
        kV = (k1V+2*k2V+2*k3V+k4V)/6;
        kK = (k1K+2*k2K+2*k3K+k4K)/6;
        kJ = (k1J+2*k2J+2*k3J+k4J)/6;
        
        V(i+1) = V(i) + kV*h;
        K(i+1) = K(i) + kK*h;
        J(i+1) = J(i) + kJ*h;
        
        AreaFrontal = ar_frontal(t,Af0,K(i+1));
        ForcaLift = lift(1,DensidadeAr,Al0,Cl,V(i+1));
        ForcaDrag = drag(1,DensidadeAr,AreaFrontal,Cd,V(i+1));
        MassaMotor = massa_motor(t(i+1),M0_motor,taxa_motor);
        MassaFoguete = massa_foguete(1,M_seco,M0_motor,MassaMotor);
        CentroGravidade = cg(1,M0_motor,M_seco,Xcg0,MassaMotor);
        Comprimento1 = dist_1(1,CentroGravidade,Xcp);
        Velocidade_ref_vent = V_ar(1,V(i+1),Vvento,K(i+1));
        Gama = ang_gama(1,K(i+1),Velocidade_ref_vent,Vvento);
        Altitude = altitude(t(i+1),V(i+1),Gama);
        Distancia = distancia(t(i+1),V(i+1),Gama);
        DensidadeAr = densid(1,1.2928,9.81,Altitude,Temp);
        PressaoAtmosferica = pressao_atm(1,1.2928,9.81,Altitude);
        Empuxo = empuxo(1,Vexaus,taxa_motor,Pbocal,PressaoAtmosferica,Abocal);
        
    end
    
endfunction
