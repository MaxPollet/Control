f = 0.9;


step = round((f-0.5)/0.1);


phi = 0;
if step == 0
    phi_e = 0;

else


    for i = 0:step-1

        phi_temp  = asin(sin(2*pi*(0.5+0.1*(i))*0.4+phi));
        if phi_temp > 0
            if abs(phi_temp - mod(2*pi*(0.5+0.1*(i))*0.4+phi,2*pi)) < 10^(-14)
                phi = phi_temp;
            else
                phi = pi - phi_temp;
            end
        else
             if abs((pi - phi_temp) - mod(2*pi*(0.5+0.1*(i))*0.4+phi,2*pi)) < 10^(-14)
                phi = pi - phi_temp;
            else
                phi = 2*pi - abs(phi_temp);

             end


        end
        
    end

    phi_e = phi;
end




