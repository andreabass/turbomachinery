function Y = y_AM_Re(Y_in,Re_out)
        
        Re_ref = 2e5;
        Y = Y_in * (Re_ref / Re_out)^0.2;
        
end

