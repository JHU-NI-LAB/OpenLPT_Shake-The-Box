function I = ParticleProjection(a, b, c, alpha)

   for x = 1:5
       for y = 1:5
            xx = (x-3)*cos(alpha) + (y-3)*sin(alpha);
            yy = -(x-3)*sin(alpha) + (y-3)*cos(alpha);
            I(y,x) = uint8(a*exp(-(b*(xx)^2 + c*(yy)^2))); % not sure if max is the right thing to use for overlapping particles
       end
   end
end

