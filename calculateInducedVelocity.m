function [coeffX, coeffY, coeffZ] = calculateInducedVelocity(xm, ym, zm, xA, yA, zA, xB, yB, zB)

    % Bound Vortex
    r0 = [xB - xA; yB - yA; zB - zA];
    r1 = [xm - xA; ym - yA; zm - zA];
    r2 = [xm - xB; ym - yB; zm - zB];

    rCross = cross(r1, r2);
    fac1 = rCross/(norm(rCross).^2);
    fac2 = dot(r0, (r1/norm(r1)) - (r2/norm(r2)));

    vBoundCoeff = (fac1 .* fac2) / (4*pi);

    % Trailing Vortex (A-side)
    v1 = [0; zm - zA; yA - ym];
    fac1 = v1/(norm(v1).^2);
    fac2 = 1.0 + (xm - xA)/norm(r1);

    vTrailACoeff = (fac1 .* fac2) / (4*pi);

    % Trailing Vortex (B-side)
    v2 = [0; zm - zB; yB - ym];
    fac1 = v2/(norm(v2).^2);
    fac2 = 1.0 + (xm - xB)/norm(r2);

    vTrailBCoeff = -(fac1 .* fac2) / (4*pi);
               
    % Sum the contributions
    vCoeff = vBoundCoeff + vTrailACoeff + vTrailBCoeff;

    coeffX = vCoeff(1);
    coeffY = vCoeff(2);
    coeffZ = vCoeff(3);
end