function powerQual = power_quality(powerRho, powerSm)

powerSmTemp = powerSm;
powerSmTemp( powerSm < 0 ) = 0;
powerSmTemp( powerRho < 0 ) = 0;
%Power ranking objective function:
powerQual = powerSmTemp.*powerRho;