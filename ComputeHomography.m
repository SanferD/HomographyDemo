function H = ComputeHomography(pointsTo, pointsFrom)
	eqns = GetMatrixOfEquations(pointsFrom, pointsTo);
	soln = SolveMatrixOfEquations(eqns);
	H = ConvertSolutionVectorToMatrix(soln);
	H = DivideMatrixByLargestSingularValue(H);
end

function eqns = GetMatrixOfEquations(u, v)
	eqns = [];
	for i = 1 : size(u,1)
		nextTwoEqs = GetNextTwoEquations(u(i, 1), u(i, 2), v(i, 1), v(i, 2));
	    eqns = [eqns; nextTwoEqs];
	end
end

function twoEqs = GetNextTwoEquations(ux, uy, vx, vy)
	firstEq = GetFirstEquation(ux, uy, vx, vy);
	secondEq = GetSecondEquation(ux, uy, vx, vy);
	twoEqs = [firstEq; secondEq];
end

function row = GetFirstEquation(ux, uy, vx, vy)
	row = [ux uy 1 0 0 0 -ux*vx -uy*vx -vx];
end

function row = GetSecondEquation(ux, uy, vx, vy)
	row = [0 0 0 ux uy 1 -ux*vy -uy*vy -vy];
end

function soln = SolveMatrixOfEquations(eqns)
	[~, ~, v] = svd(eqns);
	soln = NullSpaceFromSVD(v);
end

function nullspace = NullSpaceFromSVD(lowerTriangle)
	nullspace = lowerTriangle(:, end);
end

function H = ConvertSolutionVectorToMatrix(soln)
	H = [soln(1:3)'; soln(4:6)'; soln(7:9)']; 
end

function normalized = DivideMatrixByLargestSingularValue(M)
	normalized = M/norm(M);
end
