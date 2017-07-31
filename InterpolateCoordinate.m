function [Rset, Cset] = InterpolateCoordinate(R1, C1, R2, C2, n)

	Cx = linspace(C1(1), C2(1), n+1);
	Cy = linspace(C1(2), C2(2), n+1);
	Cz = linspace(C1(3), C2(3), n+1);

	Cset = [Cx; Cy; Cz];

	w = 0 : 1/n : 1;
	q1 = Rotation2Quaternion(R1);
	q2 = Rotation2Quaternion(R2);

	omega = acos(q1'*q2);
	for i = 1 : length(w)
	    q = sin(omega*(1-w(i)))/sin(omega) * q1 + sin(omega*w(i))/sin(omega) * q2;
	    Rset{i} = Quaternion2Rotation(q);
	end
end

function R = Quaternion2Rotation(q)
	q = q/norm(q);

	qw = q(1);
	qx = q(2);
	qy = q(3);
	qz = q(4);

	R(1,1) = 1 - 2*qy^2 - 2*qz^2;
	R(1,2) = 2*qx*qy - 2*qz*qw;
	R(1,3) = 2*qx*qz + 2*qy*qw;

	R(2,1) = 2*qx*qy + 2*qz*qw;
	R(2,2) = 1 - 2*qx^2 - 2*qz^2;
	R(2,3) = 2*qy*qz - 2*qx*qw;

	R(3,1) = 2*qx*qz - 2*qy*qw;
	R(3,2) = 2*qy*qz + 2*qx*qw;
	R(3,3) = 1 - 2*qx^2 - 2*qy^2;
end

function q = Rotation2Quaternion(R)
	m00 = R(1,1); m01 = R(1,2); m02 = R(1,3);
	m10 = R(2,1); m11 = R(2,2); m12 = R(2,3);
	m20 = R(3,1); m21 = R(3,2); m22 = R(3,3);

	qw= sqrt(1 + m00 + m11 + m22) /2
	qx = (m21 - m12)/( 4 *qw)
	qy = (m02 - m20)/( 4 *qw)
	qz = (m10 - m01)/( 4 *qw)

	q = [qw; qx; qy; qz];
end
