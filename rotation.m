function [] = rotation()
	global K
	imLeft = imread('Rot/1.jpg');
	imRight = imread('Rot/2.jpg');

	pointsLeft = [ [701 1925 1]
					[1542 1779 1]
					[1729 2101 1]
					[737 2320 1] ];

	pointsRight = [ [2600 1905 1]
					[3442 1965 1]
					[3560 2340 1]
					[2544 2246 1] ];

	f = 4500;
	px = size(imLeft, 1)/2;
	py = size(imRight, 2)/2;
	K = [f 0 px; 0 f py; 0 0 1];

	H = ComputeHomography(pointsRight, pointsLeft);
	R = GetRotationFromHomography(H);
	
	n = 25;
	RSet = InterpolateRotaton(eye(3), R, n);
	
	for i=1:n
		Ri = MakeRotationOrthonormal(RSet{i});
		
		HLeft = K*Ri*inv(K);
		HRight = HLeft*inv(H);
		
		im1 = ImageWarping(imLeft, HLeft);
		im2 = ImageWarping(imRight, HRight);

		w = (n - i + 1)/n;
		imInterp = w*double(im1) + (1-w)*double(im2);
        
        name = sprintf('output/%d.jpg', i);
        imwrite(uint8(imInterp), name);
    end
end

function R = GetRotationFromHomography(H)
	global K
	R_ = inv(K)*H*K;
	R_ = 1/det(R_)^(1/3) * R_;
	R = MakeRotationOrthonormal(R_);
end

function Rnew = MakeRotationOrthonormal(R)
	[U, ~, V] = svd(R);
	Rnew = U*transpose(V);
	lambda = det(R);
	if lambda < 0
		lambda = -lambda;
	end
	Rnew = 1/lambda^(1/3) * Rnew;
end

function Rset = InterpolateRotaton(fromR, toR, n)
	q1 = Rotation2Quaternion(fromR);
	q2 = Rotation2Quaternion(toR);
	omega = acos(q1'*q2);

	w = 0 : 1/n : 1;
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

	qw= sqrt(1 + m00 + m11 + m22) /2;
	qx = (m21 - m12)/( 4 *qw);
	qy = (m02 - m20)/( 4 *qw);
	qz = (m10 - m01)/( 4 *qw);

	q = [qw; qx; qy; qz];
end
