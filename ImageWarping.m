function im_warped = ImageWarping(im, H)
	[u_x, u_y] = GetPointsToTransform(size(im,2), size(im,1));
	[v_x, v_y] = TransformPointsUsingHomography(inv(H), u_x, u_y);
	im_warped = BuildWarpedImage(double(im), v_x, v_y);
	im_warped = uint8(im_warped);
end

function [u_x, u_y] = GetPointsToTransform(width, height)
	[u_x, u_y] = meshgrid(1:width, 1:height);
end

function [v_x, v_y] = TransformPointsUsingHomography(H, u_x, u_y)
	v_x = H(1,1)*u_x + H(1,2)*u_y + H(1,3);
	v_y = H(2,1)*u_x + H(2,2)*u_y + H(2,3);
	v_z = H(3,1)*u_x + H(3,2)*u_y + H(3,3);

	v_x = v_x./v_z;
	v_y = v_y./v_z;
end

function im_warped = BuildWarpedImage(im, v_x, v_y)
	h = size(v_x, 1); w = size(v_x,2);
	im_warped(:,:,1) = reshape(interp2(im(:,:,1), v_x(:), v_y(:)), [h, w]);
	im_warped(:,:,2) = reshape(interp2(im(:,:,2), v_x(:), v_y(:)), [h, w]);
	im_warped(:,:,3) = reshape(interp2(im(:,:,3), v_x(:), v_y(:)), [h, w]);
end

