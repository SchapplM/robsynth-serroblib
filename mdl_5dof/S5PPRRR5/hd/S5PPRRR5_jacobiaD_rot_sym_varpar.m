% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PPRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:26:01
	% EndTime: 2019-12-29 15:26:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:26:06
	% EndTime: 2019-12-29 15:26:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:26:00
	% EndTime: 2019-12-29 15:26:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:26:00
	% EndTime: 2019-12-29 15:26:01
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (46->7), mult. (124->20), div. (18->4), fcn. (140->4), ass. (0->15)
	t33 = sin(pkin(8));
	t34 = cos(pkin(8));
	t35 = sin(qJ(3));
	t36 = cos(qJ(3));
	t31 = t33 * t35 + t34 * t36;
	t28 = 0.1e1 / t31 ^ 2;
	t43 = qJD(3) * t28;
	t39 = t33 * t36 - t34 * t35;
	t27 = t39 ^ 2;
	t24 = t27 * t28 + 0.1e1;
	t40 = t31 * t43;
	t41 = t39 / t31 * t43;
	t42 = (-t27 * t41 - t39 * t40) / t24 ^ 2;
	t22 = 0.1e1 / t24;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t42 - 0.2e1 * (t22 * t40 - (-t22 * t41 - t28 * t42) * t39) * t39, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:25:53
	% EndTime: 2019-12-29 15:25:53
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (260->9), mult. (248->20), div. (36->4), fcn. (280->4), ass. (0->17)
	t50 = qJ(3) + qJ(4);
	t47 = sin(t50);
	t48 = cos(t50);
	t51 = sin(pkin(8));
	t52 = cos(pkin(8));
	t45 = t51 * t47 + t52 * t48;
	t42 = 0.1e1 / t45 ^ 2;
	t58 = t42 * (qJD(3) + qJD(4));
	t44 = t52 * t47 - t51 * t48;
	t41 = t44 ^ 2;
	t38 = t41 * t42 + 0.1e1;
	t55 = t45 * t58;
	t56 = t44 / t45 * t58;
	t57 = (t41 * t56 + t44 * t55) / t38 ^ 2;
	t36 = 0.1e1 / t38;
	t34 = -0.2e1 * t57 + 0.2e1 * (t36 * t55 + (t36 * t56 - t42 * t57) * t44) * t44;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, t34, t34, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:26:06
	% EndTime: 2019-12-29 15:26:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end