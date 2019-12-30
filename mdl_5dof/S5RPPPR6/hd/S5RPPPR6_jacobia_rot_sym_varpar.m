% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR6
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RPPPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPPPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->10), mult. (67->21), div. (19->8), fcn. (105->7), ass. (0->18)
	t22 = cos(pkin(7));
	t21 = sin(pkin(7));
	t23 = sin(qJ(1));
	t27 = t23 * t21;
	t14 = atan2(-t27, -t22);
	t12 = sin(t14);
	t13 = cos(t14);
	t11 = -t12 * t27 - t13 * t22;
	t24 = cos(qJ(1));
	t20 = t24 ^ 2;
	t30 = 0.1e1 / t11 ^ 2 * t20;
	t17 = t21 ^ 2;
	t25 = t22 ^ 2;
	t26 = t23 ^ 2;
	t15 = 0.1e1 / (0.1e1 + t26 * t17 / t25);
	t29 = t15 / t22;
	t28 = 0.1e1 / t26 * t20;
	t1 = [t24 * t21 * t29, 0, 0, 0, 0; (-0.1e1 / t11 * t27 - (-t13 * t17 * t23 * t29 + (t15 - 0.1e1) * t21 * t12) * t21 * t30) / (t17 * t30 + 0.1e1), 0, 0, 0, 0; (-0.1e1 - t28) * t22 / (t25 * t28 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (36->14), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t25 = sin(pkin(7));
	t27 = cos(pkin(7));
	t28 = sin(qJ(1));
	t32 = t28 * t27;
	t20 = atan2(-t32, t25);
	t18 = sin(t20);
	t19 = cos(t20);
	t13 = -t18 * t32 + t19 * t25;
	t29 = cos(qJ(1));
	t36 = 0.1e1 / t13 ^ 2 * t29 ^ 2;
	t23 = t27 ^ 2;
	t21 = 0.1e1 / (0.1e1 + t28 ^ 2 * t23 / t25 ^ 2);
	t35 = t21 / t25;
	t24 = sin(pkin(8));
	t34 = t28 * t24;
	t26 = cos(pkin(8));
	t33 = t28 * t26;
	t31 = t29 * t24;
	t30 = t29 * t26;
	t17 = t25 * t31 + t33;
	t16 = -t25 * t30 + t34;
	t15 = 0.1e1 / t17 ^ 2;
	t1 = [-t29 * t27 * t35, 0, 0, 0, 0; (-0.1e1 / t13 * t32 - (t19 * t23 * t28 * t35 + (t21 - 0.1e1) * t27 * t18) * t27 * t36) / (t23 * t36 + 0.1e1), 0, 0, 0, 0; ((t25 * t33 + t31) / t17 - (-t25 * t34 + t30) * t16 * t15) / (t16 ^ 2 * t15 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (77->17), mult. (231->39), div. (30->10), fcn. (341->11), ass. (0->32)
	t47 = cos(pkin(8));
	t48 = cos(pkin(7));
	t60 = t48 * t47;
	t46 = sin(pkin(7));
	t45 = sin(pkin(8));
	t52 = cos(qJ(1));
	t55 = t52 * t45;
	t50 = sin(qJ(1));
	t56 = t50 * t47;
	t40 = t46 * t56 + t55;
	t37 = atan2(t40, t60);
	t34 = sin(t37);
	t35 = cos(t37);
	t29 = t34 * t40 + t35 * t60;
	t54 = t52 * t47;
	t57 = t50 * t45;
	t38 = -t46 * t54 + t57;
	t62 = t38 ^ 2 / t29 ^ 2;
	t61 = 0.1e1 / t60;
	t59 = t48 * t50;
	t58 = t48 * t52;
	t39 = t46 * t55 + t56;
	t49 = sin(qJ(5));
	t51 = cos(qJ(5));
	t33 = t39 * t51 + t49 * t58;
	t31 = 0.1e1 / t33 ^ 2;
	t32 = t39 * t49 - t51 * t58;
	t53 = t32 ^ 2 * t31 + 0.1e1;
	t41 = -t46 * t57 + t54;
	t36 = 0.1e1 / (0.1e1 + t40 ^ 2 / t48 ^ 2 / t47 ^ 2);
	t30 = 0.1e1 / t53;
	t1 = [-t38 * t36 * t61, 0, 0, 0, 0; (t40 / t29 - (-t34 + (-t35 * t40 * t61 + t34) * t36) * t62) / (0.1e1 + t62), 0, 0, 0, 0; ((t41 * t49 + t51 * t59) / t33 - (t41 * t51 - t49 * t59) * t32 * t31) * t30, 0, 0, 0, t53 * t30;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end