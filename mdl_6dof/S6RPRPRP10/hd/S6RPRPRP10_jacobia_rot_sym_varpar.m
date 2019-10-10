% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPRPRP10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (58->14), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->27)
	t27 = sin(qJ(1));
	t23 = t27 ^ 2;
	t26 = sin(qJ(3));
	t28 = cos(qJ(3));
	t29 = cos(qJ(1));
	t32 = t29 * t28;
	t18 = atan2(-t32, t26);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t32 + t17 * t26;
	t13 = 0.1e1 / t14 ^ 2;
	t38 = t13 * t28;
	t37 = t16 * t26;
	t24 = t28 ^ 2;
	t30 = t26 ^ 2;
	t36 = 0.1e1 / t30 * t24;
	t31 = t29 ^ 2;
	t35 = t23 / t31;
	t34 = t27 * t28;
	t20 = 0.1e1 / (t31 * t36 + 0.1e1);
	t33 = t29 * t20;
	t21 = 0.1e1 / t26;
	t19 = 0.1e1 / (t30 * t35 + 0.1e1);
	t15 = (0.1e1 + t36) * t33;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t23 * t24 * t13 + 0.1e1);
	t1 = [t21 * t20 * t34, 0, t15, 0, 0, 0; (-t12 * t32 + (-t17 * t21 * t24 * t33 + (-t20 + 0.1e1) * t28 * t16) * t23 * t38) * t11, 0, (t26 * t12 + (t29 * t37 + t17 * t28 + (-t17 * t32 - t37) * t15) * t38) * t27 * t11, 0, 0, 0; (0.1e1 + t35) * t26 * t19, 0, 0.1e1 / t29 * t19 * t34, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (64->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
	t40 = cos(qJ(3));
	t37 = sin(qJ(3));
	t41 = cos(qJ(1));
	t44 = t41 * t37;
	t32 = atan2(t44, t40);
	t29 = sin(t32);
	t30 = cos(t32);
	t22 = t29 * t44 + t30 * t40;
	t21 = 0.1e1 / t22 ^ 2;
	t38 = sin(qJ(1));
	t52 = t21 * t38 ^ 2;
	t36 = sin(qJ(5));
	t39 = cos(qJ(5));
	t43 = t41 * t39;
	t47 = t38 * t40;
	t28 = -t36 * t47 + t43;
	t26 = 0.1e1 / t28 ^ 2;
	t45 = t41 * t36;
	t27 = t39 * t47 + t45;
	t51 = t26 * t27;
	t50 = t29 * t40;
	t33 = t37 ^ 2;
	t49 = t33 / t40 ^ 2;
	t48 = t37 * t38;
	t31 = 0.1e1 / (t41 ^ 2 * t49 + 0.1e1);
	t46 = t41 * t31;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t34 = 0.1e1 / t40;
	t25 = 0.1e1 / t28;
	t24 = (0.1e1 + t49) * t46;
	t23 = 0.1e1 / t42;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t52 + 0.1e1);
	t1 = [-t34 * t31 * t48, 0, t24, 0, 0, 0; (t20 * t44 - (-t30 * t33 * t34 * t46 + (t31 - 0.1e1) * t37 * t29) * t37 * t52) * t19, 0, (t40 * t20 - (t41 * t50 - t30 * t37 + (t30 * t44 - t50) * t24) * t37 * t21) * t38 * t19, 0, 0, 0; ((-t38 * t36 + t40 * t43) * t25 - (-t38 * t39 - t40 * t45) * t51) * t23, 0, (-t25 * t39 - t36 * t51) * t23 * t48, 0, t42 * t23, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (195->26), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->39)
	t50 = sin(qJ(3));
	t66 = t50 ^ 2;
	t49 = sin(qJ(5));
	t51 = sin(qJ(1));
	t53 = cos(qJ(3));
	t52 = cos(qJ(5));
	t54 = cos(qJ(1));
	t56 = t54 * t52;
	t37 = t51 * t49 - t53 * t56;
	t59 = t50 * t52;
	t32 = atan2(-t37, -t59);
	t30 = sin(t32);
	t31 = cos(t32);
	t29 = -t30 * t37 - t31 * t59;
	t27 = 0.1e1 / t29 ^ 2;
	t57 = t54 * t49;
	t58 = t51 * t53;
	t40 = t52 * t58 + t57;
	t65 = t27 * t40;
	t63 = t31 * t37;
	t62 = t40 ^ 2 * t27;
	t44 = 0.1e1 / t50;
	t47 = 0.1e1 / t52;
	t61 = t44 * t47;
	t60 = t50 * t51;
	t41 = -t49 * t58 + t56;
	t36 = 0.1e1 / t41 ^ 2;
	t55 = t51 ^ 2 * t66 * t36;
	t48 = 0.1e1 / t52 ^ 2;
	t45 = 0.1e1 / t66;
	t39 = t51 * t52 + t53 * t57;
	t35 = 0.1e1 / t41;
	t34 = 0.1e1 / (t37 ^ 2 * t45 * t48 + 0.1e1);
	t33 = 0.1e1 / (0.1e1 + t55);
	t28 = (-t37 * t45 * t47 * t53 + t54) * t34;
	t26 = 0.1e1 / t29;
	t25 = 0.1e1 / (0.1e1 + t62);
	t24 = (t37 * t48 * t49 + t39 * t47) * t44 * t34;
	t1 = [t40 * t34 * t61, 0, t28, 0, t24, 0; (-t37 * t26 - (-t30 + (-t61 * t63 + t30) * t34) * t62) * t25, 0, (t28 * t63 * t65 + (-t26 * t60 - (-t31 * t53 + (t28 - t54) * t50 * t30) * t65) * t52) * t25, 0, (t41 * t26 - (t31 * t50 * t49 - t30 * t39 + (t30 * t59 - t63) * t24) * t65) * t25, 0; (-t36 * t39 * t51 - t35 * t54) * t50 * t33, 0, (-t35 * t58 + t49 * t55) * t33, 0, -t40 * t36 * t33 * t60, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end