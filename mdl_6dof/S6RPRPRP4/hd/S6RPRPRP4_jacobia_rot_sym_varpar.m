% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP4
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
%   Wie in S6RPRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (153->17), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
	t26 = qJ(1) + pkin(9);
	t25 = cos(t26);
	t23 = t25 ^ 2;
	t31 = cos(qJ(3));
	t24 = sin(t26);
	t30 = sin(qJ(3));
	t36 = t24 * t30;
	t20 = atan2(-t36, -t31);
	t17 = sin(t20);
	t18 = cos(t20);
	t15 = -t17 * t36 - t18 * t31;
	t14 = 0.1e1 / t15 ^ 2;
	t40 = t14 * t30;
	t39 = t17 * t31;
	t32 = t24 ^ 2;
	t38 = 0.1e1 / t32 * t23;
	t27 = t30 ^ 2;
	t33 = t31 ^ 2;
	t34 = t27 / t33;
	t21 = 0.1e1 / (t32 * t34 + 0.1e1);
	t37 = t24 * t21;
	t35 = t25 * t30;
	t28 = 0.1e1 / t31;
	t19 = 0.1e1 / (t33 * t38 + 0.1e1);
	t16 = (0.1e1 + t34) * t37;
	t13 = 0.1e1 / t15;
	t12 = 0.1e1 / (t23 * t27 * t14 + 0.1e1);
	t1 = [t28 * t21 * t35, 0, t16, 0, 0, 0; (-t13 * t36 - (-t18 * t27 * t28 * t37 + (t21 - 0.1e1) * t30 * t17) * t23 * t40) * t12, 0, (t31 * t13 - (-t24 * t39 + t18 * t30 + (-t18 * t36 + t39) * t16) * t40) * t25 * t12, 0, 0, 0; (-0.1e1 - t38) * t31 * t19, 0, -0.1e1 / t24 * t19 * t35, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (196->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t37 = qJ(1) + pkin(9);
	t36 = cos(t37);
	t55 = t36 ^ 2;
	t42 = sin(qJ(3));
	t35 = sin(t37);
	t44 = cos(qJ(3));
	t50 = t35 * t44;
	t33 = atan2(-t50, t42);
	t31 = sin(t33);
	t32 = cos(t33);
	t26 = -t31 * t50 + t32 * t42;
	t23 = 0.1e1 / t26 ^ 2;
	t54 = t23 * t44;
	t43 = cos(qJ(5));
	t41 = sin(qJ(5));
	t47 = t41 * t42;
	t30 = t35 * t43 + t36 * t47;
	t28 = 0.1e1 / t30 ^ 2;
	t46 = t42 * t43;
	t29 = t35 * t41 - t36 * t46;
	t53 = t28 * t29;
	t52 = t31 * t42;
	t40 = t44 ^ 2;
	t48 = 0.1e1 / t42 ^ 2 * t40;
	t34 = 0.1e1 / (t35 ^ 2 * t48 + 0.1e1);
	t51 = t35 * t34;
	t49 = t36 * t44;
	t45 = t29 ^ 2 * t28 + 0.1e1;
	t38 = 0.1e1 / t42;
	t27 = 0.1e1 / t30;
	t25 = (0.1e1 + t48) * t51;
	t24 = 0.1e1 / t45;
	t22 = 0.1e1 / t26;
	t21 = 0.1e1 / (t55 * t40 * t23 + 0.1e1);
	t1 = [-t38 * t34 * t49, 0, t25, 0, 0, 0; (-t22 * t50 - (t32 * t38 * t40 * t51 + (t34 - 0.1e1) * t44 * t31) * t55 * t54) * t21, 0, (-t42 * t22 - (t35 * t52 + t32 * t44 + (-t32 * t50 - t52) * t25) * t54) * t36 * t21, 0, 0, 0; ((t35 * t46 + t36 * t41) * t27 - (-t35 * t47 + t36 * t43) * t53) * t24, 0, (-t27 * t43 - t41 * t53) * t24 * t49, 0, t45 * t24, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:15
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (378->26), mult. (509->71), div. (100->11), fcn. (770->9), ass. (0->39)
	t55 = cos(qJ(3));
	t66 = t55 ^ 2;
	t46 = qJ(1) + pkin(9);
	t44 = sin(t46);
	t45 = cos(t46);
	t52 = sin(qJ(5));
	t53 = sin(qJ(3));
	t54 = cos(qJ(5));
	t58 = t53 * t54;
	t40 = t44 * t58 + t45 * t52;
	t57 = t55 * t54;
	t35 = atan2(t40, t57);
	t32 = sin(t35);
	t33 = cos(t35);
	t30 = t32 * t40 + t33 * t57;
	t29 = 0.1e1 / t30 ^ 2;
	t38 = t44 * t52 - t45 * t58;
	t65 = t29 * t38;
	t63 = t33 * t40;
	t62 = t38 ^ 2 * t29;
	t61 = t45 * t55;
	t47 = 0.1e1 / t54;
	t50 = 0.1e1 / t55;
	t60 = t47 * t50;
	t59 = t52 * t53;
	t39 = t44 * t54 + t45 * t59;
	t37 = 0.1e1 / t39 ^ 2;
	t56 = t45 ^ 2 * t66 * t37;
	t51 = 0.1e1 / t66;
	t48 = 0.1e1 / t54 ^ 2;
	t41 = -t44 * t59 + t45 * t54;
	t36 = 0.1e1 / t39;
	t34 = 0.1e1 / (t40 ^ 2 * t51 * t48 + 0.1e1);
	t31 = 0.1e1 / (0.1e1 + t56);
	t28 = 0.1e1 / t30;
	t27 = (t40 * t47 * t51 * t53 + t44) * t34;
	t26 = 0.1e1 / (0.1e1 + t62);
	t25 = (t40 * t48 * t52 + t41 * t47) * t50 * t34;
	t1 = [-t38 * t34 * t60, 0, t27, 0, t25, 0; (t40 * t28 - (-t32 + (-t60 * t63 + t32) * t34) * t62) * t26, 0, (-t27 * t63 * t65 + (-t28 * t61 - (-t33 * t53 + (-t27 + t44) * t32 * t55) * t65) * t54) * t26, 0, (t39 * t28 - (-t33 * t55 * t52 + t32 * t41 + (-t32 * t57 + t63) * t25) * t65) * t26, 0; (t37 * t41 * t45 + t36 * t44) * t55 * t31, 0, (t36 * t45 * t53 + t52 * t56) * t31, 0, -t38 * t37 * t31 * t61, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end