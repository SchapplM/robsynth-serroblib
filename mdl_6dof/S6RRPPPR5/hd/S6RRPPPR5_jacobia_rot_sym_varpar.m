% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR5
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
%   Wie in S6RRPPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (99->20), mult. (197->54), div. (47->9), fcn. (297->9), ass. (0->32)
	t34 = cos(qJ(2));
	t32 = sin(qJ(2));
	t33 = sin(qJ(1));
	t39 = t33 * t32;
	t25 = atan2(-t39, -t34);
	t23 = sin(t25);
	t24 = cos(t25);
	t16 = -t23 * t39 - t24 * t34;
	t15 = 0.1e1 / t16 ^ 2;
	t35 = cos(qJ(1));
	t44 = t15 * t35 ^ 2;
	t30 = sin(pkin(9));
	t31 = cos(pkin(9));
	t36 = t35 * t31;
	t22 = t33 * t30 + t34 * t36;
	t20 = 0.1e1 / t22 ^ 2;
	t37 = t35 * t30;
	t21 = -t33 * t31 + t34 * t37;
	t43 = t20 * t21;
	t27 = t32 ^ 2;
	t42 = t27 / t34 ^ 2;
	t41 = t32 * t35;
	t26 = 0.1e1 / (t33 ^ 2 * t42 + 0.1e1);
	t40 = t33 * t26;
	t38 = t33 * t34;
	t28 = 0.1e1 / t34;
	t19 = 0.1e1 / t22;
	t18 = (0.1e1 + t42) * t40;
	t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
	t14 = 0.1e1 / t16;
	t13 = 0.1e1 / (t27 * t44 + 0.1e1);
	t1 = [t28 * t26 * t41, t18, 0, 0, 0, 0; (-t14 * t39 - (-t24 * t27 * t28 * t40 + (t26 - 0.1e1) * t32 * t23) * t32 * t44) * t13, (t34 * t14 - (-t23 * t38 + t24 * t32 + (t23 * t34 - t24 * t39) * t18) * t32 * t15) * t35 * t13, 0, 0, 0, 0; ((-t30 * t38 - t36) * t19 - (-t31 * t38 + t37) * t43) * t17, (-t19 * t30 + t31 * t43) * t17 * t41, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (100->20), mult. (313->54), div. (67->11), fcn. (484->9), ass. (0->34)
	t41 = sin(pkin(9));
	t42 = cos(pkin(9));
	t46 = cos(qJ(1));
	t47 = t46 * t42;
	t44 = sin(qJ(1));
	t45 = cos(qJ(2));
	t49 = t44 * t45;
	t31 = t41 * t49 + t47;
	t43 = sin(qJ(2));
	t50 = t43 * t41;
	t29 = atan2(-t31, t50);
	t26 = sin(t29);
	t27 = cos(t29);
	t25 = -t26 * t31 + t27 * t50;
	t24 = 0.1e1 / t25 ^ 2;
	t48 = t46 * t41;
	t33 = -t44 * t42 + t45 * t48;
	t56 = t24 * t33;
	t54 = t27 * t31;
	t53 = t33 ^ 2 * t24;
	t36 = 0.1e1 / t41;
	t37 = 0.1e1 / t43;
	t52 = t36 * t37;
	t38 = 0.1e1 / t43 ^ 2;
	t51 = t38 * t45;
	t40 = 0.1e1 / t46 ^ 2;
	t39 = 0.1e1 / t46;
	t34 = t44 * t41 + t45 * t47;
	t30 = 0.1e1 / (t34 ^ 2 * t40 * t38 + 0.1e1);
	t28 = 0.1e1 / (0.1e1 + t31 ^ 2 * t38 / t41 ^ 2);
	t23 = 0.1e1 / t25;
	t22 = (t31 * t36 * t51 + t44) * t28;
	t21 = 0.1e1 / (0.1e1 + t53);
	t1 = [-t33 * t28 * t52, t22, 0, 0, 0, 0; (-t31 * t23 - (-t26 + (t52 * t54 + t26) * t28) * t53) * t21, (t22 * t54 * t56 + (-t46 * t43 * t23 - (t27 * t45 + (-t22 + t44) * t26 * t43) * t56) * t41) * t21, 0, 0, 0, 0; ((-t42 * t49 + t48) * t39 + t44 * t34 * t40) * t37 * t30, (-t34 * t39 * t51 - t42) * t30, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (102->20), mult. (329->52), div. (61->11), fcn. (494->9), ass. (0->34)
	t46 = sin(qJ(2));
	t60 = t46 ^ 2;
	t45 = cos(pkin(9));
	t44 = sin(pkin(9));
	t49 = cos(qJ(1));
	t52 = t49 * t44;
	t47 = sin(qJ(1));
	t48 = cos(qJ(2));
	t53 = t47 * t48;
	t34 = t45 * t53 - t52;
	t54 = t46 * t45;
	t30 = atan2(-t34, t54);
	t27 = sin(t30);
	t28 = cos(t30);
	t26 = -t27 * t34 + t28 * t54;
	t25 = 0.1e1 / t26 ^ 2;
	t51 = t49 * t45;
	t37 = t47 * t44 + t48 * t51;
	t59 = t25 * t37;
	t57 = t28 * t34;
	t56 = t37 ^ 2 * t25;
	t39 = 0.1e1 / t45;
	t55 = t39 / t46;
	t36 = -t47 * t45 + t48 * t52;
	t33 = 0.1e1 / t36 ^ 2;
	t50 = t49 ^ 2 * t60 * t33;
	t42 = 0.1e1 / t60;
	t32 = 0.1e1 / t36;
	t31 = 0.1e1 / (0.1e1 + t50);
	t29 = 0.1e1 / (0.1e1 + t34 ^ 2 * t42 / t45 ^ 2);
	t24 = 0.1e1 / t26;
	t23 = (t34 * t39 * t42 * t48 + t47) * t29;
	t22 = 0.1e1 / (0.1e1 + t56);
	t1 = [-t37 * t29 * t55, t23, 0, 0, 0, 0; (-t34 * t24 - (-t27 + (t55 * t57 + t27) * t29) * t56) * t22, (t23 * t57 * t59 + (-t49 * t46 * t24 - (t28 * t48 + (-t23 + t47) * t27 * t46) * t59) * t45) * t22, 0, 0, 0, 0; (-t47 * t32 - (-t44 * t53 - t51) * t49 * t33) * t46 * t31, (t32 * t48 * t49 + t44 * t50) * t31, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (157->26), mult. (347->67), div. (52->9), fcn. (503->11), ass. (0->39)
	t50 = cos(pkin(9));
	t53 = sin(qJ(1));
	t55 = cos(qJ(2));
	t49 = sin(pkin(9));
	t56 = cos(qJ(1));
	t59 = t56 * t49;
	t40 = -t53 * t50 + t55 * t59;
	t58 = t56 * t50;
	t41 = t53 * t49 + t55 * t58;
	t51 = sin(qJ(6));
	t54 = cos(qJ(6));
	t33 = t40 * t54 + t41 * t51;
	t31 = 0.1e1 / t33 ^ 2;
	t32 = t40 * t51 - t41 * t54;
	t66 = t31 * t32;
	t52 = sin(qJ(2));
	t61 = t53 * t52;
	t44 = atan2(-t61, -t55);
	t42 = sin(t44);
	t43 = cos(t44);
	t36 = -t42 * t61 - t43 * t55;
	t35 = 0.1e1 / t36 ^ 2;
	t65 = t35 * t56 ^ 2;
	t46 = t52 ^ 2;
	t64 = t46 / t55 ^ 2;
	t63 = t52 * t56;
	t45 = 0.1e1 / (t53 ^ 2 * t64 + 0.1e1);
	t62 = t53 * t45;
	t60 = t53 * t55;
	t57 = t32 ^ 2 * t31 + 0.1e1;
	t47 = 0.1e1 / t55;
	t39 = -t50 * t60 + t59;
	t38 = -t49 * t60 - t58;
	t37 = (0.1e1 + t64) * t62;
	t34 = 0.1e1 / t36;
	t30 = 0.1e1 / t33;
	t29 = 0.1e1 / (t46 * t65 + 0.1e1);
	t28 = 0.1e1 / t57;
	t1 = [t47 * t45 * t63, t37, 0, 0, 0, 0; (-t34 * t61 - (-t43 * t46 * t47 * t62 + (t45 - 0.1e1) * t52 * t42) * t52 * t65) * t29, (t55 * t34 - (-t42 * t60 + t43 * t52 + (t42 * t55 - t43 * t61) * t37) * t52 * t35) * t56 * t29, 0, 0, 0, 0; ((t38 * t51 - t39 * t54) * t30 - (t38 * t54 + t39 * t51) * t66) * t28, ((-t49 * t51 + t50 * t54) * t30 - (-t49 * t54 - t50 * t51) * t66) * t28 * t63, 0, 0, 0, t57 * t28;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end