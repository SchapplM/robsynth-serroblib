% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR4
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
%   Wie in S6RRPPPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:23
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->16), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->27)
	t29 = cos(qJ(1));
	t25 = t29 ^ 2;
	t28 = cos(qJ(2));
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t32 = t27 * t26;
	t18 = atan2(-t32, -t28);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t32 - t17 * t28;
	t13 = 0.1e1 / t14 ^ 2;
	t38 = t13 * t26;
	t37 = t16 * t28;
	t21 = t26 ^ 2;
	t31 = t28 ^ 2;
	t36 = t21 / t31;
	t30 = t27 ^ 2;
	t35 = 0.1e1 / t30 * t25;
	t34 = t26 * t29;
	t19 = 0.1e1 / (t30 * t36 + 0.1e1);
	t33 = t27 * t19;
	t23 = 0.1e1 / t28;
	t20 = 0.1e1 / (t31 * t35 + 0.1e1);
	t15 = (0.1e1 + t36) * t33;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t25 * t21 * t13 + 0.1e1);
	t1 = [t23 * t19 * t34, t15, 0, 0, 0, 0; (-t12 * t32 - (-t17 * t21 * t23 * t33 + (t19 - 0.1e1) * t26 * t16) * t25 * t38) * t11, (t28 * t12 - (-t27 * t37 + t17 * t26 + (-t17 * t32 + t37) * t15) * t38) * t29 * t11, 0, 0, 0, 0; (-0.1e1 - t35) * t28 * t20, -0.1e1 / t27 * t20 * t34, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (77->20), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->34)
	t32 = sin(qJ(2));
	t33 = sin(qJ(1));
	t34 = cos(qJ(2));
	t39 = t33 * t34;
	t25 = atan2(-t39, t32);
	t23 = sin(t25);
	t24 = cos(t25);
	t16 = -t23 * t39 + t24 * t32;
	t15 = 0.1e1 / t16 ^ 2;
	t35 = cos(qJ(1));
	t46 = t15 * t35 ^ 2;
	t30 = sin(pkin(9));
	t37 = t35 * t30;
	t31 = cos(pkin(9));
	t40 = t33 * t31;
	t22 = t32 * t37 + t40;
	t20 = 0.1e1 / t22 ^ 2;
	t36 = t35 * t31;
	t41 = t33 * t30;
	t21 = -t32 * t36 + t41;
	t45 = t20 * t21;
	t44 = t23 * t32;
	t29 = t34 ^ 2;
	t43 = 0.1e1 / t32 ^ 2 * t29;
	t26 = 0.1e1 / (t33 ^ 2 * t43 + 0.1e1);
	t42 = t33 * t26;
	t38 = t34 * t35;
	t27 = 0.1e1 / t32;
	t19 = 0.1e1 / t22;
	t18 = (0.1e1 + t43) * t42;
	t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
	t14 = 0.1e1 / t16;
	t13 = 0.1e1 / (t29 * t46 + 0.1e1);
	t1 = [-t27 * t26 * t38, t18, 0, 0, 0, 0; (-t14 * t39 - (t24 * t27 * t29 * t42 + (t26 - 0.1e1) * t34 * t23) * t34 * t46) * t13, (-t32 * t14 - (t33 * t44 + t24 * t34 + (-t24 * t39 - t44) * t18) * t34 * t15) * t35 * t13, 0, 0, 0, 0; ((t32 * t40 + t37) * t19 - (-t32 * t41 + t36) * t45) * t17, (-t19 * t31 - t30 * t45) * t17 * t38, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (102->19), mult. (329->51), div. (61->11), fcn. (494->9), ass. (0->35)
	t47 = cos(qJ(2));
	t60 = t47 ^ 2;
	t45 = sin(qJ(2));
	t43 = sin(pkin(9));
	t48 = cos(qJ(1));
	t51 = t48 * t43;
	t44 = cos(pkin(9));
	t46 = sin(qJ(1));
	t53 = t46 * t44;
	t36 = t45 * t53 + t51;
	t52 = t47 * t44;
	t30 = atan2(t36, t52);
	t27 = sin(t30);
	t28 = cos(t30);
	t26 = t27 * t36 + t28 * t52;
	t25 = 0.1e1 / t26 ^ 2;
	t50 = t48 * t44;
	t54 = t46 * t43;
	t34 = -t45 * t50 + t54;
	t59 = t25 * t34;
	t57 = t28 * t36;
	t56 = t34 ^ 2 * t25;
	t38 = 0.1e1 / t44;
	t55 = t38 / t47;
	t35 = t45 * t51 + t53;
	t33 = 0.1e1 / t35 ^ 2;
	t49 = t48 ^ 2 * t60 * t33;
	t41 = 0.1e1 / t60;
	t32 = 0.1e1 / t35;
	t31 = 0.1e1 / (0.1e1 + t49);
	t29 = 0.1e1 / (0.1e1 + t36 ^ 2 * t41 / t44 ^ 2);
	t24 = 0.1e1 / t26;
	t23 = (t36 * t38 * t41 * t45 + t46) * t29;
	t22 = 0.1e1 / (0.1e1 + t56);
	t1 = [-t34 * t29 * t55, t23, 0, 0, 0, 0; (t36 * t24 - (-t27 + (-t55 * t57 + t27) * t29) * t56) * t22, (-t23 * t57 * t59 + (-t48 * t47 * t24 - (-t28 * t45 + (-t23 + t46) * t27 * t47) * t59) * t44) * t22, 0, 0, 0, 0; (t46 * t32 + (-t45 * t54 + t50) * t48 * t33) * t47 * t31, (t32 * t45 * t48 + t43 * t49) * t31, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (135->26), mult. (347->66), div. (52->9), fcn. (503->11), ass. (0->41)
	t52 = sin(qJ(2));
	t50 = cos(pkin(9));
	t56 = cos(qJ(1));
	t58 = t56 * t50;
	t49 = sin(pkin(9));
	t53 = sin(qJ(1));
	t63 = t53 * t49;
	t38 = -t52 * t58 + t63;
	t59 = t56 * t49;
	t62 = t53 * t50;
	t39 = t52 * t59 + t62;
	t51 = sin(qJ(6));
	t54 = cos(qJ(6));
	t33 = t38 * t51 + t39 * t54;
	t31 = 0.1e1 / t33 ^ 2;
	t32 = -t38 * t54 + t39 * t51;
	t68 = t31 * t32;
	t55 = cos(qJ(2));
	t61 = t53 * t55;
	t44 = atan2(t61, -t52);
	t42 = sin(t44);
	t43 = cos(t44);
	t36 = t42 * t61 - t43 * t52;
	t35 = 0.1e1 / t36 ^ 2;
	t67 = t35 * t56 ^ 2;
	t66 = t42 * t52;
	t48 = t55 ^ 2;
	t65 = 0.1e1 / t52 ^ 2 * t48;
	t45 = 0.1e1 / (t53 ^ 2 * t65 + 0.1e1);
	t64 = t53 * t45;
	t60 = t55 * t56;
	t57 = t32 ^ 2 * t31 + 0.1e1;
	t46 = 0.1e1 / t52;
	t41 = -t52 * t63 + t58;
	t40 = t52 * t62 + t59;
	t37 = (0.1e1 + t65) * t64;
	t34 = 0.1e1 / t36;
	t30 = 0.1e1 / t33;
	t29 = 0.1e1 / (t48 * t67 + 0.1e1);
	t28 = 0.1e1 / t57;
	t1 = [-t46 * t45 * t60, t37, 0, 0, 0, 0; (t34 * t61 + (-t43 * t46 * t48 * t64 + (-t45 + 0.1e1) * t55 * t42) * t55 * t67) * t29, (t52 * t34 + (-t53 * t66 - t43 * t55 + (t43 * t61 + t66) * t37) * t55 * t35) * t56 * t29, 0, 0, 0, 0; ((-t40 * t54 + t41 * t51) * t30 - (t40 * t51 + t41 * t54) * t68) * t28, ((t49 * t51 + t50 * t54) * t30 - (t49 * t54 - t50 * t51) * t68) * t28 * t60, 0, 0, 0, t57 * t28;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end