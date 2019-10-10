% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP5
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
%   Wie in S6RRPPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:34
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:34:00
	% EndTime: 2019-10-10 09:34:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:34:00
	% EndTime: 2019-10-10 09:34:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:34:00
	% EndTime: 2019-10-10 09:34:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:34:00
	% EndTime: 2019-10-10 09:34:01
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
	% StartTime: 2019-10-10 09:34:01
	% EndTime: 2019-10-10 09:34:01
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
	% StartTime: 2019-10-10 09:34:01
	% EndTime: 2019-10-10 09:34:01
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (135->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t45 = sin(qJ(2));
	t46 = sin(qJ(1));
	t47 = cos(qJ(2));
	t53 = t46 * t47;
	t37 = atan2(-t53, t45);
	t35 = sin(t37);
	t36 = cos(t37);
	t29 = -t35 * t53 + t36 * t45;
	t28 = 0.1e1 / t29 ^ 2;
	t48 = cos(qJ(1));
	t60 = t28 * t48 ^ 2;
	t41 = pkin(9) + qJ(5);
	t39 = sin(t41);
	t51 = t48 * t39;
	t40 = cos(t41);
	t54 = t46 * t40;
	t34 = t45 * t51 + t54;
	t32 = 0.1e1 / t34 ^ 2;
	t50 = t48 * t40;
	t55 = t46 * t39;
	t33 = -t45 * t50 + t55;
	t59 = t32 * t33;
	t58 = t35 * t45;
	t44 = t47 ^ 2;
	t57 = 0.1e1 / t45 ^ 2 * t44;
	t38 = 0.1e1 / (t46 ^ 2 * t57 + 0.1e1);
	t56 = t46 * t38;
	t52 = t47 * t48;
	t49 = t33 ^ 2 * t32 + 0.1e1;
	t42 = 0.1e1 / t45;
	t31 = 0.1e1 / t34;
	t30 = (0.1e1 + t57) * t56;
	t27 = 0.1e1 / t29;
	t26 = 0.1e1 / t49;
	t25 = 0.1e1 / (t44 * t60 + 0.1e1);
	t1 = [-t42 * t38 * t52, t30, 0, 0, 0, 0; (-t27 * t53 - (t36 * t42 * t44 * t56 + (t38 - 0.1e1) * t47 * t35) * t47 * t60) * t25, (-t45 * t27 - (t46 * t58 + t36 * t47 + (-t36 * t53 - t58) * t30) * t47 * t28) * t48 * t25, 0, 0, 0, 0; ((t45 * t54 + t51) * t31 - (-t45 * t55 + t50) * t59) * t26, (-t31 * t40 - t39 * t59) * t26 * t52, 0, 0, t49 * t26, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:34:01
	% EndTime: 2019-10-10 09:34:01
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (446->26), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
	t57 = cos(qJ(2));
	t71 = t57 ^ 2;
	t55 = sin(qJ(2));
	t50 = pkin(9) + qJ(5);
	t48 = sin(t50);
	t58 = cos(qJ(1));
	t61 = t58 * t48;
	t49 = cos(t50);
	t56 = sin(qJ(1));
	t64 = t56 * t49;
	t43 = t55 * t64 + t61;
	t63 = t57 * t49;
	t37 = atan2(t43, t63);
	t34 = sin(t37);
	t35 = cos(t37);
	t33 = t34 * t43 + t35 * t63;
	t32 = 0.1e1 / t33 ^ 2;
	t60 = t58 * t49;
	t65 = t56 * t48;
	t41 = -t55 * t60 + t65;
	t70 = t32 * t41;
	t68 = t35 * t43;
	t67 = t41 ^ 2 * t32;
	t46 = 0.1e1 / t49;
	t52 = 0.1e1 / t57;
	t66 = t46 * t52;
	t62 = t57 * t58;
	t42 = t55 * t61 + t64;
	t40 = 0.1e1 / t42 ^ 2;
	t59 = t58 ^ 2 * t71 * t40;
	t53 = 0.1e1 / t71;
	t47 = 0.1e1 / t49 ^ 2;
	t44 = -t55 * t65 + t60;
	t39 = 0.1e1 / t42;
	t38 = 0.1e1 / (0.1e1 + t59);
	t36 = 0.1e1 / (t43 ^ 2 * t53 * t47 + 0.1e1);
	t31 = 0.1e1 / t33;
	t30 = (t43 * t46 * t53 * t55 + t56) * t36;
	t29 = 0.1e1 / (0.1e1 + t67);
	t28 = (t43 * t47 * t48 + t44 * t46) * t52 * t36;
	t1 = [-t41 * t36 * t66, t30, 0, 0, t28, 0; (t43 * t31 - (-t34 + (-t66 * t68 + t34) * t36) * t67) * t29, (-t30 * t68 * t70 + (-t31 * t62 - (-t35 * t55 + (-t30 + t56) * t57 * t34) * t70) * t49) * t29, 0, 0, (t42 * t31 - (-t35 * t57 * t48 + t34 * t44 + (-t34 * t63 + t68) * t28) * t70) * t29, 0; (t40 * t44 * t58 + t39 * t56) * t57 * t38, (t39 * t55 * t58 + t48 * t59) * t38, 0, 0, -t41 * t40 * t38 * t62, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end