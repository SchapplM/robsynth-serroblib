% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPP3
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
%   Wie in S6RRPRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
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
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (157->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
	t46 = cos(qJ(2));
	t44 = sin(qJ(2));
	t45 = sin(qJ(1));
	t52 = t45 * t44;
	t36 = atan2(-t52, -t46);
	t34 = sin(t36);
	t35 = cos(t36);
	t28 = -t34 * t52 - t35 * t46;
	t27 = 0.1e1 / t28 ^ 2;
	t47 = cos(qJ(1));
	t57 = t27 * t47 ^ 2;
	t40 = pkin(9) + qJ(4);
	t38 = sin(t40);
	t39 = cos(t40);
	t49 = t47 * t39;
	t33 = t45 * t38 + t46 * t49;
	t31 = 0.1e1 / t33 ^ 2;
	t50 = t47 * t38;
	t32 = -t45 * t39 + t46 * t50;
	t56 = t31 * t32;
	t41 = t44 ^ 2;
	t55 = t41 / t46 ^ 2;
	t54 = t44 * t47;
	t37 = 0.1e1 / (t45 ^ 2 * t55 + 0.1e1);
	t53 = t45 * t37;
	t51 = t45 * t46;
	t48 = t32 ^ 2 * t31 + 0.1e1;
	t42 = 0.1e1 / t46;
	t30 = 0.1e1 / t33;
	t29 = (0.1e1 + t55) * t53;
	t26 = 0.1e1 / t28;
	t25 = 0.1e1 / t48;
	t24 = 0.1e1 / (t41 * t57 + 0.1e1);
	t1 = [t42 * t37 * t54, t29, 0, 0, 0, 0; (-t26 * t52 - (-t35 * t41 * t42 * t53 + (t37 - 0.1e1) * t44 * t34) * t44 * t57) * t24, (t46 * t26 - (-t34 * t51 + t35 * t44 + (t34 * t46 - t35 * t52) * t29) * t44 * t27) * t47 * t24, 0, 0, 0, 0; ((-t38 * t51 - t49) * t30 - (-t39 * t51 + t50) * t56) * t25, (-t30 * t38 + t39 * t56) * t25 * t54, 0, t48 * t25, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (437->27), mult. (486->70), div. (108->11), fcn. (755->9), ass. (0->39)
	t50 = pkin(9) + qJ(4);
	t48 = sin(t50);
	t49 = cos(t50);
	t58 = cos(qJ(1));
	t59 = t58 * t49;
	t56 = sin(qJ(1));
	t57 = cos(qJ(2));
	t61 = t56 * t57;
	t39 = t48 * t61 + t59;
	t55 = sin(qJ(2));
	t62 = t55 * t48;
	t37 = atan2(-t39, t62);
	t34 = sin(t37);
	t35 = cos(t37);
	t33 = -t34 * t39 + t35 * t62;
	t32 = 0.1e1 / t33 ^ 2;
	t60 = t58 * t48;
	t42 = -t56 * t49 + t57 * t60;
	t69 = t32 * t42;
	t67 = t35 * t39;
	t43 = t56 * t48 + t57 * t59;
	t52 = 0.1e1 / t55 ^ 2;
	t54 = 0.1e1 / t58 ^ 2;
	t38 = 0.1e1 / (t43 ^ 2 * t54 * t52 + 0.1e1);
	t51 = 0.1e1 / t55;
	t66 = t38 * t51;
	t65 = t42 ^ 2 * t32;
	t46 = 0.1e1 / t48;
	t64 = t46 * t51;
	t63 = t52 * t57;
	t53 = 0.1e1 / t58;
	t47 = 0.1e1 / t48 ^ 2;
	t41 = t49 * t61 - t60;
	t36 = 0.1e1 / (t39 ^ 2 * t52 * t47 + 0.1e1);
	t31 = 0.1e1 / t33;
	t30 = (t39 * t46 * t63 + t56) * t36;
	t29 = 0.1e1 / (0.1e1 + t65);
	t28 = (t39 * t47 * t49 - t41 * t46) * t51 * t36;
	t1 = [-t42 * t36 * t64, t30, 0, t28, 0, 0; (-t39 * t31 - (-t34 + (t64 * t67 + t34) * t36) * t65) * t29, (t30 * t67 * t69 + (-t58 * t55 * t31 - (t35 * t57 + (-t30 + t56) * t55 * t34) * t69) * t48) * t29, 0, (t43 * t31 - (t35 * t55 * t49 - t34 * t41 + (-t34 * t62 - t67) * t28) * t69) * t29, 0, 0; (t43 * t54 * t56 - t41 * t53) * t66, (-t43 * t53 * t63 - t49) * t38, 0, -t42 * t53 * t66, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (437->27), mult. (486->70), div. (108->11), fcn. (755->9), ass. (0->39)
	t49 = pkin(9) + qJ(4);
	t48 = cos(t49);
	t47 = sin(t49);
	t57 = cos(qJ(1));
	t59 = t57 * t47;
	t55 = sin(qJ(1));
	t56 = cos(qJ(2));
	t60 = t55 * t56;
	t40 = t48 * t60 - t59;
	t54 = sin(qJ(2));
	t61 = t54 * t48;
	t37 = atan2(-t40, t61);
	t34 = sin(t37);
	t35 = cos(t37);
	t33 = -t34 * t40 + t35 * t61;
	t32 = 0.1e1 / t33 ^ 2;
	t58 = t57 * t48;
	t43 = t55 * t47 + t56 * t58;
	t68 = t32 * t43;
	t66 = t35 * t40;
	t42 = t55 * t48 - t56 * t59;
	t51 = 0.1e1 / t54 ^ 2;
	t53 = 0.1e1 / t57 ^ 2;
	t38 = 0.1e1 / (t42 ^ 2 * t53 * t51 + 0.1e1);
	t50 = 0.1e1 / t54;
	t65 = t38 * t50;
	t64 = t43 ^ 2 * t32;
	t45 = 0.1e1 / t48;
	t63 = t45 * t50;
	t62 = t51 * t56;
	t52 = 0.1e1 / t57;
	t46 = 0.1e1 / t48 ^ 2;
	t39 = t47 * t60 + t58;
	t36 = 0.1e1 / (t40 ^ 2 * t51 * t46 + 0.1e1);
	t31 = 0.1e1 / t33;
	t30 = (t40 * t45 * t62 + t55) * t36;
	t29 = 0.1e1 / (0.1e1 + t64);
	t28 = (-t40 * t46 * t47 + t39 * t45) * t50 * t36;
	t1 = [-t43 * t36 * t63, t30, 0, t28, 0, 0; (-t40 * t31 - (-t34 + (t63 * t66 + t34) * t36) * t64) * t29, (t30 * t66 * t68 + (-t57 * t54 * t31 - (t35 * t56 + (-t30 + t55) * t54 * t34) * t68) * t48) * t29, 0, (t42 * t31 - (-t35 * t54 * t47 + t34 * t39 + (-t34 * t61 - t66) * t28) * t68) * t29, 0, 0; (t42 * t53 * t55 + t39 * t52) * t65, (-t42 * t52 * t62 + t47) * t38, 0, -t43 * t52 * t65, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end