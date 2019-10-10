% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR8
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
%   Wie in S6RRPPRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:48
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
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
	t30 = sin(pkin(10));
	t31 = cos(pkin(10));
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
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (102->20), mult. (329->52), div. (61->11), fcn. (494->9), ass. (0->34)
	t45 = sin(qJ(2));
	t59 = t45 ^ 2;
	t43 = sin(pkin(10));
	t44 = cos(pkin(10));
	t48 = cos(qJ(1));
	t50 = t48 * t44;
	t46 = sin(qJ(1));
	t47 = cos(qJ(2));
	t52 = t46 * t47;
	t33 = t43 * t52 + t50;
	t53 = t45 * t43;
	t29 = atan2(-t33, t53);
	t26 = sin(t29);
	t27 = cos(t29);
	t25 = -t26 * t33 + t27 * t53;
	t24 = 0.1e1 / t25 ^ 2;
	t51 = t48 * t43;
	t35 = -t46 * t44 + t47 * t51;
	t58 = t24 * t35;
	t56 = t27 * t33;
	t55 = t35 ^ 2 * t24;
	t38 = 0.1e1 / t43;
	t54 = t38 / t45;
	t36 = t46 * t43 + t47 * t50;
	t32 = 0.1e1 / t36 ^ 2;
	t49 = t48 ^ 2 * t59 * t32;
	t41 = 0.1e1 / t59;
	t31 = 0.1e1 / t36;
	t30 = 0.1e1 / (0.1e1 + t49);
	t28 = 0.1e1 / (0.1e1 + t33 ^ 2 * t41 / t43 ^ 2);
	t23 = 0.1e1 / t25;
	t22 = (t33 * t38 * t41 * t47 + t46) * t28;
	t21 = 0.1e1 / (0.1e1 + t55);
	t1 = [-t35 * t28 * t54, t22, 0, 0, 0, 0; (-t33 * t23 - (-t26 + (t54 * t56 + t26) * t28) * t55) * t21, (t22 * t56 * t58 + (-t48 * t45 * t23 - (t27 * t47 + (-t22 + t46) * t45 * t26) * t58) * t43) * t21, 0, 0, 0, 0; (t46 * t31 + (-t44 * t52 + t51) * t48 * t32) * t45 * t30, (-t31 * t47 * t48 - t44 * t49) * t30, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (111->24), mult. (347->67), div. (52->9), fcn. (503->11), ass. (0->39)
	t49 = cos(pkin(10));
	t52 = sin(qJ(1));
	t54 = cos(qJ(2));
	t48 = sin(pkin(10));
	t55 = cos(qJ(1));
	t58 = t55 * t48;
	t39 = -t52 * t49 + t54 * t58;
	t57 = t55 * t49;
	t40 = t52 * t48 + t54 * t57;
	t50 = sin(qJ(5));
	t53 = cos(qJ(5));
	t32 = t39 * t50 + t40 * t53;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = -t39 * t53 + t40 * t50;
	t65 = t30 * t31;
	t51 = sin(qJ(2));
	t60 = t52 * t51;
	t44 = atan2(t60, t54);
	t41 = sin(t44);
	t42 = cos(t44);
	t35 = t41 * t60 + t42 * t54;
	t34 = 0.1e1 / t35 ^ 2;
	t64 = t34 * t55 ^ 2;
	t45 = t51 ^ 2;
	t63 = t45 / t54 ^ 2;
	t62 = t51 * t55;
	t43 = 0.1e1 / (t52 ^ 2 * t63 + 0.1e1);
	t61 = t52 * t43;
	t59 = t52 * t54;
	t56 = t31 ^ 2 * t30 + 0.1e1;
	t46 = 0.1e1 / t54;
	t38 = -t49 * t59 + t58;
	t37 = -t48 * t59 - t57;
	t36 = (0.1e1 + t63) * t61;
	t33 = 0.1e1 / t35;
	t29 = 0.1e1 / t32;
	t28 = 0.1e1 / (t45 * t64 + 0.1e1);
	t27 = 0.1e1 / t56;
	t1 = [t46 * t43 * t62, t36, 0, 0, 0, 0; (t33 * t60 + (t42 * t45 * t46 * t61 + (-t43 + 0.1e1) * t51 * t41) * t51 * t64) * t28, (-t54 * t33 + (t41 * t59 - t42 * t51 + (-t41 * t54 + t42 * t60) * t36) * t51 * t34) * t55 * t28, 0, 0, 0, 0; ((-t37 * t53 + t38 * t50) * t29 - (t37 * t50 + t38 * t53) * t65) * t27, ((t48 * t53 - t49 * t50) * t29 - (-t48 * t50 - t49 * t53) * t65) * t27 * t62, 0, 0, t56 * t27, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (199->25), mult. (409->67), div. (57->9), fcn. (587->11), ass. (0->41)
	t69 = cos(pkin(10));
	t71 = sin(qJ(1));
	t72 = cos(qJ(2));
	t68 = sin(pkin(10));
	t73 = cos(qJ(1));
	t76 = t73 * t68;
	t56 = -t71 * t69 + t72 * t76;
	t75 = t73 * t69;
	t57 = t71 * t68 + t72 * t75;
	t67 = qJ(5) + qJ(6);
	t62 = sin(t67);
	t63 = cos(t67);
	t48 = t56 * t62 + t57 * t63;
	t46 = 0.1e1 / t48 ^ 2;
	t47 = -t56 * t63 + t57 * t62;
	t83 = t46 * t47;
	t70 = sin(qJ(2));
	t78 = t71 * t70;
	t61 = atan2(t78, t72);
	t58 = sin(t61);
	t59 = cos(t61);
	t52 = t58 * t78 + t59 * t72;
	t51 = 0.1e1 / t52 ^ 2;
	t82 = t51 * t73 ^ 2;
	t64 = t70 ^ 2;
	t81 = t64 / t72 ^ 2;
	t80 = t70 * t73;
	t60 = 0.1e1 / (t71 ^ 2 * t81 + 0.1e1);
	t79 = t71 * t60;
	t77 = t71 * t72;
	t74 = t47 ^ 2 * t46 + 0.1e1;
	t65 = 0.1e1 / t72;
	t55 = -t69 * t77 + t76;
	t54 = -t68 * t77 - t75;
	t53 = (0.1e1 + t81) * t79;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (t64 * t82 + 0.1e1);
	t45 = 0.1e1 / t48;
	t44 = 0.1e1 / t74;
	t43 = t74 * t44;
	t1 = [t65 * t60 * t80, t53, 0, 0, 0, 0; (t50 * t78 + (t59 * t64 * t65 * t79 + (-t60 + 0.1e1) * t70 * t58) * t70 * t82) * t49, (-t72 * t50 + (t58 * t77 - t59 * t70 + (-t58 * t72 + t59 * t78) * t53) * t70 * t51) * t73 * t49, 0, 0, 0, 0; ((-t54 * t63 + t55 * t62) * t45 - (t54 * t62 + t55 * t63) * t83) * t44, ((-t62 * t69 + t63 * t68) * t45 - (-t62 * t68 - t63 * t69) * t83) * t44 * t80, 0, 0, t43, t43;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end