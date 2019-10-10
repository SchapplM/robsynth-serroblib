% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14V3
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
%   Wie in S6RRPRRR14V3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14V3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:22
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->16), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
	t27 = cos(qJ(1));
	t28 = t27 ^ 2;
	t26 = cos(qJ(2));
	t24 = sin(qJ(2));
	t25 = sin(qJ(1));
	t29 = t25 * t24;
	t16 = atan2(-t29, -t26);
	t14 = sin(t16);
	t15 = cos(t16);
	t12 = -t14 * t29 - t15 * t26;
	t11 = 0.1e1 / t12 ^ 2;
	t34 = t11 * t24;
	t33 = t14 * t26;
	t19 = t24 ^ 2;
	t22 = 0.1e1 / t26 ^ 2;
	t32 = t19 * t22;
	t20 = t25 ^ 2;
	t31 = t20 / t28;
	t17 = 0.1e1 / (t20 * t32 + 0.1e1);
	t30 = t25 * t17;
	t21 = 0.1e1 / t26;
	t18 = 0.1e1 / (t22 * t31 + 0.1e1);
	t13 = (0.1e1 + t32) * t30;
	t10 = 0.1e1 / t12;
	t9 = 0.1e1 / (t28 * t19 * t11 + 0.1e1);
	t1 = [t27 * t24 * t21 * t17, t13, 0, 0, 0, 0; (-t10 * t29 - (-t15 * t19 * t21 * t30 + (t17 - 0.1e1) * t24 * t14) * t28 * t34) * t9, (t26 * t10 - (-t25 * t33 + t15 * t24 + (-t15 * t29 + t33) * t13) * t34) * t9 * t27, 0, 0, 0, 0; (-0.1e1 - t31) * t21 * t18, -0.1e1 / t27 * t22 * t18 * t29, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t39 = cos(qJ(2));
	t36 = sin(qJ(2));
	t37 = sin(qJ(1));
	t45 = t37 * t36;
	t30 = atan2(-t45, -t39);
	t28 = sin(t30);
	t29 = cos(t30);
	t21 = -t28 * t45 - t29 * t39;
	t20 = 0.1e1 / t21 ^ 2;
	t40 = cos(qJ(1));
	t50 = t20 * t40 ^ 2;
	t35 = sin(qJ(4));
	t38 = cos(qJ(4));
	t42 = t40 * t38;
	t27 = t37 * t35 + t39 * t42;
	t25 = 0.1e1 / t27 ^ 2;
	t43 = t40 * t35;
	t26 = -t37 * t38 + t39 * t43;
	t49 = t25 * t26;
	t32 = t36 ^ 2;
	t48 = t32 / t39 ^ 2;
	t47 = t36 * t40;
	t31 = 0.1e1 / (t37 ^ 2 * t48 + 0.1e1);
	t46 = t37 * t31;
	t44 = t37 * t39;
	t41 = t26 ^ 2 * t25 + 0.1e1;
	t33 = 0.1e1 / t39;
	t24 = 0.1e1 / t27;
	t23 = (0.1e1 + t48) * t46;
	t22 = 0.1e1 / t41;
	t19 = 0.1e1 / t21;
	t18 = 0.1e1 / (t32 * t50 + 0.1e1);
	t1 = [t33 * t31 * t47, t23, 0, 0, 0, 0; (-t19 * t45 - (-t29 * t32 * t33 * t46 + (t31 - 0.1e1) * t36 * t28) * t36 * t50) * t18, (t39 * t19 - (-t28 * t44 + t29 * t36 + (t28 * t39 - t29 * t45) * t23) * t36 * t20) * t40 * t18, 0, 0, 0, 0; ((-t35 * t44 - t42) * t24 - (-t38 * t44 + t43) * t49) * t22, (-t24 * t35 + t38 * t49) * t22 * t47, 0, t41 * t22, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (216->32), mult. (662->85), div. (108->11), fcn. (985->11), ass. (0->45)
	t59 = sin(qJ(4));
	t63 = cos(qJ(4));
	t65 = cos(qJ(1));
	t67 = t65 * t63;
	t61 = sin(qJ(1));
	t64 = cos(qJ(2));
	t69 = t61 * t64;
	t47 = t59 * t69 + t67;
	t60 = sin(qJ(2));
	t73 = t60 * t59;
	t46 = atan2(-t47, t73);
	t43 = sin(t46);
	t44 = cos(t46);
	t37 = -t43 * t47 + t44 * t73;
	t36 = 0.1e1 / t37 ^ 2;
	t68 = t65 * t59;
	t50 = -t61 * t63 + t64 * t68;
	t78 = t36 * t50;
	t51 = t61 * t59 + t64 * t67;
	t58 = sin(qJ(5));
	t62 = cos(qJ(5));
	t70 = t60 * t65;
	t42 = t51 * t62 + t58 * t70;
	t40 = 0.1e1 / t42 ^ 2;
	t41 = t51 * t58 - t62 * t70;
	t77 = t40 * t41;
	t76 = t44 * t47;
	t75 = t50 ^ 2 * t36;
	t54 = 0.1e1 / t59;
	t56 = 0.1e1 / t60;
	t74 = t54 * t56;
	t72 = t60 * t61;
	t71 = t60 * t63;
	t66 = t41 ^ 2 * t40 + 0.1e1;
	t57 = 0.1e1 / t60 ^ 2;
	t55 = 0.1e1 / t59 ^ 2;
	t49 = t63 * t69 - t68;
	t45 = 0.1e1 / (t47 ^ 2 * t57 * t55 + 0.1e1);
	t39 = 0.1e1 / t42;
	t38 = 0.1e1 / t66;
	t35 = 0.1e1 / t37;
	t34 = (t47 * t54 * t57 * t64 + t61) * t45;
	t33 = 0.1e1 / (0.1e1 + t75);
	t32 = (t47 * t55 * t63 - t49 * t54) * t56 * t45;
	t1 = [-t50 * t45 * t74, t34, 0, t32, 0, 0; (-t47 * t35 - (-t43 + (t74 * t76 + t43) * t45) * t75) * t33, (t34 * t76 * t78 + (-t35 * t70 - (t44 * t64 + (-t34 * t60 + t72) * t43) * t78) * t59) * t33, 0, (t51 * t35 - (t44 * t71 - t43 * t49 + (-t43 * t73 - t76) * t32) * t78) * t33, 0, 0; ((-t49 * t58 + t62 * t72) * t39 - (-t49 * t62 - t58 * t72) * t77) * t38, ((-t58 * t71 - t62 * t64) * t39 - (t58 * t64 - t62 * t71) * t77) * t38 * t65, 0, (-t58 * t39 + t62 * t77) * t50 * t38, t66 * t38, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (609->45), mult. (1737->115), div. (115->9), fcn. (2479->13), ass. (0->58)
	t95 = sin(qJ(2));
	t98 = cos(qJ(5));
	t110 = t95 * t98;
	t101 = cos(qJ(1));
	t100 = cos(qJ(2));
	t99 = cos(qJ(4));
	t107 = t100 * t99;
	t94 = sin(qJ(4));
	t96 = sin(qJ(1));
	t84 = -t101 * t94 + t96 * t107;
	t93 = sin(qJ(5));
	t71 = -t96 * t110 + t84 * t93;
	t109 = t95 * t99;
	t81 = t100 * t98 + t93 * t109;
	t70 = atan2(-t71, t81);
	t67 = sin(t70);
	t68 = cos(t70);
	t61 = -t67 * t71 + t68 * t81;
	t60 = 0.1e1 / t61 ^ 2;
	t106 = t101 * t95;
	t105 = t100 * t101;
	t108 = t96 * t94;
	t87 = t99 * t105 + t108;
	t75 = -t98 * t106 + t87 * t93;
	t117 = t60 * t75;
	t116 = t60 * t75 ^ 2;
	t76 = t93 * t106 + t87 * t98;
	t86 = t94 * t105 - t96 * t99;
	t92 = sin(qJ(6));
	t97 = cos(qJ(6));
	t66 = t76 * t97 + t86 * t92;
	t64 = 0.1e1 / t66 ^ 2;
	t65 = t76 * t92 - t86 * t97;
	t115 = t64 * t65;
	t114 = t68 * t71;
	t80 = 0.1e1 / t81 ^ 2;
	t113 = t71 * t80;
	t112 = t86 * t98;
	t111 = t94 * t95;
	t104 = t94 * t106;
	t103 = t64 * t65 ^ 2 + 0.1e1;
	t102 = -t67 * t81 - t114;
	t73 = t96 * t95 * t93 + t84 * t98;
	t82 = -t100 * t93 + t98 * t109;
	t85 = t93 * t107 - t110;
	t83 = -t100 * t108 - t101 * t99;
	t79 = 0.1e1 / t81;
	t78 = t82 * t101;
	t77 = t81 * t96;
	t69 = 0.1e1 / (t71 ^ 2 * t80 + 0.1e1);
	t63 = 0.1e1 / t66;
	t62 = 0.1e1 / t103;
	t59 = 0.1e1 / t61;
	t58 = 0.1e1 / (0.1e1 + t116);
	t57 = (-t111 * t113 - t79 * t83) * t93 * t69;
	t56 = (t85 * t113 + t77 * t79) * t69;
	t55 = (t82 * t113 - t73 * t79) * t69;
	t1 = [-t75 * t79 * t69, t56, 0, t57, t55, 0; (-t71 * t59 - (-t67 + (t79 * t114 + t67) * t69) * t116) * t58, (-(t102 * t56 + t67 * t77 + t68 * t85) * t117 - t81 * t59 * t101) * t58, 0, (-t86 * t93 * t59 - ((-t68 * t111 - t67 * t83) * t93 + t102 * t57) * t117) * t58, (t76 * t59 - (t102 * t55 - t67 * t73 + t68 * t82) * t117) * t58, 0; ((-t73 * t92 - t83 * t97) * t63 - (-t73 * t97 + t83 * t92) * t115) * t62, ((t97 * t104 - t78 * t92) * t63 - (-t92 * t104 - t78 * t97) * t115) * t62, 0, ((-t92 * t112 - t87 * t97) * t63 - (-t97 * t112 + t87 * t92) * t115) * t62, (t97 * t115 - t92 * t63) * t75 * t62, t103 * t62;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end