% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPP1
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
%   Wie in S6RRRPPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t40 = cos(qJ(2));
	t37 = sin(qJ(2));
	t38 = sin(qJ(1));
	t46 = t38 * t37;
	t31 = atan2(-t46, -t40);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t46 - t30 * t40;
	t21 = 0.1e1 / t22 ^ 2;
	t41 = cos(qJ(1));
	t51 = t21 * t41 ^ 2;
	t36 = sin(qJ(3));
	t39 = cos(qJ(3));
	t43 = t41 * t39;
	t28 = t38 * t36 + t40 * t43;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t41 * t36;
	t27 = -t38 * t39 + t40 * t44;
	t50 = t26 * t27;
	t33 = t37 ^ 2;
	t49 = t33 / t40 ^ 2;
	t48 = t37 * t41;
	t32 = 0.1e1 / (t38 ^ 2 * t49 + 0.1e1);
	t47 = t38 * t32;
	t45 = t38 * t40;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t34 = 0.1e1 / t40;
	t25 = 0.1e1 / t28;
	t24 = (0.1e1 + t49) * t47;
	t23 = 0.1e1 / t42;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t51 + 0.1e1);
	t1 = [t34 * t32 * t48, t24, 0, 0, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19, (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23, (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (423->38), mult. (1241->99), div. (80->9), fcn. (1765->13), ass. (0->55)
	t68 = sin(qJ(3));
	t67 = cos(pkin(6));
	t69 = sin(qJ(2));
	t84 = t69 * t67;
	t65 = sin(pkin(6));
	t72 = cos(qJ(2));
	t85 = t65 * t72;
	t92 = t68 * t84 + t85;
	t71 = cos(qJ(3));
	t73 = cos(qJ(1));
	t78 = t73 * t71;
	t70 = sin(qJ(1));
	t80 = t70 * t72;
	t58 = t68 * t80 + t78;
	t83 = t69 * t70;
	t50 = -t58 * t65 - t67 * t83;
	t56 = t69 * t68 * t65 - t72 * t67;
	t49 = atan2(t50, t56);
	t46 = sin(t49);
	t47 = cos(t49);
	t40 = t46 * t50 + t47 * t56;
	t39 = 0.1e1 / t40 ^ 2;
	t79 = t73 * t68;
	t60 = t70 * t71 - t72 * t79;
	t81 = t69 * t73;
	t51 = t60 * t65 - t67 * t81;
	t91 = t39 * t51;
	t61 = t70 * t68 + t72 * t78;
	t64 = sin(pkin(10));
	t66 = cos(pkin(10));
	t74 = t60 * t67 + t65 * t81;
	t45 = t61 * t66 + t64 * t74;
	t43 = 0.1e1 / t45 ^ 2;
	t44 = t61 * t64 - t66 * t74;
	t90 = t43 * t44;
	t89 = t47 * t50;
	t55 = 0.1e1 / t56 ^ 2;
	t88 = t50 * t55;
	t87 = t51 ^ 2 * t39;
	t86 = t61 * t67;
	t82 = t69 * t71;
	t76 = -t46 * t56 + t89;
	t75 = -t58 * t67 + t65 * t83;
	t59 = -t71 * t80 + t79;
	t57 = t68 * t85 + t84;
	t54 = 0.1e1 / t56;
	t53 = t56 * t70;
	t48 = 0.1e1 / (t50 ^ 2 * t55 + 0.1e1);
	t42 = 0.1e1 / t45;
	t41 = 0.1e1 / (t44 ^ 2 * t43 + 0.1e1);
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / (0.1e1 + t87);
	t36 = (t54 * t59 - t82 * t88) * t65 * t48;
	t35 = (t53 * t54 - t57 * t88) * t48;
	t1 = [t51 * t54 * t48, t35, t36, 0, 0, 0; (t50 * t38 + (t46 + (t54 * t89 - t46) * t48) * t87) * t37, ((t76 * t35 + t46 * t53 + t47 * t57) * t91 - t56 * t38 * t73) * t37, (t61 * t65 * t38 + ((t46 * t59 + t47 * t82) * t65 + t76 * t36) * t91) * t37, 0, 0, 0; ((t59 * t64 + t66 * t75) * t42 - (t59 * t66 - t64 * t75) * t90) * t41, ((-t64 * t82 - t92 * t66) * t42 - (t92 * t64 - t66 * t82) * t90) * t41 * t73, ((t60 * t64 + t66 * t86) * t42 - (t60 * t66 - t64 * t86) * t90) * t41, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:13
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (665->45), mult. (2053->108), div. (80->9), fcn. (2819->13), ass. (0->60)
	t105 = sin(pkin(6));
	t107 = cos(pkin(6));
	t112 = cos(qJ(2));
	t108 = sin(qJ(3));
	t109 = sin(qJ(2));
	t126 = t108 * t109;
	t114 = t105 * t112 + t107 * t126;
	t104 = sin(pkin(10));
	t106 = cos(pkin(10));
	t111 = cos(qJ(3));
	t124 = t109 * t111;
	t87 = t104 * t124 + t114 * t106;
	t110 = sin(qJ(1));
	t130 = t105 * t109;
	t117 = t106 * t130;
	t128 = t106 * t107;
	t113 = cos(qJ(1));
	t119 = t113 * t111;
	t122 = t110 * t108;
	t95 = t112 * t122 + t119;
	t120 = t113 * t108;
	t121 = t111 * t112;
	t96 = t110 * t121 - t120;
	t79 = t96 * t104 - t110 * t117 + t95 * t128;
	t77 = atan2(-t79, t87);
	t74 = sin(t77);
	t75 = cos(t77);
	t73 = -t74 * t79 + t75 * t87;
	t72 = 0.1e1 / t73 ^ 2;
	t98 = t112 * t119 + t122;
	t131 = t98 * t104;
	t97 = t110 * t111 - t112 * t120;
	t81 = -t113 * t117 - t97 * t128 + t131;
	t133 = t81 ^ 2 * t72;
	t70 = 0.1e1 / (0.1e1 + t133);
	t71 = 0.1e1 / t73;
	t137 = t70 * t71;
	t136 = t72 * t81;
	t135 = t75 * t79;
	t86 = 0.1e1 / t87 ^ 2;
	t134 = t79 * t86;
	t123 = t109 * t113;
	t82 = t98 * t106 + (t105 * t123 + t107 * t97) * t104;
	t91 = -t97 * t105 + t107 * t123;
	t90 = 0.1e1 / t91 ^ 2;
	t132 = t82 * t90;
	t127 = t107 * t112;
	t125 = t109 * t110;
	t115 = -t74 * t87 - t135;
	t94 = (-t104 * t108 + t111 * t128) * t109;
	t89 = 0.1e1 / t91;
	t88 = t104 * t121 + (t108 * t127 - t130) * t106;
	t85 = 0.1e1 / t87;
	t84 = t87 * t110;
	t83 = -t95 * t104 + t96 * t128;
	t78 = 0.1e1 / (t82 ^ 2 * t90 + 0.1e1);
	t76 = 0.1e1 / (t79 ^ 2 * t86 + 0.1e1);
	t69 = (t94 * t134 - t83 * t85) * t76;
	t68 = (t88 * t134 + t84 * t85) * t76;
	t1 = [-t81 * t85 * t76, t68, t69, 0, 0, 0; -t79 * t137 - (-t74 + (t85 * t135 + t74) * t76) * t70 * t133, -(t115 * t68 + t74 * t84 + t75 * t88) * t70 * t136 - t87 * t113 * t137, ((t97 * t104 + t98 * t128) * t71 - (t115 * t69 - t74 * t83 + t75 * t94) * t136) * t70, 0, 0, 0; ((-t96 * t106 + (-t105 * t125 + t107 * t95) * t104) * t89 - (-t95 * t105 - t107 * t125) * t132) * t78, ((t114 * t104 - t106 * t124) * t89 - (-t105 * t126 + t127) * t132) * t78 * t113, ((t97 * t106 - t107 * t131) * t89 - t98 * t105 * t132) * t78, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:13
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (665->43), mult. (2053->105), div. (80->9), fcn. (2819->13), ass. (0->60)
	t118 = sin(pkin(6));
	t120 = cos(pkin(6));
	t125 = cos(qJ(2));
	t121 = sin(qJ(3));
	t122 = sin(qJ(2));
	t140 = t121 * t122;
	t127 = t118 * t125 + t120 * t140;
	t117 = sin(pkin(10));
	t119 = cos(pkin(10));
	t124 = cos(qJ(3));
	t138 = t122 * t124;
	t151 = t127 * t117 - t119 * t138;
	t123 = sin(qJ(1));
	t126 = cos(qJ(1));
	t134 = t126 * t121;
	t114 = t123 * t124 - t125 * t134;
	t137 = t122 * t126;
	t128 = t114 * t120 + t118 * t137;
	t133 = t126 * t124;
	t136 = t123 * t121;
	t115 = t125 * t133 + t136;
	t144 = t115 * t119;
	t101 = t128 * t117 + t144;
	t135 = t124 * t125;
	t113 = t123 * t135 - t134;
	t112 = t125 * t136 + t133;
	t139 = t122 * t123;
	t129 = t112 * t120 - t118 * t139;
	t98 = t113 * t119 - t129 * t117;
	t96 = atan2(-t98, -t151);
	t93 = sin(t96);
	t94 = cos(t96);
	t92 = -t151 * t94 - t93 * t98;
	t91 = 0.1e1 / t92 ^ 2;
	t147 = t101 ^ 2 * t91;
	t89 = 0.1e1 / (0.1e1 + t147);
	t90 = 0.1e1 / t92;
	t150 = t89 * t90;
	t149 = t94 * t98;
	t148 = t101 * t91;
	t105 = 0.1e1 / t151 ^ 2;
	t146 = t105 * t98;
	t100 = -t115 * t117 + t128 * t119;
	t110 = -t114 * t118 + t120 * t137;
	t109 = 0.1e1 / t110 ^ 2;
	t145 = t100 * t109;
	t143 = t117 * t120;
	t141 = t120 * t125;
	t130 = t151 * t93 - t149;
	t111 = (-t119 * t121 - t124 * t143) * t122;
	t108 = 0.1e1 / t110;
	t107 = t119 * t135 + (t118 * t122 - t121 * t141) * t117;
	t104 = 0.1e1 / t151;
	t103 = t151 * t123;
	t102 = -t112 * t119 - t113 * t143;
	t97 = 0.1e1 / (t100 ^ 2 * t109 + 0.1e1);
	t95 = 0.1e1 / (t98 ^ 2 * t105 + 0.1e1);
	t88 = (t102 * t104 + t111 * t146) * t95;
	t87 = (t103 * t104 + t107 * t146) * t95;
	t1 = [t101 * t104 * t95, t87, t88, 0, 0, 0; -t98 * t150 - (-t93 + (-t104 * t149 + t93) * t95) * t89 * t147, -(-t93 * t103 + t94 * t107 + t130 * t87) * t89 * t148 + t151 * t126 * t150, ((t114 * t119 - t115 * t143) * t90 - (-t93 * t102 + t94 * t111 + t130 * t88) * t148) * t89, 0, 0, 0; ((t113 * t117 + t129 * t119) * t108 - (-t112 * t118 - t120 * t139) * t145) * t97, ((t117 * t138 + t127 * t119) * t108 - (-t118 * t140 + t141) * t145) * t97 * t126, ((-t114 * t117 - t120 * t144) * t108 - t115 * t118 * t145) * t97, 0, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end