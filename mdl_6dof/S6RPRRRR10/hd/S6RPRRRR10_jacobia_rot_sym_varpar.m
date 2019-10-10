% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR10
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
%   Wie in S6RPRRRR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (25->13), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(6));
	t26 = sin(pkin(6));
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t22 = atan2(t32, t28);
	t19 = sin(t22);
	t20 = cos(t22);
	t14 = t19 * t32 + t20 * t28;
	t29 = sin(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t29 ^ 2;
	t23 = t26 ^ 2;
	t21 = 0.1e1 / (0.1e1 + t30 ^ 2 * t23 / t28 ^ 2);
	t36 = t21 / t28;
	t25 = sin(pkin(13));
	t35 = t29 * t25;
	t27 = cos(pkin(13));
	t34 = t29 * t27;
	t33 = t30 * t25;
	t31 = t30 * t27;
	t18 = -t28 * t35 + t31;
	t17 = t28 * t34 + t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [-t29 * t26 * t36, 0, 0, 0, 0, 0; (0.1e1 / t14 * t32 - (-t20 * t23 * t30 * t36 + (t21 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0, 0, 0; ((t28 * t31 - t35) / t18 - (-t28 * t33 - t34) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (156->21), mult. (450->48), div. (25->9), fcn. (637->13), ass. (0->38)
	t61 = cos(pkin(6));
	t59 = cos(pkin(13));
	t65 = cos(qJ(1));
	t69 = t65 * t59;
	t56 = sin(pkin(13));
	t63 = sin(qJ(1));
	t72 = t63 * t56;
	t51 = -t61 * t69 + t72;
	t57 = sin(pkin(7));
	t60 = cos(pkin(7));
	t58 = sin(pkin(6));
	t73 = t58 * t65;
	t46 = -t51 * t57 + t60 * t73;
	t50 = -t58 * t59 * t57 + t61 * t60;
	t45 = atan2(t46, t50);
	t42 = sin(t45);
	t43 = cos(t45);
	t37 = t42 * t46 + t43 * t50;
	t70 = t65 * t56;
	t71 = t63 * t59;
	t53 = -t61 * t71 - t70;
	t74 = t58 * t63;
	t47 = t53 * t57 - t60 * t74;
	t75 = t47 ^ 2 / t37 ^ 2;
	t54 = -t61 * t72 + t69;
	t62 = sin(qJ(3));
	t64 = cos(qJ(3));
	t66 = t53 * t60 + t57 * t74;
	t41 = t54 * t64 + t66 * t62;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t54 * t62 - t66 * t64;
	t68 = t40 ^ 2 * t39 + 0.1e1;
	t67 = t51 * t60 + t57 * t73;
	t52 = -t61 * t70 - t71;
	t49 = 0.1e1 / t50;
	t44 = 0.1e1 / (0.1e1 + t46 ^ 2 / t50 ^ 2);
	t38 = 0.1e1 / t68;
	t1 = [t47 * t49 * t44, 0, 0, 0, 0, 0; (t46 / t37 + (t42 + (t43 * t46 * t49 - t42) * t44) * t75) / (0.1e1 + t75), 0, 0, 0, 0, 0; ((t52 * t62 - t67 * t64) / t41 - (t52 * t64 + t67 * t62) * t40 * t39) * t38, 0, t68 * t38, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (559->39), mult. (1669->87), div. (55->9), fcn. (2293->15), ass. (0->57)
	t115 = sin(qJ(1));
	t88 = sin(pkin(6));
	t102 = t88 * t115;
	t87 = sin(pkin(7));
	t90 = cos(pkin(7));
	t89 = cos(pkin(13));
	t100 = t115 * t89;
	t86 = sin(pkin(13));
	t96 = cos(qJ(1));
	t106 = t96 * t86;
	t91 = cos(pkin(6));
	t98 = t91 * t100 + t106;
	t116 = -t87 * t102 + t98 * t90;
	t101 = t115 * t86;
	t105 = t96 * t89;
	t83 = -t91 * t101 + t105;
	t93 = sin(qJ(3));
	t95 = cos(qJ(3));
	t72 = -t116 * t93 + t83 * t95;
	t78 = t90 * t102 + t98 * t87;
	t92 = sin(qJ(4));
	t94 = cos(qJ(4));
	t62 = t72 * t94 + t78 * t92;
	t60 = 0.1e1 / t62 ^ 2;
	t61 = t72 * t92 - t78 * t94;
	t114 = t60 * t61;
	t109 = t88 * t96;
	t104 = t87 * t109;
	t107 = t90 * t95;
	t82 = t91 * t106 + t100;
	t111 = t82 * t93;
	t81 = -t91 * t105 + t101;
	t67 = t95 * t104 + t81 * t107 + t111;
	t110 = t87 * t91;
	t75 = -t95 * t110 + (-t89 * t107 + t86 * t93) * t88;
	t66 = atan2(-t67, t75);
	t64 = cos(t66);
	t113 = t64 * t67;
	t63 = sin(t66);
	t57 = -t63 * t67 + t64 * t75;
	t56 = 0.1e1 / t57 ^ 2;
	t71 = t116 * t95 + t83 * t93;
	t112 = t71 ^ 2 * t56;
	t108 = t90 * t93;
	t103 = t61 ^ 2 * t60 + 0.1e1;
	t70 = t93 * t104 + t81 * t108 - t82 * t95;
	t77 = t90 * t109 - t81 * t87;
	t76 = t93 * t110 + (t89 * t108 + t86 * t95) * t88;
	t74 = 0.1e1 / t75 ^ 2;
	t73 = 0.1e1 / t75;
	t65 = 0.1e1 / (t67 ^ 2 * t74 + 0.1e1);
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / t103;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (0.1e1 + t112);
	t53 = (t67 * t74 * t76 + t70 * t73) * t65;
	t1 = [-t71 * t73 * t65, 0, t53, 0, 0, 0; ((-t111 + (-t81 * t90 - t104) * t95) * t55 - (-t63 + (t73 * t113 + t63) * t65) * t112) * t54, 0, (t72 * t55 - (t63 * t70 + t64 * t76 + (-t63 * t75 - t113) * t53) * t71 * t56) * t54, 0, 0, 0; ((t70 * t92 - t77 * t94) * t59 - (t70 * t94 + t77 * t92) * t114) * t58, 0, (t94 * t114 - t92 * t59) * t71 * t58, t103 * t58, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (673->40), mult. (1815->87), div. (60->9), fcn. (2489->15), ass. (0->59)
	t105 = sin(pkin(7));
	t108 = cos(pkin(7));
	t109 = cos(pkin(6));
	t107 = cos(pkin(13));
	t131 = sin(qJ(1));
	t116 = t131 * t107;
	t104 = sin(pkin(13));
	t112 = cos(qJ(1));
	t122 = t112 * t104;
	t114 = t109 * t116 + t122;
	t106 = sin(pkin(6));
	t118 = t106 * t131;
	t132 = -t105 * t118 + t114 * t108;
	t111 = cos(qJ(3));
	t125 = t106 * t112;
	t119 = t105 * t125;
	t123 = t108 * t111;
	t110 = sin(qJ(3));
	t97 = t109 * t122 + t116;
	t127 = t97 * t110;
	t117 = t131 * t104;
	t121 = t112 * t107;
	t96 = -t109 * t121 + t117;
	t82 = t111 * t119 + t96 * t123 + t127;
	t126 = t105 * t109;
	t90 = -t111 * t126 + (t104 * t110 - t107 * t123) * t106;
	t81 = atan2(-t82, t90);
	t78 = sin(t81);
	t79 = cos(t81);
	t72 = -t78 * t82 + t79 * t90;
	t71 = 0.1e1 / t72 ^ 2;
	t98 = -t109 * t117 + t121;
	t86 = t98 * t110 + t132 * t111;
	t130 = t71 * t86 ^ 2;
	t103 = qJ(4) + qJ(5);
	t101 = sin(t103);
	t102 = cos(t103);
	t87 = -t132 * t110 + t98 * t111;
	t93 = t114 * t105 + t108 * t118;
	t77 = t93 * t101 + t87 * t102;
	t75 = 0.1e1 / t77 ^ 2;
	t76 = t87 * t101 - t93 * t102;
	t129 = t75 * t76;
	t128 = t79 * t82;
	t124 = t108 * t110;
	t120 = t75 * t76 ^ 2 + 0.1e1;
	t85 = t110 * t119 - t97 * t111 + t96 * t124;
	t92 = -t96 * t105 + t108 * t125;
	t91 = t110 * t126 + (t104 * t111 + t107 * t124) * t106;
	t89 = 0.1e1 / t90 ^ 2;
	t88 = 0.1e1 / t90;
	t80 = 0.1e1 / (t82 ^ 2 * t89 + 0.1e1);
	t74 = 0.1e1 / t77;
	t73 = 0.1e1 / t120;
	t70 = 0.1e1 / t72;
	t69 = 0.1e1 / (0.1e1 + t130);
	t68 = (t82 * t89 * t91 + t85 * t88) * t80;
	t67 = t120 * t73;
	t1 = [-t86 * t88 * t80, 0, t68, 0, 0, 0; ((-t127 + (-t108 * t96 - t119) * t111) * t70 - (-t78 + (t88 * t128 + t78) * t80) * t130) * t69, 0, (t87 * t70 - (t78 * t85 + t79 * t91 + (-t78 * t90 - t128) * t68) * t86 * t71) * t69, 0, 0, 0; ((t85 * t101 - t92 * t102) * t74 - (t92 * t101 + t85 * t102) * t129) * t73, 0, (-t101 * t74 + t102 * t129) * t86 * t73, t67, t67, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (2501->51), mult. (5521->119), div. (115->9), fcn. (7584->17), ass. (0->71)
	t145 = cos(pkin(6));
	t143 = cos(pkin(13));
	t175 = sin(qJ(1));
	t158 = t175 * t143;
	t140 = sin(pkin(13));
	t150 = cos(qJ(1));
	t162 = t150 * t140;
	t134 = t145 * t162 + t158;
	t147 = sin(qJ(3));
	t149 = cos(qJ(3));
	t159 = t175 * t140;
	t161 = t150 * t143;
	t133 = -t145 * t161 + t159;
	t141 = sin(pkin(7));
	t144 = cos(pkin(7));
	t142 = sin(pkin(6));
	t164 = t142 * t150;
	t155 = t133 * t144 + t141 * t164;
	t123 = -t134 * t149 + t155 * t147;
	t130 = -t133 * t141 + t144 * t164;
	t139 = qJ(4) + qJ(5);
	t137 = sin(t139);
	t138 = cos(t139);
	t177 = t123 * t137 - t130 * t138;
	t113 = t123 * t138 + t130 * t137;
	t154 = t145 * t158 + t162;
	t160 = t142 * t175;
	t176 = -t141 * t160 + t154 * t144;
	t135 = -t145 * t159 + t161;
	t125 = t135 * t149 - t176 * t147;
	t151 = t154 * t141 + t144 * t160;
	t114 = t125 * t137 - t151 * t138;
	t163 = t143 * t144;
	t165 = t141 * t145;
	t129 = t147 * t165 + (t140 * t149 + t147 * t163) * t142;
	t132 = -t142 * t143 * t141 + t145 * t144;
	t118 = t129 * t137 - t132 * t138;
	t109 = atan2(t177, t118);
	t102 = sin(t109);
	t103 = cos(t109);
	t100 = t102 * t177 + t103 * t118;
	t99 = 0.1e1 / t100 ^ 2;
	t174 = t114 * t99;
	t173 = t114 ^ 2 * t99;
	t115 = t125 * t138 + t151 * t137;
	t148 = cos(qJ(6));
	t124 = t135 * t147 + t176 * t149;
	t146 = sin(qJ(6));
	t170 = t124 * t146;
	t108 = t115 * t148 + t170;
	t105 = 0.1e1 / t108 ^ 2;
	t169 = t124 * t148;
	t107 = t115 * t146 - t169;
	t172 = t105 * t107;
	t117 = 0.1e1 / t118 ^ 2;
	t171 = t177 * t117;
	t157 = t107 ^ 2 * t105 + 0.1e1;
	t152 = -t134 * t147 - t155 * t149;
	t128 = t149 * t165 + (-t140 * t147 + t149 * t163) * t142;
	t119 = t129 * t138 + t132 * t137;
	t116 = 0.1e1 / t118;
	t106 = 0.1e1 / (t117 * t177 ^ 2 + 0.1e1);
	t104 = 0.1e1 / t108;
	t101 = 0.1e1 / t157;
	t98 = 0.1e1 / t100;
	t97 = 0.1e1 / (0.1e1 + t173);
	t96 = (-t116 * t152 - t128 * t171) * t137 * t106;
	t95 = (t113 * t116 - t119 * t171) * t106;
	t94 = (-t146 * t104 + t148 * t172) * t114 * t101;
	t93 = (t115 * t98 - ((t177 * t95 + t119) * t103 + (-t118 * t95 + t113) * t102) * t174) * t97;
	t1 = [-t114 * t116 * t106, 0, t96, t95, t95, 0; (t177 * t98 - (-t102 + (-t103 * t116 * t177 + t102) * t106) * t173) * t97, 0, (-t124 * t137 * t98 - ((t128 * t137 + t177 * t96) * t103 + (-t118 * t96 - t137 * t152) * t102) * t174) * t97, t93, t93, 0; ((t113 * t146 - t148 * t152) * t104 - (t113 * t148 + t146 * t152) * t172) * t101, 0, ((-t125 * t148 - t138 * t170) * t104 - (t125 * t146 - t138 * t169) * t172) * t101, t94, t94, t157 * t101;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end