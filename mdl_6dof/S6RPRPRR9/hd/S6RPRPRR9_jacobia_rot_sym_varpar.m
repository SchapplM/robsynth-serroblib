% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR9
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
%   Wie in S6RPRPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
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
	t25 = sin(pkin(12));
	t35 = t29 * t25;
	t27 = cos(pkin(12));
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (156->21), mult. (450->48), div. (25->9), fcn. (637->13), ass. (0->38)
	t61 = cos(pkin(6));
	t59 = cos(pkin(12));
	t65 = cos(qJ(1));
	t69 = t65 * t59;
	t56 = sin(pkin(12));
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (218->26), mult. (618->57), div. (26->9), fcn. (869->15), ass. (0->45)
	t74 = sin(pkin(7));
	t72 = sin(pkin(13));
	t76 = cos(pkin(13));
	t80 = sin(qJ(3));
	t82 = cos(qJ(3));
	t84 = t82 * t72 + t80 * t76;
	t59 = t84 * t74;
	t78 = cos(pkin(7));
	t61 = t84 * t78;
	t79 = cos(pkin(6));
	t73 = sin(pkin(12));
	t83 = cos(qJ(1));
	t86 = t83 * t73;
	t77 = cos(pkin(12));
	t81 = sin(qJ(1));
	t87 = t81 * t77;
	t66 = -t79 * t87 - t86;
	t85 = t83 * t77;
	t88 = t81 * t73;
	t67 = -t79 * t88 + t85;
	t68 = t80 * t72 - t82 * t76;
	t75 = sin(pkin(6));
	t90 = t75 * t81;
	t50 = t59 * t90 + t66 * t61 - t67 * t68;
	t47 = 0.1e1 / t50 ^ 2;
	t58 = t68 * t74;
	t60 = t68 * t78;
	t48 = -t58 * t90 - t66 * t60 - t67 * t84;
	t92 = t48 ^ 2 * t47;
	t64 = -t79 * t85 + t88;
	t89 = t75 * t83;
	t55 = -t64 * t74 + t78 * t89;
	t63 = -t75 * t77 * t74 + t79 * t78;
	t54 = atan2(t55, t63);
	t51 = sin(t54);
	t52 = cos(t54);
	t45 = t51 * t55 + t52 * t63;
	t56 = t66 * t74 - t78 * t90;
	t91 = t56 ^ 2 / t45 ^ 2;
	t65 = -t79 * t86 - t87;
	t62 = 0.1e1 / t63;
	t53 = 0.1e1 / (0.1e1 + t55 ^ 2 / t63 ^ 2);
	t46 = 0.1e1 / t50;
	t42 = 0.1e1 / (0.1e1 + t92);
	t1 = [t56 * t62 * t53, 0, 0, 0, 0, 0; (t55 / t45 + (t51 + (t52 * t55 * t62 - t51) * t53) * t91) / (0.1e1 + t91), 0, 0, 0, 0, 0; ((t58 * t89 + t64 * t60 + t65 * t84) * t46 + (t59 * t89 + t64 * t61 - t65 * t68) * t48 * t47) * t42, 0, (t50 * t46 + t92) * t42, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (976->38), mult. (2781->89), div. (55->9), fcn. (3822->17), ass. (0->59)
	t103 = sin(pkin(13));
	t107 = cos(pkin(13));
	t112 = sin(qJ(3));
	t115 = cos(qJ(3));
	t118 = t115 * t103 + t112 * t107;
	t106 = sin(pkin(6));
	t116 = cos(qJ(1));
	t124 = t106 * t116;
	t105 = sin(pkin(7));
	t99 = t112 * t103 - t115 * t107;
	t91 = t99 * t105;
	t109 = cos(pkin(7));
	t93 = t99 * t109;
	t110 = cos(pkin(6));
	t108 = cos(pkin(12));
	t120 = t116 * t108;
	t104 = sin(pkin(12));
	t113 = sin(qJ(1));
	t123 = t113 * t104;
	t95 = -t110 * t120 + t123;
	t121 = t116 * t104;
	t122 = t113 * t108;
	t96 = t110 * t121 + t122;
	t79 = -t118 * t96 + t91 * t124 + t95 * t93;
	t87 = t110 * t91 + (t104 * t118 + t108 * t93) * t106;
	t73 = atan2(t79, t87);
	t71 = cos(t73);
	t128 = t71 * t79;
	t111 = sin(qJ(5));
	t114 = cos(qJ(5));
	t125 = t106 * t113;
	t92 = t118 * t105;
	t94 = t118 * t109;
	t97 = -t110 * t122 - t121;
	t98 = -t110 * t123 + t120;
	t83 = t92 * t125 + t97 * t94 - t98 * t99;
	t89 = -t97 * t105 + t109 * t125;
	t77 = t89 * t111 + t83 * t114;
	t75 = 0.1e1 / t77 ^ 2;
	t76 = t83 * t111 - t89 * t114;
	t127 = t75 * t76;
	t70 = sin(t73);
	t68 = t70 * t79 + t71 * t87;
	t67 = 0.1e1 / t68 ^ 2;
	t81 = -t118 * t98 - t91 * t125 - t97 * t93;
	t126 = t81 ^ 2 * t67;
	t119 = t76 ^ 2 * t75 + 0.1e1;
	t117 = t92 * t124 + t95 * t94 + t96 * t99;
	t88 = -t95 * t105 + t109 * t124;
	t86 = t110 * t92 + (-t104 * t99 + t108 * t94) * t106;
	t85 = 0.1e1 / t87 ^ 2;
	t84 = 0.1e1 / t87;
	t74 = 0.1e1 / t77;
	t72 = 0.1e1 / (t79 ^ 2 * t85 + 0.1e1);
	t69 = 0.1e1 / t119;
	t66 = 0.1e1 / t68;
	t65 = 0.1e1 / (0.1e1 + t126);
	t64 = (-t79 * t85 * t86 + t117 * t84) * t72;
	t1 = [t81 * t84 * t72, 0, t64, 0, 0, 0; (t79 * t66 + (t70 + (t84 * t128 - t70) * t72) * t126) * t65, 0, (t83 * t66 + (t70 * t117 + t71 * t86 + (-t70 * t87 + t128) * t64) * t81 * t67) * t65, 0, 0, 0; ((t111 * t117 - t88 * t114) * t74 - (t88 * t111 + t114 * t117) * t127) * t69, 0, (t111 * t74 - t114 * t127) * t81 * t69, 0, t119 * t69, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (2181->53), mult. (6097->124), div. (85->9), fcn. (8377->19), ass. (0->73)
	t130 = sin(pkin(7));
	t128 = sin(pkin(13));
	t132 = cos(pkin(13));
	t138 = sin(qJ(3));
	t142 = cos(qJ(3));
	t148 = t142 * t128 + t138 * t132;
	t118 = t148 * t130;
	t134 = cos(pkin(7));
	t120 = t148 * t134;
	t135 = cos(pkin(6));
	t133 = cos(pkin(12));
	t143 = cos(qJ(1));
	t151 = t143 * t133;
	t129 = sin(pkin(12));
	t139 = sin(qJ(1));
	t154 = t139 * t129;
	t122 = -t135 * t151 + t154;
	t152 = t143 * t129;
	t153 = t139 * t133;
	t123 = t135 * t152 + t153;
	t125 = t138 * t128 - t142 * t132;
	t131 = sin(pkin(6));
	t155 = t131 * t143;
	t108 = t118 * t155 + t122 * t120 + t123 * t125;
	t116 = -t122 * t130 + t134 * t155;
	t137 = sin(qJ(5));
	t141 = cos(qJ(5));
	t165 = t108 * t137 - t116 * t141;
	t98 = t108 * t141 + t116 * t137;
	t113 = t135 * t118 + (t120 * t133 - t125 * t129) * t131;
	t121 = -t131 * t133 * t130 + t135 * t134;
	t103 = t113 * t137 - t121 * t141;
	t94 = atan2(t165, t103);
	t91 = sin(t94);
	t92 = cos(t94);
	t85 = t92 * t103 + t165 * t91;
	t84 = 0.1e1 / t85 ^ 2;
	t124 = -t135 * t154 + t151;
	t147 = -t135 * t153 - t152;
	t156 = t131 * t139;
	t144 = t118 * t156 + t147 * t120 - t124 * t125;
	t145 = -t147 * t130 + t134 * t156;
	t99 = t137 * t144 - t145 * t141;
	t164 = t84 * t99;
	t100 = t145 * t137 + t141 * t144;
	t117 = t125 * t130;
	t119 = t125 * t134;
	t110 = -t117 * t156 - t147 * t119 - t124 * t148;
	t136 = sin(qJ(6));
	t140 = cos(qJ(6));
	t90 = t100 * t140 - t110 * t136;
	t88 = 0.1e1 / t90 ^ 2;
	t89 = t100 * t136 + t110 * t140;
	t163 = t88 * t89;
	t162 = t92 * t165;
	t161 = t99 ^ 2 * t84;
	t102 = 0.1e1 / t103 ^ 2;
	t160 = t102 * t165;
	t159 = t110 * t141;
	t150 = t89 ^ 2 * t88 + 0.1e1;
	t149 = -t103 * t91 + t162;
	t146 = t117 * t155 + t122 * t119 - t123 * t148;
	t112 = -t135 * t117 + (-t119 * t133 - t129 * t148) * t131;
	t104 = t113 * t141 + t121 * t137;
	t101 = 0.1e1 / t103;
	t93 = 0.1e1 / (t102 * t165 ^ 2 + 0.1e1);
	t87 = 0.1e1 / t90;
	t86 = 0.1e1 / t150;
	t83 = 0.1e1 / t85;
	t82 = 0.1e1 / (0.1e1 + t161);
	t81 = (-t101 * t146 - t112 * t160) * t93 * t137;
	t80 = (t101 * t98 - t104 * t160) * t93;
	t1 = [-t99 * t101 * t93, 0, t81, 0, t80, 0; (t165 * t83 - (-t91 + (-t101 * t162 + t91) * t93) * t161) * t82, 0, (t110 * t137 * t83 - (t149 * t81 + (t112 * t92 - t146 * t91) * t137) * t164) * t82, 0, (t100 * t83 - (t92 * t104 + t149 * t80 + t91 * t98) * t164) * t82, 0; ((t98 * t136 - t140 * t146) * t87 - (t136 * t146 + t98 * t140) * t163) * t86, 0, ((t136 * t159 - t140 * t144) * t87 - (t136 * t144 + t140 * t159) * t163) * t86, 0, (-t136 * t87 + t140 * t163) * t99 * t86, t150 * t86;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end