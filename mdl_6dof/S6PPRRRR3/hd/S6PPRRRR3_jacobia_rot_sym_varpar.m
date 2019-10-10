% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR3
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
%   Wie in S6PPRRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (90->0), div. (5->0), fcn. (119->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (558->35), mult. (1619->83), div. (35->9), fcn. (2209->17), ass. (0->50)
	t70 = sin(pkin(14));
	t75 = cos(pkin(14));
	t76 = cos(pkin(13));
	t71 = sin(pkin(13));
	t79 = cos(pkin(6));
	t93 = t71 * t79;
	t69 = -t70 * t93 + t76 * t75;
	t81 = sin(qJ(3));
	t83 = cos(qJ(3));
	t68 = -t76 * t70 - t75 * t93;
	t78 = cos(pkin(7));
	t73 = sin(pkin(7));
	t74 = sin(pkin(6));
	t92 = t73 * t74;
	t84 = t68 * t78 + t71 * t92;
	t61 = t69 * t83 + t84 * t81;
	t80 = sin(qJ(4));
	t95 = t61 * t80;
	t82 = cos(qJ(4));
	t94 = t61 * t82;
	t91 = t73 * t79;
	t90 = t74 * t78;
	t89 = t75 * t78;
	t88 = t76 * t79;
	t60 = -t69 * t81 + t84 * t83;
	t65 = -t68 * t73 + t71 * t90;
	t72 = sin(pkin(8));
	t77 = cos(pkin(8));
	t86 = t60 * t77 + t65 * t72;
	t51 = t86 * t80 + t94;
	t49 = 0.1e1 / t51 ^ 2;
	t50 = -t86 * t82 + t95;
	t87 = t50 ^ 2 * t49 + 0.1e1;
	t66 = -t71 * t70 + t75 * t88;
	t85 = -t66 * t78 + t76 * t92;
	t67 = t70 * t88 + t71 * t75;
	t64 = -t81 * t91 + (-t70 * t83 - t81 * t89) * t74;
	t62 = -(t83 * t91 + (-t70 * t81 + t83 * t89) * t74) * t72 + (-t75 * t92 + t79 * t78) * t77;
	t59 = -t67 * t83 + t85 * t81;
	t58 = 0.1e1 / t62 ^ 2;
	t57 = -t60 * t72 + t65 * t77;
	t56 = (-t67 * t81 - t85 * t83) * t72 - (-t66 * t73 - t76 * t90) * t77;
	t55 = atan2(t56, t62);
	t53 = cos(t55);
	t52 = sin(t55);
	t48 = 0.1e1 / t87;
	t47 = t52 * t56 + t53 * t62;
	t46 = 0.1e1 / t47 ^ 2;
	t44 = (t59 / t62 + t64 * t56 * t58) * t72 / (t56 ^ 2 * t58 + 0.1e1);
	t1 = [0, 0, t44, 0, 0, 0; 0, 0, (t61 * t72 / t47 - ((t52 * t59 - t53 * t64) * t72 + (-t52 * t62 + t53 * t56) * t44) * t57 * t46) / (t57 ^ 2 * t46 + 0.1e1), 0, 0, 0; 0, 0, ((t60 * t80 + t77 * t94) / t51 - (t60 * t82 - t77 * t95) * t50 * t49) * t48, t87 * t48, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (1774->52), mult. (5219->122), div. (65->9), fcn. (7059->19), ass. (0->74)
	t108 = sin(pkin(8));
	t109 = sin(pkin(7));
	t113 = cos(pkin(8));
	t114 = cos(pkin(7));
	t118 = sin(qJ(3));
	t121 = cos(qJ(3));
	t106 = sin(pkin(14));
	t107 = sin(pkin(13));
	t111 = cos(pkin(14));
	t112 = cos(pkin(13));
	t115 = cos(pkin(6));
	t133 = t112 * t115;
	t127 = -t106 * t107 + t111 * t133;
	t110 = sin(pkin(6));
	t134 = t112 * t110;
	t125 = -t109 * t134 + t114 * t127;
	t126 = t106 * t133 + t107 * t111;
	t123 = -t118 * t126 + t121 * t125;
	t145 = t123 * t113 + (-t109 * t127 - t114 * t134) * t108;
	t117 = sin(qJ(4));
	t120 = cos(qJ(4));
	t95 = t118 * t125 + t121 * t126;
	t80 = t95 * t117 - t120 * t145;
	t135 = t111 * t114;
	t137 = t109 * t115;
	t101 = t118 * t137 + (t106 * t121 + t118 * t135) * t110;
	t100 = t121 * t137 + (-t106 * t118 + t121 * t135) * t110;
	t136 = t110 * t109;
	t129 = t100 * t113 + (-t111 * t136 + t114 * t115) * t108;
	t89 = t101 * t117 - t120 * t129;
	t79 = atan2(-t80, t89);
	t76 = sin(t79);
	t77 = cos(t79);
	t70 = -t76 * t80 + t77 * t89;
	t69 = 0.1e1 / t70 ^ 2;
	t132 = t113 * t120;
	t138 = t107 * t115;
	t104 = -t106 * t112 - t111 * t138;
	t102 = t107 * t110 * t114 - t104 * t109;
	t139 = t102 * t108;
	t105 = -t106 * t138 + t111 * t112;
	t128 = t104 * t114 + t107 * t136;
	t97 = t105 * t121 + t118 * t128;
	t140 = t97 * t117;
	t96 = -t105 * t118 + t121 * t128;
	t83 = -t120 * t139 - t132 * t96 + t140;
	t144 = t69 * t83;
	t116 = sin(qJ(5));
	t119 = cos(qJ(5));
	t84 = t97 * t120 + (t113 * t96 + t139) * t117;
	t91 = t102 * t113 - t108 * t96;
	t75 = t116 * t91 + t119 * t84;
	t73 = 0.1e1 / t75 ^ 2;
	t74 = t116 * t84 - t119 * t91;
	t143 = t73 * t74;
	t88 = 0.1e1 / t89 ^ 2;
	t142 = t80 * t88;
	t141 = t108 * t97;
	t131 = t73 * t74 ^ 2 + 0.1e1;
	t130 = -t76 * t89 - t77 * t80;
	t92 = t100 * t117 + t101 * t132;
	t90 = t101 * t120 + t117 * t129;
	t87 = 0.1e1 / t89;
	t86 = -t113 * t140 + t120 * t96;
	t85 = t117 * t123 + t132 * t95;
	t82 = t117 * t145 + t95 * t120;
	t78 = 0.1e1 / (t80 ^ 2 * t88 + 0.1e1);
	t72 = 0.1e1 / t75;
	t71 = 0.1e1 / t131;
	t68 = 0.1e1 / t70;
	t67 = 0.1e1 / (t69 * t83 ^ 2 + 0.1e1);
	t66 = (t142 * t92 - t85 * t87) * t78;
	t65 = (t142 * t90 - t82 * t87) * t78;
	t1 = [0, 0, t66, t65, 0, 0; 0, 0, ((t96 * t117 + t132 * t97) * t68 - (t130 * t66 - t76 * t85 + t77 * t92) * t144) * t67, (t84 * t68 - (t130 * t65 - t76 * t82 + t77 * t90) * t144) * t67, 0, 0; 0, 0, ((t116 * t86 - t119 * t141) * t72 - (t116 * t141 + t119 * t86) * t143) * t71, (-t116 * t72 + t119 * t143) * t83 * t71, t131 * t71, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:22
	% DurationCPUTime: 0.99s
	% Computational Cost: add. (4253->74), mult. (12284->171), div. (95->9), fcn. (16637->21), ass. (0->89)
	t147 = sin(pkin(14));
	t148 = sin(pkin(13));
	t152 = cos(pkin(14));
	t153 = cos(pkin(13));
	t156 = cos(pkin(6));
	t179 = t153 * t156;
	t142 = -t147 * t148 + t152 * t179;
	t143 = t147 * t179 + t148 * t152;
	t160 = sin(qJ(3));
	t164 = cos(qJ(3));
	t150 = sin(pkin(7));
	t151 = sin(pkin(6));
	t182 = t150 * t151;
	t175 = t153 * t182;
	t155 = cos(pkin(7));
	t176 = t155 * t160;
	t134 = t142 * t176 + t143 * t164 - t160 * t175;
	t159 = sin(qJ(4));
	t163 = cos(qJ(4));
	t149 = sin(pkin(8));
	t154 = cos(pkin(8));
	t167 = -t143 * t160 + (t142 * t155 - t175) * t164;
	t180 = t151 * t155;
	t172 = -t142 * t150 - t153 * t180;
	t165 = t149 * t172 + t154 * t167;
	t121 = t134 * t163 + t159 * t165;
	t158 = sin(qJ(5));
	t162 = cos(qJ(5));
	t166 = -t149 * t167 + t154 * t172;
	t111 = t121 * t158 - t162 * t166;
	t181 = t150 * t156;
	t140 = t160 * t181 + (t147 * t164 + t152 * t176) * t151;
	t139 = t164 * t181 + (t152 * t155 * t164 - t147 * t160) * t151;
	t141 = -t152 * t182 + t155 * t156;
	t173 = t139 * t154 + t141 * t149;
	t132 = t140 * t163 + t159 * t173;
	t137 = -t139 * t149 + t141 * t154;
	t124 = t132 * t158 - t137 * t162;
	t110 = atan2(-t111, t124);
	t107 = sin(t110);
	t108 = cos(t110);
	t101 = -t107 * t111 + t108 * t124;
	t100 = 0.1e1 / t101 ^ 2;
	t184 = t148 * t156;
	t145 = -t147 * t184 + t152 * t153;
	t144 = -t147 * t153 - t152 * t184;
	t170 = t144 * t155 + t148 * t182;
	t135 = -t145 * t160 + t164 * t170;
	t136 = t145 * t164 + t160 * t170;
	t171 = -t144 * t150 + t148 * t180;
	t169 = t171 * t149;
	t123 = t136 * t163 + (t135 * t154 + t169) * t159;
	t168 = -t135 * t149 + t154 * t171;
	t114 = t123 * t158 - t162 * t168;
	t188 = t100 * t114;
	t115 = t123 * t162 + t158 * t168;
	t177 = t154 * t163;
	t122 = -t135 * t177 + t136 * t159 - t163 * t169;
	t157 = sin(qJ(6));
	t161 = cos(qJ(6));
	t106 = t115 * t161 + t122 * t157;
	t104 = 0.1e1 / t106 ^ 2;
	t105 = t115 * t157 - t122 * t161;
	t187 = t104 * t105;
	t119 = 0.1e1 / t124 ^ 2;
	t186 = t111 * t119;
	t185 = t122 * t162;
	t183 = t149 * t162;
	t178 = t154 * t159;
	t174 = t104 * t105 ^ 2 + 0.1e1;
	t131 = -t140 * t159 + t163 * t173;
	t128 = (t139 * t163 - t140 * t178) * t158 - t140 * t183;
	t127 = t135 * t163 - t136 * t178;
	t126 = t135 * t159 + t136 * t177;
	t125 = t132 * t162 + t137 * t158;
	t120 = -t134 * t159 + t163 * t165;
	t118 = 0.1e1 / t124;
	t117 = t136 * t149 * t158 + t127 * t162;
	t116 = (-t134 * t178 + t163 * t167) * t158 - t134 * t183;
	t113 = t121 * t162 + t158 * t166;
	t109 = 0.1e1 / (t111 ^ 2 * t119 + 0.1e1);
	t103 = 0.1e1 / t106;
	t102 = 0.1e1 / t174;
	t99 = 0.1e1 / t101;
	t98 = 0.1e1 / (t100 * t114 ^ 2 + 0.1e1);
	t97 = (-t118 * t120 + t131 * t186) * t158 * t109;
	t96 = (-t116 * t118 + t128 * t186) * t109;
	t95 = (-t113 * t118 + t125 * t186) * t109;
	t1 = [0, 0, t96, t97, t95, 0; 0, 0, ((t127 * t158 - t136 * t183) * t99 - ((-t111 * t96 + t128) * t108 + (-t124 * t96 - t116) * t107) * t188) * t98, (-t122 * t158 * t99 - ((-t111 * t97 + t131 * t158) * t108 + (-t120 * t158 - t124 * t97) * t107) * t188) * t98, (t115 * t99 - ((-t111 * t95 + t125) * t108 + (-t124 * t95 - t113) * t107) * t188) * t98, 0; 0, 0, ((t117 * t157 - t126 * t161) * t103 - (t117 * t161 + t126 * t157) * t187) * t102, ((-t123 * t161 - t157 * t185) * t103 - (t123 * t157 - t161 * t185) * t187) * t102, (-t103 * t157 + t161 * t187) * t114 * t102, t174 * t102;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end