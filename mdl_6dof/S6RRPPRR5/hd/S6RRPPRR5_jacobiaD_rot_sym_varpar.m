% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR5
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
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRPPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:43
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(6));
	t93 = t99 ^ 2;
	t100 = cos(pkin(6));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (1133->73), mult. (3290->173), div. (633->14), fcn. (4296->9), ass. (0->78)
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t157 = cos(pkin(6));
	t141 = t126 * t157;
	t134 = -t123 * t141 - t124 * t125;
	t170 = t134 * qJD(1);
	t139 = t125 * t141;
	t153 = t123 * t124;
	t105 = -t139 + t153;
	t122 = sin(pkin(6));
	t116 = 0.1e1 / t122;
	t119 = 0.1e1 / t125;
	t144 = t105 * t116 * t119;
	t154 = t122 * t125;
	t93 = atan2(-t105, -t154);
	t91 = sin(t93);
	t92 = cos(t93);
	t100 = t105 ^ 2;
	t117 = 0.1e1 / t122 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t99 = t100 * t117 * t120 + 0.1e1;
	t96 = 0.1e1 / t99;
	t169 = (t92 * t144 - t91) * t96 + t91;
	t86 = -t105 * t91 - t92 * t154;
	t83 = 0.1e1 / t86;
	t142 = t124 * t157;
	t140 = t123 * t142;
	t152 = t126 * t125;
	t109 = -t140 + t152;
	t102 = 0.1e1 / t109;
	t103 = 0.1e1 / t109 ^ 2;
	t84 = 0.1e1 / t86 ^ 2;
	t150 = qJD(2) * t123;
	t159 = t125 * t91;
	t165 = t105 * t92;
	t143 = t120 * t150;
	t162 = t116 * t96;
	t135 = -t126 * t123 - t125 * t142;
	t89 = -t135 * qJD(1) - t134 * qJD(2);
	t76 = (t105 * t143 + t119 * t89) * t162;
	t73 = -t76 * t165 - t91 * t89 + (t92 * t150 + t76 * t159) * t122;
	t168 = t73 * t83 * t84;
	t155 = t120 * t123;
	t136 = t105 * t155 - t119 * t134;
	t77 = t136 * t162;
	t167 = t76 * t77;
	t88 = t135 * qJD(2) + t170;
	t166 = t102 * t103 * t88;
	t164 = t135 * t84;
	t163 = t135 * t92;
	t161 = t119 * t96;
	t115 = t122 ^ 2;
	t118 = t124 ^ 2;
	t98 = t103 * t115 * t118 + 0.1e1;
	t94 = 0.1e1 / t98;
	t160 = t124 * t94;
	t158 = t91 * t135;
	t156 = t103 * t124;
	t151 = qJD(1) * t126;
	t101 = t135 ^ 2;
	t80 = t101 * t84 + 0.1e1;
	t138 = qJD(2) * t157 + qJD(1);
	t87 = -qJD(1) * t139 - qJD(2) * t152 + t138 * t153;
	t149 = 0.2e1 * (-t101 * t168 + t87 * t164) / t80 ^ 2;
	t148 = 0.2e1 * t168;
	t147 = 0.2e1 * (-t118 * t166 + t151 * t156) * t115 / t98 ^ 2;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t100 * t121 * t150 + t105 * t120 * t89) * t117 / t99 ^ 2;
	t145 = 0.2e1 * t166;
	t133 = t119 * t146 + t96 * t143;
	t90 = -qJD(1) * t140 - t124 * t150 + t138 * t152;
	t78 = 0.1e1 / t80;
	t75 = t169 * t135;
	t74 = -t77 * t165 + t91 * t134 + (t123 * t92 + t77 * t159) * t122;
	t72 = (t136 * t146 + (t89 * t155 + t119 * t90 + (-t134 * t155 + (0.2e1 * t121 * t123 ^ 2 + t119) * t105) * qJD(2)) * t96) * t116;
	t1 = [(-t133 * t135 - t87 * t161) * t116, t72, 0, 0, 0, 0; t105 * t83 * t149 + (-t89 * t83 + (t105 * t73 + t75 * t87) * t84) * t78 - (t75 * t148 * t78 + (t75 * t149 + ((t76 * t96 * t144 + t146) * t158 + ((t96 - 0.1e1) * t76 + (-t133 * t105 - t89 * t161) * t116) * t163 - t169 * t87) * t78) * t84) * t135, (-t109 * t83 - t74 * t164) * t149 + (-t74 * t135 * t148 + t88 * t83 + (-t109 * t73 + t74 * t87 + (t105 * t167 - t90) * t158 + (-t105 * t72 + t134 * t76 - t77 * t89) * t163) * t84 + ((-qJD(2) * t77 - t76) * t91 * t123 + (t72 * t91 + (qJD(2) + t167) * t92) * t125) * t122 * t164) * t78, 0, 0, 0, 0; ((t102 * t126 - t134 * t156) * t147 + ((qJD(1) * t102 - t134 * t145) * t124 + (-t124 * t90 + (t88 + t170) * t126) * t103) * t94) * t122, (-t135 * t145 * t160 + (t87 * t160 - (t124 * t147 - t94 * t151) * t135) * t103) * t122, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (351->41), mult. (876->109), div. (130->12), fcn. (1072->9), ass. (0->57)
	t104 = cos(pkin(6));
	t103 = sin(pkin(6));
	t108 = cos(qJ(1));
	t128 = t108 * t103;
	t91 = atan2(-t128, -t104);
	t89 = sin(t91);
	t90 = cos(t91);
	t75 = -t90 * t104 - t89 * t128;
	t72 = 0.1e1 / t75;
	t107 = cos(qJ(2));
	t126 = t108 * t107;
	t105 = sin(qJ(2));
	t106 = sin(qJ(1));
	t130 = t106 * t105;
	t117 = t104 * t130 - t126;
	t82 = 0.1e1 / t117;
	t98 = 0.1e1 / t104;
	t99 = 0.1e1 / t104 ^ 2;
	t73 = 0.1e1 / t75 ^ 2;
	t83 = 0.1e1 / t117 ^ 2;
	t127 = t108 * t105;
	t129 = t106 * t107;
	t119 = t104 * t129 + t127;
	t81 = t119 ^ 2;
	t138 = t81 * t83;
	t84 = t82 * t83;
	t137 = t81 * t84;
	t118 = t104 * t127 + t129;
	t136 = t118 * t119;
	t102 = t108 ^ 2;
	t97 = t103 ^ 2;
	t94 = t102 * t97 * t99 + 0.1e1;
	t93 = 0.1e1 / t94 ^ 2;
	t135 = t93 * t103 * t97;
	t134 = t97 * t98;
	t101 = t106 ^ 2;
	t133 = t101 * t73;
	t132 = t106 * t73;
	t131 = t108 * t73;
	t125 = qJD(1) * t108;
	t85 = -t104 * t126 + t130;
	t76 = t85 * qJD(1) + t117 * qJD(2);
	t123 = t119 * t83 * t76;
	t77 = t118 * qJD(1) + t119 * qJD(2);
	t80 = 0.1e1 + t138;
	t124 = 0.2e1 * (-t77 * t137 - t123) / t80 ^ 2;
	t92 = 0.1e1 / t94;
	t122 = t92 * t134;
	t121 = t90 * t122;
	t120 = (-t92 + 0.1e1) * t89 * t103;
	t68 = (t108 * t121 + t120) * t106;
	t100 = t98 * t99;
	t78 = 0.1e1 / t80;
	t74 = t72 * t73;
	t71 = t97 * t133 + 0.1e1;
	t67 = qJD(1) * t68;
	t1 = [(-0.2e1 * t100 * t101 * t135 - t103 * t92 * t98) * t125, 0, 0, 0, 0, 0; (0.2e1 * (t108 * t72 - t68 * t132) / t71 ^ 2 * (-t101 * t67 * t74 + t125 * t132) * t97 + ((-0.2e1 * t106 * t68 * t74 + t131) * t67 + (t68 * t131 + (t102 * t73 * t121 + t72 + t120 * t131 + (-t108 * t89 * t99 * t135 + (-0.2e1 * t122 + (0.2e1 * t100 * t102 * t97 ^ 2 + t134) * t93) * t90) * t133) * t106) * qJD(1)) / t71) * t103, 0, 0, 0, 0, 0; (t83 * t136 + t82 * t85) * t124 + (-(t119 * qJD(1) + t118 * qJD(2)) * t82 + 0.2e1 * t84 * t77 * t136 + (t85 * t77 + (t117 * qJD(1) + t85 * qJD(2)) * t119 + t118 * t76) * t83) * t78, (t117 * t82 + t138) * t124 + (0.2e1 * t123 + (t117 * t83 + 0.2e1 * t137 - t82) * t77) * t78, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (1334->89), mult. (4303->202), div. (668->14), fcn. (5516->11), ass. (0->93)
	t171 = sin(qJ(2));
	t172 = sin(qJ(1));
	t174 = cos(qJ(2));
	t175 = cos(qJ(1));
	t222 = cos(pkin(6));
	t191 = t175 * t222;
	t154 = t171 * t191 + t172 * t174;
	t192 = t172 * t222;
	t155 = t175 * t171 + t174 * t192;
	t135 = t155 * qJD(1) + t154 * qJD(2);
	t188 = t174 * t191;
	t207 = t172 * t171;
	t153 = -t188 + t207;
	t151 = t153 ^ 2;
	t169 = sin(pkin(6));
	t165 = 0.1e1 / t169 ^ 2;
	t167 = 0.1e1 / t174 ^ 2;
	t150 = t151 * t165 * t167 + 0.1e1;
	t166 = 0.1e1 / t174;
	t168 = t166 * t167;
	t204 = qJD(2) * t171;
	t217 = (t135 * t153 * t167 + t151 * t168 * t204) * t165 / t150 ^ 2;
	t225 = -0.2e1 * t217;
	t164 = 0.1e1 / t169;
	t224 = t153 * t164;
	t209 = t169 * t174;
	t149 = atan2(t153, t209);
	t145 = sin(t149);
	t146 = cos(t149);
	t147 = 0.1e1 / t150;
	t196 = t166 * t224;
	t223 = (t146 * t196 - t145) * t147 + t145;
	t129 = t145 * t153 + t146 * t209;
	t126 = 0.1e1 / t129;
	t189 = t171 * t192;
	t206 = t175 * t174;
	t157 = -t189 + t206;
	t170 = sin(qJ(5));
	t173 = cos(qJ(5));
	t210 = t169 * t172;
	t144 = t157 * t173 - t170 * t210;
	t138 = 0.1e1 / t144;
	t127 = 0.1e1 / t129 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t152 = t155 ^ 2;
	t122 = t152 * t127 + 0.1e1;
	t187 = qJD(2) * t222 + qJD(1);
	t203 = qJD(2) * t174;
	t133 = -qJD(1) * t188 - t175 * t203 + t187 * t207;
	t215 = t133 * t127;
	t193 = t167 * t204;
	t182 = (t135 * t166 + t153 * t193) * t164;
	t118 = t147 * t182;
	t184 = -t145 * t209 + t146 * t153;
	t197 = t146 * t169 * t171;
	t114 = -qJD(2) * t197 + t184 * t118 + t145 * t135;
	t220 = t114 * t126 * t127;
	t221 = (-t152 * t220 - t155 * t215) / t122 ^ 2;
	t211 = t167 * t171;
	t183 = t153 * t211 + t154 * t166;
	t119 = t183 * t164 * t147;
	t115 = t184 * t119 + t145 * t154 - t197;
	t219 = t115 * t155;
	t134 = t154 * qJD(1) + t155 * qJD(2);
	t205 = qJD(1) * t169;
	t194 = t175 * t205;
	t124 = t144 * qJD(5) - t134 * t170 + t173 * t194;
	t143 = t157 * t170 + t173 * t210;
	t137 = t143 ^ 2;
	t132 = t137 * t139 + 0.1e1;
	t214 = t139 * t143;
	t202 = qJD(5) * t143;
	t125 = -t134 * t173 - t170 * t194 - t202;
	t216 = t125 * t138 * t139;
	t218 = (t124 * t214 - t137 * t216) / t132 ^ 2;
	t213 = t145 * t155;
	t212 = t146 * t155;
	t208 = t169 * t175;
	t201 = -0.2e1 * t221;
	t200 = -0.2e1 * t220;
	t199 = 0.2e1 * t218;
	t198 = t143 * t216;
	t195 = t172 * t205;
	t190 = t166 * t225;
	t185 = t170 * t138 - t173 * t214;
	t141 = -t154 * t170 + t173 * t208;
	t142 = -t154 * t173 - t170 * t208;
	t136 = -qJD(1) * t189 - t172 * t204 + t187 * t206;
	t130 = 0.1e1 / t132;
	t120 = 0.1e1 / t122;
	t117 = t223 * t155;
	t113 = (t183 * t225 + (t135 * t211 + t136 * t166 + (t154 * t211 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t147) * t164;
	t1 = [(t155 * t190 + (-t133 * t166 + t155 * t193) * t147) * t164, t113, 0, 0, 0, 0; t153 * t126 * t201 + (t135 * t126 + (-t114 * t153 - t117 * t133) * t127) * t120 + ((t117 * t200 - t223 * t215) * t120 + (t117 * t201 + ((-t118 * t147 * t196 + 0.2e1 * t217) * t213 + (t190 * t224 + t118 + (-t118 + t182) * t147) * t212) * t120) * t127) * t155, 0.2e1 * (t126 * t157 - t127 * t219) * t221 + (t200 * t219 + t134 * t126 + (t157 * t114 - t115 * t133 + (-t169 * t203 + t113 * t153 + t119 * t135 + (-t119 * t209 + t154) * t118) * t212 + (-t118 * t119 * t153 + t136 + (-t113 * t174 + (qJD(2) * t119 + t118) * t171) * t169) * t213) * t127) * t120, 0, 0, 0, 0; (-t138 * t141 + t142 * t214) * t199 + ((qJD(5) * t142 - t136 * t170 - t173 * t195) * t138 + 0.2e1 * t142 * t198 + (-t141 * t125 - (-qJD(5) * t141 - t136 * t173 + t170 * t195) * t143 - t142 * t124) * t139) * t130, t185 * t155 * t199 + (t185 * t133 + ((-qJD(5) * t138 - 0.2e1 * t198) * t173 + (t124 * t173 + (t125 - t202) * t170) * t139) * t155) * t130, 0, 0, -0.2e1 * t218 + 0.2e1 * (t124 * t139 * t130 + (-t130 * t216 - t139 * t218) * t143) * t143, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:09
	% DurationCPUTime: 2.05s
	% Computational Cost: add. (4522->145), mult. (13478->293), div. (726->12), fcn. (17045->13), ass. (0->124)
	t251 = cos(pkin(6));
	t254 = sin(qJ(2));
	t257 = cos(qJ(5));
	t250 = sin(pkin(6));
	t253 = sin(qJ(5));
	t306 = t250 * t253;
	t237 = t251 * t257 + t254 * t306;
	t255 = sin(qJ(1));
	t258 = cos(qJ(2));
	t300 = t255 * t258;
	t259 = cos(qJ(1));
	t301 = t254 * t259;
	t240 = t251 * t301 + t300;
	t303 = t250 * t259;
	t272 = -t240 * t253 + t257 * t303;
	t215 = atan2(t272, t237);
	t210 = sin(t215);
	t211 = cos(t215);
	t193 = t210 * t272 + t211 * t237;
	t191 = 0.1e1 / t193 ^ 2;
	t302 = t254 * t255;
	t284 = t251 * t302;
	t299 = t258 * t259;
	t270 = t284 - t299;
	t305 = t250 * t257;
	t231 = -t253 * t270 + t255 * t305;
	t223 = t231 ^ 2;
	t189 = t191 * t223 + 0.1e1;
	t271 = t251 * t300 + t301;
	t219 = t240 * qJD(1) + t271 * qJD(2);
	t298 = qJD(1) * t250;
	t282 = t259 * t298;
	t285 = t255 * t306;
	t197 = -t219 * t253 - qJD(5) * t285 + (-qJD(5) * t270 + t282) * t257;
	t318 = t191 * t231;
	t222 = t272 ^ 2;
	t235 = 0.1e1 / t237 ^ 2;
	t214 = t222 * t235 + 0.1e1;
	t212 = 0.1e1 / t214;
	t297 = qJD(2) * t254;
	t221 = -qJD(1) * t284 - t255 * t297 + (qJD(2) * t251 + qJD(1)) * t299;
	t228 = t240 * t257 + t253 * t303;
	t283 = t255 * t298;
	t199 = t228 * qJD(5) + t221 * t253 + t257 * t283;
	t238 = -t251 * t253 + t254 * t305;
	t304 = t250 * t258;
	t281 = qJD(2) * t304;
	t224 = t238 * qJD(5) + t253 * t281;
	t234 = 0.1e1 / t237;
	t310 = t272 * t235;
	t275 = -t199 * t234 - t224 * t310;
	t181 = t275 * t212;
	t276 = -t210 * t237 + t211 * t272;
	t176 = t276 * t181 - t210 * t199 + t211 * t224;
	t190 = 0.1e1 / t193;
	t192 = t190 * t191;
	t323 = t176 * t192;
	t294 = 0.2e1 * (t197 * t318 - t223 * t323) / t189 ^ 2;
	t328 = t224 * t235;
	t239 = -t251 * t299 + t302;
	t269 = t234 * t239 - t304 * t310;
	t327 = t253 * t269;
	t200 = t272 * qJD(5) + t221 * t257 - t253 * t283;
	t232 = -t257 * t270 - t285;
	t252 = sin(qJ(6));
	t256 = cos(qJ(6));
	t209 = t232 * t256 - t252 * t271;
	t203 = 0.1e1 / t209;
	t204 = 0.1e1 / t209 ^ 2;
	t326 = 0.2e1 * t272;
	t325 = 0.2e1 * t231;
	t198 = -t231 * qJD(5) - t219 * t257 - t253 * t282;
	t218 = t239 * qJD(1) + t270 * qJD(2);
	t185 = t209 * qJD(6) + t198 * t252 - t218 * t256;
	t208 = t232 * t252 + t256 * t271;
	t202 = t208 ^ 2;
	t196 = t202 * t204 + 0.1e1;
	t316 = t204 * t208;
	t295 = qJD(6) * t208;
	t186 = t198 * t256 + t218 * t252 - t295;
	t320 = t186 * t203 * t204;
	t322 = (t185 * t316 - t202 * t320) / t196 ^ 2;
	t312 = t234 * t328;
	t321 = (-t199 * t310 - t222 * t312) / t214 ^ 2;
	t319 = t191 * t197;
	t317 = t203 * t252;
	t315 = t208 * t256;
	t314 = t210 * t231;
	t313 = t211 * t231;
	t311 = t272 * t234;
	t308 = t271 * t253;
	t307 = t271 * t257;
	t296 = qJD(5) * t257;
	t293 = -0.2e1 * t322;
	t292 = 0.2e1 * t322;
	t291 = -0.2e1 * t321;
	t290 = t192 * t325;
	t289 = t234 * t321;
	t288 = t208 * t320;
	t287 = t191 * t314;
	t286 = t191 * t313;
	t280 = 0.2e1 * t288;
	t279 = t312 * t326;
	t277 = -qJD(6) * t307 - t219;
	t207 = -t228 * t256 + t239 * t252;
	t206 = -t228 * t252 - t239 * t256;
	t274 = t204 * t315 - t317;
	t273 = -t228 * t234 - t238 * t310;
	t267 = -t210 + (-t211 * t311 + t210) * t212;
	t266 = qJD(5) * t308 + qJD(6) * t270 + t218 * t257;
	t225 = -t237 * qJD(5) + t257 * t281;
	t220 = t271 * qJD(1) + t240 * qJD(2);
	t217 = t252 * t270 - t256 * t307;
	t216 = -t252 * t307 - t256 * t270;
	t194 = 0.1e1 / t196;
	t187 = 0.1e1 / t189;
	t184 = t212 * t327;
	t183 = t273 * t212;
	t180 = t267 * t231;
	t178 = (t210 * t239 + t211 * t304) * t253 + t276 * t184;
	t177 = t276 * t183 - t210 * t228 + t211 * t238;
	t175 = t273 * t291 + (t238 * t279 - t200 * t234 + (t199 * t238 + t224 * t228 - t225 * t272) * t235) * t212;
	t173 = t291 * t327 + (t269 * t296 + (t279 * t304 + t220 * t234 + (-t224 * t239 + (t199 * t258 + t272 * t297) * t250) * t235) * t253) * t212;
	t1 = [t289 * t325 + (-t197 * t234 + t231 * t328) * t212, t173, 0, 0, t175, 0; -t272 * t190 * t294 + (-t199 * t190 + (-t176 * t272 - t180 * t197) * t191) * t187 + (t180 * t191 * t294 + (0.2e1 * t180 * t323 - (t181 * t212 * t311 + t291) * t287 - (t289 * t326 - t181 + (t181 - t275) * t212) * t286 - t267 * t319) * t187) * t231, (t178 * t318 + t190 * t308) * t294 + (-t178 * t319 + (t218 * t253 - t271 * t296) * t190 + (t178 * t290 + t191 * t308) * t176 - (t173 * t272 - t184 * t199 + (-t253 * t297 + t258 * t296) * t250 + (-t184 * t237 + t239 * t253) * t181) * t286 - (t239 * t296 - t173 * t237 - t184 * t224 + t220 * t253 + (-t184 * t272 - t253 * t304) * t181) * t287) * t187, 0, 0, (t177 * t318 - t190 * t232) * t294 + (t177 * t176 * t290 + t198 * t190 + (-t232 * t176 - t177 * t197 - (t175 * t272 - t183 * t199 + t225 + (-t183 * t237 - t228) * t181) * t313 - (-t175 * t237 - t183 * t224 - t200 + (-t183 * t272 - t238) * t181) * t314) * t191) * t187, 0; (-t203 * t206 + t207 * t316) * t292 + ((t207 * qJD(6) - t200 * t252 - t220 * t256) * t203 + t207 * t280 + (-t206 * t186 - (-t206 * qJD(6) - t200 * t256 + t220 * t252) * t208 - t207 * t185) * t204) * t194, (-t203 * t216 + t217 * t316) * t292 + (t217 * t280 + t277 * t203 * t256 + t266 * t317 + (t277 * t208 * t252 - t217 * t185 - t216 * t186 - t266 * t315) * t204) * t194, 0, 0, t274 * t231 * t293 + (t274 * t197 + ((-qJD(6) * t203 - 0.2e1 * t288) * t256 + (t185 * t256 + (t186 - t295) * t252) * t204) * t231) * t194, t293 + 0.2e1 * (t185 * t194 * t204 + (-t194 * t320 - t204 * t322) * t208) * t208;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end