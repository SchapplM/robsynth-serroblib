% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR7
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
%   Wie in S6RRPPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:41
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
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.84s
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
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (351->43), mult. (876->111), div. (130->12), fcn. (1072->9), ass. (0->58)
	t106 = cos(pkin(6));
	t101 = 0.1e1 / t106 ^ 2;
	t110 = cos(qJ(1));
	t104 = t110 ^ 2;
	t105 = sin(pkin(6));
	t99 = t105 ^ 2;
	t95 = t101 * t104 * t99 + 0.1e1;
	t94 = 0.1e1 / t95 ^ 2;
	t139 = t101 * t94;
	t128 = t110 * t105;
	t92 = atan2(-t128, -t106);
	t90 = sin(t92);
	t91 = cos(t92);
	t76 = -t106 * t91 - t128 * t90;
	t73 = 0.1e1 / t76;
	t107 = sin(qJ(2));
	t127 = t110 * t107;
	t108 = sin(qJ(1));
	t109 = cos(qJ(2));
	t129 = t108 * t109;
	t88 = t106 * t129 + t127;
	t82 = 0.1e1 / t88;
	t100 = 0.1e1 / t106;
	t74 = 0.1e1 / t76 ^ 2;
	t83 = 0.1e1 / t88 ^ 2;
	t84 = t82 * t83;
	t126 = t110 * t109;
	t130 = t108 * t107;
	t119 = t106 * t130 - t126;
	t85 = t119 ^ 2;
	t138 = t84 * t85;
	t137 = t85 * t83;
	t122 = t106 * t126;
	t86 = t122 - t130;
	t136 = t86 * t119;
	t135 = t100 * t99;
	t134 = t100 * t139;
	t103 = t108 ^ 2;
	t133 = t103 * t74;
	t132 = t108 * t74;
	t131 = t110 * t74;
	t125 = qJD(1) * t110;
	t87 = -t106 * t127 - t129;
	t78 = qJD(1) * t87 - qJD(2) * t88;
	t123 = t119 * t83 * t78;
	t77 = -qJD(1) * t122 - qJD(2) * t126 + (qJD(2) * t106 + qJD(1)) * t130;
	t81 = 0.1e1 + t137;
	t124 = 0.2e1 * (t138 * t77 - t123) / t81 ^ 2;
	t93 = 0.1e1 / t95;
	t121 = t91 * t93 * t135;
	t120 = (-t93 + 0.1e1) * t90 * t105;
	t69 = (t110 * t121 + t120) * t108;
	t98 = t105 * t99;
	t79 = 0.1e1 / t81;
	t75 = t73 * t74;
	t72 = t133 * t99 + 0.1e1;
	t68 = qJD(1) * t69;
	t1 = [(-t100 * t105 * t93 - 0.2e1 * t103 * t134 * t98) * t125, 0, 0, 0, 0, 0; (0.2e1 * (t110 * t73 - t132 * t69) / t72 ^ 2 * (-t103 * t68 * t75 + t125 * t132) * t99 + ((-0.2e1 * t108 * t69 * t75 + t131) * t68 + (t69 * t131 + (t104 * t74 * t121 + t73 + t120 * t131 + (-t110 * t90 * t98 * t139 + (0.2e1 * t104 * t99 ^ 2 * t134 + (-0.2e1 * t93 + t94) * t135) * t91) * t133) * t108) * qJD(1)) / t72) * t105, 0, 0, 0, 0, 0; (-t136 * t83 - t82 * t87) * t124 + ((qJD(1) * t119 - qJD(2) * t86) * t82 + 0.2e1 * t84 * t77 * t136 + (t87 * t77 + (-qJD(1) * t88 + qJD(2) * t87) * t119 - t86 * t78) * t83) * t79, (t82 * t88 + t137) * t124 + (0.2e1 * t123 + (-t83 * t88 - 0.2e1 * t138 + t82) * t77) * t79, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
	t170 = sin(qJ(2));
	t171 = sin(qJ(1));
	t173 = cos(qJ(2));
	t174 = cos(qJ(1));
	t221 = cos(pkin(6));
	t190 = t174 * t221;
	t153 = t170 * t190 + t171 * t173;
	t168 = sin(pkin(6));
	t209 = t168 * t170;
	t148 = atan2(-t153, t209);
	t144 = sin(t148);
	t145 = cos(t148);
	t150 = t153 ^ 2;
	t164 = 0.1e1 / t168 ^ 2;
	t166 = 0.1e1 / t170 ^ 2;
	t149 = t150 * t164 * t166 + 0.1e1;
	t146 = 0.1e1 / t149;
	t163 = 0.1e1 / t168;
	t165 = 0.1e1 / t170;
	t195 = t153 * t163 * t165;
	t222 = t146 * (t145 * t195 + t144) - t144;
	t128 = -t144 * t153 + t145 * t209;
	t125 = 0.1e1 / t128;
	t191 = t171 * t221;
	t155 = t174 * t170 + t173 * t191;
	t169 = sin(qJ(5));
	t172 = cos(qJ(5));
	t208 = t168 * t171;
	t143 = t155 * t172 - t169 * t208;
	t137 = 0.1e1 / t143;
	t126 = 0.1e1 / t128 ^ 2;
	t138 = 0.1e1 / t143 ^ 2;
	t186 = qJD(2) * t221 + qJD(1);
	t188 = t170 * t191;
	t203 = qJD(2) * t170;
	t205 = t174 * t173;
	t135 = -qJD(1) * t188 - t171 * t203 + t186 * t205;
	t202 = qJD(2) * t173;
	t192 = t166 * t202;
	t181 = -t135 * t165 + t153 * t192;
	t211 = t146 * t163;
	t117 = t181 * t211;
	t183 = -t144 * t209 - t145 * t153;
	t196 = t145 * t168 * t173;
	t113 = qJD(2) * t196 + t117 * t183 - t144 * t135;
	t220 = t113 * t125 * t126;
	t187 = t173 * t190;
	t206 = t171 * t170;
	t152 = -t187 + t206;
	t210 = t166 * t173;
	t182 = t152 * t165 + t153 * t210;
	t118 = t182 * t211;
	t114 = t118 * t183 + t144 * t152 + t196;
	t156 = -t188 + t205;
	t219 = t114 * t156;
	t132 = -qJD(1) * t187 - t174 * t202 + t186 * t206;
	t204 = qJD(1) * t168;
	t193 = t174 * t204;
	t123 = qJD(5) * t143 - t132 * t169 + t172 * t193;
	t142 = t155 * t169 + t172 * t208;
	t136 = t142 ^ 2;
	t131 = t136 * t138 + 0.1e1;
	t214 = t138 * t142;
	t201 = qJD(5) * t142;
	t124 = -t132 * t172 - t169 * t193 - t201;
	t216 = t124 * t137 * t138;
	t218 = (t123 * t214 - t136 * t216) / t131 ^ 2;
	t167 = t165 * t166;
	t217 = (t135 * t153 * t166 - t150 * t167 * t202) * t164 / t149 ^ 2;
	t133 = qJD(1) * t153 + t155 * qJD(2);
	t215 = t133 * t126;
	t213 = t144 * t156;
	t212 = t145 * t156;
	t207 = t168 * t174;
	t151 = t156 ^ 2;
	t121 = t151 * t126 + 0.1e1;
	t200 = 0.2e1 * (-t151 * t220 - t156 * t215) / t121 ^ 2;
	t199 = 0.2e1 * t220;
	t198 = 0.2e1 * t218;
	t197 = -0.2e1 * t217;
	t194 = t171 * t204;
	t189 = 0.2e1 * t142 * t216;
	t184 = -t169 * t137 + t172 * t214;
	t140 = -t152 * t169 + t172 * t207;
	t141 = -t152 * t172 - t169 * t207;
	t134 = qJD(1) * t155 + qJD(2) * t153;
	t129 = 0.1e1 / t131;
	t119 = 0.1e1 / t121;
	t116 = t222 * t156;
	t112 = (t182 * t197 + (t135 * t210 + t134 * t165 + (-t152 * t210 + (-0.2e1 * t167 * t173 ^ 2 - t165) * t153) * qJD(2)) * t146) * t163;
	t1 = [(0.2e1 * t156 * t165 * t217 + (t133 * t165 + t156 * t192) * t146) * t163, t112, 0, 0, 0, 0; t153 * t125 * t200 + (-t135 * t125 + (t113 * t153 + t116 * t133) * t126) * t119 + ((t116 * t199 + t222 * t215) * t119 + (t116 * t200 + (-(-t117 * t146 * t195 + t197) * t213 - (t195 * t197 - t117 + (-t163 * t181 + t117) * t146) * t212) * t119) * t126) * t156, (t125 * t155 + t126 * t219) * t200 + (t199 * t219 + t132 * t125 + (t155 * t113 + t114 * t133 - (-t168 * t203 - t112 * t153 - t118 * t135 + (-t118 * t209 + t152) * t117) * t212 - (t117 * t118 * t153 + t134 + (-t112 * t170 + (-qJD(2) * t118 - t117) * t173) * t168) * t213) * t126) * t119, 0, 0, 0, 0; (-t137 * t140 + t141 * t214) * t198 + ((qJD(5) * t141 - t134 * t169 - t172 * t194) * t137 + t141 * t189 + (-t140 * t124 - (-qJD(5) * t140 - t134 * t172 + t169 * t194) * t142 - t141 * t123) * t138) * t129, t184 * t156 * t198 + (t184 * t133 + ((qJD(5) * t137 + t189) * t172 + (-t123 * t172 + (-t124 + t201) * t169) * t138) * t156) * t129, 0, 0, -0.2e1 * t218 + 0.2e1 * (t123 * t138 * t129 + (-t129 * t216 - t138 * t218) * t142) * t142, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:43
	% DurationCPUTime: 2.07s
	% Computational Cost: add. (4522->147), mult. (13478->295), div. (726->12), fcn. (17045->13), ass. (0->127)
	t250 = cos(pkin(6));
	t252 = sin(qJ(5));
	t256 = cos(qJ(5));
	t249 = sin(pkin(6));
	t257 = cos(qJ(2));
	t306 = t249 * t257;
	t269 = -t250 * t256 + t252 * t306;
	t258 = cos(qJ(1));
	t300 = t258 * t257;
	t283 = t250 * t300;
	t253 = sin(qJ(2));
	t254 = sin(qJ(1));
	t303 = t254 * t253;
	t237 = -t283 + t303;
	t305 = t249 * t258;
	t270 = -t237 * t252 + t256 * t305;
	t213 = atan2(t270, -t269);
	t208 = sin(t213);
	t209 = cos(t213);
	t191 = t208 * t270 - t209 * t269;
	t189 = 0.1e1 / t191 ^ 2;
	t301 = t258 * t253;
	t302 = t254 * t257;
	t239 = t250 * t302 + t301;
	t307 = t249 * t254;
	t229 = t239 * t252 + t256 * t307;
	t221 = t229 ^ 2;
	t187 = t221 * t189 + 0.1e1;
	t276 = qJD(2) * t250 + qJD(1);
	t297 = qJD(2) * t257;
	t216 = -qJD(1) * t283 - t258 * t297 + t276 * t303;
	t299 = qJD(1) * t249;
	t281 = t258 * t299;
	t285 = t252 * t307;
	t195 = -t216 * t252 - qJD(5) * t285 + (qJD(5) * t239 + t281) * t256;
	t319 = t189 * t229;
	t220 = t270 ^ 2;
	t233 = 0.1e1 / t269 ^ 2;
	t212 = t220 * t233 + 0.1e1;
	t210 = 0.1e1 / t212;
	t238 = t250 * t301 + t302;
	t218 = t239 * qJD(1) + t238 * qJD(2);
	t226 = t237 * t256 + t252 * t305;
	t282 = t254 * t299;
	t197 = t226 * qJD(5) + t218 * t252 + t256 * t282;
	t236 = -t250 * t252 - t256 * t306;
	t298 = qJD(2) * t253;
	t280 = t249 * t298;
	t222 = t236 * qJD(5) + t252 * t280;
	t232 = 0.1e1 / t269;
	t312 = t270 * t233;
	t273 = t197 * t232 - t222 * t312;
	t179 = t273 * t210;
	t274 = t208 * t269 + t209 * t270;
	t174 = t274 * t179 - t208 * t197 + t209 * t222;
	t188 = 0.1e1 / t191;
	t190 = t188 * t189;
	t324 = t174 * t190;
	t294 = 0.2e1 * (t195 * t319 - t221 * t324) / t187 ^ 2;
	t329 = t222 * t233;
	t308 = t249 * t253;
	t268 = t232 * t238 - t308 * t312;
	t328 = t252 * t268;
	t198 = t270 * qJD(5) + t218 * t256 - t252 * t282;
	t230 = t239 * t256 - t285;
	t284 = t250 * t303;
	t240 = -t284 + t300;
	t251 = sin(qJ(6));
	t255 = cos(qJ(6));
	t207 = t230 * t255 + t240 * t251;
	t201 = 0.1e1 / t207;
	t202 = 0.1e1 / t207 ^ 2;
	t327 = 0.2e1 * t270;
	t326 = 0.2e1 * t229;
	t196 = -t229 * qJD(5) - t216 * t256 - t252 * t281;
	t217 = -t238 * qJD(1) - t239 * qJD(2);
	t183 = t207 * qJD(6) + t196 * t251 - t217 * t255;
	t206 = t230 * t251 - t240 * t255;
	t200 = t206 ^ 2;
	t194 = t200 * t202 + 0.1e1;
	t318 = t202 * t206;
	t295 = qJD(6) * t206;
	t184 = t196 * t255 + t217 * t251 - t295;
	t321 = t184 * t201 * t202;
	t323 = (t183 * t318 - t200 * t321) / t194 ^ 2;
	t314 = t232 * t329;
	t322 = (-t197 * t312 + t220 * t314) / t212 ^ 2;
	t320 = t189 * t195;
	t317 = t206 * t255;
	t316 = t208 * t229;
	t315 = t209 * t229;
	t313 = t270 * t232;
	t310 = t240 * t252;
	t309 = t240 * t256;
	t304 = t251 * t201;
	t296 = qJD(5) * t256;
	t293 = -0.2e1 * t323;
	t292 = 0.2e1 * t323;
	t291 = -0.2e1 * t322;
	t290 = t190 * t326;
	t289 = t232 * t322;
	t288 = t189 * t316;
	t287 = t189 * t315;
	t286 = t206 * t321;
	t279 = 0.2e1 * t286;
	t278 = t314 * t327;
	t275 = qJD(6) * t309 - t216;
	t205 = -t226 * t255 - t238 * t251;
	t204 = -t226 * t251 + t238 * t255;
	t272 = t202 * t317 - t304;
	t271 = t226 * t232 - t236 * t312;
	t266 = -t208 + (t209 * t313 + t208) * t210;
	t265 = -qJD(5) * t310 - qJD(6) * t239 + t217 * t256;
	t223 = t269 * qJD(5) + t256 * t280;
	t219 = -qJD(1) * t284 - t254 * t298 + t276 * t300;
	t215 = -t239 * t251 + t255 * t309;
	t214 = t239 * t255 + t251 * t309;
	t192 = 0.1e1 / t194;
	t185 = 0.1e1 / t187;
	t182 = t210 * t328;
	t181 = t271 * t210;
	t178 = t266 * t229;
	t176 = (-t208 * t238 + t209 * t308) * t252 + t274 * t182;
	t175 = t274 * t181 - t208 * t226 + t209 * t236;
	t173 = t271 * t291 + (-t236 * t278 + t198 * t232 + (t197 * t236 + t222 * t226 - t223 * t270) * t233) * t210;
	t171 = t291 * t328 + (t268 * t296 + (-t278 * t308 + t219 * t232 + (t222 * t238 + (t197 * t253 - t270 * t297) * t249) * t233) * t252) * t210;
	t1 = [-t289 * t326 + (t195 * t232 + t229 * t329) * t210, t171, 0, 0, t173, 0; -t270 * t188 * t294 + (-t197 * t188 + (-t174 * t270 - t178 * t195) * t189) * t185 + (t178 * t189 * t294 + (0.2e1 * t178 * t324 - (-t179 * t210 * t313 + t291) * t288 - (-t289 * t327 - t179 + (t179 - t273) * t210) * t287 - t266 * t320) * t185) * t229, (t176 * t319 - t188 * t310) * t294 + (-t176 * t320 + (t217 * t252 + t240 * t296) * t188 + (t176 * t290 - t189 * t310) * t174 - (t171 * t270 - t182 * t197 + (t252 * t297 + t253 * t296) * t249 + (t182 * t269 - t238 * t252) * t179) * t287 - (-t238 * t296 + t171 * t269 - t182 * t222 - t219 * t252 + (-t182 * t270 - t252 * t308) * t179) * t288) * t185, 0, 0, (t175 * t319 - t188 * t230) * t294 + (t175 * t174 * t290 + t196 * t188 + (-t230 * t174 - t175 * t195 - (t173 * t270 - t181 * t197 + t223 + (t181 * t269 - t226) * t179) * t315 - (t173 * t269 - t181 * t222 - t198 + (-t181 * t270 - t236) * t179) * t316) * t189) * t185, 0; (-t201 * t204 + t205 * t318) * t292 + ((t205 * qJD(6) - t198 * t251 + t219 * t255) * t201 + t205 * t279 + (-t204 * t184 - (-t204 * qJD(6) - t198 * t255 - t219 * t251) * t206 - t205 * t183) * t202) * t192, (-t201 * t214 + t215 * t318) * t292 + (t215 * t279 + t275 * t201 * t255 + t265 * t304 + (t275 * t206 * t251 - t215 * t183 - t214 * t184 - t265 * t317) * t202) * t192, 0, 0, t272 * t229 * t293 + (t272 * t195 + ((-qJD(6) * t201 - 0.2e1 * t286) * t255 + (t183 * t255 + (t184 - t295) * t251) * t202) * t229) * t192, t293 + 0.2e1 * (t183 * t202 * t192 + (-t192 * t321 - t202 * t323) * t206) * t206;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end