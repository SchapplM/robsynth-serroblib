% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR13
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
%   Wie in S6RRPRPR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR13_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR13_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
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
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:25
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (1114->72), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->74)
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t156 = cos(pkin(6));
	t141 = t126 * t156;
	t139 = t125 * t141;
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t152 = t124 * t123;
	t104 = -t139 + t152;
	t122 = sin(pkin(6));
	t114 = 0.1e1 / t122;
	t119 = 0.1e1 / t125;
	t144 = t104 * t114 * t119;
	t153 = t122 * t125;
	t94 = atan2(-t104, -t153);
	t92 = sin(t94);
	t93 = cos(t94);
	t101 = t104 ^ 2;
	t115 = 0.1e1 / t122 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t99 = t101 * t115 * t120 + 0.1e1;
	t95 = 0.1e1 / t99;
	t166 = (t93 * t144 - t92) * t95 + t92;
	t87 = -t92 * t104 - t93 * t153;
	t84 = 0.1e1 / t87;
	t116 = 0.1e1 / t124;
	t117 = 0.1e1 / t124 ^ 2;
	t85 = 0.1e1 / t87 ^ 2;
	t149 = qJD(2) * t123;
	t158 = t125 * t92;
	t163 = t104 * t93;
	t143 = t120 * t149;
	t160 = t114 * t95;
	t134 = -t123 * t141 - t124 * t125;
	t142 = t124 * t156;
	t135 = -t126 * t123 - t125 * t142;
	t90 = -t135 * qJD(1) - t134 * qJD(2);
	t77 = (t104 * t143 + t119 * t90) * t160;
	t74 = -t77 * t163 - t92 * t90 + (t93 * t149 + t77 * t158) * t122;
	t165 = t74 * t84 * t85;
	t154 = t120 * t123;
	t136 = t104 * t154 - t119 * t134;
	t78 = t136 * t160;
	t164 = t77 * t78;
	t162 = t135 * t85;
	t161 = t135 * t93;
	t159 = t119 * t95;
	t157 = t92 * t135;
	t155 = t117 * t126;
	t151 = t126 * t125;
	t150 = qJD(1) * t126;
	t102 = t135 ^ 2;
	t81 = t102 * t85 + 0.1e1;
	t138 = qJD(2) * t156 + qJD(1);
	t88 = -qJD(1) * t139 - qJD(2) * t151 + t138 * t152;
	t148 = 0.2e1 * (-t102 * t165 + t88 * t162) / t81 ^ 2;
	t147 = 0.2e1 * t165;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
	t140 = t123 * t142;
	t108 = -t140 + t151;
	t103 = t108 ^ 2;
	t100 = t103 * t117 * t115 + 0.1e1;
	t118 = t116 * t117;
	t89 = t134 * qJD(1) + t135 * qJD(2);
	t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
	t133 = t119 * t146 + t95 * t143;
	t97 = 0.1e1 / t100;
	t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
	t79 = 0.1e1 / t81;
	t76 = t166 * t135;
	t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t78 * t158) * t122;
	t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
	t1 = [(-t133 * t135 - t88 * t159) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t77 * t95 * t144 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t133 * t104 - t90 * t159) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t75 * t162) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t97 * t150 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:25
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
	t174 = sin(qJ(2));
	t175 = sin(qJ(1));
	t177 = cos(qJ(2));
	t178 = cos(qJ(1));
	t225 = cos(pkin(6));
	t194 = t178 * t225;
	t157 = t174 * t194 + t175 * t177;
	t172 = sin(pkin(6));
	t213 = t172 * t174;
	t151 = atan2(-t157, t213);
	t147 = sin(t151);
	t148 = cos(t151);
	t154 = t157 ^ 2;
	t168 = 0.1e1 / t172 ^ 2;
	t170 = 0.1e1 / t174 ^ 2;
	t152 = t154 * t168 * t170 + 0.1e1;
	t149 = 0.1e1 / t152;
	t167 = 0.1e1 / t172;
	t169 = 0.1e1 / t174;
	t199 = t157 * t167 * t169;
	t226 = (t148 * t199 + t147) * t149 - t147;
	t131 = -t147 * t157 + t148 * t213;
	t128 = 0.1e1 / t131;
	t195 = t175 * t225;
	t159 = t178 * t174 + t177 * t195;
	t173 = sin(qJ(4));
	t176 = cos(qJ(4));
	t212 = t172 * t175;
	t144 = t159 * t173 + t176 * t212;
	t140 = 0.1e1 / t144;
	t129 = 0.1e1 / t131 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t190 = qJD(2) * t225 + qJD(1);
	t192 = t174 * t195;
	t207 = qJD(2) * t174;
	t209 = t178 * t177;
	t138 = -qJD(1) * t192 - t175 * t207 + t190 * t209;
	t206 = qJD(2) * t177;
	t196 = t170 * t206;
	t185 = -t138 * t169 + t157 * t196;
	t215 = t149 * t167;
	t120 = t185 * t215;
	t187 = -t147 * t213 - t148 * t157;
	t200 = t148 * t172 * t177;
	t116 = qJD(2) * t200 + t187 * t120 - t147 * t138;
	t224 = t116 * t128 * t129;
	t191 = t177 * t194;
	t210 = t175 * t174;
	t156 = -t191 + t210;
	t214 = t170 * t177;
	t186 = t156 * t169 + t157 * t214;
	t121 = t186 * t215;
	t117 = t187 * t121 + t147 * t156 + t200;
	t160 = -t192 + t209;
	t223 = t117 * t160;
	t135 = -qJD(1) * t191 - t178 * t206 + t190 * t210;
	t208 = qJD(1) * t172;
	t197 = t178 * t208;
	t126 = t144 * qJD(4) + t135 * t176 + t173 * t197;
	t143 = -t159 * t176 + t173 * t212;
	t139 = t143 ^ 2;
	t134 = t139 * t141 + 0.1e1;
	t218 = t141 * t143;
	t205 = qJD(4) * t143;
	t127 = -t135 * t173 + t176 * t197 - t205;
	t220 = t127 * t140 * t141;
	t222 = (t126 * t218 - t139 * t220) / t134 ^ 2;
	t171 = t169 * t170;
	t221 = (t138 * t157 * t170 - t154 * t171 * t206) * t168 / t152 ^ 2;
	t136 = t157 * qJD(1) + t159 * qJD(2);
	t219 = t136 * t129;
	t217 = t147 * t160;
	t216 = t148 * t160;
	t211 = t172 * t178;
	t155 = t160 ^ 2;
	t124 = t155 * t129 + 0.1e1;
	t204 = 0.2e1 * (-t155 * t224 - t160 * t219) / t124 ^ 2;
	t203 = 0.2e1 * t224;
	t202 = 0.2e1 * t222;
	t201 = -0.2e1 * t221;
	t198 = t175 * t208;
	t193 = 0.2e1 * t143 * t220;
	t188 = t176 * t140 + t173 * t218;
	t146 = -t156 * t173 + t176 * t211;
	t145 = t156 * t176 + t173 * t211;
	t137 = t159 * qJD(1) + qJD(2) * t157;
	t132 = 0.1e1 / t134;
	t122 = 0.1e1 / t124;
	t119 = t226 * t160;
	t115 = (t186 * t201 + (t138 * t214 + t137 * t169 + (-t156 * t214 + (-0.2e1 * t171 * t177 ^ 2 - t169) * t157) * qJD(2)) * t149) * t167;
	t1 = [(0.2e1 * t160 * t169 * t221 + (t136 * t169 + t160 * t196) * t149) * t167, t115, 0, 0, 0, 0; t157 * t128 * t204 + (-t138 * t128 + (t116 * t157 + t119 * t136) * t129) * t122 + ((t119 * t203 + t226 * t219) * t122 + (t119 * t204 + (-(-t120 * t149 * t199 + t201) * t217 - (t199 * t201 - t120 + (-t167 * t185 + t120) * t149) * t216) * t122) * t129) * t160, (t128 * t159 + t129 * t223) * t204 + (t203 * t223 + t135 * t128 + (t159 * t116 + t117 * t136 - (-t172 * t207 - t115 * t157 - t121 * t138 + (-t121 * t213 + t156) * t120) * t216 - (t120 * t121 * t157 + t137 + (-t115 * t174 + (-qJD(2) * t121 - t120) * t177) * t172) * t217) * t129) * t122, 0, 0, 0, 0; (-t140 * t145 + t146 * t218) * t202 + ((qJD(4) * t146 + t137 * t176 - t173 * t198) * t140 + t146 * t193 + (-t145 * t127 - (-qJD(4) * t145 - t137 * t173 - t176 * t198) * t143 - t146 * t126) * t141) * t132, t188 * t160 * t202 + (t188 * t136 + ((qJD(4) * t140 + t193) * t173 + (-t126 * t173 + (t127 - t205) * t176) * t141) * t160) * t132, 0, -0.2e1 * t222 + 0.2e1 * (t126 * t141 * t132 + (-t132 * t220 - t141 * t222) * t143) * t143, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:27
	% DurationCPUTime: 2.34s
	% Computational Cost: add. (4128->138), mult. (12381->287), div. (702->12), fcn. (15714->13), ass. (0->123)
	t222 = sin(qJ(2));
	t223 = sin(qJ(1));
	t225 = cos(qJ(2));
	t226 = cos(qJ(1));
	t299 = cos(pkin(6));
	t255 = t226 * t299;
	t209 = t222 * t223 - t225 * t255;
	t221 = sin(qJ(4));
	t224 = cos(qJ(4));
	t219 = sin(pkin(6));
	t277 = t219 * t226;
	t245 = t209 * t224 + t221 * t277;
	t193 = t245 ^ 2;
	t278 = t219 * t225;
	t241 = -t299 * t221 - t224 * t278;
	t205 = 0.1e1 / t241 ^ 2;
	t184 = t193 * t205 + 0.1e1;
	t182 = 0.1e1 / t184;
	t208 = -t221 * t278 + t299 * t224;
	t272 = qJD(2) * t222;
	t258 = t219 * t272;
	t195 = t208 * qJD(4) - t224 * t258;
	t204 = 0.1e1 / t241;
	t285 = t245 * t205;
	t240 = -t222 * t255 - t223 * t225;
	t256 = t223 * t299;
	t242 = t226 * t222 + t225 * t256;
	t234 = t242 * qJD(1) - t240 * qJD(2);
	t261 = t224 * t277;
	t273 = qJD(1) * t223;
	t300 = t209 * qJD(4) + t219 * t273;
	t307 = qJD(4) * t261 - t300 * t221 + t234 * t224;
	t248 = -t195 * t285 - t204 * t307;
	t151 = t248 * t182;
	t185 = atan2(t245, -t241);
	t180 = sin(t185);
	t181 = cos(t185);
	t250 = t180 * t241 + t181 * t245;
	t146 = t250 * t151 + t180 * t307 + t181 * t195;
	t161 = t180 * t245 - t181 * t241;
	t159 = 0.1e1 / t161 ^ 2;
	t308 = t146 * t159;
	t270 = qJD(4) * t221;
	t257 = t219 * t270;
	t168 = t234 * t221 + t300 * t224 + t226 * t257;
	t158 = 0.1e1 / t161;
	t306 = t158 * t308;
	t238 = t242 * t224;
	t196 = t223 * t219 * t221 - t238;
	t254 = 0.2e1 * t196 * t306;
	t279 = t219 * t224;
	t197 = t242 * t221 + t223 * t279;
	t251 = t222 * t256;
	t275 = t226 * t225;
	t237 = (t299 * qJD(1) + qJD(2)) * t275 - qJD(2) * t251 - t222 * t273;
	t259 = qJD(1) * t277;
	t170 = t197 * qJD(4) + t221 * t259 - t237 * t224;
	t292 = t170 * t159;
	t305 = -t292 + t254;
	t171 = qJD(4) * t238 + t237 * t221 - t223 * t257 + t224 * t259;
	t190 = t240 * qJD(1) - t242 * qJD(2);
	t218 = sin(pkin(11));
	t220 = cos(pkin(11));
	t163 = t171 * t220 + t190 * t218;
	t211 = -t251 + t275;
	t177 = t197 * t220 + t211 * t218;
	t174 = 0.1e1 / t177 ^ 2;
	t303 = t163 * t174;
	t302 = t195 * t205;
	t280 = t219 * t222;
	t262 = t245 * t280;
	t243 = t204 * t240 + t205 * t262;
	t301 = t224 * t243;
	t173 = 0.1e1 / t177;
	t192 = t196 ^ 2;
	t157 = t159 * t192 + 0.1e1;
	t298 = (-t192 * t306 + t196 * t292) / t157 ^ 2;
	t162 = t171 * t218 - t190 * t220;
	t176 = t197 * t218 - t211 * t220;
	t172 = t176 ^ 2;
	t166 = t172 * t174 + 0.1e1;
	t291 = t174 * t176;
	t293 = t173 * t303;
	t297 = (t162 * t291 - t172 * t293) / t166 ^ 2;
	t287 = t204 * t302;
	t295 = (t193 * t287 + t285 * t307) / t184 ^ 2;
	t294 = t159 * t196;
	t283 = t211 * t221;
	t187 = -t242 * t218 + t220 * t283;
	t290 = t174 * t187;
	t289 = t180 * t196;
	t288 = t181 * t196;
	t286 = t245 * t204;
	t284 = t245 * t208;
	t282 = t211 * t224;
	t281 = t218 * t173;
	t276 = t220 * t176;
	t271 = qJD(2) * t225;
	t268 = 0.2e1 * t298;
	t267 = 0.2e1 * t297;
	t266 = -0.2e1 * t295;
	t265 = 0.2e1 * t295;
	t263 = t176 * t293;
	t253 = t204 * t265;
	t252 = 0.2e1 * t263;
	t246 = -t209 * t221 + t261;
	t247 = t204 * t246 + t205 * t284;
	t244 = qJD(4) * t282 + t190 * t221;
	t239 = -t180 + (t181 * t286 + t180) * t182;
	t194 = t241 * qJD(4) + t221 * t258;
	t191 = -qJD(1) * t251 - t223 * t272 + (qJD(2) * t299 + qJD(1)) * t275;
	t186 = t218 * t283 + t242 * t220;
	t179 = t218 * t240 + t220 * t246;
	t178 = t218 * t246 - t220 * t240;
	t164 = 0.1e1 / t166;
	t155 = 0.1e1 / t157;
	t154 = t182 * t301;
	t153 = t247 * t182;
	t148 = (-t180 * t240 - t181 * t280) * t224 + t250 * t154;
	t147 = -t250 * t153 + t180 * t246 + t181 * t208;
	t145 = t247 * t265 + (-0.2e1 * t284 * t287 + t168 * t204 + (-t194 * t245 - t195 * t246 - t208 * t307) * t205) * t182;
	t143 = t266 * t301 + (-t243 * t270 + (0.2e1 * t262 * t287 - t191 * t204 + (t195 * t240 + (t222 * t307 + t245 * t271) * t219) * t205) * t224) * t182;
	t1 = [-t196 * t253 + (t170 * t204 + t196 * t302) * t182, t143, 0, t145, 0, 0; -0.2e1 * t245 * t158 * t298 + (t307 * t158 - t245 * t308 - (t239 * t170 + ((-t151 * t182 * t286 + t266) * t180 + (-t245 * t253 - t151 + (t151 - t248) * t182) * t181) * t196) * t294) * t155 + (t305 * t155 + t294 * t268) * t239 * t196, (t148 * t294 + t158 * t282) * t268 + ((-t190 * t224 + t211 * t270) * t158 + t305 * t148 + (t282 * t146 - (t143 * t245 + t154 * t307 + (t222 * t270 - t224 * t271) * t219 + (t154 * t241 - t224 * t240) * t151) * t288 - (t240 * t270 + t143 * t241 - t154 * t195 + t191 * t224 + (-t154 * t245 + t222 * t279) * t151) * t289) * t159) * t155, 0, (t147 * t294 - t158 * t197) * t268 + (t147 * t254 + t171 * t158 + (-t197 * t146 - t147 * t170 - (t145 * t245 - t153 * t307 + t194 + (-t153 * t241 + t246) * t151) * t288 - (t145 * t241 + t153 * t195 - t168 + (t153 * t245 - t208) * t151) * t289) * t159) * t155, 0, 0; (-t173 * t178 + t179 * t291) * t267 + ((-t168 * t218 + t191 * t220) * t173 + t179 * t252 + (-t178 * t163 - (-t168 * t220 - t191 * t218) * t176 - t179 * t162) * t174) * t164, -0.2e1 * t186 * t173 * t297 + t176 * t267 * t290 + ((t244 * t218 + t220 * t237) * t173 - t186 * t303 - (-t218 * t237 + t244 * t220) * t291 - t162 * t290 + t187 * t252) * t164, 0, (-t174 * t276 + t281) * t196 * t267 + (-0.2e1 * t196 * t220 * t263 - t170 * t281 + (t170 * t276 + (t162 * t220 + t163 * t218) * t196) * t174) * t164, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:25
	% EndTime: 2019-10-10 10:26:27
	% DurationCPUTime: 2.69s
	% Computational Cost: add. (4943->148), mult. (13478->296), div. (726->12), fcn. (17045->13), ass. (0->124)
	t261 = sin(qJ(2));
	t262 = sin(qJ(1));
	t264 = cos(qJ(2));
	t265 = cos(qJ(1));
	t343 = cos(pkin(6));
	t295 = t265 * t343;
	t280 = -t261 * t295 - t262 * t264;
	t296 = t262 * t343;
	t282 = t265 * t261 + t264 * t296;
	t259 = sin(pkin(6));
	t323 = t259 * t265;
	t356 = t282 * qJD(1) - t280 * qJD(2) + qJD(4) * t323;
	t247 = t262 * t261 - t264 * t295;
	t260 = sin(qJ(4));
	t263 = cos(qJ(4));
	t284 = t247 * t263 + t260 * t323;
	t230 = t284 ^ 2;
	t324 = t259 * t264;
	t281 = -t343 * t260 - t263 * t324;
	t243 = 0.1e1 / t281 ^ 2;
	t221 = t230 * t243 + 0.1e1;
	t219 = 0.1e1 / t221;
	t246 = -t260 * t324 + t343 * t263;
	t317 = qJD(2) * t261;
	t300 = t259 * t317;
	t232 = t246 * qJD(4) - t263 * t300;
	t242 = 0.1e1 / t281;
	t329 = t284 * t243;
	t318 = qJD(1) * t262;
	t344 = t247 * qJD(4) + t259 * t318;
	t354 = -t344 * t260 + t356 * t263;
	t288 = -t232 * t329 - t242 * t354;
	t188 = t288 * t219;
	t222 = atan2(t284, -t281);
	t217 = sin(t222);
	t218 = cos(t222);
	t290 = t217 * t281 + t218 * t284;
	t183 = t290 * t188 + t217 * t354 + t218 * t232;
	t200 = t217 * t284 - t218 * t281;
	t198 = 0.1e1 / t200 ^ 2;
	t355 = t183 * t198;
	t205 = t356 * t260 + t344 * t263;
	t325 = t259 * t263;
	t234 = t282 * t260 + t262 * t325;
	t291 = t261 * t296;
	t320 = t265 * t264;
	t249 = -t291 + t320;
	t258 = pkin(11) + qJ(6);
	t256 = sin(t258);
	t257 = cos(t258);
	t213 = t234 * t256 - t249 * t257;
	t353 = 0.2e1 * t213;
	t197 = 0.1e1 / t200;
	t352 = t197 * t355;
	t345 = -t262 * t259 * t260 + t282 * t263;
	t292 = -0.2e1 * t345 * t352;
	t276 = (t343 * qJD(1) + qJD(2)) * t320 - qJD(2) * t291 - t261 * t318;
	t301 = qJD(1) * t323;
	t207 = t234 * qJD(4) + t260 * t301 - t276 * t263;
	t335 = t207 * t198;
	t351 = -t335 + t292;
	t349 = t232 * t243;
	t326 = t259 * t261;
	t304 = t284 * t326;
	t283 = t242 * t280 + t243 * t304;
	t348 = t263 * t283;
	t347 = -t249 * t260 * qJD(6) - t276;
	t327 = t249 * t263;
	t346 = qJD(4) * t327 - t282 * qJD(6);
	t214 = t234 * t257 + t249 * t256;
	t210 = 0.1e1 / t214;
	t211 = 0.1e1 / t214 ^ 2;
	t229 = t345 ^ 2;
	t196 = t229 * t198 + 0.1e1;
	t342 = (-t229 * t352 - t335 * t345) / t196 ^ 2;
	t208 = t345 * qJD(4) + t276 * t260 + t263 * t301;
	t227 = t280 * qJD(1) - t282 * qJD(2);
	t191 = t214 * qJD(6) + t208 * t256 - t227 * t257;
	t209 = t213 ^ 2;
	t203 = t209 * t211 + 0.1e1;
	t334 = t211 * t213;
	t314 = qJD(6) * t213;
	t192 = t208 * t257 + t227 * t256 - t314;
	t338 = t192 * t210 * t211;
	t341 = (t191 * t334 - t209 * t338) / t203 ^ 2;
	t331 = t242 * t349;
	t339 = (t230 * t331 + t329 * t354) / t221 ^ 2;
	t337 = t198 * t345;
	t201 = 0.1e1 / t203;
	t336 = t201 * t211;
	t333 = t217 * t345;
	t332 = t218 * t345;
	t330 = t284 * t242;
	t328 = t284 * t246;
	t322 = t260 * t256;
	t321 = t260 * t257;
	t316 = qJD(2) * t264;
	t315 = qJD(4) * t260;
	t312 = 0.2e1 * t342;
	t311 = -0.2e1 * t341;
	t310 = -0.2e1 * t339;
	t309 = 0.2e1 * t339;
	t307 = t211 * t341;
	t306 = t191 * t336;
	t305 = t213 * t338;
	t294 = t242 * t309;
	t293 = 0.2e1 * t305;
	t285 = -t247 * t260 + t263 * t323;
	t216 = t256 * t280 + t257 * t285;
	t215 = t256 * t285 - t257 * t280;
	t287 = -t256 * t210 + t257 * t334;
	t286 = t242 * t285 + t243 * t328;
	t279 = -t217 + (t218 * t330 + t217) * t219;
	t231 = t281 * qJD(4) + t260 * t300;
	t228 = -qJD(1) * t291 - t262 * t317 + (qJD(2) * t343 + qJD(1)) * t320;
	t224 = t249 * t321 - t282 * t256;
	t194 = 0.1e1 / t196;
	t193 = t219 * t348;
	t190 = t286 * t219;
	t185 = (-t217 * t280 - t218 * t326) * t263 + t290 * t193;
	t184 = -t290 * t190 + t217 * t285 + t218 * t246;
	t182 = t286 * t309 + (-0.2e1 * t328 * t331 + t205 * t242 + (-t231 * t284 - t232 * t285 - t246 * t354) * t243) * t219;
	t180 = t310 * t348 + (-t283 * t315 + (0.2e1 * t304 * t331 - t228 * t242 + (t232 * t280 + (t261 * t354 + t284 * t316) * t259) * t243) * t263) * t219;
	t1 = [t345 * t294 + (t207 * t242 - t345 * t349) * t219, t180, 0, t182, 0, 0; -0.2e1 * t284 * t197 * t342 + (t354 * t197 - t284 * t355 + (t279 * t207 - ((-t188 * t219 * t330 + t310) * t217 + (-t284 * t294 - t188 + (t188 - t288) * t219) * t218) * t345) * t337) * t194 - (t351 * t194 - t337 * t312) * t279 * t345, (-t185 * t337 + t197 * t327) * t312 + ((-t227 * t263 + t249 * t315) * t197 + t351 * t185 + (t327 * t183 + (t180 * t284 + t193 * t354 + (t261 * t315 - t263 * t316) * t259 + (t193 * t281 - t263 * t280) * t188) * t332 + (t280 * t315 + t180 * t281 - t193 * t232 + t228 * t263 + (-t193 * t284 + t261 * t325) * t188) * t333) * t198) * t194, 0, (-t184 * t337 - t197 * t234) * t312 + (t184 * t292 + t208 * t197 + (-t234 * t183 - t184 * t207 + (t182 * t284 - t190 * t354 + t231 + (-t190 * t281 + t285) * t188) * t332 + (t182 * t281 + t190 * t232 - t205 + (t190 * t284 - t246) * t188) * t333) * t198) * t194, 0, 0; 0.2e1 * (-t210 * t215 + t216 * t334) * t341 + ((t216 * qJD(6) - t205 * t256 + t228 * t257) * t210 + t216 * t293 + (-t215 * t192 - (-t215 * qJD(6) - t205 * t257 - t228 * t256) * t213 - t216 * t191) * t211) * t201, (t307 * t353 - t306) * t224 + (-t192 * t336 + t210 * t311) * (t249 * t322 + t282 * t257) + (t224 * t293 + (t322 * t210 - t321 * t334) * t227 + (-t347 * t210 - t346 * t334) * t257 + (t346 * t210 - t347 * t334) * t256) * t201, 0, -t287 * t345 * t311 + (t287 * t207 - ((-qJD(6) * t210 - 0.2e1 * t305) * t257 + (t191 * t257 + (t192 - t314) * t256) * t211) * t345) * t201, 0, t311 + (t306 + (-t201 * t338 - t307) * t213) * t353;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end