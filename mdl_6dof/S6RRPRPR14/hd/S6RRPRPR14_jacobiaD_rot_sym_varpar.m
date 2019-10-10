% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR14
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
%   Wie in S6RRPRPR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR14_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR14_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:16
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
	% DurationCPUTime: 0.76s
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:17
	% DurationCPUTime: 0.91s
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:18
	% DurationCPUTime: 1.80s
	% Computational Cost: add. (3617->117), mult. (10934->241), div. (683->12), fcn. (13943->11), ass. (0->106)
	t196 = cos(qJ(4));
	t192 = cos(pkin(6));
	t197 = cos(qJ(2));
	t198 = cos(qJ(1));
	t235 = t197 * t198;
	t194 = sin(qJ(2));
	t195 = sin(qJ(1));
	t238 = t194 * t195;
	t213 = t192 * t235 - t238;
	t193 = sin(qJ(4));
	t191 = sin(pkin(6));
	t239 = t191 * t198;
	t225 = t193 * t239;
	t214 = -t213 * t196 + t225;
	t167 = t214 ^ 2;
	t240 = t191 * t197;
	t184 = t192 * t193 + t196 * t240;
	t179 = 0.1e1 / t184 ^ 2;
	t158 = t167 * t179 + 0.1e1;
	t154 = 0.1e1 / t158;
	t222 = t196 * t239;
	t189 = qJD(4) * t222;
	t236 = t195 * t197;
	t237 = t194 * t198;
	t212 = t192 * t236 + t237;
	t242 = t191 * t195;
	t187 = t192 * t237 + t236;
	t264 = t187 * qJD(2);
	t148 = -(t212 * qJD(1) + t264) * t196 - t189 + (qJD(1) * t242 - qJD(4) * t213) * t193;
	t185 = t192 * t196 - t193 * t240;
	t243 = t191 * t194;
	t221 = qJD(2) * t243;
	t169 = t185 * qJD(4) - t196 * t221;
	t178 = 0.1e1 / t184;
	t246 = t214 * t179;
	t217 = -t148 * t178 - t169 * t246;
	t136 = t217 * t154;
	t159 = atan2(t214, t184);
	t152 = sin(t159);
	t153 = cos(t159);
	t218 = -t152 * t184 + t153 * t214;
	t132 = t218 * t136 - t152 * t148 + t153 * t169;
	t147 = t152 * t214 + t153 * t184;
	t145 = 0.1e1 / t147 ^ 2;
	t271 = t132 * t145;
	t144 = 0.1e1 / t147;
	t270 = t144 * t271;
	t241 = t191 * t196;
	t171 = t212 * t193 + t195 * t241;
	t211 = t192 * t238 - t235;
	t205 = t211 * qJD(2);
	t268 = -t213 * qJD(1) + t205;
	t150 = qJD(1) * t225 + t171 * qJD(4) + t268 * t196;
	t261 = -t193 * t242 + t212 * t196;
	t260 = -0.2e1 * t261;
	t220 = t260 * t270;
	t269 = -t145 * t150 + t220;
	t215 = t213 * t193 + t222;
	t149 = t171 * qJD(1) + t214 * qJD(4) + t193 * t264;
	t163 = -t187 * qJD(1) - t212 * qJD(2);
	t182 = 0.1e1 / t211 ^ 2;
	t251 = t163 * t182;
	t267 = t169 * t179;
	t226 = t214 * t243;
	t210 = t178 * t187 + t179 * t226;
	t266 = t196 * t210;
	t181 = 0.1e1 / t211;
	t165 = t261 ^ 2;
	t143 = t145 * t165 + 0.1e1;
	t254 = t145 * t261;
	t259 = (-t150 * t254 - t165 * t270) / t143 ^ 2;
	t151 = t215 * qJD(1) + qJD(4) * t261 - t193 * t205;
	t166 = t171 ^ 2;
	t160 = t166 * t182 + 0.1e1;
	t248 = t171 * t182;
	t250 = t181 * t251;
	t257 = (t151 * t248 + t166 * t250) / t160 ^ 2;
	t249 = t178 * t267;
	t256 = (-t148 * t246 - t167 * t249) / t158 ^ 2;
	t253 = t152 * t261;
	t252 = t153 * t261;
	t247 = t214 * t178;
	t245 = t214 * t185;
	t244 = t211 * t196;
	t234 = qJD(2) * t197;
	t233 = qJD(4) * t193;
	t232 = 0.2e1 * t259;
	t231 = -0.2e1 * t256;
	t230 = 0.2e1 * t256;
	t228 = -0.2e1 * t171 * t187;
	t227 = t181 * t257;
	t219 = t178 * t230;
	t216 = -t178 * t215 + t179 * t245;
	t209 = t182 * t212;
	t208 = -t152 + (-t153 * t247 + t152) * t154;
	t168 = -t184 * qJD(4) + t193 * t221;
	t164 = -t211 * qJD(1) + t213 * qJD(2);
	t156 = 0.1e1 / t160;
	t141 = 0.1e1 / t143;
	t140 = t154 * t266;
	t139 = t216 * t154;
	t134 = (t152 * t187 - t153 * t243) * t196 + t218 * t140;
	t133 = -t218 * t139 + t152 * t215 + t153 * t185;
	t131 = t216 * t230 + (0.2e1 * t245 * t249 - t149 * t178 + (t148 * t185 - t168 * t214 - t169 * t215) * t179) * t154;
	t129 = t231 * t266 + (-t210 * t233 + (-0.2e1 * t226 * t249 + t164 * t178 + (-t169 * t187 + (-t148 * t194 + t214 * t234) * t191) * t179) * t196) * t154;
	t1 = [-t261 * t219 + (-t150 * t178 - t261 * t267) * t154, t129, 0, t131, 0, 0; -0.2e1 * t214 * t144 * t259 + ((t261 * qJD(1) + t196 * t264 + t213 * t233 + t189) * t144 - t214 * t271 + (t208 * t150 - ((t136 * t154 * t247 + t231) * t152 + (t214 * t219 - t136 + (t136 - t217) * t154) * t153) * t261) * t254) * t141 - (t269 * t141 - t254 * t232) * t208 * t261, (-t134 * t254 - t144 * t244) * t232 + ((-t163 * t196 - t211 * t233) * t144 + t269 * t134 + (-t244 * t132 + (t129 * t214 - t140 * t148 + (t194 * t233 - t196 * t234) * t191 + (-t140 * t184 + t187 * t196) * t136) * t252 + (-t187 * t233 - t129 * t184 - t140 * t169 + t164 * t196 + (-t140 * t214 + t194 * t241) * t136) * t253) * t145) * t141, 0, (-t133 * t254 - t144 * t171) * t232 + (t133 * t220 + t151 * t144 + (-t171 * t132 - t133 * t150 + (t131 * t214 + t139 * t148 + t168 + (t139 * t184 + t215) * t136) * t252 + (-t131 * t184 + t139 * t169 - t149 + (t139 * t214 - t185) * t136) * t253) * t145) * t141, 0, 0; t182 * t228 * t257 + 0.2e1 * t215 * t227 + (t187 * t151 * t182 + t149 * t181 + t164 * t248 - t215 * t251 - t228 * t250) * t156, 0.2e1 * (-t171 * t209 - t193) * t257 + (qJD(4) * t196 + t151 * t209 + (-t268 * t182 + 0.2e1 * t212 * t250) * t171) * t156, 0, -t227 * t260 + (t150 * t181 - t251 * t261) * t156, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:18
	% DurationCPUTime: 2.29s
	% Computational Cost: add. (4522->147), mult. (13478->287), div. (726->12), fcn. (17045->13), ass. (0->124)
	t256 = sin(qJ(4));
	t254 = cos(pkin(6));
	t261 = cos(qJ(2));
	t262 = cos(qJ(1));
	t303 = t262 * t261;
	t250 = t254 * t303;
	t257 = sin(qJ(2));
	t258 = sin(qJ(1));
	t307 = t258 * t257;
	t282 = -t250 + t307;
	t260 = cos(qJ(4));
	t253 = sin(pkin(6));
	t309 = t253 * t262;
	t291 = t260 * t309;
	t330 = -t282 * t256 + t291;
	t227 = t330 ^ 2;
	t310 = t253 * t261;
	t273 = -t254 * t260 + t256 * t310;
	t239 = 0.1e1 / t273 ^ 2;
	t217 = t227 * t239 + 0.1e1;
	t215 = 0.1e1 / t217;
	t304 = t262 * t257;
	t306 = t258 * t261;
	t243 = t254 * t304 + t306;
	t244 = t254 * t306 + t304;
	t223 = t244 * qJD(1) + t243 * qJD(2);
	t302 = qJD(1) * t253;
	t289 = t258 * t302;
	t329 = t256 * t309 + t282 * t260;
	t202 = qJD(4) * t329 + t223 * t256 + t260 * t289;
	t241 = -t254 * t256 - t260 * t310;
	t301 = qJD(2) * t257;
	t288 = t253 * t301;
	t228 = t241 * qJD(4) + t256 * t288;
	t238 = 0.1e1 / t273;
	t315 = t330 * t239;
	t276 = t202 * t238 - t228 * t315;
	t184 = t276 * t215;
	t218 = atan2(t330, -t273);
	t213 = sin(t218);
	t214 = cos(t218);
	t279 = t213 * t273 + t214 * t330;
	t179 = t279 * t184 - t213 * t202 + t214 * t228;
	t196 = t213 * t330 - t214 * t273;
	t194 = 0.1e1 / t196 ^ 2;
	t337 = t179 * t194;
	t336 = -0.2e1 * t330;
	t193 = 0.1e1 / t196;
	t335 = t193 * t337;
	t311 = t253 * t258;
	t231 = t244 * t256 + t260 * t311;
	t284 = 0.2e1 * t231 * t335;
	t281 = qJD(2) * t254 + qJD(1);
	t300 = qJD(2) * t261;
	t221 = -qJD(1) * t250 - t262 * t300 + t281 * t307;
	t290 = t262 * t302;
	t294 = t256 * t311;
	t204 = -t221 * t256 - qJD(4) * t294 + (qJD(4) * t244 + t290) * t260;
	t321 = t204 * t194;
	t334 = -t321 + t284;
	t230 = -t244 * t260 + t294;
	t293 = t254 * t307;
	t245 = -t293 + t303;
	t255 = sin(qJ(6));
	t259 = cos(qJ(6));
	t278 = t230 * t259 - t245 * t255;
	t333 = t278 * qJD(6);
	t200 = -qJD(4) * t291 - t223 * t260 + (t282 * qJD(4) + t289) * t256;
	t332 = t228 * t239;
	t312 = t253 * t257;
	t272 = -t238 * t243 + t312 * t315;
	t331 = t256 * t272;
	t212 = t230 * t255 + t245 * t259;
	t206 = 0.1e1 / t212;
	t207 = 0.1e1 / t212 ^ 2;
	t226 = t231 ^ 2;
	t192 = t226 * t194 + 0.1e1;
	t327 = (-t226 * t335 + t231 * t321) / t192 ^ 2;
	t203 = t231 * qJD(4) + t221 * t260 + t256 * t290;
	t222 = -t243 * qJD(1) - t244 * qJD(2);
	t188 = t212 * qJD(6) - t203 * t259 + t222 * t255;
	t205 = t278 ^ 2;
	t199 = t205 * t207 + 0.1e1;
	t320 = t207 * t278;
	t189 = t203 * t255 + t222 * t259 + t333;
	t323 = t189 * t206 * t207;
	t326 = (-t188 * t320 - t205 * t323) / t199 ^ 2;
	t317 = t238 * t332;
	t324 = (-t202 * t315 + t227 * t317) / t217 ^ 2;
	t322 = t194 * t231;
	t319 = t213 * t231;
	t318 = t214 * t231;
	t316 = t330 * t238;
	t314 = t245 * t256;
	t313 = t245 * t260;
	t308 = t255 * t278;
	t305 = t259 * t206;
	t299 = qJD(4) * t260;
	t298 = 0.2e1 * t327;
	t297 = 0.2e1 * t326;
	t296 = 0.2e1 * t324;
	t287 = t238 * t296;
	t286 = -0.2e1 * t278 * t323;
	t285 = t317 * t336;
	t280 = -qJD(6) * t313 + t221;
	t277 = t243 * t255 + t259 * t329;
	t210 = -t243 * t259 + t255 * t329;
	t275 = -t207 * t308 + t305;
	t274 = -t238 * t329 + t241 * t315;
	t270 = -t213 + (t214 * t316 + t213) * t215;
	t269 = qJD(4) * t314 + qJD(6) * t244 - t222 * t260;
	t229 = t273 * qJD(4) + t260 * t288;
	t224 = -qJD(1) * t293 - t258 * t301 + t281 * t303;
	t220 = -t244 * t259 - t255 * t313;
	t219 = -t244 * t255 + t259 * t313;
	t197 = 0.1e1 / t199;
	t190 = 0.1e1 / t192;
	t187 = t215 * t331;
	t186 = t274 * t215;
	t181 = (-t213 * t243 + t214 * t312) * t256 - t279 * t187;
	t180 = -t279 * t186 - t213 * t329 + t214 * t241;
	t178 = t274 * t296 + (t241 * t285 - t200 * t238 + (t202 * t241 + t228 * t329 - t229 * t330) * t239) * t215;
	t176 = t296 * t331 + (-t272 * t299 + (t285 * t312 + t224 * t238 + (t228 * t243 + (t202 * t257 - t300 * t330) * t253) * t239) * t256) * t215;
	t1 = [-t231 * t287 + (t204 * t238 + t231 * t332) * t215, t176, 0, t178, 0, 0; t193 * t327 * t336 + (-t202 * t193 - t330 * t337 - (t270 * t204 + ((-t184 * t215 * t316 - 0.2e1 * t324) * t213 + (-t330 * t287 - t184 + (t184 - t276) * t215) * t214) * t231) * t322) * t190 + (t334 * t190 + t322 * t298) * t270 * t231, (t181 * t322 - t193 * t314) * t298 + ((t222 * t256 + t245 * t299) * t193 + t334 * t181 + (-t314 * t179 - (t176 * t330 + t187 * t202 + (t256 * t300 + t257 * t299) * t253 + (-t187 * t273 - t243 * t256) * t184) * t318 - (-t243 * t299 + t176 * t273 + t187 * t228 - t224 * t256 + (t187 * t330 - t256 * t312) * t184) * t319) * t194) * t190, 0, (t180 * t322 + t193 * t230) * t298 + (t180 * t284 - t203 * t193 + (t230 * t179 - t180 * t204 - (t178 * t330 + t186 * t202 + t229 + (-t186 * t273 - t329) * t184) * t318 - (t178 * t273 + t186 * t228 + t200 + (t186 * t330 - t241) * t184) * t319) * t194) * t190, 0, 0; (t206 * t277 - t210 * t320) * t297 + ((t210 * qJD(6) + t200 * t259 - t224 * t255) * t206 + t210 * t286 + (t277 * t189 + (t277 * qJD(6) - t200 * t255 - t224 * t259) * t278 - t210 * t188) * t207) * t197, (-t206 * t219 - t220 * t320) * t297 + (t220 * t286 + t280 * t206 * t255 - t269 * t305 + (t259 * t278 * t280 - t220 * t188 - t219 * t189 + t269 * t308) * t207) * t197, 0, t275 * t231 * t297 + (-t275 * t204 + ((qJD(6) * t206 + t286) * t255 + (-t188 * t255 + (t189 + t333) * t259) * t207) * t231) * t197, 0, -0.2e1 * t326 - 0.2e1 * (t188 * t207 * t197 - (-t197 * t323 - t207 * t326) * t278) * t278;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end