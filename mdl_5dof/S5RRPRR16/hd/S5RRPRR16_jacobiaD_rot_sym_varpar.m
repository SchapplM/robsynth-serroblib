% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPRR16_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRR16_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR16_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(5));
	t93 = t99 ^ 2;
	t100 = cos(pkin(5));
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
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (1114->72), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->74)
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t156 = cos(pkin(5));
	t141 = t126 * t156;
	t139 = t125 * t141;
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t152 = t124 * t123;
	t104 = -t139 + t152;
	t122 = sin(pkin(5));
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
	t166 = (t144 * t93 - t92) * t95 + t92;
	t87 = -t104 * t92 - t153 * t93;
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
	t135 = -t123 * t126 - t125 * t142;
	t90 = -qJD(1) * t135 - qJD(2) * t134;
	t77 = (t104 * t143 + t119 * t90) * t160;
	t74 = -t77 * t163 - t92 * t90 + (t149 * t93 + t158 * t77) * t122;
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
	t148 = 0.2e1 * (-t102 * t165 + t162 * t88) / t81 ^ 2;
	t147 = 0.2e1 * t165;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
	t140 = t123 * t142;
	t108 = -t140 + t151;
	t103 = t108 ^ 2;
	t100 = t103 * t115 * t117 + 0.1e1;
	t118 = t116 * t117;
	t89 = qJD(1) * t134 + qJD(2) * t135;
	t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
	t133 = t119 * t146 + t143 * t95;
	t97 = 0.1e1 / t100;
	t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
	t79 = 0.1e1 / t81;
	t76 = t166 * t135;
	t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t158 * t78) * t122;
	t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
	t1 = [(-t133 * t135 - t159 * t88) * t114, t73, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t144 * t77 * t95 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t104 * t133 - t159 * t90) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t162 * t75) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t150 * t97 + t116 * t145) * t135) * t114, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:29
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
	t174 = sin(qJ(2));
	t175 = sin(qJ(1));
	t177 = cos(qJ(2));
	t178 = cos(qJ(1));
	t225 = cos(pkin(5));
	t194 = t178 * t225;
	t157 = t174 * t194 + t175 * t177;
	t172 = sin(pkin(5));
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
	t1 = [(0.2e1 * t160 * t169 * t221 + (t136 * t169 + t160 * t196) * t149) * t167, t115, 0, 0, 0; t157 * t128 * t204 + (-t138 * t128 + (t116 * t157 + t119 * t136) * t129) * t122 + ((t119 * t203 + t226 * t219) * t122 + (t119 * t204 + (-(-t120 * t149 * t199 + t201) * t217 - (t199 * t201 - t120 + (-t185 * t167 + t120) * t149) * t216) * t122) * t129) * t160, (t128 * t159 + t129 * t223) * t204 + (t203 * t223 + t135 * t128 + (t159 * t116 + t117 * t136 - (-t172 * t207 - t115 * t157 - t121 * t138 + (-t121 * t213 + t156) * t120) * t216 - (t120 * t121 * t157 + t137 + (-t115 * t174 + (-qJD(2) * t121 - t120) * t177) * t172) * t217) * t129) * t122, 0, 0, 0; (-t140 * t145 + t146 * t218) * t202 + ((t146 * qJD(4) + t137 * t176 - t173 * t198) * t140 + t146 * t193 + (-t145 * t127 - (-t145 * qJD(4) - t137 * t173 - t176 * t198) * t143 - t146 * t126) * t141) * t132, t188 * t160 * t202 + (t188 * t136 + ((qJD(4) * t140 + t193) * t173 + (-t126 * t173 + (t127 - t205) * t176) * t141) * t160) * t132, 0, -0.2e1 * t222 + 0.2e1 * (t126 * t141 * t132 + (-t132 * t220 - t141 * t222) * t143) * t143, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:30
	% DurationCPUTime: 1.74s
	% Computational Cost: add. (4522->152), mult. (13478->303), div. (726->12), fcn. (17045->13), ass. (0->127)
	t255 = sin(qJ(2));
	t256 = sin(qJ(1));
	t259 = cos(qJ(2));
	t260 = cos(qJ(1));
	t338 = cos(pkin(5));
	t290 = t260 * t338;
	t243 = t255 * t256 - t259 * t290;
	t254 = sin(qJ(4));
	t258 = cos(qJ(4));
	t252 = sin(pkin(5));
	t318 = t252 * t260;
	t279 = t243 * t258 + t254 * t318;
	t226 = t279 ^ 2;
	t319 = t252 * t259;
	t276 = -t338 * t254 - t258 * t319;
	t239 = 0.1e1 / t276 ^ 2;
	t217 = t226 * t239 + 0.1e1;
	t215 = 0.1e1 / t217;
	t313 = qJD(1) * t256;
	t297 = t252 * t313;
	t275 = -t255 * t290 - t256 * t259;
	t291 = t256 * t338;
	t277 = t260 * t255 + t259 * t291;
	t268 = t277 * qJD(1) - t275 * qJD(2);
	t298 = t258 * t318;
	t346 = -qJD(4) * t298 - t268 * t258;
	t200 = (qJD(4) * t243 + t297) * t254 + t346;
	t242 = -t254 * t319 + t338 * t258;
	t312 = qJD(2) * t255;
	t295 = t252 * t312;
	t228 = t242 * qJD(4) - t258 * t295;
	t238 = 0.1e1 / t276;
	t324 = t279 * t239;
	t283 = t200 * t238 - t228 * t324;
	t184 = t283 * t215;
	t218 = atan2(t279, -t276);
	t213 = sin(t218);
	t214 = cos(t218);
	t285 = t213 * t276 + t214 * t279;
	t179 = t285 * t184 - t213 * t200 + t214 * t228;
	t196 = t213 * t279 - t214 * t276;
	t194 = 0.1e1 / t196 ^ 2;
	t347 = t194 * t179;
	t310 = qJD(4) * t254;
	t292 = t252 * t310;
	t309 = qJD(4) * t258;
	t201 = t243 * t309 + t268 * t254 + t258 * t297 + t260 * t292;
	t320 = t252 * t258;
	t230 = t277 * t254 + t256 * t320;
	t286 = t255 * t291;
	t315 = t260 * t259;
	t245 = -t286 + t315;
	t253 = sin(qJ(5));
	t257 = cos(qJ(5));
	t209 = t230 * t253 - t245 * t257;
	t345 = 0.2e1 * t209;
	t193 = 0.1e1 / t196;
	t344 = t193 * t347;
	t271 = (t338 * qJD(1) + qJD(2)) * t315 - qJD(2) * t286 - t255 * t313;
	t296 = qJD(1) * t318;
	t203 = t230 * qJD(4) + t254 * t296 - t271 * t258;
	t273 = t277 * t258;
	t229 = t256 * t252 * t254 - t273;
	t289 = 0.2e1 * t229 * t344;
	t343 = -t194 * t203 + t289;
	t342 = t228 * t239;
	t321 = t252 * t255;
	t299 = t279 * t321;
	t278 = t238 * t275 + t239 * t299;
	t341 = t258 * t278;
	t340 = -t245 * t254 * qJD(5) - t271;
	t339 = -qJD(5) * t277 + t245 * t309;
	t210 = t230 * t257 + t245 * t253;
	t206 = 0.1e1 / t210;
	t207 = 0.1e1 / t210 ^ 2;
	t225 = t229 ^ 2;
	t192 = t194 * t225 + 0.1e1;
	t332 = t194 * t229;
	t337 = (t203 * t332 - t225 * t344) / t192 ^ 2;
	t326 = t238 * t342;
	t335 = (-t200 * t324 + t226 * t326) / t217 ^ 2;
	t204 = qJD(4) * t273 + t271 * t254 - t256 * t292 + t258 * t296;
	t223 = t275 * qJD(1) - t277 * qJD(2);
	t308 = qJD(5) * t209;
	t189 = t204 * t257 + t223 * t253 - t308;
	t334 = t189 * t206 * t207;
	t205 = t209 ^ 2;
	t199 = t205 * t207 + 0.1e1;
	t197 = 0.1e1 / t199;
	t331 = t197 * t207;
	t188 = t210 * qJD(5) + t204 * t253 - t223 * t257;
	t329 = t207 * t209;
	t330 = 0.1e1 / t199 ^ 2 * (t188 * t329 - t205 * t334);
	t328 = t213 * t229;
	t327 = t214 * t229;
	t325 = t279 * t238;
	t323 = t279 * t242;
	t322 = t245 * t258;
	t317 = t253 * t254;
	t316 = t254 * t257;
	t311 = qJD(2) * t259;
	t307 = 0.2e1 * t337;
	t306 = -0.2e1 * t335;
	t305 = 0.2e1 * t335;
	t303 = -0.2e1 * t330;
	t302 = t207 * t330;
	t301 = t188 * t331;
	t300 = t209 * t334;
	t288 = t238 * t305;
	t287 = 0.2e1 * t300;
	t280 = -t243 * t254 + t298;
	t212 = t253 * t275 + t257 * t280;
	t211 = t253 * t280 - t257 * t275;
	t282 = -t206 * t253 + t257 * t329;
	t281 = t238 * t280 + t239 * t323;
	t274 = -t213 + (t214 * t325 + t213) * t215;
	t227 = t276 * qJD(4) + t254 * t295;
	t224 = -qJD(1) * t286 - t256 * t312 + (qJD(2) * t338 + qJD(1)) * t315;
	t220 = t245 * t316 - t277 * t253;
	t190 = 0.1e1 / t192;
	t187 = t215 * t341;
	t186 = t281 * t215;
	t181 = (-t213 * t275 - t214 * t321) * t258 + t285 * t187;
	t180 = -t285 * t186 + t213 * t280 + t214 * t242;
	t178 = t281 * t305 + (-0.2e1 * t323 * t326 + t201 * t238 + (t200 * t242 - t227 * t279 - t228 * t280) * t239) * t215;
	t176 = t306 * t341 + (-t278 * t310 + (0.2e1 * t299 * t326 - t224 * t238 + (t228 * t275 + (-t200 * t255 + t279 * t311) * t252) * t239) * t258) * t215;
	t1 = [-t229 * t288 + (t203 * t238 + t229 * t342) * t215, t176, 0, t178, 0; -0.2e1 * t279 * t193 * t337 + ((-t243 * t310 - t254 * t297 - t346) * t193 - t279 * t347 - (t274 * t203 + ((-t184 * t215 * t325 + t306) * t213 + (-t279 * t288 - t184 + (t184 - t283) * t215) * t214) * t229) * t332) * t190 + (t190 * t343 + t332 * t307) * t274 * t229, (t181 * t332 + t193 * t322) * t307 + ((-t223 * t258 + t245 * t310) * t193 + t343 * t181 + (t322 * t179 - (t176 * t279 - t187 * t200 + (t255 * t310 - t258 * t311) * t252 + (t187 * t276 - t258 * t275) * t184) * t327 - (t275 * t310 + t176 * t276 - t187 * t228 + t224 * t258 + (-t187 * t279 + t255 * t320) * t184) * t328) * t194) * t190, 0, (t180 * t332 - t193 * t230) * t307 + (t180 * t289 + t204 * t193 + (-t230 * t179 - t180 * t203 - (t178 * t279 + t186 * t200 + t227 + (-t186 * t276 + t280) * t184) * t327 - (t178 * t276 + t186 * t228 - t201 + (t186 * t279 - t242) * t184) * t328) * t194) * t190, 0; 0.2e1 * (-t206 * t211 + t212 * t329) * t330 + ((t212 * qJD(5) - t201 * t253 + t224 * t257) * t206 + t212 * t287 + (-t211 * t189 - (-t211 * qJD(5) - t201 * t257 - t224 * t253) * t209 - t212 * t188) * t207) * t197, (t302 * t345 - t301) * t220 + (-t189 * t331 + t206 * t303) * (t245 * t317 + t277 * t257) + (t220 * t287 + (t317 * t206 - t316 * t329) * t223 + (-t340 * t206 - t339 * t329) * t257 + (t339 * t206 - t340 * t329) * t253) * t197, 0, t282 * t229 * t303 + (t282 * t203 + ((-qJD(5) * t206 - 0.2e1 * t300) * t257 + (t188 * t257 + (t189 - t308) * t253) * t207) * t229) * t197, t303 + (t301 + (-t197 * t334 - t302) * t209) * t345;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end