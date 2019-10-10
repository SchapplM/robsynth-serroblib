% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP13
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
%   Wie in S6RRPRRP13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:48
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP13_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP13_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:31
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
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:32
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
	t1 = [(-t133 * t135 - t159 * t88) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t144 * t77 * t95 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t104 * t133 - t159 * t90) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t162 * t75) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t150 * t97 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:32
	% DurationCPUTime: 0.88s
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
	t226 = t149 * (t148 * t199 + t147) - t147;
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
	t116 = qJD(2) * t200 + t120 * t187 - t147 * t138;
	t224 = t116 * t128 * t129;
	t191 = t177 * t194;
	t210 = t175 * t174;
	t156 = -t191 + t210;
	t214 = t170 * t177;
	t186 = t156 * t169 + t157 * t214;
	t121 = t186 * t215;
	t117 = t121 * t187 + t147 * t156 + t200;
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
	t136 = qJD(1) * t157 + qJD(2) * t159;
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
	t137 = qJD(1) * t159 + qJD(2) * t157;
	t132 = 0.1e1 / t134;
	t122 = 0.1e1 / t124;
	t119 = t226 * t160;
	t115 = (t186 * t201 + (t138 * t214 + t137 * t169 + (-t156 * t214 + (-0.2e1 * t171 * t177 ^ 2 - t169) * t157) * qJD(2)) * t149) * t167;
	t1 = [(0.2e1 * t160 * t169 * t221 + (t136 * t169 + t160 * t196) * t149) * t167, t115, 0, 0, 0, 0; t157 * t128 * t204 + (-t138 * t128 + (t116 * t157 + t119 * t136) * t129) * t122 + ((t119 * t203 + t226 * t219) * t122 + (t119 * t204 + (-(-t120 * t149 * t199 + t201) * t217 - (t199 * t201 - t120 + (-t167 * t185 + t120) * t149) * t216) * t122) * t129) * t160, (t128 * t159 + t129 * t223) * t204 + (t203 * t223 + t135 * t128 + (t159 * t116 + t117 * t136 - (-t172 * t207 - t115 * t157 - t121 * t138 + (-t121 * t213 + t156) * t120) * t216 - (t120 * t121 * t157 + t137 + (-t115 * t174 + (-qJD(2) * t121 - t120) * t177) * t172) * t217) * t129) * t122, 0, 0, 0, 0; (-t140 * t145 + t146 * t218) * t202 + ((t146 * qJD(4) + t137 * t176 - t173 * t198) * t140 + t146 * t193 + (-t145 * t127 - (-t145 * qJD(4) - t137 * t173 - t176 * t198) * t143 - t146 * t126) * t141) * t132, t188 * t160 * t202 + (t188 * t136 + ((qJD(4) * t140 + t193) * t173 + (-t126 * t173 + (t127 - t205) * t176) * t141) * t160) * t132, 0, -0.2e1 * t222 + 0.2e1 * (t126 * t141 * t132 + (-t132 * t220 - t141 * t222) * t143) * t143, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:34
	% DurationCPUTime: 2.54s
	% Computational Cost: add. (4522->152), mult. (13478->303), div. (726->12), fcn. (17045->13), ass. (0->127)
	t255 = sin(qJ(2));
	t256 = sin(qJ(1));
	t259 = cos(qJ(2));
	t260 = cos(qJ(1));
	t338 = cos(pkin(6));
	t290 = t260 * t338;
	t243 = t255 * t256 - t259 * t290;
	t254 = sin(qJ(4));
	t258 = cos(qJ(4));
	t252 = sin(pkin(6));
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
	t179 = t285 * t184 - t200 * t213 + t214 * t228;
	t196 = t213 * t279 - t214 * t276;
	t194 = 0.1e1 / t196 ^ 2;
	t347 = t179 * t194;
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
	t331 = t194 * t229;
	t337 = (t203 * t331 - t225 * t344) / t192 ^ 2;
	t204 = qJD(4) * t273 + t271 * t254 - t256 * t292 + t258 * t296;
	t223 = t275 * qJD(1) - t277 * qJD(2);
	t188 = t210 * qJD(5) + t204 * t253 - t223 * t257;
	t205 = t209 ^ 2;
	t199 = t205 * t207 + 0.1e1;
	t329 = t207 * t209;
	t308 = qJD(5) * t209;
	t189 = t204 * t257 + t223 * t253 - t308;
	t333 = t189 * t206 * t207;
	t336 = (t188 * t329 - t205 * t333) / t199 ^ 2;
	t326 = t238 * t342;
	t334 = (-t200 * t324 + t226 * t326) / t217 ^ 2;
	t197 = 0.1e1 / t199;
	t330 = t197 * t207;
	t328 = t213 * t229;
	t327 = t214 * t229;
	t325 = t279 * t238;
	t323 = t279 * t242;
	t322 = t245 * t258;
	t317 = t253 * t254;
	t316 = t254 * t257;
	t311 = qJD(2) * t259;
	t307 = 0.2e1 * t337;
	t306 = -0.2e1 * t336;
	t305 = -0.2e1 * t334;
	t304 = 0.2e1 * t334;
	t302 = t207 * t336;
	t301 = t188 * t330;
	t300 = t209 * t333;
	t288 = t238 * t304;
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
	t178 = t281 * t304 + (-0.2e1 * t323 * t326 + t201 * t238 + (t200 * t242 - t227 * t279 - t228 * t280) * t239) * t215;
	t176 = t305 * t341 + (-t278 * t310 + (0.2e1 * t299 * t326 - t224 * t238 + (t228 * t275 + (-t200 * t255 + t279 * t311) * t252) * t239) * t258) * t215;
	t1 = [-t229 * t288 + (t203 * t238 + t229 * t342) * t215, t176, 0, t178, 0, 0; -0.2e1 * t279 * t193 * t337 + ((-t243 * t310 - t254 * t297 - t346) * t193 - t279 * t347 - (t274 * t203 + ((-t184 * t215 * t325 + t305) * t213 + (-t279 * t288 - t184 + (t184 - t283) * t215) * t214) * t229) * t331) * t190 + (t343 * t190 + t331 * t307) * t274 * t229, (t181 * t331 + t193 * t322) * t307 + ((-t223 * t258 + t245 * t310) * t193 + t343 * t181 + (t322 * t179 - (t176 * t279 - t187 * t200 + (t255 * t310 - t258 * t311) * t252 + (t187 * t276 - t258 * t275) * t184) * t327 - (t275 * t310 + t176 * t276 - t187 * t228 + t224 * t258 + (-t187 * t279 + t255 * t320) * t184) * t328) * t194) * t190, 0, (t180 * t331 - t193 * t230) * t307 + (t180 * t289 + t204 * t193 + (-t230 * t179 - t180 * t203 - (t178 * t279 + t186 * t200 + t227 + (-t186 * t276 + t280) * t184) * t327 - (t178 * t276 + t186 * t228 - t201 + (t186 * t279 - t242) * t184) * t328) * t194) * t190, 0, 0; 0.2e1 * (-t206 * t211 + t212 * t329) * t336 + ((t212 * qJD(5) - t201 * t253 + t224 * t257) * t206 + t212 * t287 + (-t211 * t189 - (-t211 * qJD(5) - t201 * t257 - t224 * t253) * t209 - t212 * t188) * t207) * t197, (t302 * t345 - t301) * t220 + (-t189 * t330 + t206 * t306) * (t245 * t317 + t277 * t257) + (t220 * t287 + (t317 * t206 - t316 * t329) * t223 + (-t340 * t206 - t339 * t329) * t257 + (t339 * t206 - t340 * t329) * t253) * t197, 0, t282 * t229 * t306 + (t282 * t203 + ((-qJD(5) * t206 - 0.2e1 * t300) * t257 + (t188 * t257 + (t189 - t308) * t253) * t207) * t229) * t197, t306 + (t301 + (-t197 * t333 - t302) * t209) * t345, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:34
	% DurationCPUTime: 2.80s
	% Computational Cost: add. (4522->147), mult. (13478->296), div. (726->12), fcn. (17045->13), ass. (0->123)
	t261 = sin(qJ(2));
	t262 = sin(qJ(1));
	t265 = cos(qJ(2));
	t266 = cos(qJ(1));
	t344 = cos(pkin(6));
	t296 = t266 * t344;
	t281 = -t261 * t296 - t262 * t265;
	t297 = t262 * t344;
	t283 = t266 * t261 + t265 * t297;
	t258 = sin(pkin(6));
	t324 = t258 * t266;
	t357 = t283 * qJD(1) - t281 * qJD(2) + qJD(4) * t324;
	t249 = t262 * t261 - t265 * t296;
	t260 = sin(qJ(4));
	t264 = cos(qJ(4));
	t285 = t249 * t264 + t260 * t324;
	t232 = t285 ^ 2;
	t325 = t258 * t265;
	t282 = -t344 * t260 - t264 * t325;
	t245 = 0.1e1 / t282 ^ 2;
	t223 = t232 * t245 + 0.1e1;
	t221 = 0.1e1 / t223;
	t248 = -t260 * t325 + t344 * t264;
	t318 = qJD(2) * t261;
	t301 = t258 * t318;
	t234 = t248 * qJD(4) - t264 * t301;
	t244 = 0.1e1 / t282;
	t330 = t285 * t245;
	t319 = qJD(1) * t262;
	t345 = t249 * qJD(4) + t258 * t319;
	t355 = -t345 * t260 + t264 * t357;
	t289 = -t234 * t330 - t244 * t355;
	t190 = t289 * t221;
	t224 = atan2(t285, -t282);
	t219 = sin(t224);
	t220 = cos(t224);
	t291 = t219 * t282 + t220 * t285;
	t185 = t190 * t291 + t219 * t355 + t220 * t234;
	t202 = t219 * t285 - t220 * t282;
	t200 = 0.1e1 / t202 ^ 2;
	t356 = t185 * t200;
	t207 = t260 * t357 + t345 * t264;
	t326 = t258 * t264;
	t236 = t283 * t260 + t262 * t326;
	t292 = t261 * t297;
	t321 = t266 * t265;
	t251 = -t292 + t321;
	t259 = sin(qJ(5));
	t263 = cos(qJ(5));
	t215 = t236 * t259 - t251 * t263;
	t354 = 0.2e1 * t215;
	t199 = 0.1e1 / t202;
	t353 = t199 * t356;
	t346 = -t262 * t258 * t260 + t283 * t264;
	t293 = -0.2e1 * t346 * t353;
	t277 = (t344 * qJD(1) + qJD(2)) * t321 - qJD(2) * t292 - t261 * t319;
	t302 = qJD(1) * t324;
	t209 = qJD(4) * t236 + t260 * t302 - t264 * t277;
	t336 = t209 * t200;
	t352 = -t336 + t293;
	t350 = t234 * t245;
	t327 = t258 * t261;
	t305 = t285 * t327;
	t284 = t244 * t281 + t245 * t305;
	t349 = t264 * t284;
	t348 = -t251 * t260 * qJD(5) - t277;
	t328 = t251 * t264;
	t347 = qJD(4) * t328 - t283 * qJD(5);
	t216 = t236 * t263 + t251 * t259;
	t212 = 0.1e1 / t216;
	t213 = 0.1e1 / t216 ^ 2;
	t231 = t346 ^ 2;
	t198 = t231 * t200 + 0.1e1;
	t343 = (-t231 * t353 - t336 * t346) / t198 ^ 2;
	t210 = t346 * qJD(4) + t277 * t260 + t264 * t302;
	t229 = t281 * qJD(1) - t283 * qJD(2);
	t194 = qJD(5) * t216 + t210 * t259 - t229 * t263;
	t211 = t215 ^ 2;
	t205 = t211 * t213 + 0.1e1;
	t335 = t213 * t215;
	t315 = qJD(5) * t215;
	t195 = t210 * t263 + t229 * t259 - t315;
	t339 = t195 * t212 * t213;
	t342 = (t194 * t335 - t211 * t339) / t205 ^ 2;
	t332 = t244 * t350;
	t340 = (t232 * t332 + t330 * t355) / t223 ^ 2;
	t338 = t200 * t346;
	t203 = 0.1e1 / t205;
	t337 = t203 * t213;
	t334 = t219 * t346;
	t333 = t220 * t346;
	t331 = t285 * t244;
	t329 = t285 * t248;
	t323 = t260 * t259;
	t322 = t260 * t263;
	t317 = qJD(2) * t265;
	t316 = qJD(4) * t260;
	t313 = 0.2e1 * t343;
	t312 = -0.2e1 * t342;
	t311 = -0.2e1 * t340;
	t310 = 0.2e1 * t340;
	t308 = t213 * t342;
	t307 = t194 * t337;
	t306 = t215 * t339;
	t295 = t244 * t310;
	t294 = 0.2e1 * t306;
	t286 = -t249 * t260 + t264 * t324;
	t218 = t259 * t281 + t263 * t286;
	t217 = t259 * t286 - t263 * t281;
	t288 = -t259 * t212 + t263 * t335;
	t287 = t244 * t286 + t245 * t329;
	t280 = -t219 + (t220 * t331 + t219) * t221;
	t233 = t282 * qJD(4) + t260 * t301;
	t230 = -qJD(1) * t292 - t262 * t318 + (qJD(2) * t344 + qJD(1)) * t321;
	t226 = t251 * t322 - t283 * t259;
	t196 = 0.1e1 / t198;
	t193 = t221 * t349;
	t192 = t287 * t221;
	t187 = (-t219 * t281 - t220 * t327) * t264 + t291 * t193;
	t186 = -t192 * t291 + t219 * t286 + t220 * t248;
	t184 = t287 * t310 + (-0.2e1 * t329 * t332 + t207 * t244 + (-t233 * t285 - t234 * t286 - t248 * t355) * t245) * t221;
	t182 = t311 * t349 + (-t284 * t316 + (0.2e1 * t305 * t332 - t230 * t244 + (t234 * t281 + (t261 * t355 + t285 * t317) * t258) * t245) * t264) * t221;
	t1 = [t346 * t295 + (t209 * t244 - t346 * t350) * t221, t182, 0, t184, 0, 0; -0.2e1 * t285 * t199 * t343 + (t355 * t199 - t285 * t356 + (t280 * t209 - ((-t190 * t221 * t331 + t311) * t219 + (-t285 * t295 - t190 + (t190 - t289) * t221) * t220) * t346) * t338) * t196 - (t352 * t196 - t338 * t313) * t280 * t346, (-t187 * t338 + t199 * t328) * t313 + ((-t229 * t264 + t251 * t316) * t199 + t352 * t187 + (t328 * t185 + (t182 * t285 + t193 * t355 + (t261 * t316 - t264 * t317) * t258 + (t193 * t282 - t264 * t281) * t190) * t333 + (t281 * t316 + t182 * t282 - t193 * t234 + t230 * t264 + (-t193 * t285 + t261 * t326) * t190) * t334) * t200) * t196, 0, (-t186 * t338 - t199 * t236) * t313 + (t186 * t293 + t210 * t199 + (-t236 * t185 - t186 * t209 + (t184 * t285 - t192 * t355 + t233 + (-t192 * t282 + t286) * t190) * t333 + (t184 * t282 + t192 * t234 - t207 + (t192 * t285 - t248) * t190) * t334) * t200) * t196, 0, 0; 0.2e1 * (-t212 * t217 + t218 * t335) * t342 + ((qJD(5) * t218 - t207 * t259 + t230 * t263) * t212 + t218 * t294 + (-t217 * t195 - (-qJD(5) * t217 - t207 * t263 - t230 * t259) * t215 - t218 * t194) * t213) * t203, (t308 * t354 - t307) * t226 + (-t195 * t337 + t212 * t312) * (t251 * t323 + t263 * t283) + (t226 * t294 + (t212 * t323 - t322 * t335) * t229 + (-t348 * t212 - t347 * t335) * t263 + (t347 * t212 - t348 * t335) * t259) * t203, 0, -t288 * t346 * t312 + (t288 * t209 - ((-qJD(5) * t212 - 0.2e1 * t306) * t263 + (t194 * t263 + (t195 - t315) * t259) * t213) * t346) * t203, t312 + (t307 + (-t203 * t339 - t308) * t215) * t354, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end