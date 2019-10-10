% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR3
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
%   Wie in S6RRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.42s
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
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (454->46), mult. (1493->122), div. (133->12), fcn. (1882->11), ass. (0->62)
	t133 = sin(qJ(1));
	t134 = cos(qJ(1));
	t131 = cos(pkin(6));
	t129 = sin(pkin(11));
	t132 = sin(qJ(2));
	t161 = cos(pkin(11));
	t166 = cos(qJ(2));
	t145 = -t132 * t129 + t166 * t161;
	t144 = t145 * t131;
	t146 = t129 * t166 + t132 * t161;
	t99 = -t133 * t144 - t134 * t146;
	t92 = t99 ^ 2;
	t109 = t146 * t131;
	t147 = t133 * t109 - t134 * t145;
	t94 = 0.1e1 / t147 ^ 2;
	t169 = t92 * t94;
	t125 = 0.1e1 / t131 ^ 2;
	t130 = sin(pkin(6));
	t123 = t130 ^ 2;
	t128 = t134 ^ 2;
	t119 = t128 * t123 * t125 + 0.1e1;
	t127 = t133 ^ 2;
	t159 = 0.1e1 / t119 ^ 2 * t127;
	t168 = t125 * t159;
	t167 = qJD(1) * t144 + t145 * qJD(2);
	t93 = 0.1e1 / t147;
	t155 = t134 * t130;
	t118 = atan2(t155, t131);
	t114 = sin(t118);
	t115 = cos(t118);
	t105 = t114 * t155 + t115 * t131;
	t102 = 0.1e1 / t105;
	t124 = 0.1e1 / t131;
	t103 = 0.1e1 / t105 ^ 2;
	t163 = t94 * t99;
	t111 = t146 * qJD(2);
	t108 = t131 * t111;
	t156 = t133 * t146;
	t83 = -qJD(1) * t156 - t133 * t108 + t167 * t134;
	t151 = t83 * t163;
	t107 = qJD(2) * t144;
	t97 = -t134 * t109 - t133 * t145;
	t84 = qJD(1) * t97 - t133 * t107 - t134 * t111;
	t95 = t93 * t94;
	t162 = t95 * t84;
	t88 = 0.1e1 + t169;
	t165 = (t162 * t92 - t151) / t88 ^ 2;
	t160 = t103 * t133;
	t158 = t123 * t124;
	t154 = qJD(1) * t134;
	t152 = -0.2e1 * t97 * t99;
	t116 = 0.1e1 / t119;
	t150 = (t116 - 0.1e1) * t130;
	t149 = -0.2e1 * t124 * t168;
	t85 = (-t115 * t116 * t134 * t158 + t114 * t150) * t133;
	t122 = t130 * t123;
	t104 = t102 * t103;
	t96 = t134 * t144 - t156;
	t91 = t127 * t123 * t103 + 0.1e1;
	t86 = 0.1e1 / t88;
	t82 = qJD(1) * t85;
	t1 = [(-t116 * t124 * t130 + t122 * t149) * t154, 0, 0, 0, 0, 0; (0.2e1 * (-t102 * t134 + t160 * t85) / t91 ^ 2 * (-t104 * t127 * t82 + t154 * t160) * t123 + ((0.2e1 * t104 * t133 * t85 - t103 * t134) * t82 + (-t133 * t102 + ((-t85 + (-t122 * t168 - t150) * t133 * t114) * t134 - (t128 * t123 ^ 2 * t149 + (-t159 + (0.2e1 * t127 - t128) * t116) * t158) * t133 * t115) * t103) * qJD(1)) / t91) * t130, 0, 0, 0, 0, 0; (t94 * t152 + 0.2e1 * t96 * t93) * t165 + ((-t97 * t83 - t96 * t84) * t94 - (-t134 * t108 - t167 * t133 - t146 * t154) * t93 + (qJD(1) * t147 - t134 * t107 + t133 * t111) * t163 - t152 * t162) * t86, 0.2e1 * (-t147 * t93 - t169) * t165 + (-0.2e1 * t151 + (t147 * t94 + 0.2e1 * t92 * t95 - t93) * t84) * t86, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:20
	% DurationCPUTime: 1.18s
	% Computational Cost: add. (2974->87), mult. (8937->192), div. (424->12), fcn. (11581->13), ass. (0->95)
	t184 = sin(pkin(6));
	t237 = sin(pkin(11));
	t238 = cos(pkin(11));
	t240 = sin(qJ(2));
	t241 = cos(qJ(2));
	t198 = t241 * t237 + t240 * t238;
	t174 = t198 * t184;
	t165 = qJD(2) * t174;
	t177 = t240 * t237 - t241 * t238;
	t173 = t177 * t184;
	t170 = 0.1e1 / t173 ^ 2;
	t225 = t165 * t170;
	t239 = cos(pkin(6));
	t206 = t239 * t237;
	t207 = t239 * t238;
	t243 = t240 * t206 - t241 * t207;
	t197 = t177 * qJD(2);
	t186 = sin(qJ(1));
	t187 = cos(qJ(1));
	t157 = -t186 * t198 - t187 * t243;
	t145 = atan2(t157, t173);
	t140 = sin(t145);
	t141 = cos(t145);
	t154 = t157 ^ 2;
	t144 = t154 * t170 + 0.1e1;
	t142 = 0.1e1 / t144;
	t169 = 0.1e1 / t173;
	t227 = t157 * t169;
	t242 = (t141 * t227 - t140) * t142 + t140;
	t129 = t140 * t157 + t141 * t173;
	t126 = 0.1e1 / t129;
	t183 = sin(pkin(12));
	t185 = cos(pkin(12));
	t175 = t241 * t206 + t240 * t207;
	t203 = t186 * t175 + t187 * t177;
	t222 = t184 * t186;
	t153 = t183 * t222 - t185 * t203;
	t147 = 0.1e1 / t153;
	t127 = 0.1e1 / t129 ^ 2;
	t148 = 0.1e1 / t153 ^ 2;
	t194 = t186 * t243;
	t160 = -t187 * t198 + t194;
	t155 = t160 ^ 2;
	t125 = t155 * t127 + 0.1e1;
	t168 = t175 * qJD(2);
	t135 = t157 * qJD(1) - t186 * t168 - t187 * t197;
	t231 = t135 * t127;
	t219 = qJD(1) * t187;
	t138 = qJD(1) * t194 - t187 * t168 + t186 * t197 - t198 * t219;
	t202 = t138 * t169 - t157 * t225;
	t120 = t202 * t142;
	t205 = -t140 * t173 + t141 * t157;
	t116 = t205 * t120 + t140 * t138 + t141 * t165;
	t235 = t116 * t126 * t127;
	t236 = (-t155 * t235 - t160 * t231) / t125 ^ 2;
	t204 = -t187 * t175 + t186 * t177;
	t226 = t157 * t174;
	t201 = -t169 * t204 + t170 * t226;
	t121 = t201 * t142;
	t117 = -t205 * t121 + t140 * t204 + t141 * t174;
	t234 = t117 * t160;
	t224 = t169 * t225;
	t233 = (t157 * t170 * t138 - t154 * t224) / t144 ^ 2;
	t167 = t243 * qJD(2);
	t176 = t198 * qJD(2);
	t136 = t204 * qJD(1) + t186 * t167 - t187 * t176;
	t213 = t184 * t219;
	t134 = t136 * t185 + t183 * t213;
	t232 = t134 * t147 * t148;
	t230 = t140 * t160;
	t229 = t141 * t160;
	t152 = -t183 * t203 - t185 * t222;
	t228 = t148 * t152;
	t223 = t183 * t147;
	t221 = t184 * t187;
	t220 = t185 * t152;
	t218 = -0.2e1 * t236;
	t217 = -0.2e1 * t235;
	t146 = t152 ^ 2;
	t132 = t146 * t148 + 0.1e1;
	t133 = t136 * t183 - t185 * t213;
	t216 = 0.2e1 * (t133 * t228 - t146 * t232) / t132 ^ 2;
	t215 = 0.2e1 * t233;
	t214 = qJD(1) * t222;
	t212 = -0.2e1 * t169 * t233;
	t211 = 0.2e1 * t152 * t232;
	t196 = t203 * qJD(1) + t187 * t167 + t186 * t176;
	t166 = t184 * t197;
	t151 = t183 * t221 + t185 * t204;
	t150 = t183 * t204 - t185 * t221;
	t130 = 0.1e1 / t132;
	t123 = 0.1e1 / t125;
	t118 = t242 * t160;
	t115 = t201 * t215 + (0.2e1 * t224 * t226 + t196 * t169 + (-t138 * t174 + t157 * t166 - t165 * t204) * t170) * t142;
	t1 = [t160 * t212 + (-t135 * t169 - t160 * t225) * t142, t115, 0, 0, 0, 0; t157 * t126 * t218 + (t138 * t126 + (-t116 * t157 - t118 * t135) * t127) * t123 + ((t118 * t217 - t242 * t231) * t123 + (t118 * t218 + ((-t120 * t142 * t227 + t215) * t230 + (t157 * t212 + t120 + (-t120 + t202) * t142) * t229) * t123) * t127) * t160, 0.2e1 * (t126 * t203 - t127 * t234) * t236 + (t136 * t126 + t217 * t234 + (t203 * t116 - t117 * t135 + (t115 * t157 - t121 * t138 - t166 + (t121 * t173 + t204) * t120) * t229 + (-t115 * t173 + t121 * t165 + t196 + (t121 * t157 - t174) * t120) * t230) * t127) * t123, 0, 0, 0, 0; (-t147 * t150 + t151 * t228) * t216 + ((t183 * t196 + t185 * t214) * t147 + t151 * t211 + (-t150 * t134 - (-t183 * t214 + t185 * t196) * t152 - t151 * t133) * t148) * t130, (t148 * t220 - t223) * t160 * t216 + (t160 * t185 * t211 - t135 * t223 + (t135 * t220 + (-t133 * t185 - t134 * t183) * t160) * t148) * t130, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:20
	% DurationCPUTime: 1.28s
	% Computational Cost: add. (3615->96), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->97)
	t226 = pkin(12) + qJ(5);
	t224 = sin(t226);
	t225 = cos(t226);
	t281 = sin(pkin(11));
	t283 = cos(pkin(6));
	t251 = t283 * t281;
	t282 = cos(pkin(11));
	t252 = t283 * t282;
	t284 = sin(qJ(2));
	t285 = cos(qJ(2));
	t215 = t285 * t251 + t284 * t252;
	t217 = t284 * t281 - t285 * t282;
	t228 = sin(qJ(1));
	t229 = cos(qJ(1));
	t248 = t228 * t215 + t229 * t217;
	t227 = sin(pkin(6));
	t266 = t227 * t228;
	t243 = t224 * t248 + t225 * t266;
	t288 = t243 * qJD(5);
	t240 = t285 * t281 + t284 * t282;
	t214 = t240 * t227;
	t205 = qJD(2) * t214;
	t213 = t217 * t227;
	t210 = 0.1e1 / t213 ^ 2;
	t268 = t205 * t210;
	t287 = t284 * t251 - t285 * t252;
	t239 = t217 * qJD(2);
	t197 = -t228 * t240 - t229 * t287;
	t185 = atan2(t197, t213);
	t180 = sin(t185);
	t181 = cos(t185);
	t194 = t197 ^ 2;
	t184 = t194 * t210 + 0.1e1;
	t182 = 0.1e1 / t184;
	t209 = 0.1e1 / t213;
	t270 = t197 * t209;
	t286 = (t181 * t270 - t180) * t182 + t180;
	t169 = t180 * t197 + t181 * t213;
	t166 = 0.1e1 / t169;
	t193 = t224 * t266 - t225 * t248;
	t187 = 0.1e1 / t193;
	t167 = 0.1e1 / t169 ^ 2;
	t188 = 0.1e1 / t193 ^ 2;
	t236 = t228 * t287;
	t200 = -t229 * t240 + t236;
	t195 = t200 ^ 2;
	t165 = t195 * t167 + 0.1e1;
	t208 = t215 * qJD(2);
	t175 = t197 * qJD(1) - t228 * t208 - t229 * t239;
	t274 = t175 * t167;
	t264 = qJD(1) * t229;
	t178 = qJD(1) * t236 - t229 * t208 + t228 * t239 - t240 * t264;
	t247 = t178 * t209 - t197 * t268;
	t160 = t247 * t182;
	t250 = -t180 * t213 + t181 * t197;
	t156 = t250 * t160 + t180 * t178 + t181 * t205;
	t279 = t156 * t166 * t167;
	t280 = (-t195 * t279 - t200 * t274) / t165 ^ 2;
	t249 = -t229 * t215 + t228 * t217;
	t269 = t197 * t214;
	t245 = -t209 * t249 + t210 * t269;
	t161 = t245 * t182;
	t157 = -t250 * t161 + t180 * t249 + t181 * t214;
	t278 = t157 * t200;
	t207 = t287 * qJD(2);
	t216 = t240 * qJD(2);
	t176 = t249 * qJD(1) + t228 * t207 - t229 * t216;
	t258 = t227 * t264;
	t170 = t193 * qJD(5) + t176 * t224 - t225 * t258;
	t186 = t243 ^ 2;
	t174 = t186 * t188 + 0.1e1;
	t271 = t188 * t243;
	t171 = t176 * t225 + t224 * t258 + t288;
	t275 = t171 * t187 * t188;
	t277 = (-t170 * t271 - t186 * t275) / t174 ^ 2;
	t267 = t209 * t268;
	t276 = (t197 * t210 * t178 - t194 * t267) / t184 ^ 2;
	t273 = t180 * t200;
	t272 = t181 * t200;
	t265 = t227 * t229;
	t263 = -0.2e1 * t280;
	t262 = -0.2e1 * t279;
	t261 = 0.2e1 * t277;
	t260 = 0.2e1 * t276;
	t259 = qJD(1) * t266;
	t257 = -0.2e1 * t209 * t276;
	t256 = -0.2e1 * t243 * t275;
	t246 = -t224 * t187 - t225 * t271;
	t244 = -t224 * t249 + t225 * t265;
	t191 = t224 * t265 + t225 * t249;
	t238 = t248 * qJD(1) + t229 * t207 + t228 * t216;
	t206 = t227 * t239;
	t172 = 0.1e1 / t174;
	t163 = 0.1e1 / t165;
	t159 = t286 * t200;
	t155 = t245 * t260 + (0.2e1 * t267 * t269 + t238 * t209 + (-t178 * t214 + t197 * t206 - t205 * t249) * t210) * t182;
	t1 = [t200 * t257 + (-t175 * t209 - t200 * t268) * t182, t155, 0, 0, 0, 0; t197 * t166 * t263 + (t178 * t166 + (-t156 * t197 - t159 * t175) * t167) * t163 + ((t159 * t262 - t286 * t274) * t163 + (t159 * t263 + ((-t160 * t182 * t270 + t260) * t273 + (t197 * t257 + t160 + (-t160 + t247) * t182) * t272) * t163) * t167) * t200, 0.2e1 * (t166 * t248 - t167 * t278) * t280 + (t176 * t166 + t262 * t278 + (t248 * t156 - t157 * t175 + (t155 * t197 - t161 * t178 - t206 + (t161 * t213 + t249) * t160) * t272 + (-t155 * t213 + t161 * t205 + t238 + (t161 * t197 - t214) * t160) * t273) * t167) * t163, 0, 0, 0, 0; (t187 * t244 - t191 * t271) * t261 + ((t191 * qJD(5) + t224 * t238 + t225 * t259) * t187 + t191 * t256 + (t244 * t171 + (t244 * qJD(5) - t224 * t259 + t225 * t238) * t243 - t191 * t170) * t188) * t172, t246 * t200 * t261 + (t246 * t175 + ((qJD(5) * t187 + t256) * t225 + (-t170 * t225 + (-t171 - t288) * t224) * t188) * t200) * t172, 0, 0, -0.2e1 * t277 - 0.2e1 * (t170 * t188 * t172 - (-t172 * t275 - t188 * t277) * t243) * t243, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:20
	% EndTime: 2019-10-10 09:39:22
	% DurationCPUTime: 2.74s
	% Computational Cost: add. (12327->154), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->132)
	t326 = sin(pkin(11));
	t327 = cos(pkin(11));
	t330 = sin(qJ(2));
	t332 = cos(qJ(2));
	t313 = t326 * t330 - t332 * t327;
	t328 = cos(pkin(6));
	t310 = t313 * t328;
	t305 = qJD(2) * t310;
	t354 = t326 * t332 + t327 * t330;
	t312 = t354 * qJD(2);
	t333 = cos(qJ(1));
	t399 = sin(pkin(6));
	t363 = t333 * t399;
	t311 = t354 * t328;
	t400 = sin(qJ(1));
	t366 = t400 * t311;
	t405 = -t400 * t312 - qJD(1) * t366 + (-qJD(1) * t313 - t305) * t333 - qJD(5) * t363;
	t325 = pkin(12) + qJ(5);
	t323 = sin(t325);
	t324 = cos(t325);
	t347 = -t333 * t311 + t400 * t313;
	t283 = -t323 * t347 + t324 * t363;
	t364 = t332 * t399;
	t365 = t330 * t399;
	t341 = -t326 * t364 - t327 * t365;
	t299 = -t323 * t341 - t324 * t328;
	t271 = atan2(-t283, t299);
	t266 = sin(t271);
	t267 = cos(t271);
	t249 = -t266 * t283 + t267 * t299;
	t247 = 0.1e1 / t249 ^ 2;
	t348 = -t333 * t313 - t366;
	t359 = t400 * t399;
	t344 = -t323 * t348 + t324 * t359;
	t282 = t344 ^ 2;
	t245 = t247 * t282 + 0.1e1;
	t289 = t323 * t359 + t324 * t348;
	t340 = t347 * qJD(1) + t400 * t305 - t333 * t312;
	t358 = qJD(1) * t363;
	t253 = t289 * qJD(5) + t323 * t340 - t324 * t358;
	t392 = t253 * t247;
	t281 = t283 ^ 2;
	t297 = 0.1e1 / t299 ^ 2;
	t270 = t281 * t297 + 0.1e1;
	t268 = 0.1e1 / t270;
	t353 = qJD(1) * t359;
	t377 = qJD(5) * t324;
	t255 = t405 * t323 - t324 * t353 - t347 * t377;
	t300 = t323 * t328 - t324 * t341;
	t308 = -t326 * t365 + t327 * t364;
	t304 = t308 * qJD(2);
	t279 = t300 * qJD(5) + t304 * t323;
	t296 = 0.1e1 / t299;
	t386 = t283 * t297;
	t352 = -t255 * t296 + t279 * t386;
	t237 = t352 * t268;
	t355 = -t266 * t299 - t267 * t283;
	t232 = t355 * t237 - t255 * t266 + t267 * t279;
	t246 = 0.1e1 / t249;
	t248 = t246 * t247;
	t397 = t232 * t248;
	t375 = 0.2e1 * (-t282 * t397 - t344 * t392) / t245 ^ 2;
	t404 = t279 * t297;
	t291 = -t333 * t310 - t354 * t400;
	t349 = -t291 * t296 + t308 * t386;
	t403 = t323 * t349;
	t256 = (qJD(5) * t347 + t353) * t323 + t405 * t324;
	t331 = cos(qJ(6));
	t294 = t400 * t310 - t333 * t354;
	t329 = sin(qJ(6));
	t384 = t294 * t329;
	t265 = t289 * t331 - t384;
	t259 = 0.1e1 / t265;
	t260 = 0.1e1 / t265 ^ 2;
	t402 = -0.2e1 * t283;
	t401 = -0.2e1 * t344;
	t254 = t344 * qJD(5) + t323 * t358 + t324 * t340;
	t306 = t328 * t312;
	t346 = t313 * qJD(2);
	t275 = t291 * qJD(1) - t400 * t306 - t333 * t346;
	t241 = t265 * qJD(6) + t254 * t329 - t275 * t331;
	t383 = t294 * t331;
	t264 = t289 * t329 + t383;
	t258 = t264 ^ 2;
	t252 = t258 * t260 + 0.1e1;
	t391 = t260 * t264;
	t376 = qJD(6) * t264;
	t242 = t254 * t331 + t275 * t329 - t376;
	t394 = t242 * t259 * t260;
	t396 = (t241 * t391 - t258 * t394) / t252 ^ 2;
	t388 = t296 * t404;
	t395 = (t255 * t386 - t281 * t388) / t270 ^ 2;
	t393 = t247 * t344;
	t390 = t266 * t344;
	t389 = t267 * t344;
	t387 = t283 * t296;
	t385 = t294 * t323;
	t381 = t329 * t259;
	t380 = t331 * t264;
	t374 = -0.2e1 * t396;
	t373 = 0.2e1 * t396;
	t372 = -0.2e1 * t395;
	t371 = t248 * t401;
	t370 = t296 * t395;
	t369 = t247 * t390;
	t368 = t247 * t389;
	t367 = t264 * t394;
	t362 = 0.2e1 * t367;
	t361 = t388 * t402;
	t285 = -t323 * t363 - t324 * t347;
	t356 = qJD(6) * t294 * t324 - t340;
	t263 = -t285 * t331 + t291 * t329;
	t262 = -t285 * t329 - t291 * t331;
	t351 = t260 * t380 - t381;
	t350 = -t285 * t296 + t300 * t386;
	t345 = -t266 + (t267 * t387 + t266) * t268;
	t342 = -qJD(5) * t385 + qJD(6) * t348 - t275 * t324;
	t303 = t341 * qJD(2);
	t280 = -t299 * qJD(5) + t304 * t324;
	t277 = t294 * qJD(1) - t333 * t306 + t400 * t346;
	t273 = t324 * t383 + t329 * t348;
	t272 = t324 * t384 - t331 * t348;
	t250 = 0.1e1 / t252;
	t243 = 0.1e1 / t245;
	t240 = t268 * t403;
	t239 = t350 * t268;
	t236 = t345 * t344;
	t234 = (-t266 * t291 + t267 * t308) * t323 + t355 * t240;
	t233 = t355 * t239 - t266 * t285 + t267 * t300;
	t230 = t350 * t372 + (t300 * t361 - t256 * t296 + (t255 * t300 + t279 * t285 + t280 * t283) * t297) * t268;
	t229 = t372 * t403 + (t349 * t377 + (t308 * t361 - t277 * t296 + (t255 * t308 + t279 * t291 + t283 * t303) * t297) * t323) * t268;
	t1 = [t370 * t401 + (-t253 * t296 - t344 * t404) * t268, t229, 0, 0, t230, 0; t283 * t246 * t375 + (-t255 * t246 + (t232 * t283 + t236 * t253) * t247) * t243 - (-t236 * t247 * t375 + (-0.2e1 * t236 * t397 + (-t237 * t268 * t387 + t372) * t369 + (t370 * t402 - t237 + (t237 - t352) * t268) * t368 - t345 * t392) * t243) * t344, (-t234 * t393 - t246 * t385) * t375 + (-t234 * t392 + (-t275 * t323 + t294 * t377) * t246 + (t234 * t371 - t247 * t385) * t232 + (t308 * t377 - t229 * t283 - t240 * t255 + t303 * t323 + (-t240 * t299 - t291 * t323) * t237) * t368 + (-t291 * t377 - t229 * t299 - t240 * t279 - t277 * t323 + (t240 * t283 - t308 * t323) * t237) * t369) * t243, 0, 0, (-t233 * t393 - t246 * t289) * t375 + (t233 * t232 * t371 + t254 * t246 + (-t289 * t232 - t233 * t253 + (-t230 * t283 - t239 * t255 + t280 + (-t239 * t299 - t285) * t237) * t389 + (-t230 * t299 - t239 * t279 - t256 + (t239 * t283 - t300) * t237) * t390) * t247) * t243, 0; (-t259 * t262 + t263 * t391) * t373 + ((t263 * qJD(6) - t256 * t329 - t277 * t331) * t259 + t263 * t362 + (-t262 * t242 - (-t262 * qJD(6) - t256 * t331 + t277 * t329) * t264 - t263 * t241) * t260) * t250, (-t259 * t272 + t273 * t391) * t373 + (t273 * t362 + t356 * t259 * t331 + t342 * t381 + (t356 * t264 * t329 - t273 * t241 - t272 * t242 - t342 * t380) * t260) * t250, 0, 0, -t351 * t344 * t374 + (t351 * t253 - ((-qJD(6) * t259 - 0.2e1 * t367) * t331 + (t241 * t331 + (t242 - t376) * t329) * t260) * t344) * t250, t374 + 0.2e1 * (t241 * t260 * t250 + (-t250 * t394 - t260 * t396) * t264) * t264;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end