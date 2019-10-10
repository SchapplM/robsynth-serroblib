% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR4
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
%   Wie in S6RRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.41s
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.46s
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:43
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (3306->95), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->96)
	t216 = sin(qJ(4));
	t218 = cos(qJ(4));
	t271 = sin(pkin(11));
	t273 = cos(pkin(6));
	t241 = t273 * t271;
	t272 = cos(pkin(11));
	t242 = t273 * t272;
	t274 = sin(qJ(2));
	t275 = cos(qJ(2));
	t206 = t241 * t275 + t242 * t274;
	t208 = t271 * t274 - t272 * t275;
	t217 = sin(qJ(1));
	t219 = cos(qJ(1));
	t238 = t217 * t206 + t219 * t208;
	t215 = sin(pkin(6));
	t256 = t215 * t217;
	t233 = t216 * t238 + t218 * t256;
	t278 = qJD(4) * t233;
	t230 = t275 * t271 + t274 * t272;
	t205 = t230 * t215;
	t196 = qJD(2) * t205;
	t204 = t208 * t215;
	t201 = 0.1e1 / t204 ^ 2;
	t258 = t196 * t201;
	t277 = t274 * t241 - t275 * t242;
	t229 = t208 * qJD(2);
	t188 = -t217 * t230 - t219 * t277;
	t176 = atan2(t188, t204);
	t171 = sin(t176);
	t172 = cos(t176);
	t185 = t188 ^ 2;
	t175 = t185 * t201 + 0.1e1;
	t173 = 0.1e1 / t175;
	t200 = 0.1e1 / t204;
	t260 = t188 * t200;
	t276 = t173 * (t172 * t260 - t171) + t171;
	t160 = t171 * t188 + t172 * t204;
	t157 = 0.1e1 / t160;
	t184 = t216 * t256 - t218 * t238;
	t178 = 0.1e1 / t184;
	t158 = 0.1e1 / t160 ^ 2;
	t179 = 0.1e1 / t184 ^ 2;
	t226 = t217 * t277;
	t191 = -t219 * t230 + t226;
	t186 = t191 ^ 2;
	t156 = t186 * t158 + 0.1e1;
	t199 = t206 * qJD(2);
	t166 = qJD(1) * t188 - t217 * t199 - t219 * t229;
	t264 = t166 * t158;
	t254 = qJD(1) * t219;
	t169 = qJD(1) * t226 - t219 * t199 + t217 * t229 - t230 * t254;
	t237 = t169 * t200 - t188 * t258;
	t151 = t237 * t173;
	t240 = -t171 * t204 + t172 * t188;
	t147 = t151 * t240 + t171 * t169 + t172 * t196;
	t269 = t147 * t157 * t158;
	t270 = (-t186 * t269 - t191 * t264) / t156 ^ 2;
	t239 = -t219 * t206 + t217 * t208;
	t259 = t188 * t205;
	t235 = -t200 * t239 + t201 * t259;
	t152 = t235 * t173;
	t148 = -t152 * t240 + t171 * t239 + t172 * t205;
	t268 = t148 * t191;
	t198 = t277 * qJD(2);
	t207 = t230 * qJD(2);
	t167 = qJD(1) * t239 + t217 * t198 - t219 * t207;
	t248 = t215 * t254;
	t161 = qJD(4) * t184 + t167 * t216 - t218 * t248;
	t177 = t233 ^ 2;
	t165 = t177 * t179 + 0.1e1;
	t261 = t179 * t233;
	t162 = t167 * t218 + t216 * t248 + t278;
	t265 = t162 * t178 * t179;
	t267 = (-t161 * t261 - t177 * t265) / t165 ^ 2;
	t257 = t200 * t258;
	t266 = (t188 * t201 * t169 - t185 * t257) / t175 ^ 2;
	t263 = t171 * t191;
	t262 = t172 * t191;
	t255 = t215 * t219;
	t253 = -0.2e1 * t270;
	t252 = -0.2e1 * t269;
	t251 = 0.2e1 * t267;
	t250 = 0.2e1 * t266;
	t249 = qJD(1) * t256;
	t247 = -0.2e1 * t200 * t266;
	t246 = -0.2e1 * t233 * t265;
	t236 = -t216 * t178 - t218 * t261;
	t234 = -t216 * t239 + t218 * t255;
	t182 = t216 * t255 + t218 * t239;
	t228 = qJD(1) * t238 + t219 * t198 + t217 * t207;
	t197 = t215 * t229;
	t163 = 0.1e1 / t165;
	t154 = 0.1e1 / t156;
	t150 = t276 * t191;
	t146 = t235 * t250 + (0.2e1 * t257 * t259 + t228 * t200 + (-t169 * t205 + t188 * t197 - t196 * t239) * t201) * t173;
	t1 = [t191 * t247 + (-t166 * t200 - t191 * t258) * t173, t146, 0, 0, 0, 0; t188 * t157 * t253 + (t169 * t157 + (-t147 * t188 - t150 * t166) * t158) * t154 + ((t150 * t252 - t276 * t264) * t154 + (t150 * t253 + ((-t151 * t173 * t260 + t250) * t263 + (t188 * t247 + t151 + (-t151 + t237) * t173) * t262) * t154) * t158) * t191, 0.2e1 * (t157 * t238 - t158 * t268) * t270 + (t167 * t157 + t252 * t268 + (t238 * t147 - t148 * t166 + (t146 * t188 - t152 * t169 - t197 + (t152 * t204 + t239) * t151) * t262 + (-t146 * t204 + t152 * t196 + t228 + (t152 * t188 - t205) * t151) * t263) * t158) * t154, 0, 0, 0, 0; (t178 * t234 - t182 * t261) * t251 + ((qJD(4) * t182 + t216 * t228 + t218 * t249) * t178 + t182 * t246 + (t234 * t162 + (qJD(4) * t234 - t216 * t249 + t218 * t228) * t233 - t182 * t161) * t179) * t163, t236 * t191 * t251 + (t236 * t166 + ((qJD(4) * t178 + t246) * t218 + (-t161 * t218 + (-t162 - t278) * t216) * t179) * t191) * t163, 0, -0.2e1 * t267 - 0.2e1 * (t161 * t179 * t163 - (-t163 * t265 - t179 * t267) * t233) * t233, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:43
	% DurationCPUTime: 1.28s
	% Computational Cost: add. (3615->96), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->96)
	t227 = qJ(4) + pkin(12);
	t225 = sin(t227);
	t226 = cos(t227);
	t281 = sin(pkin(11));
	t283 = cos(pkin(6));
	t252 = t283 * t281;
	t282 = cos(pkin(11));
	t253 = t283 * t282;
	t284 = sin(qJ(2));
	t285 = cos(qJ(2));
	t216 = t285 * t252 + t284 * t253;
	t218 = t284 * t281 - t285 * t282;
	t229 = sin(qJ(1));
	t230 = cos(qJ(1));
	t249 = t229 * t216 + t218 * t230;
	t228 = sin(pkin(6));
	t267 = t228 * t229;
	t244 = t225 * t249 + t226 * t267;
	t288 = t244 * qJD(4);
	t241 = t285 * t281 + t284 * t282;
	t215 = t241 * t228;
	t206 = qJD(2) * t215;
	t214 = t218 * t228;
	t211 = 0.1e1 / t214 ^ 2;
	t269 = t206 * t211;
	t287 = t284 * t252 - t285 * t253;
	t240 = qJD(2) * t218;
	t198 = -t229 * t241 - t230 * t287;
	t186 = atan2(t198, t214);
	t181 = sin(t186);
	t182 = cos(t186);
	t195 = t198 ^ 2;
	t185 = t195 * t211 + 0.1e1;
	t183 = 0.1e1 / t185;
	t210 = 0.1e1 / t214;
	t271 = t198 * t210;
	t286 = (t182 * t271 - t181) * t183 + t181;
	t170 = t181 * t198 + t182 * t214;
	t167 = 0.1e1 / t170;
	t194 = t225 * t267 - t226 * t249;
	t188 = 0.1e1 / t194;
	t168 = 0.1e1 / t170 ^ 2;
	t189 = 0.1e1 / t194 ^ 2;
	t237 = t229 * t287;
	t201 = -t230 * t241 + t237;
	t196 = t201 ^ 2;
	t166 = t168 * t196 + 0.1e1;
	t209 = t216 * qJD(2);
	t176 = t198 * qJD(1) - t229 * t209 - t230 * t240;
	t276 = t168 * t201;
	t265 = qJD(1) * t230;
	t179 = qJD(1) * t237 - t230 * t209 + t229 * t240 - t241 * t265;
	t248 = t179 * t210 - t198 * t269;
	t161 = t248 * t183;
	t251 = -t181 * t214 + t182 * t198;
	t157 = t251 * t161 + t181 * t179 + t182 * t206;
	t279 = t157 * t167 * t168;
	t280 = (-t176 * t276 - t196 * t279) / t166 ^ 2;
	t208 = t287 * qJD(2);
	t217 = t241 * qJD(2);
	t250 = -t216 * t230 + t229 * t218;
	t177 = t250 * qJD(1) + t229 * t208 - t217 * t230;
	t259 = t228 * t265;
	t171 = t194 * qJD(4) + t177 * t225 - t226 * t259;
	t187 = t244 ^ 2;
	t175 = t187 * t189 + 0.1e1;
	t272 = t189 * t244;
	t172 = t177 * t226 + t225 * t259 + t288;
	t275 = t172 * t188 * t189;
	t278 = (-t171 * t272 - t187 * t275) / t175 ^ 2;
	t268 = t210 * t269;
	t277 = (t179 * t198 * t211 - t195 * t268) / t185 ^ 2;
	t274 = t181 * t201;
	t273 = t182 * t201;
	t270 = t198 * t215;
	t266 = t228 * t230;
	t264 = -0.2e1 * t280;
	t263 = -0.2e1 * t279;
	t262 = 0.2e1 * t278;
	t261 = 0.2e1 * t277;
	t260 = qJD(1) * t267;
	t258 = -0.2e1 * t210 * t277;
	t257 = -0.2e1 * t244 * t275;
	t247 = -t225 * t188 - t226 * t272;
	t246 = -t210 * t250 + t211 * t270;
	t245 = -t225 * t250 + t226 * t266;
	t192 = t225 * t266 + t226 * t250;
	t239 = t249 * qJD(1) + t208 * t230 + t229 * t217;
	t207 = t228 * t240;
	t173 = 0.1e1 / t175;
	t164 = 0.1e1 / t166;
	t162 = t246 * t183;
	t160 = t286 * t201;
	t158 = -t251 * t162 + t181 * t250 + t182 * t215;
	t156 = t246 * t261 + (0.2e1 * t268 * t270 + t239 * t210 + (-t179 * t215 + t198 * t207 - t206 * t250) * t211) * t183;
	t1 = [t201 * t258 + (-t176 * t210 - t201 * t269) * t183, t156, 0, 0, 0, 0; t198 * t167 * t264 + (t179 * t167 + (-t157 * t198 - t160 * t176) * t168) * t164 + (t160 * t263 * t164 + (t160 * t264 + ((-t161 * t183 * t271 + t261) * t274 + (t198 * t258 + t161 + (-t161 + t248) * t183) * t273 - t286 * t176) * t164) * t168) * t201, 0.2e1 * (-t158 * t276 + t167 * t249) * t280 + (t177 * t167 + t158 * t201 * t263 + (t249 * t157 - t158 * t176 + (t156 * t198 - t162 * t179 - t207 + (t162 * t214 + t250) * t161) * t273 + (-t156 * t214 + t162 * t206 + t239 + (t162 * t198 - t215) * t161) * t274) * t168) * t164, 0, 0, 0, 0; (t188 * t245 - t192 * t272) * t262 + ((t192 * qJD(4) + t225 * t239 + t226 * t260) * t188 + t192 * t257 + (t245 * t172 + (t245 * qJD(4) - t225 * t260 + t226 * t239) * t244 - t192 * t171) * t189) * t173, t247 * t201 * t262 + (t247 * t176 + ((qJD(4) * t188 + t257) * t226 + (-t171 * t226 + (-t172 - t288) * t225) * t189) * t201) * t173, 0, -0.2e1 * t278 - 0.2e1 * (t171 * t189 * t173 - (-t173 * t275 - t189 * t278) * t244) * t244, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:45
	% DurationCPUTime: 2.87s
	% Computational Cost: add. (12327->154), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->132)
	t327 = sin(pkin(11));
	t328 = cos(pkin(11));
	t331 = sin(qJ(2));
	t333 = cos(qJ(2));
	t314 = t331 * t327 - t333 * t328;
	t329 = cos(pkin(6));
	t311 = t314 * t329;
	t306 = qJD(2) * t311;
	t355 = t333 * t327 + t331 * t328;
	t313 = t355 * qJD(2);
	t334 = cos(qJ(1));
	t400 = sin(pkin(6));
	t364 = t334 * t400;
	t312 = t355 * t329;
	t401 = sin(qJ(1));
	t367 = t401 * t312;
	t406 = -t401 * t313 - qJD(1) * t367 + (-qJD(1) * t314 - t306) * t334 - qJD(4) * t364;
	t326 = qJ(4) + pkin(12);
	t324 = sin(t326);
	t325 = cos(t326);
	t348 = -t334 * t312 + t401 * t314;
	t284 = -t324 * t348 + t325 * t364;
	t365 = t333 * t400;
	t366 = t331 * t400;
	t342 = -t327 * t365 - t328 * t366;
	t300 = -t324 * t342 - t329 * t325;
	t272 = atan2(-t284, t300);
	t267 = sin(t272);
	t268 = cos(t272);
	t250 = -t267 * t284 + t268 * t300;
	t248 = 0.1e1 / t250 ^ 2;
	t349 = -t334 * t314 - t367;
	t360 = t401 * t400;
	t345 = -t324 * t349 + t325 * t360;
	t283 = t345 ^ 2;
	t246 = t283 * t248 + 0.1e1;
	t290 = t324 * t360 + t325 * t349;
	t341 = t348 * qJD(1) + t401 * t306 - t334 * t313;
	t358 = qJD(1) * t364;
	t254 = t290 * qJD(4) + t324 * t341 - t325 * t358;
	t393 = t254 * t248;
	t282 = t284 ^ 2;
	t298 = 0.1e1 / t300 ^ 2;
	t271 = t282 * t298 + 0.1e1;
	t269 = 0.1e1 / t271;
	t354 = qJD(1) * t360;
	t378 = qJD(4) * t325;
	t256 = t406 * t324 - t325 * t354 - t348 * t378;
	t301 = t329 * t324 - t325 * t342;
	t309 = -t327 * t366 + t328 * t365;
	t305 = t309 * qJD(2);
	t280 = t301 * qJD(4) + t305 * t324;
	t297 = 0.1e1 / t300;
	t387 = t284 * t298;
	t353 = -t256 * t297 + t280 * t387;
	t238 = t353 * t269;
	t356 = -t267 * t300 - t268 * t284;
	t233 = t356 * t238 - t267 * t256 + t268 * t280;
	t247 = 0.1e1 / t250;
	t249 = t247 * t248;
	t398 = t233 * t249;
	t376 = 0.2e1 * (-t283 * t398 - t345 * t393) / t246 ^ 2;
	t405 = t280 * t298;
	t292 = -t334 * t311 - t355 * t401;
	t350 = -t292 * t297 + t309 * t387;
	t404 = t324 * t350;
	t257 = (qJD(4) * t348 + t354) * t324 + t406 * t325;
	t332 = cos(qJ(6));
	t295 = t401 * t311 - t334 * t355;
	t330 = sin(qJ(6));
	t385 = t295 * t330;
	t266 = t290 * t332 - t385;
	t260 = 0.1e1 / t266;
	t261 = 0.1e1 / t266 ^ 2;
	t403 = -0.2e1 * t284;
	t402 = -0.2e1 * t345;
	t255 = t345 * qJD(4) + t324 * t358 + t325 * t341;
	t307 = t329 * t313;
	t347 = t314 * qJD(2);
	t276 = t292 * qJD(1) - t401 * t307 - t334 * t347;
	t242 = t266 * qJD(6) + t255 * t330 - t276 * t332;
	t384 = t295 * t332;
	t265 = t290 * t330 + t384;
	t259 = t265 ^ 2;
	t253 = t259 * t261 + 0.1e1;
	t392 = t261 * t265;
	t377 = qJD(6) * t265;
	t243 = t255 * t332 + t276 * t330 - t377;
	t395 = t243 * t260 * t261;
	t397 = (t242 * t392 - t259 * t395) / t253 ^ 2;
	t389 = t297 * t405;
	t396 = (t256 * t387 - t282 * t389) / t271 ^ 2;
	t394 = t248 * t345;
	t391 = t267 * t345;
	t390 = t268 * t345;
	t388 = t284 * t297;
	t386 = t295 * t324;
	t383 = t330 * t260;
	t381 = t332 * t265;
	t375 = -0.2e1 * t397;
	t374 = 0.2e1 * t397;
	t373 = -0.2e1 * t396;
	t372 = t249 * t402;
	t371 = t297 * t396;
	t370 = t248 * t391;
	t369 = t248 * t390;
	t368 = t265 * t395;
	t363 = 0.2e1 * t368;
	t362 = t389 * t403;
	t286 = -t324 * t364 - t325 * t348;
	t357 = qJD(6) * t295 * t325 - t341;
	t264 = -t286 * t332 + t292 * t330;
	t263 = -t286 * t330 - t292 * t332;
	t352 = t261 * t381 - t383;
	t351 = -t286 * t297 + t301 * t387;
	t346 = -t267 + (t268 * t388 + t267) * t269;
	t343 = -qJD(4) * t386 + qJD(6) * t349 - t276 * t325;
	t304 = t342 * qJD(2);
	t281 = -t300 * qJD(4) + t305 * t325;
	t278 = t295 * qJD(1) - t334 * t307 + t401 * t347;
	t274 = t325 * t384 + t330 * t349;
	t273 = t325 * t385 - t332 * t349;
	t251 = 0.1e1 / t253;
	t244 = 0.1e1 / t246;
	t241 = t269 * t404;
	t240 = t351 * t269;
	t237 = t346 * t345;
	t235 = (-t267 * t292 + t268 * t309) * t324 + t356 * t241;
	t234 = t356 * t240 - t267 * t286 + t268 * t301;
	t231 = t351 * t373 + (t301 * t362 - t257 * t297 + (t256 * t301 + t280 * t286 + t281 * t284) * t298) * t269;
	t230 = t373 * t404 + (t350 * t378 + (t309 * t362 - t278 * t297 + (t256 * t309 + t280 * t292 + t284 * t304) * t298) * t324) * t269;
	t1 = [t371 * t402 + (-t254 * t297 - t345 * t405) * t269, t230, 0, t231, 0, 0; t284 * t247 * t376 + (-t256 * t247 + (t233 * t284 + t237 * t254) * t248) * t244 - (-t237 * t248 * t376 + (-0.2e1 * t237 * t398 + (-t238 * t269 * t388 + t373) * t370 + (t371 * t403 - t238 + (t238 - t353) * t269) * t369 - t346 * t393) * t244) * t345, (-t235 * t394 - t247 * t386) * t376 + (-t235 * t393 + (-t276 * t324 + t295 * t378) * t247 + (t235 * t372 - t248 * t386) * t233 + (t309 * t378 - t230 * t284 - t241 * t256 + t304 * t324 + (-t241 * t300 - t292 * t324) * t238) * t369 + (-t292 * t378 - t230 * t300 - t241 * t280 - t278 * t324 + (t241 * t284 - t309 * t324) * t238) * t370) * t244, 0, (-t234 * t394 - t247 * t290) * t376 + (t234 * t233 * t372 + t255 * t247 + (-t290 * t233 - t234 * t254 + (-t231 * t284 - t240 * t256 + t281 + (-t240 * t300 - t286) * t238) * t390 + (-t231 * t300 - t240 * t280 - t257 + (t240 * t284 - t301) * t238) * t391) * t248) * t244, 0, 0; (-t260 * t263 + t264 * t392) * t374 + ((t264 * qJD(6) - t257 * t330 - t278 * t332) * t260 + t264 * t363 + (-t263 * t243 - (-t263 * qJD(6) - t257 * t332 + t278 * t330) * t265 - t264 * t242) * t261) * t251, (-t260 * t273 + t274 * t392) * t374 + (t274 * t363 + t357 * t260 * t332 + t343 * t383 + (t357 * t265 * t330 - t274 * t242 - t273 * t243 - t343 * t381) * t261) * t251, 0, -t352 * t345 * t375 + (t352 * t254 - ((-qJD(6) * t260 - 0.2e1 * t368) * t332 + (t242 * t332 + (t243 - t377) * t330) * t261) * t345) * t251, 0, t375 + 0.2e1 * (t242 * t261 * t251 + (-t251 * t395 - t261 * t397) * t265) * t265;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end