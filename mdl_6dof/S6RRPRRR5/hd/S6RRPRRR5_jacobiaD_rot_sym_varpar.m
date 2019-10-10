% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR5
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
%   Wie in S6RRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (454->46), mult. (1493->122), div. (133->12), fcn. (1882->11), ass. (0->62)
	t133 = sin(qJ(1));
	t134 = cos(qJ(1));
	t131 = cos(pkin(6));
	t129 = sin(pkin(12));
	t132 = sin(qJ(2));
	t161 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:46
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (3306->95), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->96)
	t216 = sin(qJ(4));
	t218 = cos(qJ(4));
	t271 = sin(pkin(12));
	t273 = cos(pkin(6));
	t241 = t273 * t271;
	t272 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:57:46
	% EndTime: 2019-10-10 10:57:48
	% DurationCPUTime: 2.79s
	% Computational Cost: add. (8421->153), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->130)
	t315 = sin(pkin(12));
	t316 = cos(pkin(12));
	t320 = sin(qJ(2));
	t323 = cos(qJ(2));
	t305 = t320 * t315 - t323 * t316;
	t317 = cos(pkin(6));
	t302 = t305 * t317;
	t297 = qJD(2) * t302;
	t345 = t323 * t315 + t320 * t316;
	t304 = t345 * qJD(2);
	t324 = cos(qJ(1));
	t389 = sin(pkin(6));
	t354 = t324 * t389;
	t303 = t345 * t317;
	t390 = sin(qJ(1));
	t357 = t390 * t303;
	t395 = -t390 * t304 - qJD(1) * t357 + (-qJD(1) * t305 - t297) * t324 - qJD(4) * t354;
	t319 = sin(qJ(4));
	t322 = cos(qJ(4));
	t338 = -t324 * t303 + t390 * t305;
	t275 = -t319 * t338 + t322 * t354;
	t355 = t323 * t389;
	t356 = t320 * t389;
	t332 = -t315 * t355 - t316 * t356;
	t291 = -t317 * t322 - t319 * t332;
	t270 = atan2(-t275, t291);
	t265 = sin(t270);
	t266 = cos(t270);
	t241 = -t265 * t275 + t266 * t291;
	t239 = 0.1e1 / t241 ^ 2;
	t339 = -t324 * t305 - t357;
	t350 = t390 * t389;
	t335 = -t319 * t339 + t322 * t350;
	t274 = t335 ^ 2;
	t237 = t274 * t239 + 0.1e1;
	t281 = t319 * t350 + t322 * t339;
	t331 = t338 * qJD(1) + t390 * t297 - t324 * t304;
	t348 = qJD(1) * t354;
	t245 = t281 * qJD(4) + t319 * t331 - t322 * t348;
	t382 = t245 * t239;
	t273 = t275 ^ 2;
	t289 = 0.1e1 / t291 ^ 2;
	t269 = t273 * t289 + 0.1e1;
	t267 = 0.1e1 / t269;
	t344 = qJD(1) * t350;
	t368 = qJD(4) * t322;
	t247 = t395 * t319 - t322 * t344 - t338 * t368;
	t292 = t317 * t319 - t322 * t332;
	t300 = -t315 * t356 + t316 * t355;
	t296 = t300 * qJD(2);
	t271 = t292 * qJD(4) + t296 * t319;
	t288 = 0.1e1 / t291;
	t376 = t275 * t289;
	t343 = -t247 * t288 + t271 * t376;
	t229 = t343 * t267;
	t346 = -t265 * t291 - t266 * t275;
	t224 = t346 * t229 - t265 * t247 + t266 * t271;
	t238 = 0.1e1 / t241;
	t240 = t238 * t239;
	t387 = t224 * t240;
	t366 = 0.2e1 * (-t274 * t387 - t335 * t382) / t237 ^ 2;
	t394 = t271 * t289;
	t283 = -t324 * t302 - t345 * t390;
	t340 = -t283 * t288 + t300 * t376;
	t393 = t319 * t340;
	t248 = t319 * (qJD(4) * t338 + t344) + t395 * t322;
	t286 = t390 * t302 - t324 * t345;
	t318 = sin(qJ(5));
	t321 = cos(qJ(5));
	t257 = t281 * t321 - t286 * t318;
	t251 = 0.1e1 / t257;
	t252 = 0.1e1 / t257 ^ 2;
	t392 = -0.2e1 * t275;
	t391 = -0.2e1 * t335;
	t246 = t335 * qJD(4) + t319 * t348 + t322 * t331;
	t298 = t317 * t304;
	t337 = t305 * qJD(2);
	t261 = t283 * qJD(1) - t390 * t298 - t324 * t337;
	t232 = t257 * qJD(5) + t246 * t318 - t261 * t321;
	t256 = t281 * t318 + t286 * t321;
	t250 = t256 ^ 2;
	t244 = t250 * t252 + 0.1e1;
	t381 = t252 * t256;
	t367 = qJD(5) * t256;
	t233 = t246 * t321 + t261 * t318 - t367;
	t384 = t233 * t251 * t252;
	t386 = (t232 * t381 - t250 * t384) / t244 ^ 2;
	t378 = t288 * t394;
	t385 = (t247 * t376 - t273 * t378) / t269 ^ 2;
	t383 = t239 * t335;
	t380 = t265 * t335;
	t379 = t266 * t335;
	t377 = t275 * t288;
	t375 = t286 * t319;
	t374 = t286 * t322;
	t373 = t318 * t251;
	t371 = t321 * t256;
	t365 = -0.2e1 * t386;
	t364 = 0.2e1 * t386;
	t363 = -0.2e1 * t385;
	t362 = t240 * t391;
	t361 = t288 * t385;
	t360 = t239 * t380;
	t359 = t239 * t379;
	t358 = t256 * t384;
	t353 = 0.2e1 * t358;
	t352 = t378 * t392;
	t277 = -t319 * t354 - t322 * t338;
	t347 = qJD(5) * t374 - t331;
	t255 = -t277 * t321 + t283 * t318;
	t254 = -t277 * t318 - t283 * t321;
	t342 = t252 * t371 - t373;
	t341 = -t277 * t288 + t292 * t376;
	t336 = -t265 + (t266 * t377 + t265) * t267;
	t333 = -qJD(4) * t375 + qJD(5) * t339 - t261 * t322;
	t295 = t332 * qJD(2);
	t272 = -t291 * qJD(4) + t296 * t322;
	t263 = t286 * qJD(1) - t324 * t298 + t390 * t337;
	t259 = t318 * t339 + t321 * t374;
	t258 = t318 * t374 - t321 * t339;
	t242 = 0.1e1 / t244;
	t235 = 0.1e1 / t237;
	t234 = t267 * t393;
	t231 = t341 * t267;
	t228 = t336 * t335;
	t226 = (-t265 * t283 + t266 * t300) * t319 + t346 * t234;
	t225 = t346 * t231 - t265 * t277 + t266 * t292;
	t223 = t341 * t363 + (t292 * t352 - t248 * t288 + (t247 * t292 + t271 * t277 + t272 * t275) * t289) * t267;
	t221 = t363 * t393 + (t340 * t368 + (t300 * t352 - t263 * t288 + (t247 * t300 + t271 * t283 + t275 * t295) * t289) * t319) * t267;
	t1 = [t361 * t391 + (-t245 * t288 - t335 * t394) * t267, t221, 0, t223, 0, 0; t275 * t238 * t366 + (-t247 * t238 + (t224 * t275 + t228 * t245) * t239) * t235 - (-t228 * t239 * t366 + (-0.2e1 * t228 * t387 + (-t229 * t267 * t377 + t363) * t360 + (t361 * t392 - t229 + (t229 - t343) * t267) * t359 - t336 * t382) * t235) * t335, (-t226 * t383 - t238 * t375) * t366 + (-t226 * t382 + (-t261 * t319 + t286 * t368) * t238 + (t226 * t362 - t239 * t375) * t224 + (t300 * t368 - t221 * t275 - t234 * t247 + t295 * t319 + (-t234 * t291 - t283 * t319) * t229) * t359 + (-t283 * t368 - t221 * t291 - t234 * t271 - t263 * t319 + (t234 * t275 - t300 * t319) * t229) * t360) * t235, 0, (-t225 * t383 - t238 * t281) * t366 + (t225 * t224 * t362 + t246 * t238 + (-t281 * t224 - t225 * t245 + (-t223 * t275 - t231 * t247 + t272 + (-t231 * t291 - t277) * t229) * t379 + (-t223 * t291 - t231 * t271 - t248 + (t231 * t275 - t292) * t229) * t380) * t239) * t235, 0, 0; (-t251 * t254 + t255 * t381) * t364 + ((t255 * qJD(5) - t248 * t318 - t263 * t321) * t251 + t255 * t353 + (-t254 * t233 - (-t254 * qJD(5) - t248 * t321 + t263 * t318) * t256 - t255 * t232) * t252) * t242, (-t251 * t258 + t259 * t381) * t364 + (t259 * t353 + t347 * t251 * t321 + t333 * t373 + (t347 * t256 * t318 - t259 * t232 - t258 * t233 - t333 * t371) * t252) * t242, 0, -t342 * t335 * t365 + (t342 * t245 - ((-qJD(5) * t251 - 0.2e1 * t358) * t321 + (t232 * t321 + (t233 - t367) * t318) * t252) * t335) * t242, t365 + 0.2e1 * (t232 * t252 * t242 + (-t242 * t384 - t252 * t386) * t256) * t256, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:46
	% EndTime: 2019-10-10 10:57:48
	% DurationCPUTime: 2.83s
	% Computational Cost: add. (9451->155), mult. (25311->308), div. (744->12), fcn. (32679->15), ass. (0->136)
	t341 = sin(pkin(12));
	t342 = cos(pkin(12));
	t345 = sin(qJ(2));
	t347 = cos(qJ(2));
	t327 = t345 * t341 - t347 * t342;
	t343 = cos(pkin(6));
	t324 = t327 * t343;
	t319 = qJD(2) * t324;
	t369 = t347 * t341 + t345 * t342;
	t326 = t369 * qJD(2);
	t348 = cos(qJ(1));
	t416 = sin(pkin(6));
	t382 = t348 * t416;
	t325 = t369 * t343;
	t417 = sin(qJ(1));
	t385 = t417 * t325;
	t422 = -t417 * t326 - qJD(1) * t385 + (-qJD(1) * t327 - t319) * t348 - qJD(4) * t382;
	t344 = sin(qJ(4));
	t346 = cos(qJ(4));
	t362 = -t348 * t325 + t417 * t327;
	t297 = -t344 * t362 + t346 * t382;
	t383 = t347 * t416;
	t384 = t345 * t416;
	t356 = -t341 * t383 - t342 * t384;
	t313 = -t343 * t346 - t344 * t356;
	t292 = atan2(-t297, t313);
	t287 = sin(t292);
	t288 = cos(t292);
	t263 = -t287 * t297 + t288 * t313;
	t261 = 0.1e1 / t263 ^ 2;
	t363 = -t348 * t327 - t385;
	t374 = t417 * t416;
	t359 = -t344 * t363 + t346 * t374;
	t296 = t359 ^ 2;
	t259 = t296 * t261 + 0.1e1;
	t303 = t344 * t374 + t346 * t363;
	t355 = t362 * qJD(1) + t417 * t319 - t348 * t326;
	t372 = qJD(1) * t382;
	t267 = t303 * qJD(4) + t344 * t355 - t346 * t372;
	t409 = t267 * t261;
	t295 = t297 ^ 2;
	t311 = 0.1e1 / t313 ^ 2;
	t291 = t295 * t311 + 0.1e1;
	t289 = 0.1e1 / t291;
	t368 = qJD(1) * t374;
	t395 = qJD(4) * t346;
	t269 = t422 * t344 - t346 * t368 - t362 * t395;
	t314 = t343 * t344 - t346 * t356;
	t322 = -t341 * t384 + t342 * t383;
	t318 = t322 * qJD(2);
	t293 = t314 * qJD(4) + t318 * t344;
	t310 = 0.1e1 / t313;
	t403 = t297 * t311;
	t367 = -t269 * t310 + t293 * t403;
	t251 = t367 * t289;
	t370 = -t287 * t313 - t288 * t297;
	t246 = t370 * t251 - t287 * t269 + t288 * t293;
	t260 = 0.1e1 / t263;
	t262 = t260 * t261;
	t414 = t246 * t262;
	t394 = 0.2e1 * (-t296 * t414 - t359 * t409) / t259 ^ 2;
	t421 = t293 * t311;
	t305 = -t348 * t324 - t369 * t417;
	t364 = -t305 * t310 + t322 * t403;
	t420 = t344 * t364;
	t270 = (qJD(4) * t362 + t368) * t344 + t422 * t346;
	t308 = t417 * t324 - t348 * t369;
	t340 = qJ(5) + qJ(6);
	t337 = sin(t340);
	t338 = cos(t340);
	t279 = t303 * t338 - t308 * t337;
	t273 = 0.1e1 / t279;
	t274 = 0.1e1 / t279 ^ 2;
	t419 = -0.2e1 * t297;
	t418 = -0.2e1 * t359;
	t320 = t343 * t326;
	t361 = t327 * qJD(2);
	t283 = t305 * qJD(1) - t417 * t320 - t348 * t361;
	t339 = qJD(5) + qJD(6);
	t377 = t303 * t339 - t283;
	t268 = t359 * qJD(4) + t344 * t372 + t346 * t355;
	t379 = -t308 * t339 + t268;
	t253 = t379 * t337 + t377 * t338;
	t278 = t303 * t337 + t308 * t338;
	t272 = t278 ^ 2;
	t266 = t272 * t274 + 0.1e1;
	t408 = t274 * t278;
	t254 = -t377 * t337 + t379 * t338;
	t411 = t254 * t273 * t274;
	t413 = (t253 * t408 - t272 * t411) / t266 ^ 2;
	t405 = t310 * t421;
	t412 = (t269 * t403 - t295 * t405) / t291 ^ 2;
	t410 = t261 * t359;
	t407 = t287 * t359;
	t406 = t288 * t359;
	t404 = t297 * t310;
	t402 = t308 * t344;
	t401 = t308 * t346;
	t400 = t337 * t273;
	t399 = t338 * t278;
	t393 = -0.2e1 * t413;
	t392 = 0.2e1 * t413;
	t391 = -0.2e1 * t412;
	t390 = t262 * t418;
	t389 = t310 * t412;
	t388 = t261 * t407;
	t387 = t261 * t406;
	t386 = t278 * t411;
	t381 = 0.2e1 * t386;
	t380 = t405 * t419;
	t378 = t305 * t339 - t270;
	t285 = t308 * qJD(1) - t348 * t320 + t417 * t361;
	t299 = -t344 * t382 - t346 * t362;
	t376 = -t299 * t339 - t285;
	t371 = t339 * t401 - t355;
	t366 = t274 * t399 - t400;
	t365 = -t299 * t310 + t314 * t403;
	t360 = -t287 + (t288 * t404 + t287) * t289;
	t358 = -qJD(4) * t402 - t283 * t346 + t339 * t363;
	t317 = t356 * qJD(2);
	t294 = -t313 * qJD(4) + t318 * t346;
	t281 = t337 * t363 + t338 * t401;
	t280 = t337 * t401 - t338 * t363;
	t277 = -t299 * t338 + t305 * t337;
	t276 = -t299 * t337 - t305 * t338;
	t264 = 0.1e1 / t266;
	t257 = 0.1e1 / t259;
	t256 = t289 * t420;
	t255 = t365 * t289;
	t250 = t360 * t359;
	t248 = (-t287 * t305 + t288 * t322) * t344 + t370 * t256;
	t247 = t370 * t255 - t287 * t299 + t288 * t314;
	t245 = t365 * t391 + (t314 * t380 - t270 * t310 + (t269 * t314 + t293 * t299 + t294 * t297) * t311) * t289;
	t243 = t391 * t420 + (t364 * t395 + (t322 * t380 - t285 * t310 + (t269 * t322 + t293 * t305 + t297 * t317) * t311) * t344) * t289;
	t242 = t393 + 0.2e1 * (t253 * t274 * t264 + (-t264 * t411 - t274 * t413) * t278) * t278;
	t1 = [t389 * t418 + (-t267 * t310 - t359 * t421) * t289, t243, 0, t245, 0, 0; t297 * t260 * t394 + (-t269 * t260 + (t246 * t297 + t250 * t267) * t261) * t257 - (-t250 * t261 * t394 + (-0.2e1 * t250 * t414 + (-t251 * t289 * t404 + t391) * t388 + (t389 * t419 - t251 + (t251 - t367) * t289) * t387 - t360 * t409) * t257) * t359, (-t248 * t410 - t260 * t402) * t394 + (-t248 * t409 + (-t283 * t344 + t308 * t395) * t260 + (t248 * t390 - t261 * t402) * t246 + (t322 * t395 - t243 * t297 - t256 * t269 + t317 * t344 + (-t256 * t313 - t305 * t344) * t251) * t387 + (-t305 * t395 - t243 * t313 - t256 * t293 - t285 * t344 + (t256 * t297 - t322 * t344) * t251) * t388) * t257, 0, (-t247 * t410 - t260 * t303) * t394 + (t247 * t246 * t390 + t268 * t260 + (-t303 * t246 - t247 * t267 + (-t245 * t297 - t255 * t269 + t294 + (-t255 * t313 - t299) * t251) * t406 + (-t245 * t313 - t255 * t293 - t270 + (t255 * t297 - t314) * t251) * t407) * t261) * t257, 0, 0; (-t273 * t276 + t277 * t408) * t392 + ((t378 * t337 + t376 * t338) * t273 + t277 * t381 + (-t276 * t254 - (-t376 * t337 + t378 * t338) * t278 - t277 * t253) * t274) * t264, (-t273 * t280 + t281 * t408) * t392 + (t281 * t381 + t371 * t273 * t338 + t358 * t400 + (t278 * t337 * t371 - t281 * t253 - t280 * t254 - t358 * t399) * t274) * t264, 0, -t366 * t359 * t393 + (t366 * t267 - ((-t273 * t339 - 0.2e1 * t386) * t338 + (t253 * t338 + (-t278 * t339 + t254) * t337) * t274) * t359) * t264, t242, t242;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end