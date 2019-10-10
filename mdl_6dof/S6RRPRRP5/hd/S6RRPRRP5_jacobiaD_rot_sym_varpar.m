% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP5
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
%   Wie in S6RRPRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
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
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:24
	% DurationCPUTime: 0.47s
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
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:25
	% DurationCPUTime: 1.28s
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
	% StartTime: 2019-10-10 10:35:24
	% EndTime: 2019-10-10 10:35:27
	% DurationCPUTime: 2.75s
	% Computational Cost: add. (8421->153), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->130)
	t315 = sin(pkin(11));
	t316 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:35:24
	% EndTime: 2019-10-10 10:35:27
	% DurationCPUTime: 2.64s
	% Computational Cost: add. (8421->153), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->130)
	t317 = sin(pkin(11));
	t318 = cos(pkin(11));
	t322 = sin(qJ(2));
	t325 = cos(qJ(2));
	t307 = t322 * t317 - t325 * t318;
	t319 = cos(pkin(6));
	t304 = t307 * t319;
	t299 = qJD(2) * t304;
	t347 = t325 * t317 + t322 * t318;
	t306 = t347 * qJD(2);
	t326 = cos(qJ(1));
	t391 = sin(pkin(6));
	t356 = t326 * t391;
	t305 = t347 * t319;
	t392 = sin(qJ(1));
	t359 = t392 * t305;
	t397 = -t392 * t306 - qJD(1) * t359 + (-qJD(1) * t307 - t299) * t326 - qJD(4) * t356;
	t321 = sin(qJ(4));
	t324 = cos(qJ(4));
	t340 = -t326 * t305 + t307 * t392;
	t277 = -t321 * t340 + t324 * t356;
	t357 = t325 * t391;
	t358 = t322 * t391;
	t334 = -t317 * t357 - t318 * t358;
	t293 = -t319 * t324 - t321 * t334;
	t272 = atan2(-t277, t293);
	t267 = sin(t272);
	t268 = cos(t272);
	t243 = -t267 * t277 + t268 * t293;
	t241 = 0.1e1 / t243 ^ 2;
	t341 = -t326 * t307 - t359;
	t352 = t392 * t391;
	t337 = -t321 * t341 + t324 * t352;
	t276 = t337 ^ 2;
	t239 = t276 * t241 + 0.1e1;
	t283 = t321 * t352 + t324 * t341;
	t333 = qJD(1) * t340 + t299 * t392 - t326 * t306;
	t350 = qJD(1) * t356;
	t247 = qJD(4) * t283 + t321 * t333 - t324 * t350;
	t384 = t247 * t241;
	t275 = t277 ^ 2;
	t291 = 0.1e1 / t293 ^ 2;
	t271 = t275 * t291 + 0.1e1;
	t269 = 0.1e1 / t271;
	t346 = qJD(1) * t352;
	t370 = qJD(4) * t324;
	t249 = t397 * t321 - t324 * t346 - t340 * t370;
	t294 = t319 * t321 - t324 * t334;
	t302 = -t317 * t358 + t318 * t357;
	t298 = t302 * qJD(2);
	t273 = qJD(4) * t294 + t298 * t321;
	t290 = 0.1e1 / t293;
	t378 = t277 * t291;
	t345 = -t249 * t290 + t273 * t378;
	t231 = t345 * t269;
	t348 = -t267 * t293 - t268 * t277;
	t226 = t231 * t348 - t267 * t249 + t268 * t273;
	t240 = 0.1e1 / t243;
	t242 = t240 * t241;
	t389 = t226 * t242;
	t368 = 0.2e1 * (-t276 * t389 - t337 * t384) / t239 ^ 2;
	t396 = t273 * t291;
	t285 = -t326 * t304 - t347 * t392;
	t342 = -t285 * t290 + t302 * t378;
	t395 = t321 * t342;
	t250 = t321 * (qJD(4) * t340 + t346) + t397 * t324;
	t288 = t304 * t392 - t326 * t347;
	t320 = sin(qJ(5));
	t323 = cos(qJ(5));
	t259 = t283 * t323 - t288 * t320;
	t253 = 0.1e1 / t259;
	t254 = 0.1e1 / t259 ^ 2;
	t394 = -0.2e1 * t277;
	t393 = -0.2e1 * t337;
	t248 = qJD(4) * t337 + t321 * t350 + t324 * t333;
	t300 = t319 * t306;
	t339 = t307 * qJD(2);
	t263 = qJD(1) * t285 - t300 * t392 - t326 * t339;
	t234 = qJD(5) * t259 + t248 * t320 - t263 * t323;
	t258 = t283 * t320 + t288 * t323;
	t252 = t258 ^ 2;
	t246 = t252 * t254 + 0.1e1;
	t383 = t254 * t258;
	t369 = qJD(5) * t258;
	t235 = t248 * t323 + t263 * t320 - t369;
	t386 = t235 * t253 * t254;
	t388 = (t234 * t383 - t252 * t386) / t246 ^ 2;
	t380 = t290 * t396;
	t387 = (t249 * t378 - t275 * t380) / t271 ^ 2;
	t385 = t241 * t337;
	t382 = t267 * t337;
	t381 = t268 * t337;
	t379 = t277 * t290;
	t377 = t288 * t321;
	t376 = t288 * t324;
	t375 = t320 * t253;
	t373 = t323 * t258;
	t367 = -0.2e1 * t388;
	t366 = 0.2e1 * t388;
	t365 = -0.2e1 * t387;
	t364 = t242 * t393;
	t363 = t290 * t387;
	t362 = t241 * t382;
	t361 = t241 * t381;
	t360 = t258 * t386;
	t355 = 0.2e1 * t360;
	t354 = t380 * t394;
	t279 = -t321 * t356 - t324 * t340;
	t349 = qJD(5) * t376 - t333;
	t257 = -t279 * t323 + t285 * t320;
	t256 = -t279 * t320 - t285 * t323;
	t344 = t373 * t254 - t375;
	t343 = -t279 * t290 + t294 * t378;
	t338 = -t267 + (t268 * t379 + t267) * t269;
	t335 = -qJD(4) * t377 + qJD(5) * t341 - t263 * t324;
	t297 = t334 * qJD(2);
	t274 = -qJD(4) * t293 + t298 * t324;
	t265 = qJD(1) * t288 - t326 * t300 + t339 * t392;
	t261 = t320 * t341 + t323 * t376;
	t260 = t320 * t376 - t323 * t341;
	t244 = 0.1e1 / t246;
	t237 = 0.1e1 / t239;
	t236 = t269 * t395;
	t233 = t343 * t269;
	t230 = t338 * t337;
	t228 = (-t267 * t285 + t268 * t302) * t321 + t348 * t236;
	t227 = t233 * t348 - t267 * t279 + t268 * t294;
	t225 = t343 * t365 + (t294 * t354 - t250 * t290 + (t249 * t294 + t273 * t279 + t274 * t277) * t291) * t269;
	t223 = t365 * t395 + (t342 * t370 + (t302 * t354 - t265 * t290 + (t249 * t302 + t273 * t285 + t277 * t297) * t291) * t321) * t269;
	t1 = [t363 * t393 + (-t247 * t290 - t337 * t396) * t269, t223, 0, t225, 0, 0; t277 * t240 * t368 + (-t249 * t240 + (t226 * t277 + t230 * t247) * t241) * t237 - (-t230 * t241 * t368 + (-0.2e1 * t230 * t389 + (-t231 * t269 * t379 + t365) * t362 + (t363 * t394 - t231 + (t231 - t345) * t269) * t361 - t338 * t384) * t237) * t337, (-t228 * t385 - t240 * t377) * t368 + (-t228 * t384 + (-t263 * t321 + t288 * t370) * t240 + (t228 * t364 - t241 * t377) * t226 + (t302 * t370 - t223 * t277 - t236 * t249 + t297 * t321 + (-t236 * t293 - t285 * t321) * t231) * t361 + (-t285 * t370 - t223 * t293 - t236 * t273 - t265 * t321 + (t236 * t277 - t302 * t321) * t231) * t362) * t237, 0, (-t227 * t385 - t240 * t283) * t368 + (t227 * t226 * t364 + t248 * t240 + (-t283 * t226 - t227 * t247 + (-t225 * t277 - t233 * t249 + t274 + (-t233 * t293 - t279) * t231) * t381 + (-t225 * t293 - t233 * t273 - t250 + (t233 * t277 - t294) * t231) * t382) * t241) * t237, 0, 0; (-t253 * t256 + t257 * t383) * t366 + ((qJD(5) * t257 - t250 * t320 - t265 * t323) * t253 + t257 * t355 + (-t256 * t235 - (-qJD(5) * t256 - t250 * t323 + t265 * t320) * t258 - t257 * t234) * t254) * t244, (-t253 * t260 + t261 * t383) * t366 + (t261 * t355 + t349 * t253 * t323 + t335 * t375 + (t258 * t320 * t349 - t261 * t234 - t260 * t235 - t335 * t373) * t254) * t244, 0, -t344 * t337 * t367 + (t344 * t247 - ((-qJD(5) * t253 - 0.2e1 * t360) * t323 + (t234 * t323 + (t235 - t369) * t320) * t254) * t337) * t244, t367 + 0.2e1 * (t234 * t254 * t244 + (-t244 * t386 - t254 * t388) * t258) * t258, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end