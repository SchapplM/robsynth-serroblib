% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR15
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
%   Wie in S6RRRPRR15_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR15_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:25
	% DurationCPUTime: 1.23s
	% Computational Cost: add. (2555->103), mult. (7918->233), div. (442->12), fcn. (10062->13), ass. (0->106)
	t218 = sin(pkin(6));
	t219 = cos(pkin(7));
	t220 = cos(pkin(6));
	t217 = sin(pkin(7));
	t225 = cos(qJ(2));
	t268 = t217 * t225;
	t206 = -t218 * t268 + t220 * t219;
	t203 = 0.1e1 / t206;
	t222 = sin(qJ(2));
	t226 = cos(qJ(1));
	t260 = t226 * t222;
	t223 = sin(qJ(1));
	t261 = t223 * t225;
	t238 = t220 * t260 + t261;
	t267 = t218 * t222;
	t204 = 0.1e1 / t206 ^ 2;
	t259 = t226 * t225;
	t262 = t223 * t222;
	t207 = -t220 * t259 + t262;
	t265 = t218 * t226;
	t242 = -t207 * t217 + t219 * t265;
	t272 = t242 * t204;
	t284 = t217 * (t203 * t238 + t267 * t272);
	t189 = atan2(t242, t206);
	t184 = sin(t189);
	t185 = cos(t189);
	t170 = t184 * t242 + t185 * t206;
	t167 = 0.1e1 / t170;
	t221 = sin(qJ(3));
	t224 = cos(qJ(3));
	t237 = t220 * t262 - t259;
	t239 = t220 * t261 + t260;
	t266 = t218 * t223;
	t252 = t217 * t266;
	t240 = -t219 * t239 + t252;
	t181 = t221 * t240 - t224 * t237;
	t175 = 0.1e1 / t181;
	t168 = 0.1e1 / t170 ^ 2;
	t176 = 0.1e1 / t181 ^ 2;
	t200 = -t217 * t239 - t219 * t266;
	t197 = t200 ^ 2;
	t163 = t197 * t168 + 0.1e1;
	t192 = qJD(1) * t207 + qJD(2) * t237;
	t258 = qJD(1) * t218;
	t249 = t226 * t258;
	t182 = t192 * t217 - t219 * t249;
	t276 = t182 * t168;
	t196 = t242 ^ 2;
	t188 = t196 * t204 + 0.1e1;
	t186 = 0.1e1 / t188;
	t194 = qJD(1) * t239 + qJD(2) * t238;
	t250 = t223 * t258;
	t183 = -t194 * t217 - t219 * t250;
	t257 = qJD(2) * t218;
	t269 = t217 * t222;
	t245 = t257 * t269;
	t244 = t204 * t245;
	t233 = t183 * t203 - t242 * t244;
	t159 = t233 * t186;
	t243 = -t184 * t206 + t185 * t242;
	t155 = t159 * t243 + t184 * t183 + t185 * t245;
	t282 = t155 * t167 * t168;
	t283 = (-t197 * t282 + t200 * t276) / t163 ^ 2;
	t193 = qJD(1) * t238 + qJD(2) * t239;
	t235 = t192 * t219 + t217 * t249;
	t165 = qJD(3) * t181 - t193 * t221 - t224 * t235;
	t263 = t219 * t224;
	t270 = t237 * t221;
	t180 = -t224 * t252 + t239 * t263 - t270;
	t174 = t180 ^ 2;
	t173 = t174 * t176 + 0.1e1;
	t277 = t176 * t180;
	t166 = -t193 * t224 + t235 * t221 + (t224 * t240 + t270) * qJD(3);
	t279 = t166 * t175 * t176;
	t281 = (t165 * t277 - t174 * t279) / t173 ^ 2;
	t205 = t203 * t204;
	t280 = (-t196 * t205 * t245 + t183 * t272) / t188 ^ 2;
	t278 = t168 * t200;
	t275 = t184 * t200;
	t274 = t185 * t200;
	t273 = t242 * t203;
	t271 = t238 * t221;
	t264 = t219 * t221;
	t256 = -0.2e1 * t283;
	t255 = -0.2e1 * t282;
	t254 = 0.2e1 * t281;
	t253 = 0.2e1 * t280;
	t251 = t217 * t265;
	t248 = -0.2e1 * t203 * t280;
	t247 = 0.2e1 * t180 * t279;
	t246 = t217 * t250;
	t241 = t207 * t219 + t251;
	t190 = -t221 * t239 - t237 * t263;
	t191 = -t224 * t239 + t237 * t264;
	t234 = t184 + (t185 * t273 - t184) * t186;
	t179 = t221 * t241 - t224 * t238;
	t216 = t217 ^ 2;
	t195 = qJD(1) * t237 + qJD(2) * t207;
	t178 = -t224 * t241 - t271;
	t171 = 0.1e1 / t173;
	t161 = 0.1e1 / t163;
	t160 = t186 * t284;
	t158 = t234 * t200;
	t156 = (-t184 * t238 + t185 * t267) * t217 - t243 * t160;
	t154 = t253 * t284 + (t195 * t203 * t217 + (-t183 * t204 * t269 + (t204 * t238 * t216 * t222 + (0.2e1 * t205 * t216 * t218 * t222 ^ 2 - t204 * t268) * t242) * qJD(2)) * t218) * t186;
	t1 = [t200 * t248 + (t182 * t203 - t200 * t244) * t186, t154, 0, 0, 0, 0; t242 * t167 * t256 + (t183 * t167 + (-t155 * t242 + t158 * t182) * t168) * t161 + ((t158 * t255 + t234 * t276) * t161 + (t158 * t256 + ((-t159 * t186 * t273 + t253) * t275 + (t242 * t248 + t159 + (-t159 + t233) * t186) * t274) * t161) * t168) * t200, 0.2e1 * (t167 * t217 * t237 - t156 * t278) * t283 + ((t243 * t154 - (-t159 * t170 + t183 * t185) * t160) * t278 + (t200 * t255 + t276) * t156 + (-t193 * t167 + (t237 * t155 + (-t159 * t238 + t225 * t257) * t274 + (t195 + (qJD(2) * t160 - t159) * t267) * t275) * t168) * t217) * t161, 0, 0, 0, 0; (-t175 * t178 + t179 * t277) * t254 + ((-t194 * t263 + t195 * t221 + t224 * t246) * t175 + t179 * t247 + (-t178 * t166 - (t194 * t264 + t195 * t224 - t221 * t246) * t180 - t179 * t165) * t176 + (t179 * t175 - (t207 * t263 + t224 * t251 + t271) * t277) * qJD(3)) * t171, (-t175 * t190 + t191 * t277) * t254 + ((qJD(3) * t191 + t192 * t221 - t193 * t263) * t175 + t191 * t247 + (-t190 * t166 - (-qJD(3) * t190 + t192 * t224 + t193 * t264) * t180 - t191 * t165) * t176) * t171, -0.2e1 * t281 + 0.2e1 * (t165 * t176 * t171 + (-t171 * t279 - t176 * t281) * t180) * t180, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:26
	% DurationCPUTime: 2.42s
	% Computational Cost: add. (6571->146), mult. (20921->276), div. (667->12), fcn. (26154->13), ass. (0->139)
	t265 = cos(pkin(6));
	t268 = cos(qJ(2));
	t355 = sin(qJ(1));
	t315 = t355 * t268;
	t267 = sin(qJ(2));
	t269 = cos(qJ(1));
	t334 = t269 * t267;
	t250 = t265 * t334 + t315;
	t266 = sin(qJ(3));
	t262 = sin(pkin(7));
	t356 = cos(qJ(3));
	t320 = t262 * t356;
	t263 = sin(pkin(6));
	t338 = t263 * t269;
	t299 = t320 * t338;
	t281 = -t250 * t266 - t299;
	t264 = cos(pkin(7));
	t261 = t355 * t267;
	t333 = t269 * t268;
	t364 = -t265 * t333 + t261;
	t291 = t364 * t356;
	t282 = t264 * t291;
	t218 = t282 - t281;
	t215 = t218 ^ 2;
	t316 = t356 * t268;
	t336 = t266 * t267;
	t287 = t264 * t316 - t336;
	t309 = t265 * t320;
	t239 = -t287 * t263 - t309;
	t237 = 0.1e1 / t239 ^ 2;
	t208 = t215 * t237 + 0.1e1;
	t206 = 0.1e1 / t208;
	t296 = t364 * t266;
	t321 = t250 * t356;
	t276 = -t264 * t296 + t321;
	t285 = t265 * t315 + t334;
	t234 = t285 * qJD(1) + t250 * qJD(2);
	t307 = t265 * t261;
	t235 = -qJD(1) * t307 - qJD(2) * t261 + (qJD(2) * t265 + qJD(1)) * t333;
	t339 = t262 * t266;
	t257 = t338 * t339;
	t319 = t263 * t355;
	t308 = t262 * t319;
	t292 = t356 * t308;
	t318 = t264 * t356;
	t280 = -qJD(1) * t292 - qJD(3) * t257 + t234 * t318 + t235 * t266;
	t201 = t276 * qJD(3) + t280;
	t317 = t356 * t267;
	t335 = t266 * t268;
	t288 = t264 * t317 + t335;
	t289 = t264 * t335 + t317;
	t323 = t265 * t339;
	t213 = qJD(3) * t323 + (t288 * qJD(2) + t289 * qJD(3)) * t263;
	t236 = 0.1e1 / t239;
	t343 = t218 * t237;
	t295 = -t201 * t236 + t213 * t343;
	t187 = t295 * t206;
	t209 = atan2(-t218, t239);
	t204 = sin(t209);
	t205 = cos(t209);
	t298 = -t204 * t239 - t205 * t218;
	t183 = t298 * t187 - t204 * t201 + t205 * t213;
	t198 = -t204 * t218 + t205 * t239;
	t196 = 0.1e1 / t198 ^ 2;
	t366 = t183 * t196;
	t365 = t213 * t237;
	t286 = t307 - t333;
	t232 = t364 * qJD(1) + t286 * qJD(2);
	t245 = t285 * t262 + t264 * t319;
	t242 = 0.1e1 / t245 ^ 2;
	t314 = qJD(1) * t338;
	t341 = (-t232 * t262 + t264 * t314) * t242;
	t224 = -t286 * t356 + (-t285 * t264 + t308) * t266;
	t217 = t224 ^ 2;
	t212 = t217 * t242 + 0.1e1;
	t210 = 0.1e1 / t212;
	t241 = 0.1e1 / t245;
	t340 = t241 * t341;
	t342 = t224 * t242;
	t279 = t285 * t356;
	t223 = t264 * t279 - t266 * t286 - t292;
	t233 = t250 * qJD(1) + t285 * qJD(2);
	t200 = -t233 * t356 + (t232 * t264 + t262 * t314) * t266 - t223 * qJD(3);
	t351 = (t200 * t342 - t217 * t340) / t212 ^ 2;
	t363 = 0.2e1 * t210 * t224 * t340 + 0.2e1 * t342 * t351;
	t328 = t241 * t351;
	t362 = -t210 * t341 - 0.2e1 * t328;
	t216 = t223 ^ 2;
	t194 = t216 * t196 + 0.1e1;
	t192 = 0.1e1 / t194;
	t195 = 0.1e1 / t198;
	t199 = -qJD(1) * t299 + t224 * qJD(3) - t232 * t318 - t233 * t266;
	t348 = t199 * t196;
	t353 = t195 * t366;
	t354 = (-t216 * t353 + t223 * t348) / t194 ^ 2;
	t361 = -t192 * t366 - 0.2e1 * t195 * t354;
	t357 = 0.2e1 * t223;
	t312 = t353 * t357;
	t332 = 0.2e1 * t354;
	t349 = t196 * t223;
	t360 = t332 * t349 + (t312 - t348) * t192;
	t306 = qJD(1) * t319;
	t359 = (-t250 * qJD(3) - t234 * t264 + t262 * t306) * t266 - qJD(3) * t299 + t235 * t356;
	t358 = -0.2e1 * t218;
	t345 = t236 * t365;
	t352 = (t201 * t343 - t215 * t345) / t208 ^ 2;
	t350 = t192 * t195;
	t347 = t210 * t241;
	t346 = t210 * t242;
	t344 = t218 * t236;
	t337 = t264 * t266;
	t331 = -0.2e1 * t352;
	t329 = t236 * t352;
	t326 = t192 * t349;
	t324 = t262 * t346;
	t310 = t345 * t358;
	t297 = t264 * t364;
	t220 = -t257 + t276;
	t240 = t289 * t263 + t323;
	t294 = -t220 * t236 + t240 * t343;
	t229 = t250 * t318 - t296;
	t249 = t288 * t263;
	t293 = -t229 * t236 + t249 * t343;
	t290 = -t264 * t336 + t316;
	t284 = -t204 + (t205 * t344 + t204) * t206;
	t283 = t356 * t297;
	t277 = t266 * t297 - t321;
	t230 = -t285 * t266 - t286 * t318;
	t231 = t286 * t337 - t279;
	t226 = (t287 * qJD(2) + t290 * qJD(3)) * t263;
	t214 = qJD(3) * t309 + (t290 * qJD(2) + t287 * qJD(3)) * t263;
	t203 = t235 * t318 - t234 * t266 + (-t250 * t337 - t291) * qJD(3);
	t202 = -qJD(3) * t282 + t359;
	t191 = t293 * t206;
	t189 = t294 * t206;
	t184 = t298 * t189 - t204 * t220 + t205 * t240;
	t182 = t293 * t331 + (t249 * t310 - t203 * t236 + (t201 * t249 + t213 * t229 + t218 * t226) * t237) * t206;
	t181 = t294 * t331 + (t240 * t310 - t202 * t236 + (t201 * t240 + t213 * t220 + t214 * t218) * t237) * t206;
	t1 = [t329 * t357 + (-t199 * t236 + t223 * t365) * t206, t182, t181, 0, 0, 0; (t277 * qJD(3) - t280) * t350 - (t284 * t199 + ((-t187 * t206 * t344 + t331) * t204 + (t329 * t358 - t187 + (t187 - t295) * t206) * t205) * t223) * t326 + t361 * (-t283 + t281) + t360 * t284 * t223, (t231 * qJD(3) + t232 * t266 - t233 * t318) * t350 - ((-t182 * t218 - t191 * t201 + t226 + (-t191 * t239 - t229) * t187) * t205 + (-t182 * t239 - t191 * t213 - t203 + (t191 * t218 - t249) * t187) * t204) * t326 + t361 * t230 + t360 * (t298 * t191 - t204 * t229 + t205 * t249), (t184 * t349 - t195 * t224) * t332 + (t184 * t312 + t200 * t195 + (-t224 * t183 - t184 * t199 + (-(-t181 * t218 - t189 * t201 + t214 + (-t189 * t239 - t220) * t187) * t205 - (-t181 * t239 - t189 * t213 - t202 + (t189 * t218 - t240) * t187) * t204) * t223) * t196) * t192, 0, 0, 0; (qJD(3) * t283 - t359) * t347 - (-t234 * t262 - t264 * t306) * t210 * t342 + t362 * (t257 + t277) + (-t200 * t346 + t363) * (-t262 * t364 + t264 * t338), (-t230 * qJD(3) + t232 * t356 + t233 * t337) * t347 + t233 * t224 * t324 + t362 * t231 - (-t200 * t324 + t363 * t262) * t286, t328 * t357 + (-t199 * t241 + t223 * t341) * t210, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:24
	% EndTime: 2019-10-10 12:18:27
	% DurationCPUTime: 2.93s
	% Computational Cost: add. (8022->178), mult. (25165->336), div. (705->12), fcn. (31370->15), ass. (0->146)
	t342 = sin(qJ(3));
	t346 = cos(qJ(3));
	t340 = cos(pkin(6));
	t343 = sin(qJ(2));
	t344 = sin(qJ(1));
	t396 = t344 * t343;
	t381 = t340 * t396;
	t347 = cos(qJ(2));
	t348 = cos(qJ(1));
	t392 = t348 * t347;
	t358 = t381 - t392;
	t393 = t348 * t343;
	t395 = t344 * t347;
	t359 = t340 * t395 + t393;
	t337 = sin(pkin(7));
	t338 = sin(pkin(6));
	t403 = t338 * t344;
	t383 = t337 * t403;
	t339 = cos(pkin(7));
	t400 = t339 * t346;
	t300 = -t342 * t358 - t346 * t383 + t359 * t400;
	t318 = t337 * t359 + t339 * t403;
	t341 = sin(qJ(5));
	t345 = cos(qJ(5));
	t370 = t300 * t345 - t318 * t341;
	t423 = t370 * qJD(5);
	t394 = t346 * t347;
	t399 = t342 * t343;
	t360 = -t339 * t399 + t394;
	t363 = t339 * t394 - t399;
	t406 = t337 * t340;
	t378 = qJD(3) * t406;
	t291 = t346 * t378 + (t360 * qJD(2) + t363 * qJD(3)) * t338;
	t397 = t343 * t346;
	t398 = t342 * t347;
	t361 = t339 * t398 + t397;
	t316 = t361 * t338 + t342 * t406;
	t313 = 0.1e1 / t316 ^ 2;
	t422 = t291 * t313;
	t325 = t340 * t393 + t395;
	t310 = t359 * qJD(1) + t325 * qJD(2);
	t311 = -qJD(1) * t381 - qJD(2) * t396 + (qJD(2) * t340 + qJD(1)) * t392;
	t324 = -t340 * t392 + t396;
	t401 = t339 * t342;
	t365 = t324 * t401 - t325 * t346;
	t402 = t338 * t348;
	t382 = t337 * t402;
	t373 = qJD(3) * t382;
	t391 = qJD(1) * t338;
	t380 = t344 * t391;
	t374 = t337 * t380;
	t271 = t365 * qJD(3) - t310 * t400 + t346 * t374 + (-t311 + t373) * t342;
	t299 = t342 * t382 + t365;
	t287 = atan2(t299, t316);
	t282 = sin(t287);
	t283 = cos(t287);
	t263 = t282 * t299 + t283 * t316;
	t260 = 0.1e1 / t263;
	t281 = t300 * t341 + t318 * t345;
	t275 = 0.1e1 / t281;
	t312 = 0.1e1 / t316;
	t261 = 0.1e1 / t263 ^ 2;
	t276 = 0.1e1 / t281 ^ 2;
	t421 = 0.2e1 * t299;
	t364 = -t339 * t359 + t383;
	t301 = t364 * t342 - t346 * t358;
	t294 = t301 ^ 2;
	t257 = t294 * t261 + 0.1e1;
	t309 = t325 * qJD(1) + t359 * qJD(2);
	t308 = t324 * qJD(1) + t358 * qJD(2);
	t379 = t348 * t391;
	t357 = t308 * t339 + t337 * t379;
	t389 = qJD(3) * t342;
	t268 = t358 * t389 + t357 * t342 + (t364 * qJD(3) - t309) * t346;
	t414 = t268 * t261;
	t293 = t299 ^ 2;
	t286 = t293 * t313 + 0.1e1;
	t284 = 0.1e1 / t286;
	t329 = t346 * t373;
	t390 = qJD(3) * t324;
	t270 = t342 * t374 + t311 * t346 - t325 * t389 - t329 + (-t310 * t342 - t346 * t390) * t339;
	t408 = t299 * t313;
	t369 = -t270 * t312 - t291 * t408;
	t251 = t369 * t284;
	t372 = -t282 * t316 + t283 * t299;
	t246 = t372 * t251 - t282 * t270 + t283 * t291;
	t419 = t246 * t260 * t261;
	t420 = (-t294 * t419 + t301 * t414) / t257 ^ 2;
	t267 = t301 * qJD(3) - t309 * t342 - t357 * t346;
	t302 = -t308 * t337 + t339 * t379;
	t258 = t281 * qJD(5) - t267 * t345 + t302 * t341;
	t274 = t370 ^ 2;
	t266 = t274 * t276 + 0.1e1;
	t413 = t276 * t370;
	t259 = t267 * t341 + t302 * t345 + t423;
	t416 = t259 * t275 * t276;
	t418 = (-t258 * t413 - t274 * t416) / t266 ^ 2;
	t410 = t312 * t422;
	t417 = (-t270 * t408 - t293 * t410) / t286 ^ 2;
	t415 = t261 * t301;
	t412 = t282 * t301;
	t411 = t283 * t301;
	t409 = t299 * t312;
	t407 = t325 * t342;
	t405 = t337 * t341;
	t404 = t337 * t345;
	t388 = 0.2e1 * t420;
	t387 = 0.2e1 * t419;
	t386 = 0.2e1 * t418;
	t385 = -0.2e1 * t417;
	t384 = t312 * t417;
	t377 = -0.2e1 * t370 * t416;
	t376 = t410 * t421;
	t375 = t301 * t387;
	t298 = -t407 + (-t324 * t339 - t382) * t346;
	t317 = -t324 * t337 + t339 * t402;
	t371 = t298 * t345 - t317 * t341;
	t279 = t298 * t341 + t317 * t345;
	t368 = t345 * t275 - t341 * t413;
	t295 = t324 * t400 + t346 * t382 + t407;
	t315 = t363 * t338 + t346 * t406;
	t367 = t295 * t312 - t315 * t408;
	t305 = -t324 * t346 - t325 * t401;
	t321 = t360 * t338;
	t366 = -t305 * t312 - t321 * t408;
	t306 = -t342 * t359 - t358 * t400;
	t289 = t306 * t341 - t358 * t404;
	t288 = -t306 * t345 - t358 * t405;
	t307 = -t346 * t359 + t358 * t401;
	t362 = -t339 * t397 - t398;
	t356 = -t282 + (-t283 * t409 + t282) * t284;
	t304 = (-t361 * qJD(2) + t362 * qJD(3)) * t338;
	t303 = -t310 * t337 - t339 * t380;
	t290 = -t342 * t378 + (t362 * qJD(2) - t361 * qJD(3)) * t338;
	t273 = -t311 * t401 - t310 * t346 + (t324 * t342 - t325 * t400) * qJD(3);
	t272 = t307 * qJD(3) + t308 * t342 - t309 * t400;
	t264 = 0.1e1 / t266;
	t255 = 0.1e1 / t257;
	t254 = t366 * t284;
	t253 = t367 * t284;
	t250 = t356 * t301;
	t248 = t372 * t254 - t282 * t305 + t283 * t321;
	t247 = t372 * t253 + t282 * t295 + t283 * t315;
	t245 = t366 * t385 + (t321 * t376 - t273 * t312 + (t270 * t321 + t291 * t305 - t299 * t304) * t313) * t284;
	t244 = t367 * t385 + (t315 * t376 - t271 * t312 + (t270 * t315 - t290 * t299 - t291 * t295) * t313) * t284;
	t1 = [0.2e1 * t301 * t384 + (-t268 * t312 + t301 * t422) * t284, t245, t244, 0, 0, 0; -0.2e1 * t299 * t260 * t420 + ((t329 + (t339 * t390 - t311) * t346 + (qJD(3) * t325 + t310 * t339 - t374) * t342) * t260 + (-t299 * t246 - t250 * t268) * t261) * t255 + ((t250 * t387 - t356 * t414) * t255 + (t250 * t388 + (-(t251 * t284 * t409 + t385) * t412 - (t384 * t421 - t251 + (t251 - t369) * t284) * t411) * t255) * t261) * t301, (t248 * t415 - t260 * t307) * t388 + ((-t306 * qJD(3) + t308 * t346 + t309 * t401) * t260 + t248 * t375 + (-t307 * t246 - t248 * t268 - (t245 * t299 - t254 * t270 + t304 + (-t254 * t316 - t305) * t251) * t411 - (-t245 * t316 - t254 * t291 - t273 + (-t254 * t299 - t321) * t251) * t412) * t261) * t255, (t247 * t415 + t260 * t300) * t388 + (t247 * t375 - t267 * t260 + (t300 * t246 - t247 * t268 - (t244 * t299 - t253 * t270 + t290 + (-t253 * t316 + t295) * t251) * t411 - (-t244 * t316 - t253 * t291 - t271 + (-t253 * t299 - t315) * t251) * t412) * t261) * t255, 0, 0, 0; (t275 * t371 - t279 * t413) * t386 + ((t279 * qJD(5) - t271 * t345 + t303 * t341) * t275 + t279 * t377 + (t371 * t259 + (t371 * qJD(5) + t271 * t341 + t303 * t345) * t370 - t279 * t258) * t276) * t264, (-t275 * t288 - t289 * t413) * t386 + ((t289 * qJD(5) - t272 * t345 - t309 * t405) * t275 + t289 * t377 + (-t288 * t259 + (-t288 * qJD(5) + t272 * t341 - t309 * t404) * t370 - t289 * t258) * t276) * t264, t368 * t301 * t386 + (-t368 * t268 + ((qJD(5) * t275 + t377) * t341 + (-t258 * t341 + (t259 + t423) * t345) * t276) * t301) * t264, 0, -0.2e1 * t418 - 0.2e1 * (t258 * t276 * t264 - (-t264 * t416 - t276 * t418) * t370) * t370, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:25
	% EndTime: 2019-10-10 12:18:33
	% DurationCPUTime: 7.74s
	% Computational Cost: add. (19869->249), mult. (59163->466), div. (983->12), fcn. (74815->17), ass. (0->203)
	t462 = sin(qJ(3));
	t572 = cos(pkin(6));
	t573 = sin(qJ(2));
	t510 = t572 * t573;
	t465 = cos(qJ(2));
	t574 = sin(qJ(1));
	t527 = t574 * t465;
	t576 = cos(qJ(1));
	t485 = t510 * t576 + t527;
	t575 = cos(qJ(3));
	t478 = t485 * t575;
	t459 = cos(pkin(7));
	t456 = t574 * t573;
	t511 = t572 * t576;
	t494 = -t465 * t511 + t456;
	t487 = t494 * t459;
	t457 = sin(pkin(7));
	t458 = sin(pkin(6));
	t531 = t458 * t576;
	t518 = t457 * t531;
	t428 = -t478 + (t487 + t518) * t462;
	t484 = t527 * t572 + t573 * t576;
	t440 = qJD(1) * t484 + qJD(2) * t485;
	t495 = t574 * t510;
	t524 = t576 * qJD(1);
	t441 = -qJD(1) * t495 - qJD(2) * t456 + (qJD(2) * t511 + t524) * t465;
	t530 = t458 * t574;
	t516 = t457 * t530;
	t452 = t575 * t516;
	t529 = t459 * t575;
	t379 = qJD(1) * t452 + qJD(3) * t428 - t440 * t529 - t441 * t462;
	t513 = qJD(1) * t530;
	t432 = t440 * t457 + t459 * t513;
	t461 = sin(qJ(5));
	t464 = cos(qJ(5));
	t591 = -t379 * t461 + t432 * t464;
	t590 = -t379 * t464 - t432 * t461;
	t445 = t457 * t494 - t459 * t531;
	t589 = t445 * t461;
	t588 = t445 * t464;
	t446 = t457 * t484 + t459 * t530;
	t479 = t484 * t575;
	t483 = -t465 * t576 + t495;
	t508 = t459 * t479 - t462 * t483 - t452;
	t412 = t446 * t464 + t461 * t508;
	t460 = sin(qJ(6));
	t429 = t462 * (-t459 * t484 + t516) - t483 * t575;
	t463 = cos(qJ(6));
	t556 = t429 * t463;
	t387 = t412 * t460 - t556;
	t587 = 0.2e1 * t387;
	t481 = t485 * t462;
	t496 = t575 * t518;
	t584 = t481 + t496;
	t486 = t494 * t575;
	t473 = t459 * t486 + t584;
	t408 = -t464 * t473 + t589;
	t526 = t573 * t462;
	t528 = t575 * t465;
	t491 = t459 * t528 - t526;
	t523 = t457 * t572;
	t503 = t575 * t523;
	t442 = -t458 * t491 - t503;
	t554 = t457 * t458;
	t534 = t465 * t554;
	t449 = t459 * t572 - t534;
	t505 = t442 * t464 - t449 * t461;
	t394 = atan2(-t408, -t505);
	t389 = sin(t394);
	t390 = cos(t394);
	t364 = -t389 * t408 - t390 * t505;
	t362 = 0.1e1 / t364 ^ 2;
	t411 = t446 * t461 - t464 * t508;
	t405 = t411 ^ 2;
	t360 = t405 * t362 + 0.1e1;
	t438 = t494 * qJD(1) + t483 * qJD(2);
	t514 = t458 * t524;
	t430 = -t438 * t457 + t459 * t514;
	t439 = qJD(1) * t485 + qJD(2) * t484;
	t474 = -qJD(1) * t496 + t429 * qJD(3) - t438 * t529 - t439 * t462;
	t368 = qJD(5) * t412 + t430 * t461 - t464 * t474;
	t565 = t362 * t411;
	t404 = t408 ^ 2;
	t422 = 0.1e1 / t505 ^ 2;
	t393 = t404 * t422 + 0.1e1;
	t391 = 0.1e1 / t393;
	t410 = t461 * t473 + t588;
	t371 = t410 * qJD(5) - t590;
	t515 = t573 * t575;
	t550 = t462 * t465;
	t489 = t459 * t515 + t550;
	t492 = t459 * t550 + t515;
	t512 = t462 * t523;
	t416 = qJD(3) * t512 + (qJD(2) * t489 + qJD(3) * t492) * t458;
	t425 = t442 * t461 + t449 * t464;
	t519 = t573 * t554;
	t504 = qJD(2) * t519;
	t395 = qJD(5) * t425 - t416 * t464 + t461 * t504;
	t421 = 0.1e1 / t505;
	t558 = t408 * t422;
	t501 = t371 * t421 + t395 * t558;
	t351 = t501 * t391;
	t507 = t389 * t505 - t390 * t408;
	t345 = t351 * t507 - t389 * t371 + t390 * t395;
	t361 = 0.1e1 / t364;
	t363 = t361 * t362;
	t570 = t345 * t363;
	t546 = 0.2e1 * (t368 * t565 - t405 * t570) / t360 ^ 2;
	t582 = t395 * t422;
	t443 = t458 * t492 + t512;
	t497 = -t421 * t428 + t443 * t558;
	t581 = t464 * t497;
	t557 = t429 * t461;
	t579 = qJD(6) * t557 + t474;
	t578 = -t462 * (-t440 * t459 + t457 * t513) - t441 * t575;
	t388 = t412 * t463 + t429 * t460;
	t382 = 0.1e1 / t388;
	t383 = 0.1e1 / t388 ^ 2;
	t577 = 0.2e1 * t411;
	t369 = -t411 * qJD(5) + t430 * t464 + t474 * t461;
	t377 = -t439 * t575 + (t438 * t459 + t457 * t514) * t462 - t508 * qJD(3);
	t354 = qJD(6) * t388 + t369 * t460 - t377 * t463;
	t381 = t387 ^ 2;
	t367 = t381 * t383 + 0.1e1;
	t563 = t383 * t387;
	t547 = qJD(6) * t387;
	t355 = t369 * t463 + t377 * t460 - t547;
	t567 = t355 * t382 * t383;
	t569 = (t354 * t563 - t381 * t567) / t367 ^ 2;
	t560 = t421 * t582;
	t568 = (t371 * t558 + t404 * t560) / t393 ^ 2;
	t566 = t362 * t368;
	t365 = 0.1e1 / t367;
	t564 = t365 * t383;
	t562 = t389 * t411;
	t561 = t390 * t411;
	t559 = t408 * t421;
	t555 = t429 * t464;
	t553 = t457 * t461;
	t552 = t457 * t464;
	t551 = t459 * t462;
	t549 = qJD(5) * t461;
	t548 = qJD(5) * t464;
	t545 = -0.2e1 * t569;
	t544 = 0.2e1 * t569;
	t543 = -0.2e1 * t568;
	t542 = t363 * t577;
	t541 = t383 * t569;
	t540 = t421 * t568;
	t539 = t354 * t564;
	t538 = t362 * t562;
	t537 = t362 * t561;
	t536 = t387 * t567;
	t535 = t408 * t560;
	t522 = t345 * t542;
	t521 = 0.2e1 * t536;
	t520 = 0.2e1 * t535;
	t427 = -t487 * t575 - t584;
	t407 = t427 * t461 - t588;
	t386 = t407 * t463 + t428 * t460;
	t385 = t407 * t460 - t428 * t463;
	t435 = -t462 * t484 - t483 * t529;
	t415 = t435 * t461 - t483 * t552;
	t436 = t483 * t551 - t479;
	t402 = t415 * t463 + t436 * t460;
	t401 = t415 * t460 - t436 * t463;
	t506 = t427 * t464 + t589;
	t500 = -t460 * t382 + t463 * t563;
	t499 = t410 * t421 + t425 * t558;
	t434 = t459 * t478 - t462 * t494;
	t482 = t457 * t485;
	t413 = t434 * t464 - t461 * t482;
	t448 = t489 * t458;
	t437 = -t448 * t464 + t461 * t519;
	t498 = -t413 * t421 + t437 * t558;
	t414 = -t435 * t464 - t483 * t553;
	t493 = -t389 + (-t390 * t559 + t389) * t391;
	t490 = -t459 * t526 + t528;
	t480 = -qJD(6) * t508 + t377 * t461 + t429 * t548;
	t417 = qJD(3) * t503 + (qJD(2) * t490 + qJD(3) * t491) * t458;
	t403 = qJD(2) * t461 * t534 + t519 * t548 - (qJD(2) * t491 + qJD(3) * t490) * t458 * t464 + t448 * t549;
	t400 = -t460 * t508 + t461 * t556;
	t398 = -qJD(3) * t435 + t438 * t575 + t439 * t551;
	t397 = qJD(3) * t436 + t438 * t462 - t439 * t529;
	t396 = qJD(5) * t505 + t416 * t461 + t464 * t504;
	t380 = -qJD(3) * t427 + t578;
	t378 = -qJD(3) * t473 - t578;
	t374 = -t441 * t553 - t482 * t548 + (t441 * t529 - t440 * t462 + (-t459 * t481 - t486) * qJD(3)) * t464 - t434 * t549;
	t373 = -qJD(5) * t414 + t397 * t461 - t439 * t552;
	t372 = -t408 * qJD(5) + t591;
	t370 = qJD(5) * t506 - t591;
	t358 = 0.1e1 / t360;
	t357 = t391 * t581;
	t356 = t498 * t391;
	t353 = t499 * t391;
	t350 = t493 * t411;
	t348 = (-t389 * t428 - t390 * t443) * t464 - t507 * t357;
	t347 = t356 * t507 + t389 * t413 + t390 * t437;
	t346 = t353 * t507 - t389 * t410 + t390 * t425;
	t344 = t498 * t543 + (t437 * t520 - t374 * t421 + (t371 * t437 - t395 * t413 + t403 * t408) * t422) * t391;
	t342 = t499 * t543 + (t425 * t520 + t372 * t421 + (t371 * t425 + t395 * t410 + t396 * t408) * t422) * t391;
	t341 = 0.2e1 * t568 * t581 + (t497 * t549 + (-0.2e1 * t443 * t535 - t378 * t421 + (-t371 * t443 + t395 * t428 - t408 * t417) * t422) * t464) * t391;
	t1 = [-t540 * t577 + (t368 * t421 + t411 * t582) * t391, t344, t341, 0, t342, 0; t506 * t361 * t546 + ((qJD(5) * t407 + t590) * t361 + (t345 * t506 - t350 * t368) * t362) * t358 + (t350 * t362 * t546 + (0.2e1 * t350 * t570 - (t351 * t391 * t559 + t543) * t538 - (0.2e1 * t408 * t540 - t351 + (t351 - t501) * t391) * t537 - t493 * t566) * t358) * t411, (t347 * t565 - t361 * t414) * t546 + ((qJD(5) * t415 - t397 * t464 - t439 * t553) * t361 + t347 * t522 + (-t414 * t345 - t347 * t368 - (-t344 * t408 - t356 * t371 + t403 + (t356 * t505 + t413) * t351) * t561 - (t344 * t505 - t356 * t395 + t374 + (t356 * t408 - t437) * t351) * t562) * t362) * t358, (t348 * t565 + t361 * t555) * t546 + (-t348 * t566 + (-t377 * t464 + t429 * t549) * t361 + (t348 * t542 + t362 * t555) * t345 - (t443 * t549 - t341 * t408 + t357 * t371 - t417 * t464 + (-t357 * t505 - t428 * t464) * t351) * t537 - (t428 * t549 + t341 * t505 + t357 * t395 + t378 * t464 + (-t357 * t408 + t443 * t464) * t351) * t538) * t358, 0, (t346 * t565 - t361 * t412) * t546 + (t346 * t522 + t369 * t361 + (-t412 * t345 - t346 * t368 - (-t342 * t408 - t353 * t371 + t396 + (t353 * t505 - t410) * t351) * t561 - (t342 * t505 - t353 * t395 - t372 + (t353 * t408 - t425) * t351) * t562) * t362) * t358, 0; (-t382 * t385 + t386 * t563) * t544 + ((qJD(6) * t386 + t370 * t460 - t380 * t463) * t382 + t386 * t521 + (-t385 * t355 - (-qJD(6) * t385 + t370 * t463 + t380 * t460) * t387 - t386 * t354) * t383) * t365, (-t382 * t401 + t402 * t563) * t544 + ((qJD(6) * t402 + t373 * t460 - t398 * t463) * t382 + t402 * t521 + (-t401 * t355 - (-qJD(6) * t401 + t373 * t463 + t398 * t460) * t387 - t402 * t354) * t383) * t365, (t541 * t587 - t539) * t400 + (-t355 * t564 + t382 * t545) * (t460 * t557 + t463 * t508) + ((t480 * t460 + t579 * t463) * t382 - (-t579 * t460 + t480 * t463) * t563 + t400 * t521) * t365, 0, t500 * t411 * t545 + (t500 * t368 + ((-qJD(6) * t382 - 0.2e1 * t536) * t463 + (t354 * t463 + (t355 - t547) * t460) * t383) * t411) * t365, t545 + (t539 + (-t365 * t567 - t541) * t387) * t587;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end