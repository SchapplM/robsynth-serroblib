% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR10
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
%   Wie in S6RPRRRR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (137->30), mult. (614->95), div. (108->12), fcn. (792->9), ass. (0->49)
	t86 = sin(pkin(6));
	t79 = t86 ^ 2;
	t88 = cos(pkin(6));
	t81 = 0.1e1 / t88 ^ 2;
	t90 = cos(qJ(1));
	t84 = t90 ^ 2;
	t77 = t79 * t81 * t84 + 0.1e1;
	t89 = sin(qJ(1));
	t83 = t89 ^ 2;
	t108 = 0.1e1 / t77 ^ 2 * t83;
	t112 = t108 * t81;
	t103 = t90 * t86;
	t76 = atan2(t103, t88);
	t72 = sin(t76);
	t73 = cos(t76);
	t58 = t72 * t103 + t73 * t88;
	t55 = 0.1e1 / t58;
	t85 = sin(pkin(13));
	t105 = t89 * t85;
	t87 = cos(pkin(13));
	t99 = t88 * t105 - t87 * t90;
	t65 = 0.1e1 / t99;
	t80 = 0.1e1 / t88;
	t56 = 0.1e1 / t58 ^ 2;
	t66 = 0.1e1 / t99 ^ 2;
	t111 = t56 * t89;
	t104 = t89 * t87;
	t70 = t88 * t104 + t85 * t90;
	t110 = t66 * t70;
	t106 = t88 * t90;
	t69 = -t85 * t106 - t104;
	t109 = t69 * t70;
	t107 = t79 * t80;
	t102 = qJD(1) * t90;
	t74 = 0.1e1 / t77;
	t101 = (t74 - 0.1e1) * t86;
	t100 = -0.2e1 * t80 * t112;
	t68 = t87 * t106 - t105;
	t51 = (-t73 * t74 * t90 * t107 + t72 * t101) * t89;
	t78 = t86 * t79;
	t67 = t65 * t66;
	t64 = t70 ^ 2;
	t63 = t69 * qJD(1);
	t62 = t68 * qJD(1);
	t61 = t64 * t66 + 0.1e1;
	t57 = t55 * t56;
	t54 = t56 * t79 * t83 + 0.1e1;
	t50 = qJD(1) * t51;
	t1 = [(-t74 * t80 * t86 + t78 * t100) * t102, 0, 0, 0, 0, 0; (0.2e1 * (t51 * t111 - t55 * t90) / t54 ^ 2 * (-t50 * t57 * t83 + t102 * t111) * t79 + ((0.2e1 * t51 * t57 * t89 - t56 * t90) * t50 + (-t89 * t55 + ((-t51 + (-t78 * t112 - t101) * t89 * t72) * t90 - (t84 * t79 ^ 2 * t100 + (-t108 + (0.2e1 * t83 - t84) * t74) * t107) * t89 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0, 0, 0; 0.2e1 * (t66 * t109 + t65 * t68) / t61 ^ 2 * (t63 * t64 * t67 + t62 * t110) + (-t69 * t62 * t66 + (-0.2e1 * t67 * t109 - t68 * t66) * t63 + (-t99 * t110 + t70 * t65) * qJD(1)) / t61, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (994->58), mult. (3107->139), div. (132->12), fcn. (4021->13), ass. (0->80)
	t183 = sin(pkin(7));
	t184 = sin(pkin(6));
	t185 = cos(pkin(13));
	t186 = cos(pkin(7));
	t187 = cos(pkin(6));
	t173 = -t184 * t185 * t183 + t187 * t186;
	t191 = cos(qJ(1));
	t211 = t191 * t185;
	t182 = sin(pkin(13));
	t189 = sin(qJ(1));
	t214 = t189 * t182;
	t174 = -t187 * t211 + t214;
	t216 = t184 * t191;
	t204 = -t174 * t183 + t186 * t216;
	t157 = atan2(t204, t173);
	t152 = sin(t157);
	t153 = cos(t157);
	t138 = t152 * t204 + t153 * t173;
	t135 = 0.1e1 / t138;
	t188 = sin(qJ(3));
	t190 = cos(qJ(3));
	t200 = t187 * t214 - t211;
	t212 = t191 * t182;
	t213 = t189 * t185;
	t201 = t187 * t213 + t212;
	t217 = t184 * t189;
	t209 = t183 * t217;
	t202 = -t186 * t201 + t209;
	t151 = t202 * t188 - t190 * t200;
	t145 = 0.1e1 / t151;
	t230 = t204 ^ 2;
	t170 = 0.1e1 / t173;
	t136 = 0.1e1 / t138 ^ 2;
	t146 = 0.1e1 / t151 ^ 2;
	t171 = 0.1e1 / t173 ^ 2;
	t229 = -0.2e1 * t170 * t171;
	t175 = -t187 * t212 - t213;
	t167 = t175 * qJD(1);
	t166 = t174 * qJD(1);
	t210 = qJD(1) * t184;
	t206 = t191 * t210;
	t199 = t166 * t186 + t183 * t206;
	t139 = t151 * qJD(3) + t167 * t188 - t199 * t190;
	t215 = t186 * t190;
	t218 = t200 * t188;
	t150 = -t190 * t209 + t201 * t215 - t218;
	t144 = t150 ^ 2;
	t143 = t144 * t146 + 0.1e1;
	t225 = t146 * t150;
	t140 = t167 * t190 + t199 * t188 + (t202 * t190 + t218) * qJD(3);
	t226 = t140 * t145 * t146;
	t228 = (t139 * t225 - t144 * t226) / t143 ^ 2;
	t168 = t201 * qJD(1);
	t207 = t189 * t210;
	t159 = -t168 * t183 - t186 * t207;
	t156 = t230 * t171 + 0.1e1;
	t154 = 0.1e1 / t156;
	t198 = t152 + (t153 * t170 * t204 - t152) * t154;
	t130 = t198 * t159;
	t227 = t130 * t135 * t136;
	t208 = t183 * t216;
	t203 = t174 * t186 + t208;
	t149 = t175 * t190 + t203 * t188;
	t224 = t149 * t150;
	t223 = t154 * t170;
	t155 = 0.1e1 / t156 ^ 2;
	t222 = t155 * t204;
	t158 = t166 * t183 - t186 * t206;
	t221 = t158 * t136;
	t163 = -t183 * t201 - t186 * t217;
	t220 = t159 * t163;
	t219 = t175 * t188;
	t205 = t183 * t207;
	t169 = t200 * qJD(1);
	t160 = t163 ^ 2;
	t148 = -t203 * t190 + t219;
	t141 = 0.1e1 / t143;
	t134 = t160 * t136 + 0.1e1;
	t131 = t198 * t163;
	t1 = [t220 * t222 * t229 + t158 * t223, 0, 0, 0, 0, 0; 0.2e1 * (-t131 * t136 * t163 - t135 * t204) / t134 ^ 2 * (-t160 * t227 + t163 * t221) + (t159 * t135 + (-t130 * t204 + t131 * t158) * t136 + (-0.2e1 * t131 * t227 + t198 * t221 + (t152 * t171 * t222 + (0.2e1 * t223 + (t230 * t229 - t170) * t155) * t153) * t136 * t220) * t163) / t134, 0, 0, 0, 0, 0; 0.2e1 * (-t145 * t148 + t146 * t224) * t228 + ((-t168 * t215 + t169 * t188 + t190 * t205) * t145 + 0.2e1 * t224 * t226 + (-t148 * t140 - (t169 * t190 + (t168 * t186 - t205) * t188) * t150 - t149 * t139) * t146 + (t149 * t145 - (t174 * t215 + t190 * t208 - t219) * t225) * qJD(3)) * t141, 0, -0.2e1 * t228 + 0.2e1 * (t139 * t146 * t141 + (-t141 * t226 - t146 * t228) * t150) * t150, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:32
	% DurationCPUTime: 1.97s
	% Computational Cost: add. (4837->111), mult. (15324->224), div. (448->12), fcn. (19446->15), ass. (0->108)
	t347 = sin(pkin(13));
	t352 = cos(pkin(6));
	t320 = t352 * t347;
	t350 = cos(pkin(13));
	t353 = sin(qJ(1));
	t355 = cos(qJ(1));
	t276 = t320 * t355 + t350 * t353;
	t287 = sin(qJ(3));
	t348 = sin(pkin(7));
	t349 = sin(pkin(6));
	t318 = t349 * t348;
	t311 = t355 * t318;
	t283 = t287 * t311;
	t322 = t352 * t350;
	t275 = -t322 * t355 + t347 * t353;
	t351 = cos(pkin(7));
	t326 = t275 * t351;
	t354 = cos(qJ(3));
	t362 = -t276 * t354 + t287 * t326 + t283;
	t304 = t354 * t311;
	t324 = t351 * t354;
	t361 = -t275 * t324 - t304;
	t249 = t276 * t287 - t361;
	t317 = t349 * t347;
	t319 = t351 * t349;
	t359 = t350 * t319 + t352 * t348;
	t295 = -t287 * t317 + t359 * t354;
	t246 = atan2(-t249, -t295);
	t241 = sin(t246);
	t242 = cos(t246);
	t222 = -t241 * t249 - t242 * t295;
	t220 = 0.1e1 / t222 ^ 2;
	t277 = -t353 * t320 + t355 * t350;
	t300 = t322 * t353 + t347 * t355;
	t298 = t300 * t351;
	t308 = t353 * t318;
	t303 = t354 * t308;
	t358 = -t277 * t287 - t354 * t298 + t303;
	t248 = t358 ^ 2;
	t218 = t248 * t220 + 0.1e1;
	t255 = t277 * t354 + (-t298 + t308) * t287;
	t272 = t276 * qJD(1);
	t297 = qJD(1) * t326;
	t228 = -qJD(1) * t304 + qJD(3) * t255 - t272 * t287 - t297 * t354;
	t339 = t228 * t220;
	t219 = 0.1e1 / t222;
	t247 = t249 ^ 2;
	t262 = 0.1e1 / t295 ^ 2;
	t245 = t247 * t262 + 0.1e1;
	t243 = 0.1e1 / t245;
	t261 = 0.1e1 / t295;
	t265 = t359 * t287 + t354 * t317;
	t257 = t265 * qJD(3);
	t335 = t257 * t262;
	t273 = t300 * qJD(1);
	t274 = t277 * qJD(1);
	t357 = qJD(1) * t303 + t362 * qJD(3) - t273 * t324 - t274 * t287;
	t314 = t249 * t335 - t261 * t357;
	t213 = t314 * t243;
	t316 = t241 * t295 - t242 * t249;
	t209 = t213 * t316 + t241 * t357 + t242 * t257;
	t360 = t209 * t220;
	t345 = t219 * t360;
	t331 = 0.2e1 * (-t248 * t345 - t339 * t358) / t218 ^ 2;
	t231 = t287 * (qJD(1) * t308 - qJD(3) * t276 - t273 * t351) + t274 * t354 + t361 * qJD(3);
	t309 = t353 * t319;
	t267 = t300 * t348 + t309;
	t286 = sin(qJ(4));
	t288 = cos(qJ(4));
	t240 = t255 * t288 + t267 * t286;
	t234 = 0.1e1 / t240;
	t235 = 0.1e1 / t240 ^ 2;
	t356 = -0.2e1 * t358;
	t229 = qJD(1) * t283 + t358 * qJD(3) - t272 * t354 + t287 * t297;
	t266 = -t275 * t348 + t319 * t355;
	t258 = t266 * qJD(1);
	t223 = qJD(4) * t240 + t229 * t286 - t258 * t288;
	t239 = t255 * t286 - t267 * t288;
	t233 = t239 ^ 2;
	t227 = t233 * t235 + 0.1e1;
	t338 = t235 * t239;
	t332 = qJD(4) * t239;
	t224 = t229 * t288 + t258 * t286 - t332;
	t340 = t224 * t234 * t235;
	t344 = (t223 * t338 - t233 * t340) / t227 ^ 2;
	t334 = t261 * t335;
	t342 = (-t249 * t262 * t357 + t247 * t334) / t245 ^ 2;
	t341 = t220 * t358;
	t337 = t249 * t261;
	t336 = t249 * t265;
	t330 = -0.2e1 * t344;
	t329 = -0.2e1 * t342;
	t328 = t261 * t342;
	t327 = t239 * t340;
	t325 = t345 * t356;
	t238 = t266 * t286 + t288 * t362;
	t237 = -t266 * t288 + t286 * t362;
	t313 = -t286 * t234 + t288 * t338;
	t312 = -t261 * t362 + t262 * t336;
	t306 = -t241 + (-t242 * t337 + t241) * t243;
	t259 = -qJD(1) * t309 - t273 * t348;
	t256 = t295 * qJD(3);
	t225 = 0.1e1 / t227;
	t216 = 0.1e1 / t218;
	t214 = t312 * t243;
	t210 = t214 * t316 + t241 * t362 + t242 * t265;
	t208 = t312 * t329 + (0.2e1 * t334 * t336 + t231 * t261 + (t249 * t256 - t257 * t362 - t265 * t357) * t262) * t243;
	t1 = [-t328 * t356 + (t228 * t261 - t335 * t358) * t243, 0, t208, 0, 0, 0; t249 * t219 * t331 + (t357 * t219 + t249 * t360 + (t306 * t228 - ((t213 * t243 * t337 + t329) * t241 + (0.2e1 * t249 * t328 - t213 + (t213 - t314) * t243) * t242) * t358) * t341) * t216 - (-t341 * t331 + (-t339 + t325) * t216) * t306 * t358, 0, (-t210 * t341 - t219 * t255) * t331 + (t210 * t325 + t229 * t219 + (-t255 * t209 - t210 * t228 - (-(-t208 * t249 + t214 * t357 + t256 + (t214 * t295 + t362) * t213) * t242 - (t208 * t295 - t214 * t257 - t231 + (t214 * t249 - t265) * t213) * t241) * t358) * t220) * t216, 0, 0, 0; 0.2e1 * (-t234 * t237 + t238 * t338) * t344 + ((qJD(4) * t238 - t231 * t286 - t259 * t288) * t234 + 0.2e1 * t238 * t327 + (-t237 * t224 - (-qJD(4) * t237 - t231 * t288 + t259 * t286) * t239 - t238 * t223) * t235) * t225, 0, -t313 * t358 * t330 + (t313 * t228 - ((-qJD(4) * t234 - 0.2e1 * t327) * t288 + (t223 * t288 + (t224 - t332) * t286) * t235) * t358) * t225, t330 + 0.2e1 * (t223 * t235 * t225 + (-t225 * t340 - t235 * t344) * t239) * t239, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:32
	% DurationCPUTime: 2.14s
	% Computational Cost: add. (5639->113), mult. (16379->225), div. (466->12), fcn. (20741->15), ass. (0->112)
	t378 = sin(pkin(13));
	t383 = cos(pkin(6));
	t348 = t383 * t378;
	t381 = cos(pkin(13));
	t384 = sin(qJ(1));
	t386 = cos(qJ(1));
	t302 = t386 * t348 + t384 * t381;
	t316 = sin(qJ(3));
	t379 = sin(pkin(7));
	t380 = sin(pkin(6));
	t346 = t380 * t379;
	t339 = t386 * t346;
	t309 = t316 * t339;
	t350 = t383 * t381;
	t301 = -t386 * t350 + t384 * t378;
	t382 = cos(pkin(7));
	t358 = t301 * t382;
	t385 = cos(qJ(3));
	t393 = -t302 * t385 + t316 * t358 + t309;
	t332 = t385 * t339;
	t352 = t382 * t385;
	t392 = -t301 * t352 - t332;
	t275 = t302 * t316 - t392;
	t345 = t380 * t378;
	t347 = t382 * t380;
	t390 = t381 * t347 + t383 * t379;
	t323 = -t316 * t345 + t390 * t385;
	t272 = atan2(-t275, -t323);
	t267 = sin(t272);
	t268 = cos(t272);
	t250 = -t267 * t275 - t268 * t323;
	t246 = 0.1e1 / t250 ^ 2;
	t303 = -t384 * t348 + t386 * t381;
	t328 = t384 * t350 + t386 * t378;
	t326 = t328 * t382;
	t336 = t384 * t346;
	t331 = t385 * t336;
	t389 = -t303 * t316 - t385 * t326 + t331;
	t274 = t389 ^ 2;
	t244 = t274 * t246 + 0.1e1;
	t281 = t303 * t385 + (-t326 + t336) * t316;
	t298 = t302 * qJD(1);
	t325 = qJD(1) * t358;
	t254 = -qJD(1) * t332 + t281 * qJD(3) - t298 * t316 - t385 * t325;
	t370 = t254 * t246;
	t245 = 0.1e1 / t250;
	t273 = t275 ^ 2;
	t288 = 0.1e1 / t323 ^ 2;
	t271 = t273 * t288 + 0.1e1;
	t269 = 0.1e1 / t271;
	t287 = 0.1e1 / t323;
	t291 = t390 * t316 + t385 * t345;
	t283 = t291 * qJD(3);
	t366 = t283 * t288;
	t299 = t328 * qJD(1);
	t300 = t303 * qJD(1);
	t388 = qJD(1) * t331 + t393 * qJD(3) - t299 * t352 - t300 * t316;
	t342 = t275 * t366 - t287 * t388;
	t239 = t342 * t269;
	t344 = t267 * t323 - t268 * t275;
	t235 = t344 * t239 + t267 * t388 + t268 * t283;
	t391 = t235 * t246;
	t376 = t245 * t391;
	t363 = 0.2e1 * (-t274 * t376 - t370 * t389) / t244 ^ 2;
	t257 = (qJD(1) * t336 - qJD(3) * t302 - t382 * t299) * t316 + t300 * t385 + t392 * qJD(3);
	t337 = t384 * t347;
	t293 = t328 * t379 + t337;
	t315 = qJ(4) + qJ(5);
	t312 = sin(t315);
	t313 = cos(t315);
	t266 = t281 * t313 + t293 * t312;
	t260 = 0.1e1 / t266;
	t261 = 0.1e1 / t266 ^ 2;
	t387 = -0.2e1 * t389;
	t292 = -t301 * t379 + t386 * t347;
	t314 = qJD(4) + qJD(5);
	t353 = -t292 * qJD(1) + t281 * t314;
	t255 = qJD(1) * t309 + t389 * qJD(3) - t298 * t385 + t316 * t325;
	t356 = t293 * t314 + t255;
	t248 = t356 * t312 + t353 * t313;
	t265 = t281 * t312 - t293 * t313;
	t259 = t265 ^ 2;
	t253 = t259 * t261 + 0.1e1;
	t369 = t261 * t265;
	t249 = -t353 * t312 + t356 * t313;
	t371 = t249 * t260 * t261;
	t375 = (t248 * t369 - t259 * t371) / t253 ^ 2;
	t365 = t287 * t366;
	t373 = (-t275 * t288 * t388 + t273 * t365) / t271 ^ 2;
	t372 = t246 * t389;
	t368 = t275 * t287;
	t367 = t275 * t291;
	t362 = -0.2e1 * t375;
	t361 = -0.2e1 * t373;
	t360 = t287 * t373;
	t359 = t265 * t371;
	t357 = t376 * t387;
	t355 = t292 * t314 - t257;
	t354 = qJD(1) * t337 + t299 * t379 + t314 * t393;
	t341 = -t312 * t260 + t313 * t369;
	t340 = -t287 * t393 + t288 * t367;
	t334 = -t267 + (-t268 * t368 + t267) * t269;
	t282 = t323 * qJD(3);
	t264 = t292 * t312 + t313 * t393;
	t263 = -t292 * t313 + t312 * t393;
	t251 = 0.1e1 / t253;
	t242 = 0.1e1 / t244;
	t240 = t340 * t269;
	t236 = t344 * t240 + t267 * t393 + t268 * t291;
	t234 = t340 * t361 + (0.2e1 * t365 * t367 + t257 * t287 + (t275 * t282 - t283 * t393 - t291 * t388) * t288) * t269;
	t232 = t362 + 0.2e1 * (t248 * t261 * t251 + (-t251 * t371 - t261 * t375) * t265) * t265;
	t1 = [-t360 * t387 + (t254 * t287 - t366 * t389) * t269, 0, t234, 0, 0, 0; t275 * t245 * t363 + (t388 * t245 + t275 * t391 + (t334 * t254 - ((t239 * t269 * t368 + t361) * t267 + (0.2e1 * t275 * t360 - t239 + (t239 - t342) * t269) * t268) * t389) * t372) * t242 - (-t372 * t363 + (-t370 + t357) * t242) * t334 * t389, 0, (-t236 * t372 - t245 * t281) * t363 + (t236 * t357 + t255 * t245 + (-t281 * t235 - t236 * t254 - (-(-t234 * t275 + t240 * t388 + t282 + (t240 * t323 + t393) * t239) * t268 - (t234 * t323 - t240 * t283 - t257 + (t240 * t275 - t291) * t239) * t267) * t389) * t246) * t242, 0, 0, 0; 0.2e1 * (-t260 * t263 + t264 * t369) * t375 + ((t355 * t312 + t354 * t313) * t260 + 0.2e1 * t264 * t359 + (-t263 * t249 - (-t354 * t312 + t355 * t313) * t265 - t264 * t248) * t261) * t251, 0, -t341 * t389 * t362 + (t341 * t254 - ((-t260 * t314 - 0.2e1 * t359) * t313 + (t248 * t313 + (-t265 * t314 + t249) * t312) * t261) * t389) * t251, t232, t232, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:31
	% EndTime: 2019-10-10 09:11:38
	% DurationCPUTime: 6.47s
	% Computational Cost: add. (24403->171), mult. (55387->326), div. (989->12), fcn. (71220->17), ass. (0->154)
	t526 = sin(pkin(13));
	t531 = cos(pkin(6));
	t484 = t531 * t526;
	t529 = cos(pkin(13));
	t532 = sin(qJ(1));
	t533 = cos(qJ(1));
	t424 = -t532 * t484 + t533 * t529;
	t421 = t424 * qJD(1);
	t434 = sin(qJ(3));
	t436 = cos(qJ(3));
	t423 = t533 * t484 + t532 * t529;
	t527 = sin(pkin(7));
	t528 = sin(pkin(6));
	t482 = t528 * t527;
	t472 = t533 * t482;
	t530 = cos(pkin(7));
	t485 = t531 * t529;
	t540 = -t533 * t485 + t532 * t526;
	t541 = t540 * t530;
	t447 = t541 + t472;
	t536 = t423 * t434 + t447 * t436;
	t456 = t532 * t485 + t533 * t526;
	t420 = t456 * qJD(1);
	t469 = t532 * t482;
	t546 = qJD(1) * t469 - t420 * t530;
	t383 = t536 * qJD(3) - t421 * t436 - t546 * t434;
	t507 = t423 * t436;
	t404 = t447 * t434 - t507;
	t505 = qJ(4) + qJ(5);
	t431 = sin(t505);
	t432 = qJD(4) + qJD(5);
	t483 = t530 * t528;
	t448 = -t533 * t483 + t527 * t540;
	t470 = t532 * t483;
	t457 = qJD(1) * t470 + t420 * t527;
	t493 = cos(t505);
	t486 = t432 * t493;
	t559 = t404 * t486 + (-t448 * t432 + t383) * t431 + t457 * t493;
	t557 = -t383 * t493 + (t404 * t432 + t457) * t431;
	t445 = t448 * t493;
	t391 = t404 * t431 + t445;
	t392 = t404 * t493 - t448 * t431;
	t386 = t391 ^ 2;
	t453 = t529 * t483 + t531 * t527;
	t481 = t528 * t526;
	t415 = t453 * t434 + t436 * t481;
	t422 = -t529 * t482 + t531 * t530;
	t466 = -t415 * t431 + t422 * t493;
	t396 = 0.1e1 / t466 ^ 2;
	t374 = t386 * t396 + 0.1e1;
	t368 = 0.1e1 / t374;
	t414 = -t434 * t481 + t453 * t436;
	t409 = t414 * qJD(3);
	t384 = t415 * t486 + (t422 * t432 + t409) * t431;
	t395 = 0.1e1 / t466;
	t512 = t391 * t396;
	t478 = -t384 * t512 - t395 * t559;
	t341 = t478 * t368;
	t375 = atan2(t391, -t466);
	t362 = sin(t375);
	t363 = cos(t375);
	t480 = t362 * t466 + t363 * t391;
	t336 = t480 * t341 + t362 * t559 + t363 * t384;
	t353 = t362 * t391 - t363 * t466;
	t351 = 0.1e1 / t353 ^ 2;
	t553 = t336 * t351;
	t550 = -t456 * t530 + t469;
	t406 = t424 * t436 + t550 * t434;
	t446 = t456 * t527 + t470;
	t394 = t406 * t493 + t446 * t431;
	t538 = t550 * t436;
	t405 = t424 * t434 - t538;
	t433 = sin(qJ(6));
	t435 = cos(qJ(6));
	t372 = t394 * t433 - t405 * t435;
	t549 = 0.2e1 * t372;
	t350 = 0.1e1 / t353;
	t548 = t350 * t553;
	t443 = t446 * t493;
	t393 = t406 * t431 - t443;
	t534 = 0.2e1 * t393;
	t489 = t534 * t548;
	t419 = t423 * qJD(1);
	t504 = qJD(3) * t434;
	t539 = t447 * qJD(1);
	t379 = t538 * qJD(3) - t419 * t436 - t424 * t504 + t539 * t434;
	t444 = t448 * qJD(1);
	t357 = t406 * t486 + t493 * t444 + (t432 * t446 + t379) * t431;
	t518 = t357 * t351;
	t545 = -t518 + t489;
	t544 = (t434 * t541 - t507) * qJD(3) - t421 * t434 + t546 * t436 + t472 * t504;
	t543 = t384 * t396;
	t475 = -t395 * t536 - t414 * t512;
	t542 = t431 * t475;
	t373 = t394 * t435 + t405 * t433;
	t365 = 0.1e1 / t373;
	t366 = 0.1e1 / t373 ^ 2;
	t535 = 0.2e1 * t391;
	t387 = t393 ^ 2;
	t347 = t387 * t351 + 0.1e1;
	t525 = (-t387 * t548 + t393 * t518) / t347 ^ 2;
	t506 = t431 * t432;
	t358 = t379 * t493 - t406 * t506 - t431 * t444 + t432 * t443;
	t378 = t406 * qJD(3) - t419 * t434 - t539 * t436;
	t348 = t373 * qJD(6) + t358 * t433 - t378 * t435;
	t364 = t372 ^ 2;
	t356 = t364 * t366 + 0.1e1;
	t515 = t366 * t372;
	t503 = qJD(6) * t372;
	t349 = t358 * t435 + t378 * t433 - t503;
	t521 = t349 * t365 * t366;
	t524 = (t348 * t515 - t364 * t521) / t356 ^ 2;
	t514 = t395 * t543;
	t522 = (t386 * t514 + t512 * t559) / t374 ^ 2;
	t520 = t351 * t393;
	t354 = 0.1e1 / t356;
	t519 = t354 * t366;
	t517 = t362 * t393;
	t516 = t363 * t393;
	t513 = t391 * t395;
	t510 = t405 * t431;
	t502 = 0.2e1 * t525;
	t501 = -0.2e1 * t524;
	t500 = -0.2e1 * t522;
	t498 = t366 * t524;
	t497 = t395 * t522;
	t496 = t348 * t519;
	t495 = t372 * t521;
	t491 = 0.2e1 * t495;
	t490 = t514 * t535;
	t487 = t405 * t493;
	t371 = t392 * t435 - t433 * t536;
	t370 = t392 * t433 + t435 * t536;
	t477 = -t433 * t365 + t435 * t515;
	t399 = t415 * t493 + t422 * t431;
	t476 = -t392 * t395 - t399 * t512;
	t468 = qJD(6) * t487 + t379;
	t467 = -t362 + (t363 * t513 + t362) * t368;
	t458 = qJD(6) * t406 - t493 * t378 + t405 * t506;
	t410 = t415 * qJD(3);
	t385 = t409 * t493 + t466 * t432;
	t377 = t406 * t433 - t435 * t487;
	t361 = -t448 * t486 - t557;
	t360 = t432 * t445 + t557;
	t345 = 0.1e1 / t347;
	t344 = t368 * t542;
	t343 = t476 * t368;
	t338 = (t362 * t536 + t363 * t414) * t431 + t480 * t344;
	t337 = t480 * t343 + t362 * t392 + t363 * t399;
	t334 = t476 * t500 + (-t399 * t490 + t360 * t395 + (-t384 * t392 - t385 * t391 - t399 * t559) * t396) * t368;
	t333 = t500 * t542 + ((-t414 * t490 + t544 * t395 + (-t384 * t536 + t391 * t410 - t414 * t559) * t396) * t431 + t475 * t486) * t368;
	t332 = t477 * t393 * t501 + (t477 * t357 + ((-qJD(6) * t365 - 0.2e1 * t495) * t435 + (t348 * t435 + (t349 - t503) * t433) * t366) * t393) * t354;
	t331 = (t337 * t520 - t350 * t394) * t502 + (t337 * t489 + t358 * t350 + (-t394 * t336 - t337 * t357 - (t334 * t391 + t343 * t559 + t385 + (t343 * t466 + t392) * t341) * t516 - (t334 * t466 - t343 * t384 - t360 + (-t343 * t391 - t399) * t341) * t517) * t351) * t345;
	t1 = [-t497 * t534 + (t357 * t395 + t393 * t543) * t368, 0, t333, t334, t334, 0; -0.2e1 * t391 * t350 * t525 + (t559 * t350 - t391 * t553 - (t467 * t357 + ((-t341 * t368 * t513 + t500) * t362 + (-t497 * t535 - t341 + (t341 - t478) * t368) * t363) * t393) * t520) * t345 + (t545 * t345 + t520 * t502) * t467 * t393, 0, (t338 * t520 + t350 * t510) * t502 + ((-t378 * t431 - t405 * t486) * t350 + t545 * t338 + (t510 * t336 - (t414 * t486 + t333 * t391 + t344 * t559 - t410 * t431 + (t344 * t466 + t431 * t536) * t341) * t516 - (t536 * t486 + t333 * t466 - t344 * t384 - t544 * t431 + (-t344 * t391 - t414 * t431) * t341) * t517) * t351) * t345, t331, t331, 0; 0.2e1 * (-t365 * t370 + t371 * t515) * t524 + ((t371 * qJD(6) + t361 * t433 - t435 * t544) * t365 + t371 * t491 + (-t370 * t349 - (-t370 * qJD(6) + t361 * t435 + t433 * t544) * t372 - t371 * t348) * t366) * t354, 0, (t498 * t549 - t496) * t377 + (-t349 * t519 + t365 * t501) * (-t406 * t435 - t433 * t487) + ((t458 * t433 - t468 * t435) * t365 - (t468 * t433 + t458 * t435) * t515 + t377 * t491) * t354, t332, t332, t501 + (t496 + (-t354 * t521 - t498) * t372) * t549;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end