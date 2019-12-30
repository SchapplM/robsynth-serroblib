% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR14
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
%   Wie in S5RPRRR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRR14_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR14_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (137->30), mult. (614->95), div. (108->12), fcn. (792->9), ass. (0->49)
	t86 = sin(pkin(5));
	t79 = t86 ^ 2;
	t88 = cos(pkin(5));
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
	t85 = sin(pkin(11));
	t105 = t89 * t85;
	t87 = cos(pkin(11));
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
	t1 = [(-t74 * t80 * t86 + t78 * t100) * t102, 0, 0, 0, 0; (0.2e1 * (t51 * t111 - t55 * t90) / t54 ^ 2 * (-t50 * t57 * t83 + t102 * t111) * t79 + ((0.2e1 * t51 * t57 * t89 - t56 * t90) * t50 + (-t89 * t55 + ((-t51 + (-t78 * t112 - t101) * t89 * t72) * t90 - (t84 * t79 ^ 2 * t100 + (-t108 + (0.2e1 * t83 - t84) * t74) * t107) * t89 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0, 0; 0.2e1 * (t66 * t109 + t65 * t68) / t61 ^ 2 * (t63 * t64 * t67 + t62 * t110) + (-t69 * t62 * t66 + (-0.2e1 * t67 * t109 - t68 * t66) * t63 + (-t99 * t110 + t70 * t65) * qJD(1)) / t61, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:38
	% EndTime: 2019-12-29 18:05:39
	% DurationCPUTime: 0.89s
	% Computational Cost: add. (994->58), mult. (3107->139), div. (132->12), fcn. (4021->13), ass. (0->80)
	t183 = sin(pkin(6));
	t184 = sin(pkin(5));
	t185 = cos(pkin(11));
	t186 = cos(pkin(6));
	t187 = cos(pkin(5));
	t173 = -t184 * t185 * t183 + t187 * t186;
	t191 = cos(qJ(1));
	t211 = t191 * t185;
	t182 = sin(pkin(11));
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
	t151 = t188 * t202 - t190 * t200;
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
	t139 = qJD(3) * t151 + t167 * t188 - t190 * t199;
	t215 = t186 * t190;
	t218 = t200 * t188;
	t150 = -t190 * t209 + t201 * t215 - t218;
	t144 = t150 ^ 2;
	t143 = t144 * t146 + 0.1e1;
	t225 = t146 * t150;
	t140 = t167 * t190 + t199 * t188 + (t190 * t202 + t218) * qJD(3);
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
	t149 = t175 * t190 + t188 * t203;
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
	t148 = -t190 * t203 + t219;
	t141 = 0.1e1 / t143;
	t134 = t160 * t136 + 0.1e1;
	t131 = t198 * t163;
	t1 = [t220 * t222 * t229 + t158 * t223, 0, 0, 0, 0; 0.2e1 * (-t131 * t136 * t163 - t135 * t204) / t134 ^ 2 * (-t160 * t227 + t163 * t221) + (t159 * t135 + (-t130 * t204 + t131 * t158) * t136 + (-0.2e1 * t131 * t227 + t198 * t221 + (t152 * t171 * t222 + (0.2e1 * t223 + (t230 * t229 - t170) * t155) * t153) * t136 * t220) * t163) / t134, 0, 0, 0, 0; 0.2e1 * (-t145 * t148 + t146 * t224) * t228 + ((-t168 * t215 + t169 * t188 + t190 * t205) * t145 + 0.2e1 * t224 * t226 + (-t148 * t140 - (t169 * t190 + (t168 * t186 - t205) * t188) * t150 - t149 * t139) * t146 + (t149 * t145 - (t174 * t215 + t190 * t208 - t219) * t225) * qJD(3)) * t141, 0, -0.2e1 * t228 + 0.2e1 * (t139 * t146 * t141 + (-t141 * t226 - t146 * t228) * t150) * t150, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:38
	% EndTime: 2019-12-29 18:05:41
	% DurationCPUTime: 3.00s
	% Computational Cost: add. (4837->111), mult. (15324->224), div. (448->12), fcn. (19446->15), ass. (0->108)
	t347 = sin(pkin(11));
	t352 = cos(pkin(5));
	t320 = t352 * t347;
	t350 = cos(pkin(11));
	t353 = sin(qJ(1));
	t355 = cos(qJ(1));
	t276 = t320 * t355 + t350 * t353;
	t287 = sin(qJ(3));
	t348 = sin(pkin(6));
	t349 = sin(pkin(5));
	t318 = t349 * t348;
	t311 = t355 * t318;
	t283 = t287 * t311;
	t322 = t352 * t350;
	t275 = -t322 * t355 + t347 * t353;
	t351 = cos(pkin(6));
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
	t1 = [-t328 * t356 + (t228 * t261 - t335 * t358) * t243, 0, t208, 0, 0; t249 * t219 * t331 + (t357 * t219 + t249 * t360 + (t306 * t228 - ((t213 * t243 * t337 + t329) * t241 + (0.2e1 * t249 * t328 - t213 + (t213 - t314) * t243) * t242) * t358) * t341) * t216 - (-t341 * t331 + (-t339 + t325) * t216) * t306 * t358, 0, (-t210 * t341 - t219 * t255) * t331 + (t210 * t325 + t229 * t219 + (-t255 * t209 - t210 * t228 - (-(-t208 * t249 + t214 * t357 + t256 + (t214 * t295 + t362) * t213) * t242 - (t208 * t295 - t214 * t257 - t231 + (t214 * t249 - t265) * t213) * t241) * t358) * t220) * t216, 0, 0; 0.2e1 * (-t234 * t237 + t238 * t338) * t344 + ((qJD(4) * t238 - t231 * t286 - t259 * t288) * t234 + 0.2e1 * t238 * t327 + (-t237 * t224 - (-qJD(4) * t237 - t231 * t288 + t259 * t286) * t239 - t238 * t223) * t235) * t225, 0, -t313 * t358 * t330 + (t313 * t228 - ((-qJD(4) * t234 - 0.2e1 * t327) * t288 + (t223 * t288 + (t224 - t332) * t286) * t235) * t358) * t225, t330 + 0.2e1 * (t223 * t235 * t225 + (-t225 * t340 - t235 * t344) * t239) * t239, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:39
	% EndTime: 2019-12-29 18:05:47
	% DurationCPUTime: 8.10s
	% Computational Cost: add. (13786->165), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->145)
	t482 = sin(pkin(11));
	t487 = cos(pkin(5));
	t442 = t487 * t482;
	t485 = cos(pkin(11));
	t488 = sin(qJ(1));
	t490 = cos(qJ(1));
	t383 = -t488 * t442 + t490 * t485;
	t380 = t383 * qJD(1);
	t392 = sin(qJ(3));
	t394 = cos(qJ(3));
	t382 = t490 * t442 + t488 * t485;
	t483 = sin(pkin(6));
	t484 = sin(pkin(5));
	t440 = t484 * t483;
	t427 = t490 * t440;
	t486 = cos(pkin(6));
	t443 = t487 * t485;
	t497 = -t490 * t443 + t488 * t482;
	t498 = t497 * t486;
	t404 = t498 + t427;
	t493 = t382 * t392 + t404 * t394;
	t413 = t488 * t443 + t490 * t482;
	t379 = t413 * qJD(1);
	t424 = t488 * t440;
	t503 = qJD(1) * t424 - t379 * t486;
	t340 = t493 * qJD(3) - t380 * t394 - t503 * t392;
	t391 = sin(qJ(4));
	t441 = t486 * t484;
	t425 = t488 * t441;
	t414 = qJD(1) * t425 + t379 * t483;
	t489 = cos(qJ(4));
	t463 = t382 * t394;
	t363 = t404 * t392 - t463;
	t405 = -t490 * t441 + t483 * t497;
	t512 = -t363 * t391 - t405 * t489;
	t320 = qJD(4) * t512 + t340 * t489 - t414 * t391;
	t351 = t363 * t489 - t405 * t391;
	t514 = t351 * qJD(4) + t340 * t391 + t414 * t489;
	t345 = t512 ^ 2;
	t410 = t485 * t441 + t487 * t483;
	t439 = t484 * t482;
	t374 = t410 * t392 + t394 * t439;
	t381 = -t485 * t440 + t487 * t486;
	t430 = -t374 * t391 + t381 * t489;
	t356 = 0.1e1 / t430 ^ 2;
	t333 = t345 * t356 + 0.1e1;
	t331 = 0.1e1 / t333;
	t359 = t374 * t489 + t381 * t391;
	t373 = -t392 * t439 + t410 * t394;
	t368 = t373 * qJD(3);
	t343 = t359 * qJD(4) + t368 * t391;
	t355 = 0.1e1 / t430;
	t468 = t512 * t356;
	t435 = t343 * t468 - t355 * t514;
	t300 = t435 * t331;
	t334 = atan2(-t512, -t430);
	t329 = sin(t334);
	t330 = cos(t334);
	t438 = t329 * t430 - t330 * t512;
	t295 = t438 * t300 + t329 * t514 + t330 * t343;
	t312 = -t329 * t512 - t330 * t430;
	t310 = 0.1e1 / t312 ^ 2;
	t511 = t295 * t310;
	t508 = -t413 * t486 + t424;
	t365 = t383 * t394 + t508 * t392;
	t403 = t413 * t483 + t425;
	t353 = t365 * t489 + t403 * t391;
	t495 = t508 * t394;
	t364 = t383 * t392 - t495;
	t390 = sin(qJ(5));
	t393 = cos(qJ(5));
	t327 = t353 * t390 - t364 * t393;
	t507 = 0.2e1 * t327;
	t309 = 0.1e1 / t312;
	t506 = t309 * t511;
	t378 = t382 * qJD(1);
	t462 = qJD(3) * t392;
	t496 = t404 * qJD(1);
	t336 = t495 * qJD(3) - t378 * t394 - t383 * t462 + t496 * t392;
	t402 = qJD(1) * t405;
	t316 = t353 * qJD(4) + t336 * t391 + t489 * t402;
	t401 = t403 * t489;
	t352 = t365 * t391 - t401;
	t491 = 0.2e1 * t352;
	t447 = t491 * t506;
	t502 = -t310 * t316 + t447;
	t501 = (t392 * t498 - t463) * qJD(3) - t380 * t392 + t503 * t394 + t427 * t462;
	t346 = t352 ^ 2;
	t306 = t346 * t310 + 0.1e1;
	t475 = t310 * t352;
	t459 = 0.2e1 * (t316 * t475 - t346 * t506) / t306 ^ 2;
	t500 = t343 * t356;
	t432 = -t355 * t493 + t373 * t468;
	t499 = t391 * t432;
	t449 = qJD(4) * t489;
	t328 = t353 * t393 + t364 * t390;
	t322 = 0.1e1 / t328;
	t323 = 0.1e1 / t328 ^ 2;
	t492 = -0.2e1 * t512;
	t461 = qJD(4) * t391;
	t317 = qJD(4) * t401 + t336 * t489 - t365 * t461 - t391 * t402;
	t335 = t365 * qJD(3) - t378 * t392 - t496 * t394;
	t307 = t328 * qJD(5) + t317 * t390 - t335 * t393;
	t321 = t327 ^ 2;
	t315 = t321 * t323 + 0.1e1;
	t473 = t323 * t327;
	t460 = qJD(5) * t327;
	t308 = t317 * t393 + t335 * t390 - t460;
	t477 = t308 * t322 * t323;
	t480 = (t307 * t473 - t321 * t477) / t315 ^ 2;
	t470 = t355 * t500;
	t478 = (t345 * t470 - t468 * t514) / t333 ^ 2;
	t313 = 0.1e1 / t315;
	t474 = t313 * t323;
	t472 = t329 * t352;
	t471 = t330 * t352;
	t469 = t512 * t355;
	t466 = t364 * t391;
	t458 = -0.2e1 * t480;
	t457 = -0.2e1 * t478;
	t455 = t323 * t480;
	t454 = t355 * t478;
	t453 = t307 * t474;
	t452 = t327 * t477;
	t450 = t364 * t489;
	t446 = 0.2e1 * t452;
	t445 = t470 * t492;
	t326 = t351 * t393 - t390 * t493;
	t325 = t351 * t390 + t393 * t493;
	t436 = qJD(5) * t450 + t336;
	t434 = -t390 * t322 + t393 * t473;
	t433 = -t351 * t355 + t359 * t468;
	t423 = -t329 + (-t330 * t469 + t329) * t331;
	t417 = qJD(5) * t365 - t489 * t335 + t364 * t461;
	t369 = t374 * qJD(3);
	t344 = t430 * qJD(4) + t368 * t489;
	t342 = t365 * t390 - t393 * t450;
	t304 = 0.1e1 / t306;
	t303 = t331 * t499;
	t301 = t433 * t331;
	t297 = (t329 * t493 + t330 * t373) * t391 + t438 * t303;
	t296 = t438 * t301 + t329 * t351 + t330 * t359;
	t293 = t433 * t457 + (-t359 * t445 - t320 * t355 + (-t343 * t351 + t344 * t512 - t359 * t514) * t356) * t331;
	t292 = t457 * t499 + ((-t373 * t445 + t501 * t355 + (-t343 * t493 - t369 * t512 - t373 * t514) * t356) * t391 + t432 * t449) * t331;
	t1 = [-t454 * t491 + (t316 * t355 + t352 * t500) * t331, 0, t292, t293, 0; t512 * t309 * t459 + (t514 * t309 + t512 * t511 - (t423 * t316 + ((t300 * t331 * t469 + t457) * t329 + (-t454 * t492 - t300 + (t300 - t435) * t331) * t330) * t352) * t475) * t304 + (t502 * t304 + t475 * t459) * t423 * t352, 0, (t297 * t475 + t309 * t466) * t459 + ((-t335 * t391 - t364 * t449) * t309 + t502 * t297 + (t466 * t295 - (t373 * t449 - t292 * t512 + t303 * t514 - t369 * t391 + (t303 * t430 + t391 * t493) * t300) * t471 - (t493 * t449 + t292 * t430 - t303 * t343 - t501 * t391 + (t303 * t512 - t373 * t391) * t300) * t472) * t310) * t304, (t296 * t475 - t309 * t353) * t459 + (t296 * t447 + t317 * t309 + (-t353 * t295 - t296 * t316 - (-t293 * t512 + t301 * t514 + t344 + (t301 * t430 + t351) * t300) * t471 - (t293 * t430 - t301 * t343 + t320 + (t301 * t512 - t359) * t300) * t472) * t310) * t304, 0; 0.2e1 * (-t322 * t325 + t326 * t473) * t480 + ((t326 * qJD(5) + t320 * t390 - t393 * t501) * t322 + t326 * t446 + (-t325 * t308 - (-t325 * qJD(5) + t320 * t393 + t390 * t501) * t327 - t326 * t307) * t323) * t313, 0, (t455 * t507 - t453) * t342 + (-t308 * t474 + t322 * t458) * (-t365 * t393 - t390 * t450) + ((t417 * t390 - t436 * t393) * t322 - (t436 * t390 + t417 * t393) * t473 + t342 * t446) * t313, t434 * t352 * t458 + (t434 * t316 + ((-qJD(5) * t322 - 0.2e1 * t452) * t393 + (t307 * t393 + (t308 - t460) * t390) * t323) * t352) * t313, t458 + (t453 + (-t313 * t477 - t455) * t327) * t507;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end