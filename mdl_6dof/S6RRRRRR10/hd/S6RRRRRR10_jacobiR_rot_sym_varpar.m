% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (58->25), mult. (178->58), div. (0->0), fcn. (255->10), ass. (0->32)
	t115 = sin(qJ(2));
	t116 = sin(qJ(1));
	t118 = cos(qJ(2));
	t119 = cos(qJ(1));
	t136 = cos(pkin(6));
	t124 = t119 * t136;
	t104 = t116 * t115 - t118 * t124;
	t105 = t115 * t124 + t116 * t118;
	t113 = cos(pkin(7));
	t114 = sin(qJ(3));
	t117 = cos(qJ(3));
	t111 = sin(pkin(7));
	t112 = sin(pkin(6));
	t133 = t112 * t119;
	t126 = t111 * t133;
	t137 = (t104 * t113 + t126) * t117 + t105 * t114;
	t134 = t112 * t116;
	t132 = t113 * t114;
	t131 = t113 * t117;
	t130 = t114 * t115;
	t129 = t114 * t118;
	t128 = t115 * t117;
	t127 = t117 * t118;
	t125 = t116 * t136;
	t123 = t136 * t111;
	t106 = -t119 * t115 - t118 * t125;
	t121 = t106 * t113 + t111 * t134;
	t120 = t104 * t132 - t105 * t117 + t114 * t126;
	t107 = -t115 * t125 + t119 * t118;
	t103 = t107 * t117 + t121 * t114;
	t102 = -t107 * t114 + t121 * t117;
	t1 = [t120, t106 * t117 - t107 * t132, t102, 0, 0, 0; t103, -t104 * t117 - t105 * t132, -t137, 0, 0, 0; 0, (-t113 * t130 + t127) * t112, t117 * t123 + (t113 * t127 - t130) * t112, 0, 0, 0; t137, -t106 * t114 - t107 * t131, -t103, 0, 0, 0; t102, t104 * t114 - t105 * t131, t120, 0, 0, 0; 0, (-t113 * t128 - t129) * t112, -t114 * t123 + (-t113 * t129 - t128) * t112, 0, 0, 0; -t104 * t111 + t113 * t133, t107 * t111, 0, 0, 0, 0; -t106 * t111 + t113 * t134, t105 * t111, 0, 0, 0, 0; 0, t112 * t115 * t111, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:54
	% EndTime: 2019-10-10 13:34:54
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (239->56), mult. (721->123), div. (0->0), fcn. (987->14), ass. (0->62)
	t225 = sin(qJ(2));
	t226 = sin(qJ(1));
	t229 = cos(qJ(2));
	t230 = cos(qJ(1));
	t259 = cos(pkin(6));
	t241 = t230 * t259;
	t211 = t226 * t225 - t229 * t241;
	t212 = t225 * t241 + t226 * t229;
	t224 = sin(qJ(3));
	t228 = cos(qJ(3));
	t219 = sin(pkin(7));
	t220 = sin(pkin(6));
	t253 = t220 * t230;
	t243 = t219 * t253;
	t222 = cos(pkin(7));
	t250 = t222 * t224;
	t196 = t211 * t250 - t212 * t228 + t224 * t243;
	t223 = sin(qJ(4));
	t227 = cos(qJ(4));
	t195 = (t211 * t222 + t243) * t228 + t212 * t224;
	t205 = -t211 * t219 + t222 * t253;
	t218 = sin(pkin(8));
	t221 = cos(pkin(8));
	t238 = t195 * t221 + t205 * t218;
	t264 = t196 * t227 + t238 * t223;
	t263 = -t196 * t223 + t238 * t227;
	t257 = t218 * t219;
	t256 = t219 * t220;
	t255 = t219 * t221;
	t254 = t220 * t226;
	t252 = t221 * t223;
	t251 = t221 * t227;
	t249 = t222 * t228;
	t248 = t224 * t225;
	t247 = t224 * t229;
	t246 = t225 * t228;
	t245 = t228 * t229;
	t244 = t225 * t256;
	t242 = t226 * t259;
	t240 = t259 * t219;
	t214 = -t225 * t242 + t230 * t229;
	t213 = -t230 * t225 - t229 * t242;
	t232 = t213 * t222 + t219 * t254;
	t197 = -t214 * t224 + t232 * t228;
	t207 = -t213 * t219 + t222 * t254;
	t237 = t197 * t221 + t207 * t218;
	t203 = t228 * t240 + (t222 * t245 - t248) * t220;
	t236 = t203 * t221 + (t259 * t222 - t229 * t256) * t218;
	t199 = t211 * t224 - t212 * t249;
	t235 = t199 * t221 + t212 * t257;
	t201 = -t213 * t224 - t214 * t249;
	t234 = t201 * t221 + t214 * t257;
	t208 = (-t222 * t246 - t247) * t220;
	t231 = t208 * t221 + t218 * t244;
	t209 = (-t222 * t248 + t245) * t220;
	t204 = t224 * t240 + (t222 * t247 + t246) * t220;
	t202 = t213 * t228 - t214 * t250;
	t200 = -t211 * t228 - t212 * t250;
	t198 = t214 * t228 + t232 * t224;
	t192 = t198 * t227 + t237 * t223;
	t191 = -t198 * t223 + t237 * t227;
	t1 = [t264, t202 * t227 + t234 * t223, t197 * t227 - t198 * t252, t191, 0, 0; t192, t200 * t227 + t235 * t223, -t195 * t227 + t196 * t252, -t263, 0, 0; 0, t209 * t227 + t231 * t223, t203 * t227 - t204 * t252, -t204 * t223 + t236 * t227, 0, 0; t263, -t202 * t223 + t234 * t227, -t197 * t223 - t198 * t251, -t192, 0, 0; t191, -t200 * t223 + t235 * t227, t195 * t223 + t196 * t251, t264, 0, 0; 0, -t209 * t223 + t231 * t227, -t203 * t223 - t204 * t251, -t204 * t227 - t236 * t223, 0, 0; -t195 * t218 + t205 * t221, -t201 * t218 + t214 * t255, t198 * t218, 0, 0, 0; -t197 * t218 + t207 * t221, -t199 * t218 + t212 * t255, -t196 * t218, 0, 0, 0; 0, -t208 * t218 + t221 * t244, t204 * t218, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:55
	% EndTime: 2019-10-10 13:34:56
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (495->79), mult. (1476->167), div. (0->0), fcn. (2007->16), ass. (0->85)
	t321 = sin(qJ(2));
	t322 = sin(qJ(1));
	t326 = cos(qJ(2));
	t327 = cos(qJ(1));
	t358 = cos(pkin(6));
	t336 = t327 * t358;
	t306 = t321 * t322 - t326 * t336;
	t307 = t321 * t336 + t322 * t326;
	t320 = sin(qJ(3));
	t325 = cos(qJ(3));
	t314 = sin(pkin(7));
	t315 = sin(pkin(6));
	t349 = t315 * t327;
	t339 = t314 * t349;
	t317 = cos(pkin(7));
	t346 = t317 * t320;
	t289 = t306 * t346 - t307 * t325 + t320 * t339;
	t319 = sin(qJ(4));
	t324 = cos(qJ(4));
	t288 = (t306 * t317 + t339) * t325 + t307 * t320;
	t300 = -t306 * t314 + t317 * t349;
	t313 = sin(pkin(8));
	t316 = cos(pkin(8));
	t334 = t288 * t316 + t300 * t313;
	t268 = t289 * t324 + t334 * t319;
	t277 = t288 * t313 - t300 * t316;
	t318 = sin(qJ(5));
	t323 = cos(qJ(5));
	t366 = t268 * t318 + t277 * t323;
	t365 = t268 * t323 - t277 * t318;
	t266 = t289 * t319 - t334 * t324;
	t337 = t322 * t358;
	t308 = -t327 * t321 - t326 * t337;
	t350 = t315 * t322;
	t302 = -t308 * t314 + t317 * t350;
	t357 = t302 * t313;
	t355 = t313 * t314;
	t354 = t313 * t318;
	t353 = t313 * t323;
	t352 = t314 * t315;
	t351 = t314 * t316;
	t348 = t316 * t319;
	t347 = t316 * t324;
	t345 = t317 * t325;
	t344 = t320 * t321;
	t343 = t320 * t326;
	t342 = t321 * t325;
	t341 = t325 * t326;
	t340 = t321 * t352;
	t338 = t314 * t358;
	t298 = t325 * t338 + (t317 * t341 - t344) * t315;
	t305 = t358 * t317 - t326 * t352;
	t333 = t298 * t316 + t305 * t313;
	t292 = t306 * t320 - t307 * t345;
	t332 = t292 * t316 + t307 * t355;
	t309 = -t321 * t337 + t326 * t327;
	t294 = -t308 * t320 - t309 * t345;
	t331 = t294 * t316 + t309 * t355;
	t329 = t308 * t317 + t314 * t350;
	t303 = (-t317 * t342 - t343) * t315;
	t328 = t303 * t316 + t313 * t340;
	t304 = (-t317 * t344 + t341) * t315;
	t299 = t320 * t338 + (t317 * t343 + t342) * t315;
	t296 = -t303 * t313 + t316 * t340;
	t295 = t308 * t325 - t309 * t346;
	t293 = -t306 * t325 - t307 * t346;
	t291 = t309 * t325 + t329 * t320;
	t290 = -t309 * t320 + t329 * t325;
	t285 = -t298 * t313 + t305 * t316;
	t283 = -t294 * t313 + t309 * t351;
	t282 = -t292 * t313 + t307 * t351;
	t281 = t304 * t324 + t328 * t319;
	t280 = t298 * t324 - t299 * t348;
	t279 = -t290 * t313 + t302 * t316;
	t276 = t299 * t324 + t333 * t319;
	t275 = -t299 * t319 + t333 * t324;
	t274 = t290 * t324 - t291 * t348;
	t273 = -t288 * t324 + t289 * t348;
	t272 = t295 * t324 + t331 * t319;
	t271 = t293 * t324 + t332 * t319;
	t270 = t291 * t324 + (t290 * t316 + t357) * t319;
	t269 = -t290 * t347 + t291 * t319 - t324 * t357;
	t265 = t270 * t323 + t279 * t318;
	t264 = -t270 * t318 + t279 * t323;
	t1 = [t365, t272 * t323 + t283 * t318, t274 * t323 + t291 * t354, -t269 * t323, t264, 0; t265, t271 * t323 + t282 * t318, t273 * t323 - t289 * t354, t266 * t323, t366, 0; 0, t281 * t323 + t296 * t318, t280 * t323 + t299 * t354, t275 * t323, -t276 * t318 + t285 * t323, 0; -t366, -t272 * t318 + t283 * t323, -t274 * t318 + t291 * t353, t269 * t318, -t265, 0; t264, -t271 * t318 + t282 * t323, -t273 * t318 - t289 * t353, -t266 * t318, t365, 0; 0, -t281 * t318 + t296 * t323, -t280 * t318 + t299 * t353, -t275 * t318, -t276 * t323 - t285 * t318, 0; t266, t295 * t319 - t331 * t324, t290 * t319 + t291 * t347, t270, 0, 0; t269, t293 * t319 - t332 * t324, -t288 * t319 - t289 * t347, -t268, 0, 0; 0, t304 * t319 - t328 * t324, t298 * t319 + t299 * t347, t276, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:57
	% EndTime: 2019-10-10 13:34:59
	% DurationCPUTime: 1.50s
	% Computational Cost: add. (985->110), mult. (2912->225), div. (0->0), fcn. (3942->18), ass. (0->110)
	t415 = sin(qJ(2));
	t416 = sin(qJ(1));
	t421 = cos(qJ(2));
	t422 = cos(qJ(1));
	t460 = cos(pkin(6));
	t436 = t422 * t460;
	t399 = t416 * t415 - t421 * t436;
	t400 = t415 * t436 + t416 * t421;
	t414 = sin(qJ(3));
	t420 = cos(qJ(3));
	t407 = sin(pkin(7));
	t408 = sin(pkin(6));
	t450 = t408 * t422;
	t438 = t407 * t450;
	t410 = cos(pkin(7));
	t447 = t410 * t414;
	t382 = t399 * t447 - t400 * t420 + t414 * t438;
	t413 = sin(qJ(4));
	t419 = cos(qJ(4));
	t381 = (t399 * t410 + t438) * t420 + t400 * t414;
	t409 = cos(pkin(8));
	t395 = -t399 * t407 + t410 * t450;
	t406 = sin(pkin(8));
	t459 = t395 * t406;
	t434 = t381 * t409 + t459;
	t352 = t382 * t419 + t434 * t413;
	t368 = t381 * t406 - t395 * t409;
	t412 = sin(qJ(5));
	t418 = cos(qJ(5));
	t340 = t352 * t418 - t368 * t412;
	t411 = sin(qJ(6));
	t467 = t340 * t411;
	t417 = cos(qJ(6));
	t466 = t340 * t417;
	t338 = t352 * t412 + t368 * t418;
	t462 = t382 * t413;
	t456 = t406 * t407;
	t455 = t406 * t412;
	t454 = t406 * t418;
	t453 = t407 * t408;
	t452 = t407 * t409;
	t451 = t408 * t416;
	t449 = t409 * t413;
	t448 = t409 * t419;
	t446 = t410 * t420;
	t445 = t411 * t418;
	t444 = t414 * t415;
	t443 = t414 * t421;
	t442 = t415 * t420;
	t441 = t417 * t418;
	t440 = t420 * t421;
	t439 = t415 * t453;
	t437 = t416 * t460;
	t435 = t460 * t407;
	t385 = t399 * t414 - t400 * t446;
	t433 = t385 * t409 + t400 * t456;
	t401 = -t422 * t415 - t421 * t437;
	t402 = -t415 * t437 + t422 * t421;
	t387 = -t401 * t414 - t402 * t446;
	t432 = t387 * t409 + t402 * t456;
	t430 = -t401 * t407 + t410 * t451;
	t429 = t401 * t410 + t407 * t451;
	t397 = (-t410 * t442 - t443) * t408;
	t428 = t397 * t409 + t406 * t439;
	t427 = t460 * t410 - t421 * t453;
	t425 = t430 * t406;
	t424 = t427 * t406;
	t383 = -t402 * t414 + t429 * t420;
	t423 = -t383 * t406 + t430 * t409;
	t398 = (-t410 * t444 + t440) * t408;
	t394 = t414 * t435 + (t410 * t443 + t442) * t408;
	t393 = t420 * t435 + (t410 * t440 - t444) * t408;
	t389 = -t397 * t406 + t409 * t439;
	t388 = t401 * t420 - t402 * t447;
	t386 = -t399 * t420 - t400 * t447;
	t384 = t402 * t420 + t429 * t414;
	t378 = -t393 * t406 + t427 * t409;
	t375 = -t387 * t406 + t402 * t452;
	t374 = -t385 * t406 + t400 * t452;
	t373 = t398 * t419 + t428 * t413;
	t372 = t398 * t413 - t428 * t419;
	t371 = t393 * t419 - t394 * t449;
	t370 = t393 * t413 + t394 * t448;
	t366 = t394 * t419 + (t393 * t409 + t424) * t413;
	t365 = -t393 * t448 + t394 * t413 - t419 * t424;
	t364 = t371 * t418 + t394 * t455;
	t363 = t373 * t418 + t389 * t412;
	t362 = t383 * t419 - t384 * t449;
	t361 = t383 * t413 + t384 * t448;
	t360 = -t381 * t419 + t382 * t449;
	t359 = -t381 * t413 - t382 * t448;
	t358 = t388 * t419 + t432 * t413;
	t357 = t388 * t413 - t432 * t419;
	t356 = t386 * t419 + t433 * t413;
	t355 = t386 * t413 - t433 * t419;
	t354 = t384 * t419 + (t383 * t409 + t425) * t413;
	t353 = -t383 * t448 + t384 * t413 - t419 * t425;
	t351 = -t434 * t419 + t462;
	t349 = t381 * t448 + t419 * t459 - t462;
	t348 = t366 * t418 + t378 * t412;
	t347 = -t366 * t412 + t378 * t418;
	t346 = t362 * t418 + t384 * t455;
	t345 = t360 * t418 - t382 * t455;
	t344 = t358 * t418 + t375 * t412;
	t343 = t356 * t418 + t374 * t412;
	t342 = t354 * t418 + t423 * t412;
	t341 = t354 * t412 - t423 * t418;
	t337 = t342 * t417 + t353 * t411;
	t336 = -t342 * t411 + t353 * t417;
	t1 = [t351 * t411 + t466, t344 * t417 + t357 * t411, t346 * t417 + t361 * t411, -t353 * t441 + t354 * t411, -t341 * t417, t336; t337, t343 * t417 + t355 * t411, t345 * t417 + t359 * t411, -t349 * t441 - t352 * t411, t338 * t417, t349 * t417 + t467; 0, t363 * t417 + t372 * t411, t364 * t417 + t370 * t411, -t365 * t441 + t366 * t411, t347 * t417, -t348 * t411 + t365 * t417; t351 * t417 - t467, -t344 * t411 + t357 * t417, -t346 * t411 + t361 * t417, t353 * t445 + t354 * t417, t341 * t411, -t337; t336, -t343 * t411 + t355 * t417, -t345 * t411 + t359 * t417, t349 * t445 - t352 * t417, -t338 * t411, -t349 * t411 + t466; 0, -t363 * t411 + t372 * t417, -t364 * t411 + t370 * t417, t365 * t445 + t366 * t417, -t347 * t411, -t348 * t417 - t365 * t411; t338, t358 * t412 - t375 * t418, t362 * t412 - t384 * t454, -t353 * t412, t342, 0; t341, t356 * t412 - t374 * t418, t360 * t412 + t382 * t454, -t349 * t412, -t340, 0; 0, t373 * t412 - t389 * t418, t371 * t412 - t394 * t454, -t365 * t412, t348, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end