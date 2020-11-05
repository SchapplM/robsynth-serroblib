% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPPRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:52
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t225 = cos(pkin(12));
	t224 = sin(pkin(12));
	t1 = [t225, -t224, 0, 0; t224, t225, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t227 = sin(pkin(12));
	t228 = sin(pkin(6));
	t235 = t227 * t228;
	t231 = cos(pkin(6));
	t234 = t227 * t231;
	t230 = cos(pkin(12));
	t233 = t230 * t228;
	t232 = t230 * t231;
	t229 = cos(pkin(13));
	t226 = sin(pkin(13));
	t1 = [-t226 * t234 + t230 * t229, -t230 * t226 - t229 * t234, t235, t230 * pkin(1) + qJ(2) * t235 + 0; t226 * t232 + t227 * t229, -t227 * t226 + t229 * t232, -t233, t227 * pkin(1) - qJ(2) * t233 + 0; t228 * t226, t228 * t229, t231, t231 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->33), mult. (101->60), div. (0->0), fcn. (137->10), ass. (0->27)
	t242 = sin(pkin(13));
	t244 = sin(pkin(7));
	t262 = t242 * t244;
	t249 = cos(pkin(7));
	t261 = t242 * t249;
	t250 = cos(pkin(6));
	t260 = t242 * t250;
	t243 = sin(pkin(12));
	t259 = t243 * t242;
	t245 = sin(pkin(6));
	t258 = t245 * t244;
	t246 = cos(pkin(14));
	t257 = t246 * t249;
	t247 = cos(pkin(13));
	t248 = cos(pkin(12));
	t256 = t248 * t247;
	t255 = t249 * t245;
	t254 = t250 * t244;
	t253 = t250 * t249;
	t237 = t247 * t253 - t258;
	t241 = sin(pkin(14));
	t252 = t237 * t241 + t246 * t260;
	t238 = t247 * t254 + t255;
	t251 = pkin(2) * t260 - t245 * qJ(2) - t238 * qJ(3);
	t239 = t247 * pkin(2) + qJ(3) * t262 + pkin(1);
	t236 = -t241 * t261 + t246 * t247;
	t1 = [t248 * t236 - t252 * t243, (-t237 * t243 - t248 * t261) * t246 + t241 * (t250 * t259 - t256), t238 * t243 + t248 * t262, t239 * t248 - t251 * t243 + 0; t243 * t236 + t252 * t248, (t237 * t246 - t241 * t260) * t248 - t243 * (t241 * t247 + t242 * t257), (-t250 * t256 + t259) * t244 - t248 * t255, t239 * t243 + t251 * t248 + 0; (t247 * t255 + t254) * t241 + t245 * t242 * t246, (-t241 * t242 + t247 * t257) * t245 + t246 * t254, -t247 * t258 + t253, (t249 * qJ(3) + qJ(2)) * t250 + (-t247 * t244 * qJ(3) + t242 * pkin(2)) * t245 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:43
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (133->51), mult. (309->95), div. (0->0), fcn. (398->14), ass. (0->50)
	t287 = sin(pkin(14));
	t293 = cos(pkin(14));
	t290 = sin(pkin(8));
	t315 = pkin(9) * t290;
	t283 = -t287 * pkin(3) + t293 * t315;
	t296 = cos(pkin(8));
	t285 = t296 * pkin(9) + qJ(3);
	t291 = sin(pkin(7));
	t297 = cos(pkin(7));
	t273 = -t283 * t291 + t285 * t297 + qJ(2);
	t292 = sin(pkin(6));
	t298 = cos(pkin(6));
	t274 = t283 * t297 + t291 * t285;
	t282 = t293 * pkin(3) + t287 * t315 + pkin(2);
	t288 = sin(pkin(13));
	t294 = cos(pkin(13));
	t305 = t274 * t294 - t288 * t282;
	t316 = t292 * t273 + t305 * t298;
	t314 = t287 * t288;
	t313 = t291 * t293;
	t312 = t291 * t296;
	t311 = t291 * t298;
	t309 = t293 * t288;
	t308 = t294 * t296;
	t307 = t294 * t297;
	t306 = t297 * t296;
	t302 = -t290 * t291 + t293 * t306;
	t269 = t302 * t294 - t296 * t314;
	t280 = t290 * t297 + t293 * t312;
	t304 = t269 * t298 - t292 * t280;
	t279 = t293 * t307 - t314;
	t303 = -t279 * t298 + t292 * t313;
	t300 = cos(qJ(4));
	t299 = sin(qJ(4));
	t295 = cos(pkin(12));
	t289 = sin(pkin(12));
	t281 = t297 * t292 + t294 * t311;
	t278 = t287 * t307 + t309;
	t277 = t287 * t294 + t297 * t309;
	t276 = -t293 * t294 + t297 * t314;
	t272 = -t292 * t291 * t287 + t278 * t298;
	t271 = t278 * t292 + t287 * t311;
	t270 = t287 * t308 + t302 * t288;
	t268 = t274 * t288 + t282 * t294 + pkin(1);
	t267 = t272 * t295 - t289 * t276;
	t266 = t272 * t289 + t295 * t276;
	t265 = t269 * t292 + t298 * t280;
	t264 = -t295 * t270 - t304 * t289;
	t263 = -t289 * t270 + t304 * t295;
	t1 = [t264 * t299 - t266 * t300, t264 * t300 + t266 * t299, (t295 * t277 - t303 * t289) * t290 + (t291 * t288 * t295 + t281 * t289) * t296, t268 * t295 + t316 * t289 + 0; t263 * t299 + t267 * t300, t263 * t300 - t267 * t299, (-t296 * t281 + t303 * t290) * t295 + t289 * (t277 * t290 + t288 * t312), t268 * t289 - t316 * t295 + 0; t265 * t299 + t300 * t271, t265 * t300 - t299 * t271, (-t279 * t290 - t291 * t308) * t292 - t298 * (t290 * t313 - t306), t273 * t298 - t292 * t305 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:43
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (256->80), mult. (644->146), div. (0->0), fcn. (788->16), ass. (0->74)
	t360 = sin(pkin(14));
	t361 = sin(pkin(13));
	t366 = cos(pkin(14));
	t367 = cos(pkin(13));
	t370 = cos(pkin(7));
	t386 = t367 * t370;
	t346 = t360 * t386 + t366 * t361;
	t365 = sin(pkin(6));
	t371 = cos(pkin(6));
	t364 = sin(pkin(7));
	t397 = t360 * t364;
	t336 = t346 * t371 - t365 * t397;
	t398 = t360 * t361;
	t344 = -t366 * t367 + t370 * t398;
	t362 = sin(pkin(12));
	t368 = cos(pkin(12));
	t327 = t336 * t368 - t362 * t344;
	t373 = sin(qJ(4));
	t375 = cos(qJ(4));
	t369 = cos(pkin(8));
	t384 = t370 * t369;
	t382 = t366 * t384;
	t363 = sin(pkin(8));
	t391 = t363 * t367;
	t396 = t360 * t369;
	t331 = -t361 * t396 - t364 * t391 + t367 * t382;
	t388 = t366 * t369;
	t390 = t363 * t370;
	t347 = t364 * t388 + t390;
	t321 = t331 * t371 - t365 * t347;
	t387 = t367 * t369;
	t393 = t363 * t364;
	t334 = t360 * t387 + (t382 - t393) * t361;
	t406 = t321 * t368 - t362 * t334;
	t412 = t327 * t375 + t406 * t373;
	t354 = pkin(4) * t388 + t360 * pkin(10);
	t340 = pkin(4) * t390 + t364 * t354;
	t352 = pkin(4) * t396 - t366 * pkin(10);
	t379 = pkin(4) * t393 - t354 * t370;
	t403 = t361 * t352 + t379 * t367;
	t411 = t340 * t365 + t403 * t371;
	t353 = -t360 * pkin(4) + pkin(10) * t388;
	t339 = pkin(10) * t390 + t364 * t353;
	t351 = t366 * pkin(4) + pkin(10) * t396;
	t378 = pkin(10) * t393 - t353 * t370;
	t376 = t361 * t351 + t378 * t367;
	t410 = t339 * t365 + t376 * t371;
	t317 = t321 * t362 + t368 * t334;
	t326 = t336 * t362 + t368 * t344;
	t407 = t317 * t373 + t326 * t375;
	t392 = t363 * t366;
	t350 = -t360 * pkin(3) + pkin(9) * t392;
	t358 = t369 * pkin(9) + qJ(3);
	t337 = -t350 * t364 + t358 * t370 + qJ(2);
	t338 = t350 * t370 + t364 * t358;
	t349 = t363 * t360 * pkin(9) + t366 * pkin(3) + pkin(2);
	t381 = t338 * t367 - t361 * t349;
	t405 = t365 * t337 + t381 * t371;
	t335 = t346 * t365 + t371 * t397;
	t380 = t331 * t365 + t371 * t347;
	t404 = t375 * t335 + t380 * t373;
	t374 = cos(qJ(5));
	t372 = sin(qJ(5));
	t345 = -t364 * t392 + t384;
	t333 = t364 * t387 + (t366 * t386 - t398) * t363;
	t332 = t360 * t391 + (t364 * t369 + t366 * t390) * t361;
	t330 = t352 * t367 - t361 * t379;
	t329 = t351 * t367 - t361 * t378;
	t328 = t338 * t361 + t349 * t367 + pkin(1);
	t325 = -t333 * t365 + t371 * t345;
	t324 = t333 * t371 + t365 * t345;
	t320 = t324 * t368 - t362 * t332;
	t319 = t324 * t362 + t368 * t332;
	t1 = [t319 * t372 - t407 * t374, t374 * t319 + t407 * t372, t317 * t375 - t326 * t373, (t368 * t329 - t410 * t362) * t375 + (-t368 * t330 + t411 * t362) * t373 + t405 * t362 + t328 * t368 + 0; -t320 * t372 + t412 * t374, -t320 * t374 - t412 * t372, t327 * t373 - t375 * t406, (t362 * t329 + t410 * t368) * t375 + (-t330 * t362 - t411 * t368) * t373 - t405 * t368 + t328 * t362 + 0; t372 * t325 + t404 * t374, t374 * t325 - t404 * t372, t373 * t335 - t375 * t380, (-t339 * t371 + t376 * t365) * t375 + (t371 * t340 - t403 * t365) * t373 - t381 * t365 + t337 * t371 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:43
	% EndTime: 2020-11-04 20:52:44
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (455->101), mult. (1190->174), div. (0->0), fcn. (1515->18), ass. (0->90)
	t472 = sin(pkin(14));
	t478 = cos(pkin(14));
	t481 = cos(pkin(8));
	t507 = t478 * t481;
	t466 = pkin(4) * t507 + t472 * pkin(10);
	t476 = sin(pkin(7));
	t475 = sin(pkin(8));
	t482 = cos(pkin(7));
	t510 = t475 * t482;
	t449 = pkin(4) * t510 + t476 * t466;
	t477 = sin(pkin(6));
	t483 = cos(pkin(6));
	t512 = t472 * t481;
	t464 = pkin(4) * t512 - t478 * pkin(10);
	t473 = sin(pkin(13));
	t479 = cos(pkin(13));
	t509 = t476 * t475;
	t493 = pkin(4) * t509 - t466 * t482;
	t521 = t473 * t464 + t493 * t479;
	t531 = t449 * t477 + t521 * t483;
	t465 = -t472 * pkin(4) + pkin(10) * t507;
	t448 = pkin(10) * t510 + t476 * t465;
	t463 = t478 * pkin(4) + pkin(10) * t512;
	t492 = pkin(10) * t509 - t465 * t482;
	t490 = t473 * t463 + t492 * t479;
	t530 = t448 * t477 + t490 * t483;
	t474 = sin(pkin(12));
	t480 = cos(pkin(12));
	t503 = t482 * t481;
	t459 = t478 * t503 - t509;
	t438 = t459 * t479 - t473 * t512;
	t504 = t481 * t476;
	t457 = t478 * t504 + t510;
	t494 = t438 * t483 - t477 * t457;
	t498 = t473 * t459 + t479 * t512;
	t418 = t494 * t474 + t480 * t498;
	t462 = t478 * t475 * pkin(9) - t472 * pkin(3);
	t470 = t481 * pkin(9) + qJ(3);
	t446 = -t462 * t476 + t470 * t482 + qJ(2);
	t447 = t462 * t482 + t476 * t470;
	t514 = t472 * t475;
	t461 = t478 * pkin(3) + pkin(9) * t514 + pkin(2);
	t500 = t447 * t479 - t473 * t461;
	t527 = t477 * t446 + t500 * t483;
	t458 = t478 * t510 + t504;
	t499 = t473 * t458 + t479 * t514;
	t526 = t474 * t499;
	t525 = t480 * t499;
	t485 = sin(qJ(5));
	t488 = cos(qJ(5));
	t523 = -pkin(5) * t485 + pkin(11) * t488;
	t506 = t479 * t482;
	t456 = t472 * t506 + t478 * t473;
	t513 = t472 * t476;
	t445 = t456 * t483 - t477 * t513;
	t515 = t472 * t473;
	t454 = -t478 * t479 + t482 * t515;
	t431 = t445 * t474 + t480 * t454;
	t486 = sin(qJ(4));
	t489 = cos(qJ(4));
	t522 = t418 * t486 + t431 * t489;
	t419 = -t474 * t498 + t480 * t494;
	t440 = t479 * t504 + (t478 * t506 - t515) * t475;
	t455 = -t478 * t509 + t503;
	t430 = -t440 * t477 + t483 * t455;
	t502 = t485 * t430;
	t501 = t488 * t430;
	t497 = pkin(5) * t488 + pkin(11) * t485;
	t432 = t445 * t480 - t474 * t454;
	t496 = t419 * t486 + t432 * t489;
	t427 = t438 * t477 + t483 * t457;
	t444 = t456 * t477 + t483 * t513;
	t495 = t427 * t486 + t489 * t444;
	t487 = cos(qJ(6));
	t484 = sin(qJ(6));
	t452 = t477 * t455;
	t437 = t464 * t479 - t473 * t493;
	t436 = t463 * t479 - t473 * t492;
	t435 = t447 * t473 + t461 * t479 + pkin(1);
	t433 = (t458 * t479 - t473 * t514) * t483 + t452;
	t429 = t440 * t483 + t452;
	t422 = t427 * t489 - t486 * t444;
	t421 = t429 * t480 - t526;
	t420 = t429 * t474 + t525;
	t417 = t419 * t489 - t432 * t486;
	t416 = t418 * t489 - t431 * t486;
	t415 = t495 * t488 + t502;
	t414 = -t421 * t485 + t496 * t488;
	t413 = -t420 * t485 + t522 * t488;
	t1 = [-t413 * t487 + t416 * t484, t413 * t484 + t416 * t487, -t488 * t420 - t522 * t485, (-t497 * t418 - t480 * t437 + t531 * t474) * t486 + (-t497 * t431 + t480 * t436 - t530 * t474) * t489 + t527 * t474 + t435 * t480 + 0 - t523 * (t433 * t474 + t525); t414 * t487 - t484 * t417, -t414 * t484 - t417 * t487, t421 * t488 + t496 * t485, (t497 * t419 - t437 * t474 - t531 * t480) * t486 + (t497 * t432 + t474 * t436 + t530 * t480) * t489 - t527 * t480 + t435 * t474 + 0 + t523 * (t433 * t480 - t526); t415 * t487 - t422 * t484, -t415 * t484 - t422 * t487, t495 * t485 - t501, (t497 * t427 + t483 * t449 - t521 * t477) * t486 + (t497 * t444 - t448 * t483 + t490 * t477) * t489 - pkin(11) * t501 + pkin(5) * t502 - t500 * t477 + t446 * t483 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end