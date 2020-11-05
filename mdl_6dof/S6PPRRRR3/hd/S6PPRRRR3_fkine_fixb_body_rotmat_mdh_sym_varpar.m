% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:56
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:44
	% EndTime: 2020-11-04 20:56:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:44
	% EndTime: 2020-11-04 20:56:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t258 = cos(pkin(13));
	t257 = sin(pkin(13));
	t1 = [t258, -t257, 0, 0; t257, t258, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:44
	% EndTime: 2020-11-04 20:56:44
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t260 = sin(pkin(13));
	t261 = sin(pkin(6));
	t268 = t260 * t261;
	t264 = cos(pkin(6));
	t267 = t260 * t264;
	t263 = cos(pkin(13));
	t266 = t263 * t261;
	t265 = t263 * t264;
	t262 = cos(pkin(14));
	t259 = sin(pkin(14));
	t1 = [-t259 * t267 + t263 * t262, -t263 * t259 - t262 * t267, t268, t263 * pkin(1) + qJ(2) * t268 + 0; t259 * t265 + t260 * t262, -t260 * t259 + t262 * t265, -t266, t260 * pkin(1) - qJ(2) * t266 + 0; t261 * t259, t261 * t262, t264, t264 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:44
	% EndTime: 2020-11-04 20:56:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t280 = sin(pkin(7));
	t297 = pkin(9) * t280;
	t278 = sin(pkin(14));
	t283 = cos(pkin(13));
	t296 = t278 * t283;
	t279 = sin(pkin(13));
	t295 = t279 * t278;
	t285 = cos(pkin(6));
	t294 = t280 * t285;
	t281 = sin(pkin(6));
	t293 = t281 * t280;
	t282 = cos(pkin(14));
	t284 = cos(pkin(7));
	t292 = t282 * t284;
	t291 = t283 * t285;
	t290 = t284 * t281;
	t289 = t285 * t284;
	t275 = -t278 * pkin(2) + t282 * t297;
	t276 = t284 * pkin(9) + qJ(2);
	t288 = t275 * t285 + t281 * t276;
	t287 = cos(qJ(3));
	t286 = sin(qJ(3));
	t274 = t282 * pkin(2) + t278 * t297 + pkin(1);
	t273 = t283 * t282 - t285 * t295;
	t272 = t278 * t291 + t279 * t282;
	t271 = t282 * t289 - t293;
	t270 = -t271 * t279 - t284 * t296;
	t269 = t271 * t283 - t284 * t295;
	t1 = [t270 * t286 + t273 * t287, t270 * t287 - t273 * t286, (t282 * t294 + t290) * t279 + t280 * t296, t274 * t283 + t288 * t279 + 0; t269 * t286 + t272 * t287, t269 * t287 - t272 * t286, (-t282 * t291 + t295) * t280 - t283 * t290, t274 * t279 - t288 * t283 + 0; t286 * t294 + (t278 * t287 + t286 * t292) * t281, t287 * t294 + (-t278 * t286 + t287 * t292) * t281, -t282 * t293 + t289, -t275 * t281 + t276 * t285 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:44
	% EndTime: 2020-11-04 20:56:45
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (133->59), mult. (335->107), div. (0->0), fcn. (424->14), ass. (0->58)
	t327 = sin(pkin(8));
	t357 = pkin(10) * t327;
	t330 = cos(pkin(14));
	t356 = t330 * pkin(3);
	t325 = sin(pkin(14));
	t328 = sin(pkin(7));
	t355 = t325 * t328;
	t329 = sin(pkin(6));
	t354 = t325 * t329;
	t333 = cos(pkin(7));
	t353 = t325 * t333;
	t334 = cos(pkin(6));
	t352 = t325 * t334;
	t351 = t328 * t334;
	t350 = t329 * t328;
	t349 = t333 * t329;
	t348 = t334 * t333;
	t347 = t328 * t357;
	t346 = t325 * t357;
	t345 = t330 * t357;
	t313 = t330 * t348 - t350;
	t326 = sin(pkin(13));
	t331 = cos(pkin(13));
	t305 = t313 * t326 + t331 * t353;
	t310 = t326 * t352 - t331 * t330;
	t336 = sin(qJ(3));
	t338 = cos(qJ(3));
	t344 = t305 * t338 - t310 * t336;
	t306 = t313 * t331 - t326 * t353;
	t315 = t326 * t330 + t331 * t352;
	t343 = t306 * t338 - t315 * t336;
	t332 = cos(pkin(8));
	t322 = t332 * pkin(10) + pkin(9);
	t309 = t328 * t322 * t330 - t325 * pkin(2);
	t320 = t322 * t333 + qJ(2);
	t342 = t309 * t334 + t320 * t329;
	t319 = t333 * t356 + t346;
	t341 = pkin(3) * t350 - t319 * t334;
	t314 = t330 * t349 + t351;
	t340 = -t314 * t338 + t336 * t354;
	t317 = -t325 * pkin(3) + t333 * t345;
	t339 = -t317 * t334 + t329 * t347;
	t337 = cos(qJ(4));
	t335 = sin(qJ(4));
	t318 = pkin(3) * t353 - t345;
	t316 = t333 * t346 + t356;
	t312 = t330 * t350 - t348;
	t311 = t330 * t351 + t349;
	t308 = t330 * pkin(2) + t322 * t355 + pkin(1);
	t307 = t336 * t314 + t338 * t354;
	t304 = t311 * t331 - t326 * t355;
	t303 = t311 * t326 + t331 * t355;
	t302 = t306 * t336 + t315 * t338;
	t301 = t305 * t336 + t310 * t338;
	t300 = -t327 * t312 - t340 * t332;
	t299 = -t304 * t327 + t343 * t332;
	t298 = t303 * t327 - t344 * t332;
	t1 = [t298 * t335 - t337 * t301, t298 * t337 + t335 * t301, t303 * t332 + t344 * t327, (t331 * t316 - t339 * t326) * t338 + (-t331 * t318 + t341 * t326) * t336 + t342 * t326 + t308 * t331 + 0; t299 * t335 + t302 * t337, t299 * t337 - t302 * t335, -t304 * t332 - t343 * t327, (t326 * t316 + t339 * t331) * t338 + (-t326 * t318 - t341 * t331) * t336 - t342 * t331 + t308 * t326 + 0; t300 * t335 + t307 * t337, t300 * t337 - t307 * t335, -t332 * t312 + t340 * t327, (-t317 * t329 - t334 * t347) * t338 + (pkin(3) * t351 + t319 * t329) * t336 - t309 * t329 + t320 * t334 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:45
	% EndTime: 2020-11-04 20:56:45
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (256->114), mult. (689->208), div. (0->0), fcn. (843->16), ass. (0->88)
	t399 = sin(pkin(8));
	t454 = pkin(10) * t399;
	t402 = cos(pkin(14));
	t401 = sin(pkin(6));
	t405 = cos(pkin(7));
	t438 = t405 * t401;
	t400 = sin(pkin(7));
	t406 = cos(pkin(6));
	t445 = t400 * t406;
	t374 = t402 * t445 + t438;
	t398 = sin(pkin(13));
	t403 = cos(pkin(13));
	t397 = sin(pkin(14));
	t451 = t397 * t400;
	t367 = t374 * t398 + t403 * t451;
	t453 = t367 * t399;
	t368 = t374 * t403 - t398 * t451;
	t452 = t368 * t399;
	t450 = t397 * t401;
	t404 = cos(pkin(8));
	t449 = t397 * t404;
	t448 = t397 * t405;
	t447 = t397 * t406;
	t437 = t406 * t405;
	t444 = t401 * t400;
	t375 = t402 * t444 - t437;
	t446 = t399 * t375;
	t443 = t402 * t404;
	t442 = t402 * t405;
	t441 = t404 * t405;
	t408 = sin(qJ(4));
	t440 = t404 * t408;
	t409 = sin(qJ(3));
	t439 = t404 * t409;
	t411 = cos(qJ(4));
	t436 = t409 * t411;
	t435 = t400 * t454;
	t434 = t397 * t454;
	t433 = t402 * t454;
	t432 = t409 * t450;
	t431 = t397 * t441;
	t430 = t404 * t444;
	t429 = t404 * t445;
	t428 = t402 * t441;
	t427 = pkin(4) * t408 - pkin(11) * t411;
	t376 = t402 * t437 - t444;
	t369 = t376 * t398 + t403 * t448;
	t373 = t398 * t447 - t403 * t402;
	t412 = cos(qJ(3));
	t426 = t369 * t412 - t373 * t409;
	t370 = t376 * t403 - t398 * t448;
	t378 = t398 * t402 + t403 * t447;
	t425 = t370 * t412 - t378 * t409;
	t393 = t404 * pkin(10) + pkin(9);
	t372 = t400 * t393 * t402 - t397 * pkin(2);
	t391 = t393 * t405 + qJ(2);
	t424 = t372 * t406 + t391 * t401;
	t390 = pkin(3) * t442 + t434;
	t423 = pkin(3) * t444 - t390 * t406;
	t389 = pkin(4) * t442 + pkin(11) * t449;
	t422 = pkin(4) * t444 - t389 * t406;
	t388 = -pkin(4) * t449 + pkin(11) * t442;
	t421 = pkin(11) * t444 - t388 * t406;
	t420 = t427 * t399;
	t377 = t402 * t438 + t445;
	t419 = -t377 * t412 + t432;
	t384 = pkin(4) * t428 + t397 * pkin(11);
	t418 = pkin(4) * t430 - t384 * t406;
	t380 = -t397 * pkin(3) + t405 * t433;
	t417 = -t380 * t406 + t401 * t435;
	t383 = -t397 * pkin(4) + pkin(11) * t428;
	t416 = pkin(11) * t430 - t383 * t406;
	t415 = -(t373 * t439 + t453) * t408 + (t369 * t440 + t411 * t373) * t412 + t369 * t436;
	t414 = -(t378 * t439 + t452) * t408 + (t370 * t440 + t411 * t378) * t412 + t370 * t436;
	t413 = -(t404 * t432 + t446) * t408 + (t377 * t440 + t411 * t450) * t412 + t377 * t436;
	t410 = cos(qJ(5));
	t407 = sin(qJ(5));
	t387 = pkin(3) * t448 - t433;
	t386 = pkin(4) * t448 - pkin(11) * t443;
	t385 = pkin(4) * t443 + pkin(11) * t448;
	t382 = pkin(4) * t431 - t402 * pkin(11);
	t381 = t402 * pkin(4) + pkin(11) * t431;
	t379 = t402 * pkin(3) + t405 * t434;
	t371 = t402 * pkin(2) + t393 * t451 + pkin(1);
	t364 = -t375 * t404 + t419 * t399;
	t359 = t368 * t404 + t425 * t399;
	t358 = t367 * t404 + t426 * t399;
	t1 = [t358 * t407 - t415 * t410, t358 * t410 + t415 * t407, (t426 * t404 - t453) * t411 - t408 * (t369 * t409 + t373 * t412), ((t403 * t381 - t416 * t398) * t411 + (-t403 * t382 + t418 * t398) * t408 - t417 * t398 + t403 * t379) * t412 + ((-t403 * t386 + t422 * t398) * t411 + (-t403 * t385 + t421 * t398) * t408 + t423 * t398 - t403 * t387) * t409 + t424 * t398 + t371 * t403 + 0 + t427 * t453; -t359 * t407 + t414 * t410, -t410 * t359 - t414 * t407, (-t425 * t404 + t452) * t411 + (t370 * t409 + t378 * t412) * t408, ((t398 * t381 + t416 * t403) * t411 + (-t398 * t382 - t418 * t403) * t408 + t417 * t403 + t398 * t379) * t412 + ((-t398 * t386 - t422 * t403) * t411 + (-t398 * t385 - t421 * t403) * t408 - t423 * t403 - t398 * t387) * t409 - t424 * t403 + t371 * t398 + 0 - t368 * t420; t407 * t364 + t413 * t410, t364 * t410 - t413 * t407, (t419 * t404 + t446) * t411 + (t409 * t377 + t412 * t450) * t408, ((-pkin(11) * t429 - t383 * t401) * t411 + (pkin(4) * t429 + t384 * t401) * t408 - t380 * t401 - t406 * t435) * t412 + ((pkin(4) * t445 + t389 * t401) * t411 + (pkin(11) * t445 + t388 * t401) * t408 + t390 * t401 + pkin(3) * t445) * t409 - t372 * t401 + t391 * t406 + qJ(1) + 0 - t375 * t420; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:45
	% EndTime: 2020-11-04 20:56:46
	% DurationCPUTime: 1.19s
	% Computational Cost: add. (435->150), mult. (1204->263), div. (0->0), fcn. (1489->18), ass. (0->107)
	t504 = sin(pkin(8));
	t513 = sin(qJ(5));
	t565 = t504 * t513;
	t517 = cos(qJ(5));
	t573 = pkin(12) * t517;
	t514 = sin(qJ(4));
	t509 = cos(pkin(8));
	t536 = pkin(5) * t517 + pkin(12) * t513;
	t585 = t536 * t509;
	t586 = t514 * t585;
	t580 = -pkin(5) * t565 + t504 * t573 + t586;
	t554 = t509 * t517;
	t518 = cos(qJ(4));
	t575 = pkin(11) * t518;
	t578 = pkin(5) * t513;
	t587 = -pkin(12) * t554 + t509 * t578 + (t514 * (pkin(4) + t536) - t575) * t504;
	t506 = sin(pkin(6));
	t510 = cos(pkin(7));
	t511 = cos(pkin(6));
	t552 = t510 * t511;
	t505 = sin(pkin(7));
	t507 = cos(pkin(14));
	t562 = t505 * t507;
	t478 = -t506 * t562 + t552;
	t502 = sin(pkin(14));
	t515 = sin(qJ(3));
	t560 = t506 * t515;
	t541 = t502 * t560;
	t467 = t478 * t509 + t504 * t541;
	t558 = t507 * t510;
	t561 = t505 * t511;
	t480 = t506 * t558 + t561;
	t482 = t514 * t554 - t565;
	t519 = cos(qJ(3));
	t468 = t478 * t504 - t509 * t541;
	t550 = t515 * t518;
	t523 = t468 * t514 + t480 * t550;
	t570 = t502 * t506;
	t540 = t518 * t570;
	t583 = t523 * t517 + (t480 * t482 + t517 * t540) * t519 + t513 * t467;
	t477 = t506 * t510 + t507 * t561;
	t503 = sin(pkin(13));
	t508 = cos(pkin(13));
	t571 = t502 * t505;
	t470 = t477 * t508 - t503 * t571;
	t466 = t470 * t509;
	t563 = t505 * t506;
	t479 = t507 * t552 - t563;
	t568 = t502 * t510;
	t472 = t479 * t508 - t503 * t568;
	t567 = t502 * t511;
	t481 = t503 * t507 + t508 * t567;
	t555 = t509 * t515;
	t463 = t470 * t504 + t481 * t555;
	t524 = -t463 * t514 + t472 * t550;
	t549 = t517 * t518;
	t564 = t504 * t515;
	t582 = t524 * t517 + (t472 * t482 + t481 * t549) * t519 - (-t481 * t564 + t466) * t513;
	t469 = t477 * t503 + t508 * t571;
	t465 = t469 * t509;
	t471 = t479 * t503 + t508 * t568;
	t476 = t503 * t567 - t507 * t508;
	t462 = t469 * t504 + t476 * t555;
	t525 = -t462 * t514 + t471 * t550;
	t581 = t525 * t517 + (t471 * t482 + t476 * t549) * t519 - (-t476 * t564 + t465) * t513;
	t576 = pkin(10) * t504;
	t572 = t502 * t504;
	t569 = t502 * t509;
	t559 = t507 * t509;
	t557 = t509 * t510;
	t556 = t509 * t514;
	t553 = t509 * t518;
	t551 = t514 * t515;
	t546 = t510 * t576;
	t539 = t502 * t557;
	t538 = t509 * t563;
	t537 = t507 * t557;
	t535 = -t573 + t578;
	t498 = pkin(10) * t509 + pkin(9);
	t475 = -t502 * pkin(2) + t498 * t562;
	t496 = t498 * t510 + qJ(2);
	t531 = t475 * t511 + t496 * t506;
	t494 = pkin(3) * t558 + pkin(10) * t572;
	t530 = pkin(3) * t563 - t494 * t511;
	t493 = pkin(4) * t558 + pkin(11) * t569;
	t529 = pkin(4) * t563 - t493 * t511;
	t492 = -pkin(4) * t569 + pkin(11) * t558;
	t528 = pkin(11) * t563 - t492 * t511;
	t488 = pkin(4) * t537 + pkin(11) * t502;
	t522 = pkin(4) * t538 - t488 * t511;
	t484 = -t502 * pkin(3) + t507 * t546;
	t521 = -t484 * t511 + t563 * t576;
	t487 = -t502 * pkin(4) + pkin(11) * t537;
	t520 = pkin(11) * t538 - t487 * t511;
	t516 = cos(qJ(6));
	t512 = sin(qJ(6));
	t491 = pkin(3) * t568 - t507 * t576;
	t490 = pkin(4) * t568 - pkin(11) * t559;
	t489 = pkin(4) * t559 + pkin(11) * t568;
	t486 = pkin(4) * t539 - pkin(11) * t507;
	t485 = pkin(4) * t507 + pkin(11) * t539;
	t483 = pkin(3) * t507 + t502 * t546;
	t474 = pkin(2) * t507 + t498 * t571 + pkin(1);
	t457 = (t480 * t553 - t514 * t570) * t519 + t468 * t518 - t480 * t551;
	t456 = (t472 * t553 - t481 * t514) * t519 - t463 * t518 - t472 * t551;
	t455 = (t471 * t553 - t476 * t514) * t519 - t462 * t518 - t471 * t551;
	t1 = [t455 * t512 - t581 * t516, t516 * t455 + t581 * t512, ((-t471 * t556 - t476 * t518) * t519 - t525) * t513 - (t465 + (t471 * t519 - t476 * t515) * t504) * t517, ((-t508 * t486 + t522 * t503) * t514 + (-t536 * t476 + t508 * t485 - t520 * t503) * t518 - t521 * t503 + t508 * t483 - t580 * t471) * t519 + ((-t508 * t489 + t528 * t503) * t514 + (-t536 * t471 - t508 * t490 + t529 * t503) * t518 + t530 * t503 - t508 * t491 + t580 * t476) * t515 + t531 * t503 + t474 * t508 + 0 + t587 * t469; -t456 * t512 + t582 * t516, -t516 * t456 - t582 * t512, ((t472 * t556 + t481 * t518) * t519 + t524) * t513 + t517 * (t466 + (t472 * t519 - t481 * t515) * t504), ((-t503 * t486 - t522 * t508) * t514 + (t536 * t481 + t503 * t485 + t520 * t508) * t518 + t521 * t508 + t503 * t483 + t580 * t472) * t519 + ((-t503 * t489 - t528 * t508) * t514 + (t536 * t472 - t503 * t490 - t529 * t508) * t518 - t530 * t508 - t503 * t491 - t580 * t481) * t515 - t531 * t508 + t474 * t503 + 0 - t587 * t470; -t512 * t457 + t583 * t516, -t457 * t516 - t583 * t512, ((t480 * t556 + t540) * t519 + t523) * t513 - (-t480 * t504 * t519 + t467) * t517, ((-t576 + (pkin(4) * t514 - t575) * t509) * t561 + (t488 * t514 + (t536 * t502 - t487) * t518 - t484) * t506 + (-t535 * t504 + t586) * t480) * t519 + (pkin(11) * t515 * t561 + (-t502 * t585 + t492) * t560) * t514 + ((pkin(4) * t561 + t536 * t480) * t518 + pkin(3) * t561 + (t493 * t518 + t535 * t572 + t494) * t506) * t515 - t475 * t506 + t496 * t511 + qJ(1) + 0 + t587 * t478; 0, 0, 0, 1;];
	Tc_mdh = t1;
end