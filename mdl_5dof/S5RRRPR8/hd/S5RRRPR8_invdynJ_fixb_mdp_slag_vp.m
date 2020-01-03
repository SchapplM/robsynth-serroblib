% Calculate vector of inverse dynamics joint torques for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:36
% EndTime: 2019-12-31 21:20:42
% DurationCPUTime: 4.34s
% Computational Cost: add. (2888->405), mult. (6359->504), div. (0->0), fcn. (4393->10), ass. (0->186)
t400 = cos(qJ(2));
t498 = cos(qJ(3));
t457 = t498 * t400;
t442 = qJD(1) * t457;
t396 = sin(qJ(3));
t397 = sin(qJ(2));
t472 = qJD(1) * t397;
t456 = t396 * t472;
t337 = -t442 + t456;
t390 = qJD(2) + qJD(3);
t395 = sin(qJ(5));
t399 = cos(qJ(5));
t317 = -t399 * t337 + t390 * t395;
t348 = t396 * t400 + t498 * t397;
t339 = t348 * qJD(1);
t506 = qJD(5) + t339;
t507 = t317 * t506;
t319 = t337 * t395 + t390 * t399;
t447 = t506 * t319;
t499 = pkin(7) + pkin(6);
t355 = t499 * t400;
t351 = qJD(1) * t355;
t344 = t498 * t351;
t354 = t499 * t397;
t349 = qJD(1) * t354;
t309 = -t396 * t349 + t344;
t471 = qJD(3) * t396;
t440 = pkin(2) * t471 - t309;
t341 = t396 * t351;
t310 = -t498 * t349 - t341;
t455 = qJD(3) * t498;
t476 = -pkin(2) * t455 - qJD(4) + t310;
t394 = qJ(2) + qJ(3);
t387 = sin(t394);
t388 = cos(t394);
t474 = t388 * pkin(3) + t387 * qJ(4);
t389 = qJDD(2) + qJDD(3);
t383 = t389 * qJ(4);
t505 = -t390 * qJD(4) - t383;
t493 = qJD(2) * pkin(2);
t345 = -t349 + t493;
t305 = -t498 * t345 + t341;
t465 = qJD(4) + t305;
t479 = t396 * t397;
t436 = t390 * t479;
t452 = qJDD(1) * t498;
t461 = qJDD(1) * t400;
t443 = -t390 * t442 - t396 * t461 - t397 * t452;
t291 = qJD(1) * t436 + t443;
t287 = -qJDD(5) + t291;
t381 = -t498 * pkin(2) - pkin(3);
t373 = -pkin(8) + t381;
t495 = t337 * pkin(4);
t504 = (t495 + t440) * t506 - t373 * t287;
t458 = qJD(2) * t499;
t350 = t397 * t458;
t352 = t400 * t458;
t280 = t498 * t350 + t396 * t352 + t354 * t455 + t355 * t471;
t321 = -t396 * t354 + t498 * t355;
t398 = sin(qJ(1));
t401 = cos(qJ(1));
t437 = g(1) * t398 - g(2) * t401;
t503 = -t280 * t390 + t321 * t389 + t437 * t387;
t281 = t321 * qJD(3) - t396 * t350 + t498 * t352;
t320 = t498 * t354 + t396 * t355;
t502 = t281 * t390 + t320 * t389 - t437 * t388;
t501 = t339 ^ 2;
t500 = pkin(3) + pkin(8);
t497 = pkin(2) * t400;
t496 = pkin(3) * t389;
t377 = g(3) * t387;
t378 = g(3) * t388;
t494 = t339 * pkin(4);
t312 = t390 * t348;
t462 = qJDD(1) * t397;
t435 = t396 * t462 - t400 * t452;
t292 = t312 * qJD(1) + t435;
t467 = qJD(5) * t399;
t459 = t395 * t292 + t337 * t467 + t399 * t389;
t468 = qJD(5) * t395;
t271 = -t390 * t468 + t459;
t492 = t271 * t399;
t491 = t287 * t395;
t347 = -t457 + t479;
t382 = pkin(1) + t497;
t429 = -qJ(4) * t348 - t382;
t294 = t500 * t347 + t429;
t490 = t294 * t287;
t306 = t396 * t345 + t344;
t489 = t306 * t390;
t488 = t337 * t339;
t487 = t347 * t395;
t485 = t387 * t398;
t484 = t387 * t401;
t483 = t388 * t398;
t482 = t388 * t401;
t481 = t395 * t398;
t480 = t395 * t401;
t478 = t398 * t399;
t284 = t399 * t287;
t477 = t399 * t401;
t475 = t494 - t476;
t392 = t397 ^ 2;
t473 = -t400 ^ 2 + t392;
t353 = t382 * qJD(1);
t418 = -qJ(4) * t339 - t353;
t279 = t500 * t337 + t418;
t470 = qJD(5) * t279;
t469 = qJD(5) * t390;
t466 = t494 + t465;
t463 = qJD(1) * qJD(2);
t385 = t397 * t493;
t454 = t397 * t463;
t453 = t400 * t463;
t315 = qJDD(2) * pkin(2) + t499 * (-t453 - t462);
t316 = t499 * (-t454 + t461);
t444 = -t396 * t315 - t498 * t316 - t345 * t455 + t351 * t471;
t267 = t444 + t505;
t260 = -pkin(4) * t292 - t267;
t277 = -t500 * t390 + t466;
t265 = t277 * t395 + t279 * t399;
t451 = t260 * t399 - t265 * t337;
t450 = t395 * t506;
t449 = t399 * t506;
t302 = -qJ(4) * t390 - t306;
t283 = -t302 - t495;
t448 = t506 * t283;
t334 = pkin(2) * t454 - t382 * qJDD(1);
t410 = qJ(4) * t291 - qJD(4) * t339 + t334;
t254 = t500 * t292 + t410;
t446 = qJD(5) * t277 + t254;
t445 = -t498 * t315 + t396 * t316 + t345 * t471 + t351 * t455;
t439 = -pkin(2) * t397 - pkin(3) * t387;
t438 = g(1) * t401 + g(2) * t398;
t303 = pkin(3) * t339 + qJ(4) * t337;
t434 = -t470 + t378;
t431 = t500 * t287 - (t306 - t495) * t506;
t430 = qJDD(4) + t445;
t264 = t277 * t399 - t279 * t395;
t428 = t260 * t395 + t264 * t337 + (t339 * t399 + t467) * t283;
t299 = pkin(2) * t472 + t303;
t427 = t382 + t474;
t426 = -0.2e1 * pkin(1) * t463 - pkin(6) * qJDD(2);
t425 = t312 * t395 + t347 * t467;
t311 = -qJD(2) * t457 - t400 * t455 + t436;
t424 = qJ(4) * t311 - qJD(4) * t348 + t385;
t421 = -g(1) * t482 - g(2) * t483 - t377 - t444;
t420 = -g(1) * t484 - g(2) * t485 + t378 + t445;
t300 = t348 * pkin(4) + t320;
t419 = t260 * t347 + t283 * t312 + t300 * t287;
t417 = -t438 * t388 - t377;
t404 = qJD(2) ^ 2;
t416 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t404 + t437;
t405 = qJD(1) ^ 2;
t415 = pkin(1) * t405 - pkin(6) * qJDD(1) + t438;
t414 = -t353 * t337 - t421;
t413 = t353 * t339 - t420;
t297 = pkin(3) * t337 + t418;
t412 = t297 * t339 + qJDD(4) + t420;
t411 = -t297 * t337 + t421 - t505;
t409 = -t443 + (t337 - t456) * t390;
t332 = t339 * pkin(8);
t408 = (-qJD(5) * t373 + t299 + t332) * t506 + t417;
t407 = (qJD(5) * t500 + t303 + t332) * t506 + t417;
t286 = t399 * t292;
t272 = t319 * qJD(5) + t389 * t395 - t286;
t406 = t506 * t337 * MDP(26) + MDP(11) * t488 + ((-t272 - t447) * t399 + (-t271 + t507) * t395) * MDP(23) + (-t395 * t447 + t492) * MDP(22) + (t319 * t337 - t450 * t506 - t284) * MDP(24) + (-t317 * t337 - t449 * t506 + t491) * MDP(25) + t409 * MDP(13) - t435 * MDP(14) + (-t337 ^ 2 + t501) * MDP(12) + t389 * MDP(15);
t375 = pkin(2) * t396 + qJ(4);
t357 = qJ(4) * t482;
t356 = qJ(4) * t483;
t331 = -t387 * t481 + t477;
t330 = t387 * t478 + t480;
t329 = t387 * t480 + t478;
t328 = t387 * t477 - t481;
t304 = pkin(3) * t347 + t429;
t301 = -t347 * pkin(4) + t321;
t298 = -pkin(3) * t390 + t465;
t273 = pkin(3) * t312 + t424;
t270 = -t311 * pkin(4) + t281;
t269 = -pkin(4) * t312 - t280;
t268 = t430 - t496;
t266 = t500 * t312 + t424;
t263 = pkin(3) * t292 + t410;
t259 = -pkin(4) * t291 - t500 * t389 + t430;
t256 = t399 * t259;
t1 = [(t263 * t304 - t267 * t321 + t268 * t320 + t297 * t273 + t302 * t280 + t298 * t281 + (-g(1) * t499 - g(2) * t427) * t401 + (g(1) * t427 - g(2) * t499) * t398) * MDP(21) + (t271 * t348 - t287 * t487 - t311 * t319 + t425 * t506) * MDP(24) + (g(1) * t330 - g(2) * t328 + t265 * t311 + t269 * t319 + t301 * t271 + (-(qJD(5) * t300 + t266) * t506 + t490 - t446 * t348 + t283 * qJD(5) * t347) * t399 + (-(-qJD(5) * t294 + t270) * t506 - (t259 - t470) * t348 + t419) * t395) * MDP(28) + (-g(1) * t331 - g(2) * t329 + t256 * t348 - t264 * t311 + t269 * t317 + t301 * t272 + (-t254 * t348 - t266 * t506 + t490) * t395 + (t270 * t506 - t419) * t399 + ((-t294 * t399 - t300 * t395) * t506 - t265 * t348 + t283 * t487) * qJD(5)) * MDP(27) + (-t347 * t284 - t272 * t348 + t311 * t317 + (t312 * t399 - t347 * t468) * t506) * MDP(25) + (-t287 * t348 - t311 * t506) * MDP(26) + qJDD(1) * MDP(1) + (t426 * t397 + t416 * t400) * MDP(9) + (-t416 * t397 + t426 * t400) * MDP(10) + (t267 * t347 + t268 * t348 + t280 * t337 + t281 * t339 - t291 * t320 - t292 * t321 - t298 * t311 + t302 * t312 - t438) * MDP(18) + t438 * MDP(3) + t437 * MDP(2) + 0.2e1 * (t397 * t461 - t473 * t463) * MDP(5) + (t271 * t487 + t425 * t319) * MDP(22) + ((-t317 * t395 + t319 * t399) * t312 + (t492 - t272 * t395 + (-t317 * t399 - t319 * t395) * qJD(5)) * t347) * MDP(23) + (qJDD(1) * t392 + 0.2e1 * t397 * t453) * MDP(4) + (-t291 * t348 - t311 * t339) * MDP(11) + (t291 * t347 - t292 * t348 + t311 * t337 - t312 * t339) * MDP(12) + (-t312 * t390 - t347 * t389) * MDP(14) + (-t311 * t390 + t348 * t389) * MDP(13) + (qJDD(2) * t397 + t400 * t404) * MDP(6) + (qJDD(2) * t400 - t397 * t404) * MDP(7) + (-t292 * t382 - t312 * t353 + t334 * t347 + t337 * t385 - t502) * MDP(16) + (-t263 * t347 - t273 * t337 - t292 * t304 - t297 * t312 + t502) * MDP(19) + (t291 * t382 + t311 * t353 + t334 * t348 + t339 * t385 - t503) * MDP(17) + (-t263 * t348 - t273 * t339 + t291 * t304 + t297 * t311 + t503) * MDP(20); t406 + qJDD(2) * MDP(8) + (t299 * t337 + t440 * t390 + (-pkin(3) + t381) * t389 + t412) * MDP(19) + (t375 * t271 + t475 * t319 + t408 * t399 + (-t448 - t504) * t395 + t451) * MDP(28) + (t375 * t272 + t475 * t317 + t408 * t395 + t504 * t399 + t428) * MDP(27) + (-g(3) * t400 + t415 * t397) * MDP(9) + (-t267 * t375 + t268 * t381 - t297 * t299 - g(1) * (t439 * t401 + t357) - g(2) * (t439 * t398 + t356) - g(3) * (t474 + t497) + t476 * t302 + t440 * t298) * MDP(21) + (g(3) * t397 + t415 * t400) * MDP(10) + (t310 * t390 + (-t339 * t472 - t389 * t396 - t390 * t455) * pkin(2) + t414) * MDP(17) + (t309 * t390 + (-t337 * t472 + t498 * t389 - t390 * t471) * pkin(2) + t413) * MDP(16) + (t299 * t339 + t375 * t389 - t476 * t390 + t411) * MDP(20) + (-t291 * t381 - t292 * t375 + (-t302 + t440) * t339 + (t298 + t476) * t337) * MDP(18) + MDP(6) * t462 + MDP(7) * t461 + (-t397 * t400 * MDP(4) + t473 * MDP(5)) * t405; t406 + (qJ(4) * t272 + t466 * t317 + t407 * t395 + t431 * t399 + t428) * MDP(27) + (-t267 * qJ(4) - t268 * pkin(3) - t297 * t303 - t298 * t306 - g(1) * (-pkin(3) * t484 + t357) - g(2) * (-pkin(3) * t485 + t356) - g(3) * t474 - t465 * t302) * MDP(21) + (t413 + t489) * MDP(16) + (pkin(3) * t291 - qJ(4) * t292 + (-t302 - t306) * t339 + (t298 - t465) * t337) * MDP(18) + (t303 * t337 + t412 - t489 - 0.2e1 * t496) * MDP(19) + (t303 * t339 + t465 * t390 + t383 + t411) * MDP(20) + (qJ(4) * t271 + t466 * t319 + (-t448 - t431) * t395 + t407 * t399 + t451) * MDP(28) + (-t305 * t390 + t414) * MDP(17); t409 * MDP(18) + (t389 - t488) * MDP(19) + (-t390 ^ 2 - t501) * MDP(20) + (t302 * t390 + t412 - t496) * MDP(21) + (-t317 * t390 - t284) * MDP(27) + (-t319 * t390 + t491) * MDP(28) + (-MDP(27) * t450 - MDP(28) * t449) * t506; t319 * t317 * MDP(22) + (-t317 ^ 2 + t319 ^ 2) * MDP(23) + (t459 + t507) * MDP(24) + (t286 + t447) * MDP(25) - t287 * MDP(26) + (-g(1) * t328 - g(2) * t330 + t265 * t506 - t283 * t319 + t256) * MDP(27) + (g(1) * t329 - g(2) * t331 + t264 * t506 + t283 * t317) * MDP(28) + (-MDP(25) * t469 + t434 * MDP(27) - t446 * MDP(28)) * t399 + (-MDP(24) * t469 + (-qJD(5) * t337 - t389) * MDP(25) - t446 * MDP(27) + (-t259 - t434) * MDP(28)) * t395;];
tau = t1;
