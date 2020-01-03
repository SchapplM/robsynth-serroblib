% Calculate vector of inverse dynamics joint torques for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR13_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:38
% EndTime: 2019-12-31 19:15:44
% DurationCPUTime: 4.31s
% Computational Cost: add. (2124->401), mult. (4220->548), div. (0->0), fcn. (2824->10), ass. (0->180)
t358 = sin(qJ(1));
t362 = cos(qJ(1));
t389 = g(1) * t358 - g(2) * t362;
t363 = -pkin(1) - pkin(6);
t336 = t363 * qJDD(1) + qJDD(2);
t339 = t363 * qJD(1) + qJD(2);
t361 = cos(qJ(3));
t357 = sin(qJ(3));
t430 = qJD(3) * t357;
t278 = -qJDD(3) * pkin(3) - t336 * t361 + t339 * t430;
t435 = qJD(1) * t357;
t341 = qJD(4) + t435;
t465 = g(3) * t357;
t370 = t389 * t361 - t465;
t478 = qJD(4) * pkin(7) * t341 + t278 + t370;
t356 = sin(qJ(4));
t360 = cos(qJ(4));
t419 = t360 * qJD(3);
t434 = qJD(1) * t361;
t315 = t356 * t434 - t419;
t359 = cos(qJ(5));
t431 = qJD(3) * t356;
t317 = t360 * t434 + t431;
t355 = sin(qJ(5));
t452 = t317 * t355;
t269 = t359 * t315 + t452;
t338 = qJD(5) + t341;
t477 = t269 * t338;
t386 = t315 * t355 - t359 * t317;
t476 = t338 * t386;
t475 = qJDD(2) - t389;
t325 = pkin(3) * t357 - pkin(7) * t361 + qJ(2);
t298 = t325 * qJD(1);
t324 = t357 * t339;
t305 = qJD(3) * pkin(7) + t324;
t261 = t298 * t356 + t305 * t360;
t254 = -pkin(8) * t315 + t261;
t423 = qJD(5) * t355;
t252 = t254 * t423;
t450 = t339 * t361;
t306 = -qJD(3) * pkin(3) - t450;
t277 = pkin(4) * t315 + t306;
t354 = qJ(4) + qJ(5);
t349 = sin(t354);
t350 = cos(t354);
t447 = t357 * t358;
t286 = t349 * t362 + t350 * t447;
t445 = t357 * t362;
t288 = -t349 * t358 + t350 * t445;
t464 = g(3) * t361;
t474 = g(1) * t286 - g(2) * t288 + t269 * t277 + t350 * t464 + t252;
t285 = -t349 * t447 + t350 * t362;
t287 = t349 * t445 + t350 * t358;
t424 = qJD(4) * t361;
t375 = -t356 * t424 - t357 * t419;
t415 = qJDD(1) * t361;
t262 = t375 * qJD(1) + qJD(4) * t419 + t356 * qJDD(3) + t360 * t415;
t391 = pkin(3) * t361 + pkin(7) * t357;
t313 = t391 * qJD(3) + qJD(2);
t267 = t313 * qJD(1) + t325 * qJDD(1);
t265 = t360 * t267;
t429 = qJD(3) * t361;
t279 = qJDD(3) * pkin(7) + t336 * t357 + t339 * t429;
t418 = qJD(1) * qJD(3);
t401 = t361 * t418;
t416 = qJDD(1) * t357;
t312 = qJDD(4) + t401 + t416;
t240 = pkin(4) * t312 - pkin(8) * t262 - t261 * qJD(4) - t279 * t356 + t265;
t403 = t356 * t430;
t263 = -qJD(1) * t403 + t317 * qJD(4) - t360 * qJDD(3) + t356 * t415;
t425 = qJD(4) * t360;
t410 = -t356 * t267 - t360 * t279 - t298 * t425;
t427 = qJD(4) * t356;
t377 = -t305 * t427 - t410;
t241 = -pkin(8) * t263 + t377;
t398 = t359 * t240 - t355 * t241;
t473 = -g(1) * t285 - g(2) * t287 + t277 * t386 + t349 * t464 + t398;
t308 = qJDD(5) + t312;
t472 = t308 * MDP(25) + (-t269 ^ 2 + t386 ^ 2) * MDP(22) - t269 * t386 * MDP(21);
t320 = t355 * t360 + t356 * t359;
t290 = t320 * t361;
t463 = pkin(1) * qJDD(1);
t470 = t463 - t475;
t469 = -qJD(5) * t360 - t425;
t468 = qJD(4) + qJD(5);
t397 = t262 * t355 + t359 * t263;
t245 = -t386 * qJD(5) + t397;
t467 = pkin(7) + pkin(8);
t365 = qJD(1) ^ 2;
t462 = qJ(2) * t365;
t260 = t360 * t298 - t305 * t356;
t253 = -pkin(8) * t317 + t260;
t251 = pkin(4) * t341 + t253;
t461 = t251 * t359;
t460 = t254 * t359;
t459 = t262 * t356;
t319 = t355 * t356 - t359 * t360;
t458 = t308 * t319;
t457 = t308 * t320;
t456 = t312 * t356;
t455 = t312 * t360;
t454 = t315 * t341;
t453 = t317 * t341;
t451 = t317 * t360;
t449 = t356 * t361;
t448 = t356 * t363;
t446 = t357 * t360;
t444 = t358 * t360;
t443 = t360 * t361;
t442 = t360 * t362;
t379 = t319 * t357;
t441 = -qJD(1) * t379 - t468 * t319;
t376 = t320 * qJD(1);
t440 = t468 * t320 + t357 * t376;
t321 = t391 * qJD(1);
t439 = t356 * t321 + t339 * t443;
t337 = t363 * t446;
t438 = t356 * t325 + t337;
t353 = t361 ^ 2;
t437 = t357 ^ 2 - t353;
t364 = qJD(3) ^ 2;
t436 = -t364 - t365;
t433 = qJD(3) * t315;
t432 = qJD(3) * t317;
t428 = qJD(3) * t363;
t426 = qJD(4) * t357;
t422 = qJD(5) * t359;
t417 = qJDD(1) * qJ(2);
t414 = qJDD(3) * t357;
t413 = pkin(8) * t446;
t411 = t359 * t262 - t355 * t263 - t315 * t422;
t406 = t361 * t428;
t409 = t356 * t313 + t325 * t425 + t360 * t406;
t408 = qJD(4) * t467;
t407 = t356 * t435;
t404 = t356 * t423;
t400 = pkin(4) - t448;
t396 = t341 * t363 + t305;
t395 = -qJD(4) * t298 - t279;
t394 = qJD(5) * t251 + t241;
t393 = qJD(1) + t426;
t392 = -t324 + (t407 + t427) * pkin(4);
t390 = g(1) * t362 + g(2) * t358;
t304 = t360 * t321;
t332 = t467 * t360;
t388 = qJD(5) * t332 - t339 * t449 + t304 + (pkin(4) * t361 + t413) * qJD(1) + t360 * t408;
t331 = t467 * t356;
t387 = pkin(8) * t407 + qJD(5) * t331 + t356 * t408 + t439;
t243 = t251 * t355 + t460;
t383 = t341 * t425 + t456;
t382 = -t341 * t427 + t455;
t380 = t338 * t319;
t378 = 0.2e1 * qJ(2) * t418 + qJDD(3) * t363;
t244 = -t317 * t423 + t411;
t374 = -t360 * t424 + t403;
t373 = 0.2e1 * qJD(1) * qJD(2) - t390;
t372 = -t336 + t389 + t462;
t371 = -pkin(7) * t312 + t341 * t306;
t367 = t373 + 0.2e1 * t417;
t366 = -t363 * t364 + t367;
t347 = qJDD(3) * t361;
t344 = -pkin(4) * t360 - pkin(3);
t314 = (pkin(4) * t356 - t363) * t361;
t311 = t360 * t325;
t302 = -t356 * t358 + t357 * t442;
t301 = t356 * t445 + t444;
t300 = t356 * t362 + t357 * t444;
t299 = -t356 * t447 + t442;
t295 = t360 * t313;
t291 = t319 * t361;
t280 = -t374 * pkin(4) + t357 * t428;
t276 = -pkin(8) * t449 + t438;
t266 = -pkin(8) * t443 + t400 * t357 + t311;
t250 = -t361 * t404 + (t468 * t443 - t403) * t359 + t375 * t355;
t249 = qJD(3) * t379 - t468 * t290;
t248 = pkin(4) * t263 + t278;
t247 = t374 * pkin(8) - t426 * t448 + t409;
t246 = t295 + (-t337 + (pkin(8) * t361 - t325) * t356) * qJD(4) + (t400 * t361 + t413) * qJD(3);
t242 = -t254 * t355 + t461;
t1 = [(qJDD(1) * t353 - 0.2e1 * t357 * t401) * MDP(7) + t390 * MDP(3) + (t366 * t357 + t378 * t361) * MDP(12) + (-t378 * t357 + t366 * t361) * MDP(13) + t367 * MDP(5) + ((t315 * t360 + t317 * t356) * t430 + (-t459 - t263 * t360 + (t315 * t356 - t451) * qJD(4)) * t361) * MDP(15) + (t470 * pkin(1) + (t373 + t417) * qJ(2)) * MDP(6) + t389 * MDP(2) + (-g(1) * t302 - g(2) * t300 + t295 * t341 + t311 * t312 + (t315 * t428 - t396 * t425 + t265) * t357 + (qJD(3) * t260 - t263 * t363 + t306 * t425) * t361 + ((-qJD(4) * t325 - t406) * t341 + t278 * t361 + (-qJD(3) * t306 - t312 * t363 + t395) * t357) * t356) * MDP(19) + (-t245 * t357 - t250 * t338 - t269 * t429 - t290 * t308) * MDP(24) + ((t246 * t359 - t247 * t355) * t338 + (t266 * t359 - t276 * t355) * t308 + t398 * t357 + t242 * t429 + t280 * t269 + t314 * t245 + t248 * t290 + t277 * t250 - g(1) * t288 - g(2) * t286 + ((-t266 * t355 - t276 * t359) * t338 - t243 * t357) * qJD(5)) * MDP(26) + (t312 * t357 + t341 * t429) * MDP(18) + (t308 * t357 + t338 * t429) * MDP(25) + ((-t341 * t419 + t262) * t357 + (t382 + t432) * t361) * MDP(16) + (-t243 * t429 + g(1) * t287 - g(2) * t285 + t314 * t244 - t248 * t291 + t277 * t249 + t252 * t357 - t280 * t386 + (-(-qJD(5) * t276 + t246) * t338 - t266 * t308 - t240 * t357) * t355 + (-(qJD(5) * t266 + t247) * t338 - t276 * t308 - t394 * t357) * t359) * MDP(27) + (t244 * t357 + t249 * t338 - t291 * t308 - t386 * t429) * MDP(23) + (-t244 * t290 + t245 * t291 - t249 * t269 + t250 * t386) * MDP(22) + (-t244 * t291 - t249 * t386) * MDP(21) + (-t361 * t364 - t414) * MDP(10) + (-0.2e1 * t463 + t475) * MDP(4) + qJDD(1) * MDP(1) + (t262 * t443 + t375 * t317) * MDP(14) + ((t341 * t431 - t263) * t357 + (-t383 - t433) * t361) * MDP(17) + 0.2e1 * (-t357 * t415 + t437 * t418) * MDP(8) + (-t409 * t341 - t438 * t312 + g(1) * t301 - g(2) * t299 + (t396 * t427 + (-t306 * t360 + t317 * t363) * qJD(3) + t410) * t357 + (-qJD(3) * t261 - t262 * t363 + t278 * t360 - t306 * t427) * t361) * MDP(20) + (-t357 * t364 + t347) * MDP(9); qJDD(1) * MDP(4) - t365 * MDP(5) + (-t462 - t470) * MDP(6) + (t436 * t357 + t347) * MDP(12) + (t436 * t361 - t414) * MDP(13) + (-t361 * t263 + (t433 - t456) * t357 + (-t356 * t429 - t393 * t360) * t341) * MDP(19) + (-t262 * t361 + (t432 - t455) * t357 + (t393 * t356 - t361 * t419) * t341) * MDP(20) + (qJD(1) * t380 + (-t320 * t338 * qJD(3) - t245) * t361 + ((t355 * t427 + t469 * t359 + t404) * t338 - t457 + qJD(3) * t269) * t357) * MDP(26) + (t338 * t376 + (qJD(3) * t380 - t244) * t361 + (-(t469 * t355 - t356 * t422 - t359 * t427) * t338 + t458 - qJD(3) * t386) * t357) * MDP(27); MDP(9) * t415 - MDP(10) * t416 + qJDD(3) * MDP(11) + (-t372 * t361 + t465) * MDP(12) + (t372 * t357 + t464) * MDP(13) + (t341 * t451 + t459) * MDP(14) + ((t262 - t454) * t360 + (-t263 - t453) * t356) * MDP(15) + ((-t317 * t361 + t341 * t446) * qJD(1) + t383) * MDP(16) + ((-t341 * t356 * t357 + t315 * t361) * qJD(1) + t382) * MDP(17) + (-t315 * t324 - pkin(3) * t263 - t304 * t341 + (t341 * t450 + t371) * t356 - t478 * t360) * MDP(19) + (-pkin(3) * t262 - t317 * t324 + t439 * t341 + t356 * t478 + t371 * t360) * MDP(20) + (t244 * t320 - t386 * t441) * MDP(21) + (-t244 * t319 - t245 * t320 - t441 * t269 + t386 * t440) * MDP(22) + (t441 * t338 + t457) * MDP(23) + (-t440 * t338 - t458) * MDP(24) + ((-t331 * t359 - t332 * t355) * t308 + t344 * t245 + t248 * t319 + (t387 * t355 - t388 * t359) * t338 + t440 * t277 + t392 * t269 - t370 * t350) * MDP(26) + (-(-t331 * t355 + t332 * t359) * t308 + t344 * t244 + t248 * t320 + (t388 * t355 + t387 * t359) * t338 + t441 * t277 - t392 * t386 + t370 * t349) * MDP(27) + (-t341 * MDP(18) - t260 * MDP(19) + t261 * MDP(20) + MDP(23) * t386 + t269 * MDP(24) - t338 * MDP(25) - t242 * MDP(26) + t243 * MDP(27)) * t434 + (t361 * t357 * MDP(7) - t437 * MDP(8)) * t365; t317 * t315 * MDP(14) + (-t315 ^ 2 + t317 ^ 2) * MDP(15) + (t262 + t454) * MDP(16) + (-t263 + t453) * MDP(17) + t312 * MDP(18) + (-t305 * t425 - g(1) * t299 - g(2) * t301 + t261 * t341 - t306 * t317 + t265 + (t395 + t464) * t356) * MDP(19) + (g(1) * t300 - g(2) * t302 + g(3) * t443 + t260 * t341 + t306 * t315 - t377) * MDP(20) + (t244 + t477) * MDP(23) + (-t245 - t476) * MDP(24) + (-(-t253 * t355 - t460) * t338 - t243 * qJD(5) + (-t269 * t317 + t359 * t308 - t338 * t423) * pkin(4) + t473) * MDP(26) + ((-t254 * t338 - t240) * t355 + (t253 * t338 - t394) * t359 + (-t355 * t308 + t317 * t386 - t338 * t422) * pkin(4) + t474) * MDP(27) + t472; (t411 + t477) * MDP(23) + (-t397 - t476) * MDP(24) + (t243 * t338 + t473) * MDP(26) + (-t355 * t240 - t359 * t241 + t242 * t338 + t474) * MDP(27) + (-MDP(23) * t452 + t386 * MDP(24) - t243 * MDP(26) - MDP(27) * t461) * qJD(5) + t472;];
tau = t1;
