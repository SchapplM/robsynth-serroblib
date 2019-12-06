% Calculate vector of inverse dynamics joint torques for
% S5RRPRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:36
% EndTime: 2019-12-05 18:28:42
% DurationCPUTime: 4.39s
% Computational Cost: add. (3777->335), mult. (9156->440), div. (0->0), fcn. (7183->14), ass. (0->160)
t428 = sin(pkin(9));
t429 = cos(pkin(9));
t433 = sin(qJ(2));
t437 = cos(qJ(2));
t387 = -t428 * t433 + t429 * t437;
t377 = t387 * qJD(1);
t388 = t428 * t437 + t429 * t433;
t379 = t388 * qJD(1);
t432 = sin(qJ(4));
t436 = cos(qJ(4));
t338 = t436 * t377 - t379 * t432;
t431 = sin(qJ(5));
t435 = cos(qJ(5));
t454 = t377 * t432 + t436 * t379;
t290 = t338 * t431 + t435 * t454;
t494 = t454 * t431;
t294 = t338 * t435 - t494;
t424 = qJDD(2) + qJDD(4);
t421 = qJDD(5) + t424;
t526 = t421 * MDP(24) + (t290 ^ 2 - t294 ^ 2) * MDP(21) - t294 * t290 * MDP(20);
t378 = t388 * qJD(2);
t348 = -qJD(1) * t378 + t387 * qJDD(1);
t478 = qJD(1) * qJD(2);
t473 = t437 * t478;
t474 = t433 * t478;
t349 = t388 * qJDD(1) - t428 * t474 + t429 * t473;
t483 = qJD(4) * t436;
t484 = qJD(4) * t432;
t284 = t432 * t348 + t436 * t349 + t377 * t483 - t379 * t484;
t285 = t454 * qJD(4) - t436 * t348 + t349 * t432;
t481 = qJD(5) * t435;
t475 = t435 * t284 - t431 * t285 + t338 * t481;
t482 = qJD(5) * t431;
t265 = -t454 * t482 + t475;
t471 = t284 * t431 + t435 * t285;
t443 = -qJD(5) * t290 - t471;
t425 = qJD(2) + qJD(4);
t495 = t338 * t425;
t496 = t454 * t425;
t422 = qJD(5) + t425;
t523 = t294 * t422;
t524 = t290 * t422;
t525 = t424 * MDP(17) + (-t285 + t496) * MDP(16) + (-t338 ^ 2 + t454 ^ 2) * MDP(14) - t338 * t454 * MDP(13) + (t284 - t495) * MDP(15) + (t265 - t523) * MDP(22) + (t443 + t524) * MDP(23) + t526;
t500 = qJ(3) + pkin(6);
t472 = qJD(2) * t500;
t375 = -qJD(3) * t433 - t437 * t472;
t403 = t500 * t433;
t345 = qJDD(2) * pkin(2) + t375 * qJD(1) - qJDD(1) * t403;
t374 = qJD(3) * t437 - t433 * t472;
t404 = t500 * t437;
t352 = t374 * qJD(1) + qJDD(1) * t404;
t297 = t429 * t345 - t352 * t428;
t286 = qJDD(2) * pkin(3) - pkin(7) * t349 + t297;
t298 = t428 * t345 + t429 * t352;
t287 = pkin(7) * t348 + t298;
t393 = qJD(1) * t404;
t382 = t428 * t393;
t392 = qJD(1) * t403;
t499 = qJD(2) * pkin(2);
t386 = -t392 + t499;
t346 = t429 * t386 - t382;
t504 = pkin(7) * t379;
t312 = qJD(2) * pkin(3) + t346 - t504;
t493 = t429 * t393;
t347 = t428 * t386 + t493;
t505 = pkin(7) * t377;
t318 = t347 + t505;
t456 = -t432 * t312 - t436 * t318;
t446 = t456 * qJD(4) + t436 * t286 - t432 * t287;
t263 = pkin(4) * t424 - pkin(8) * t284 + t446;
t507 = (qJD(4) * t312 + t287) * t436 + t432 * t286 - t318 * t484;
t264 = -pkin(8) * t285 + t507;
t418 = pkin(2) * t437 + pkin(1);
t398 = -t418 * qJD(1) + qJD(3);
t355 = -pkin(3) * t377 + t398;
t303 = -pkin(4) * t338 + t355;
t423 = qJ(2) + pkin(9) + qJ(4);
t416 = qJ(5) + t423;
t410 = sin(t416);
t411 = cos(t416);
t434 = sin(qJ(1));
t438 = cos(qJ(1));
t461 = g(1) * t438 + g(2) * t434;
t517 = -g(3) * t411 + t435 * t263 - t431 * t264 - t303 * t290 + t461 * t410;
t521 = pkin(8) * t338;
t275 = -t456 + t521;
t518 = g(3) * t410 + t275 * t482 - t303 * t294 + t461 * t411;
t413 = sin(t423);
t414 = cos(t423);
t516 = g(3) * t413 - t355 * t338 + t461 * t414 - t507;
t512 = pkin(8) * t454;
t356 = -t429 * t403 - t404 * t428;
t331 = -pkin(7) * t388 + t356;
t357 = -t428 * t403 + t429 * t404;
t332 = pkin(7) * t387 + t357;
t487 = t432 * t331 + t436 * t332;
t353 = t392 * t428 - t493;
t319 = t353 - t505;
t354 = -t429 * t392 - t382;
t320 = t354 - t504;
t415 = pkin(2) * t429 + pkin(3);
t506 = pkin(2) * t428;
t462 = t436 * t415 - t432 * t506;
t510 = -t462 * qJD(4) + t432 * t319 + t436 * t320;
t373 = t415 * t432 + t436 * t506;
t509 = -t373 * qJD(4) - t436 * t319 + t320 * t432;
t508 = -g(3) * t414 - t355 * t454 + t461 * t413 + t446;
t501 = g(3) * t437;
t469 = t436 * t312 - t318 * t432;
t274 = t469 - t512;
t272 = pkin(4) * t425 + t274;
t492 = t435 * t272;
t491 = t435 * t275;
t490 = t509 + t521;
t489 = t510 - t512;
t330 = t429 * t374 + t428 * t375;
t426 = t433 ^ 2;
t485 = -t437 ^ 2 + t426;
t477 = qJDD(1) * t437;
t420 = t433 * t499;
t359 = pkin(3) * t378 + t420;
t358 = t433 * qJD(1) * pkin(2) + pkin(3) * t379;
t467 = t436 * t331 - t332 * t432;
t329 = -t374 * t428 + t429 * t375;
t464 = -qJD(5) * t272 - t264;
t460 = g(1) * t434 - g(2) * t438;
t459 = -t431 * t272 - t491;
t351 = t387 * t432 + t388 * t436;
t278 = -pkin(8) * t351 + t467;
t453 = t436 * t387 - t388 * t432;
t279 = pkin(8) * t453 + t487;
t458 = t278 * t435 - t279 * t431;
t457 = t278 * t431 + t279 * t435;
t301 = t351 * t431 - t435 * t453;
t302 = t351 * t435 + t431 * t453;
t361 = -pkin(3) * t387 - t418;
t452 = -0.2e1 * pkin(1) * t478 - pkin(6) * qJDD(2);
t381 = t387 * qJD(2);
t308 = -pkin(7) * t381 + t329;
t309 = -pkin(7) * t378 + t330;
t451 = t432 * t308 + t436 * t309 + t331 * t483 - t332 * t484;
t450 = pkin(2) * t474 - t418 * qJDD(1) + qJDD(3);
t439 = qJD(2) ^ 2;
t448 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t439 + t460;
t440 = qJD(1) ^ 2;
t447 = pkin(1) * t440 - pkin(6) * qJDD(1) + t461;
t314 = -pkin(3) * t348 + t450;
t444 = -t487 * qJD(4) + t436 * t308 - t309 * t432;
t372 = pkin(4) + t462;
t313 = -pkin(4) * t453 + t361;
t304 = pkin(4) * t454 + t358;
t300 = t351 * qJD(4) + t436 * t378 + t381 * t432;
t299 = t453 * qJD(4) - t378 * t432 + t381 * t436;
t288 = pkin(4) * t300 + t359;
t271 = pkin(4) * t285 + t314;
t270 = t302 * qJD(5) + t299 * t431 + t435 * t300;
t269 = -t301 * qJD(5) + t299 * t435 - t300 * t431;
t268 = -pkin(8) * t299 + t444;
t267 = -pkin(8) * t300 + t451;
t1 = [qJDD(1) * MDP(1) + t460 * MDP(2) + t461 * MDP(3) + (qJDD(1) * t426 + 0.2e1 * t433 * t473) * MDP(4) + 0.2e1 * (t433 * t477 - t485 * t478) * MDP(5) + (qJDD(2) * t433 + t437 * t439) * MDP(6) + (qJDD(2) * t437 - t433 * t439) * MDP(7) + (t452 * t433 + t448 * t437) * MDP(9) + (-t448 * t433 + t452 * t437) * MDP(10) + (-t297 * t388 + t298 * t387 - t329 * t379 + t330 * t377 - t346 * t381 - t347 * t378 + t348 * t357 - t349 * t356 - t461) * MDP(11) + (t298 * t357 + t347 * t330 + t297 * t356 + t346 * t329 - t450 * t418 + t398 * t420 - g(1) * (-t418 * t434 + t438 * t500) - g(2) * (t418 * t438 + t434 * t500)) * MDP(12) + (t284 * t351 + t299 * t454) * MDP(13) + (t284 * t453 - t285 * t351 + t299 * t338 - t300 * t454) * MDP(14) + (t299 * t425 + t351 * t424) * MDP(15) + (-t300 * t425 + t424 * t453) * MDP(16) + (t361 * t285 + t355 * t300 - t314 * t453 - t338 * t359 + t460 * t414 + t467 * t424 + t444 * t425) * MDP(18) + (t361 * t284 + t355 * t299 + t314 * t351 + t359 * t454 - t460 * t413 - t487 * t424 - t451 * t425) * MDP(19) + (t265 * t302 + t269 * t290) * MDP(20) + (-t265 * t301 + t269 * t294 - t270 * t290 + t302 * t443) * MDP(21) + (t269 * t422 + t302 * t421) * MDP(22) + (-t270 * t422 - t301 * t421) * MDP(23) + (-t288 * t294 - t313 * t443 + t271 * t301 + t303 * t270 + (-t457 * qJD(5) - t267 * t431 + t268 * t435) * t422 + t458 * t421 + t460 * t411) * MDP(25) + (t288 * t290 + t313 * t265 + t271 * t302 + t303 * t269 - (t458 * qJD(5) + t267 * t435 + t268 * t431) * t422 - t457 * t421 - t460 * t410) * MDP(26); (g(3) * t433 + t447 * t437) * MDP(10) + ((t347 + t353) * t379 + (t346 - t354) * t377 + (t348 * t428 - t349 * t429) * pkin(2)) * MDP(11) + (-t358 * t454 - t373 * t424 + t510 * t425 + t516) * MDP(19) + qJDD(2) * MDP(8) + (-t304 * t290 + (-t372 * t421 - t263 + (qJD(5) * t373 - t490) * t422) * t431 + (-t373 * t421 + (-qJD(5) * t372 + t489) * t422 + t464) * t435 + t518) * MDP(26) + (t338 * t358 + t462 * t424 + t509 * t425 + t508) * MDP(18) + (-t433 * t437 * MDP(4) + t485 * MDP(5)) * t440 + t433 * qJDD(1) * MDP(6) + ((t372 * t435 - t373 * t431) * t421 + t304 * t294 + (t489 * t431 + t490 * t435) * t422 + ((-t372 * t431 - t373 * t435) * t422 + t459) * qJD(5) + t517) * MDP(25) + (t447 * t433 - t501) * MDP(9) + (-t346 * t353 - t347 * t354 + (-t501 + t297 * t429 + t298 * t428 + (-qJD(1) * t398 + t461) * t433) * pkin(2)) * MDP(12) + MDP(7) * t477 + t525; (-t377 ^ 2 - t379 ^ 2) * MDP(11) + (t346 * t379 - t347 * t377 + t450 - t460) * MDP(12) + (t285 + t496) * MDP(18) + (t284 + t495) * MDP(19) + (-t443 + t524) * MDP(25) + (t265 + t523) * MDP(26); (-t456 * t425 + t508) * MDP(18) + (t469 * t425 + t516) * MDP(19) + (-(-t274 * t431 - t491) * t422 + t459 * qJD(5) + (t294 * t454 + t435 * t421 - t422 * t482) * pkin(4) + t517) * MDP(25) + ((-t275 * t422 - t263) * t431 + (t274 * t422 + t464) * t435 + (-t290 * t454 - t431 * t421 - t422 * t481) * pkin(4) + t518) * MDP(26) + t525; (t475 - t523) * MDP(22) + (-t471 + t524) * MDP(23) + (-t459 * t422 + t517) * MDP(25) + (-t435 * t264 - t431 * t263 + (-t275 * t431 + t492) * t422 + t518) * MDP(26) + (-MDP(22) * t494 - t290 * MDP(23) + t459 * MDP(25) - MDP(26) * t492) * qJD(5) + t526;];
tau = t1;
