% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:43:05
% EndTime: 2019-12-05 18:43:10
% DurationCPUTime: 2.73s
% Computational Cost: add. (2809->316), mult. (4325->399), div. (0->0), fcn. (3045->14), ass. (0->172)
t518 = qJ(4) + pkin(7);
t426 = sin(pkin(9));
t427 = cos(pkin(9));
t430 = sin(qJ(3));
t434 = cos(qJ(3));
t376 = -t426 * t430 + t427 * t434;
t422 = qJD(1) + qJD(2);
t356 = t376 * t422;
t433 = cos(qJ(5));
t343 = t433 * t356;
t377 = t426 * t434 + t427 * t430;
t357 = t377 * t422;
t429 = sin(qJ(5));
t499 = t357 * t429;
t306 = t343 - t499;
t421 = qJD(3) + qJD(5);
t500 = t306 * t421;
t425 = qJ(1) + qJ(2);
t415 = sin(t425);
t416 = cos(t425);
t517 = g(2) * t416 + g(3) * t415;
t513 = -g(2) * t415 + g(3) * t416;
t420 = qJDD(1) + qJDD(2);
t431 = sin(qJ(2));
t479 = qJDD(1) * t431;
t435 = cos(qJ(2));
t483 = qJD(2) * t435;
t361 = pkin(7) * t420 + (qJD(1) * t483 + t479) * pkin(1);
t447 = qJ(4) * t420 + qJD(4) * t422 + t361;
t511 = pkin(1) * t431;
t478 = qJD(1) * t511;
t469 = t518 * t422 + t478;
t451 = qJD(3) * t469;
t291 = qJDD(3) * pkin(3) - t430 * t447 - t434 * t451;
t294 = -t430 * t451 + t434 * t447;
t268 = t427 * t291 - t294 * t426;
t481 = qJD(3) * t434;
t471 = t422 * t481;
t482 = qJD(3) * t430;
t472 = t422 * t482;
t312 = t377 * t420 - t426 * t472 + t427 * t471;
t262 = qJDD(3) * pkin(4) - pkin(8) * t312 + t268;
t269 = t426 * t291 + t427 * t294;
t369 = t377 * qJD(3);
t311 = -t369 * t422 + t376 * t420;
t263 = pkin(8) * t311 + t269;
t345 = t469 * t430;
t339 = qJD(3) * pkin(3) - t345;
t346 = t469 * t434;
t496 = t427 * t346;
t298 = t426 * t339 + t496;
t506 = pkin(8) * t356;
t281 = t298 + t506;
t407 = pkin(3) * t434 + pkin(2);
t485 = qJD(1) * t435;
t477 = pkin(1) * t485;
t355 = -t407 * t422 + qJD(4) - t477;
t313 = -pkin(4) * t356 + t355;
t413 = qJ(3) + pkin(9) + qJ(5);
t398 = sin(t413);
t480 = qJD(5) * t429;
t399 = cos(t413);
t497 = t399 * t416;
t498 = t399 * t415;
t516 = g(1) * t398 - g(2) * t498 + g(3) * t497 - t429 * t262 - t433 * t263 + t281 * t480 - t313 * t306;
t419 = qJDD(3) + qJDD(5);
t454 = t356 * t429 + t433 * t357;
t515 = t419 * MDP(20) + (-t306 ^ 2 + t454 ^ 2) * MDP(17) - t306 * MDP(16) * t454;
t501 = t454 * t421;
t414 = t434 * qJD(4);
t470 = qJD(3) * t518;
t366 = -t430 * t470 + t414;
t367 = -qJD(4) * t430 - t434 * t470;
t492 = -t366 * t426 + t427 * t367 + t377 * t477;
t491 = t427 * t366 + t426 * t367 - t376 * t477;
t512 = -g(1) * t399 + t433 * t262 - t429 * t263 - t313 * t454 + t513 * t398;
t468 = -t433 * t311 + t312 * t429;
t267 = qJD(5) * t454 + t468;
t510 = pkin(1) * t435;
t509 = pkin(2) * t420;
t508 = pkin(2) * t422;
t507 = pkin(3) * t426;
t505 = pkin(8) * t357;
t370 = t376 * qJD(3);
t504 = pkin(8) * t370;
t503 = pkin(8) * t377;
t502 = g(1) * t434;
t333 = t426 * t346;
t495 = t430 * t434;
t297 = t427 * t339 - t333;
t279 = qJD(3) * pkin(4) + t297 - t505;
t494 = t433 * t279;
t406 = pkin(7) + t511;
t493 = -qJ(4) - t406;
t466 = qJD(3) * t493;
t476 = pkin(1) * t483;
t328 = t430 * t466 + t434 * t476 + t414;
t329 = (-qJD(4) - t476) * t430 + t434 * t466;
t293 = t427 * t328 + t426 * t329;
t300 = -t427 * t345 - t333;
t484 = qJD(2) * t431;
t412 = pkin(1) * t484;
t488 = -qJD(1) * t412 + qJDD(1) * t510;
t360 = -t488 - t509;
t388 = -t477 - t508;
t490 = t360 * t430 + t388 * t481;
t374 = t493 * t430;
t417 = t434 * qJ(4);
t375 = t406 * t434 + t417;
t320 = t426 * t374 + t427 * t375;
t390 = t518 * t430;
t391 = pkin(7) * t434 + t417;
t338 = -t426 * t390 + t427 * t391;
t423 = t430 ^ 2;
t487 = -t434 ^ 2 + t423;
t411 = pkin(3) * t482;
t475 = qJD(5) * t343 + t429 * t311 + t433 * t312;
t474 = t388 * t482 + t517 * t434;
t473 = t422 * t484;
t342 = pkin(4) * t369 + t411;
t292 = -t328 * t426 + t427 * t329;
t299 = t345 * t426 - t496;
t319 = t427 * t374 - t375 * t426;
t337 = -t427 * t390 - t391 * t426;
t467 = -t407 * t416 - t415 * t518;
t465 = t422 * t478;
t318 = pkin(3) * t472 - t407 * t420 + qJDD(4) - t488;
t280 = -pkin(4) * t311 + t318;
t322 = t376 * t429 + t377 * t433;
t285 = qJD(5) * t322 + t433 * t369 + t370 * t429;
t453 = t433 * t376 - t377 * t429;
t464 = g(2) * t497 + g(3) * t498 - t280 * t453 + t313 * t285;
t463 = -t488 - t517;
t462 = t342 - t478;
t373 = t376 * pkin(8);
t317 = t373 + t338;
t459 = qJD(5) * t317 - t492 + t504;
t316 = t337 - t503;
t365 = t369 * pkin(8);
t458 = -qJD(5) * t316 + t365 - t491;
t457 = -t429 * t279 - t433 * t281;
t301 = t319 - t503;
t302 = t373 + t320;
t456 = t301 * t433 - t302 * t429;
t455 = t301 * t429 + t302 * t433;
t452 = -t407 * t415 + t416 * t518;
t350 = -pkin(4) * t376 - t407;
t400 = pkin(3) * t427 + pkin(4);
t450 = t400 * t429 + t433 * t507;
t449 = t400 * t433 - t429 * t507;
t448 = -t268 * t377 + t269 * t376 - t297 * t370 - t298 * t369 - t513;
t266 = -t357 * t480 + t475;
t446 = -t388 * t422 - t361 + t513;
t284 = qJD(5) * t453 - t369 * t429 + t370 * t433;
t445 = t280 * t322 + t313 * t284 - t398 * t517;
t437 = qJD(3) ^ 2;
t444 = (t266 * t453 - t267 * t322 + t284 * t306 - t285 * t454) * MDP(17) + (t266 * t322 + t284 * t454) * MDP(16) + (t284 * t421 + t322 * t419) * MDP(18) + (-t285 * t421 + t419 * t453) * MDP(19) + 0.2e1 * (-qJD(3) * t422 * t487 + t420 * t495) * MDP(8) + (t420 * t423 + 0.2e1 * t430 * t471) * MDP(7) + (qJDD(3) * t434 - t430 * t437) * MDP(10) + (qJDD(3) * t430 + t434 * t437) * MDP(9) + t420 * MDP(4);
t443 = -pkin(7) * t437 + t465 + t509;
t408 = -pkin(2) - t510;
t442 = -pkin(1) * t473 - t406 * t437 - t408 * t420;
t440 = -pkin(7) * qJDD(3) + (t477 - t508) * qJD(3);
t439 = -qJDD(3) * t406 + (t408 * t422 - t476) * qJD(3);
t436 = cos(qJ(1));
t432 = sin(qJ(1));
t341 = t350 - t510;
t332 = t342 + t412;
t327 = pkin(3) * t422 * t430 + pkin(4) * t357;
t283 = t300 - t505;
t282 = t299 - t506;
t277 = -t365 + t293;
t276 = t292 - t504;
t1 = [(t332 * t454 + t341 * t266 - (qJD(5) * t456 + t276 * t429 + t277 * t433) * t421 - t455 * t419 + t445) * MDP(22) + qJDD(1) * MDP(1) + (t439 * t434 + (-t442 - t517) * t430 + t490) * MDP(13) + (t439 * t430 + (-t360 + t442) * t434 + t474) * MDP(12) + (-t332 * t306 + t341 * t267 + (-qJD(5) * t455 + t276 * t433 - t277 * t429) * t421 + t456 * t419 + t464) * MDP(21) + (((-qJDD(1) - t420) * t431 + (-qJD(1) - t422) * t483) * pkin(1) + t513) * MDP(6) + (-t292 * t357 + t293 * t356 + t311 * t320 - t312 * t319 + t448) * MDP(14) + t444 + (g(2) * t436 + g(3) * t432) * MDP(2) + (-g(2) * t432 + g(3) * t436) * MDP(3) + (t269 * t320 + t298 * t293 + t268 * t319 + t297 * t292 + t318 * (-t407 - t510) + t355 * (t412 + t411) - g(2) * (-pkin(1) * t436 + t467) - g(3) * (-pkin(1) * t432 + t452)) * MDP(15) + ((t420 * t435 - t473) * pkin(1) - t463) * MDP(5); (t269 * t338 + t268 * t337 - t318 * t407 - g(2) * t467 - g(3) * t452 + (t411 - t478) * t355 + t491 * t298 + t492 * t297) * MDP(15) + ((-t479 + (-qJD(2) + t422) * t485) * pkin(1) + t513) * MDP(6) + (t350 * t267 + (t316 * t433 - t317 * t429) * t419 + (t429 * t458 - t433 * t459) * t421 - t462 * t306 + t464) * MDP(21) + (t440 * t434 + (-t443 - t517) * t430 + t490) * MDP(13) + (t440 * t430 + (-t360 + t443) * t434 + t474) * MDP(12) + (t350 * t266 - (t316 * t429 + t317 * t433) * t419 + (t429 * t459 + t433 * t458) * t421 + t462 * t454 + t445) * MDP(22) + (t311 * t338 - t312 * t337 + t356 * t491 - t357 * t492 + t448) * MDP(14) + t444 + (-t463 + t465) * MDP(5); qJDD(3) * MDP(11) + (t430 * t446 - t502) * MDP(12) + (g(1) * t430 + t434 * t446) * MDP(13) + ((t298 + t299) * t357 + (t297 - t300) * t356 + (t311 * t426 - t312 * t427) * pkin(3)) * MDP(14) + (-t297 * t299 - t298 * t300 + (-t502 + t268 * t427 + t269 * t426 + (-t355 * t422 + t513) * t430) * pkin(3)) * MDP(15) + (t266 - t500) * MDP(18) + (-t267 + t501) * MDP(19) + (t449 * t419 + t327 * t306 - (t282 * t433 - t283 * t429) * t421 + (-t421 * t450 + t457) * qJD(5) + t512) * MDP(21) + (-t450 * t419 - t327 * t454 + (t282 * t429 + t283 * t433) * t421 + (-t421 * t449 - t494) * qJD(5) + t516) * MDP(22) + (MDP(10) * t434 + MDP(9) * t430) * t420 + (-MDP(7) * t495 + MDP(8) * t487) * t422 ^ 2 + t515; (-t356 ^ 2 - t357 ^ 2) * MDP(14) + (t297 * t357 - t298 * t356 + t318 - t517) * MDP(15) + (t267 + t501) * MDP(21) + (t266 + t500) * MDP(22); (t475 - t500) * MDP(18) + (-t468 + t501) * MDP(19) + (-t421 * t457 + t512) * MDP(21) + ((-t281 * t429 + t494) * t421 + t516) * MDP(22) + (-MDP(18) * t499 - t454 * MDP(19) + t457 * MDP(21) - MDP(22) * t494) * qJD(5) + t515;];
tau = t1;
