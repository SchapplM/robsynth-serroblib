% Calculate vector of inverse dynamics joint torques for
% S5RRRPR1
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:39:01
% EndTime: 2019-12-05 18:39:10
% DurationCPUTime: 5.28s
% Computational Cost: add. (4230->359), mult. (10328->470), div. (0->0), fcn. (7724->14), ass. (0->183)
t462 = sin(qJ(3));
t466 = cos(qJ(3));
t467 = cos(qJ(2));
t520 = qJD(1) * t467;
t463 = sin(qJ(2));
t521 = qJD(1) * t463;
t393 = -t462 * t520 - t466 * t521;
t388 = t393 * qJ(4);
t550 = pkin(6) + pkin(7);
t420 = t550 * t467;
t410 = qJD(1) * t420;
t394 = t462 * t410;
t419 = t550 * t463;
t408 = qJD(1) * t419;
t541 = qJD(2) * pkin(2);
t400 = -t408 + t541;
t503 = t466 * t400 - t394;
t353 = t388 + t503;
t455 = qJD(2) + qJD(3);
t344 = pkin(3) * t455 + t353;
t459 = sin(pkin(9));
t398 = t466 * t410;
t488 = -t400 * t462 - t398;
t509 = t466 * t520;
t510 = t462 * t521;
t391 = -t509 + t510;
t540 = qJ(4) * t391;
t354 = -t488 - t540;
t460 = cos(pkin(9));
t534 = t460 * t354;
t309 = t459 * t344 + t534;
t489 = -t391 * t460 + t393 * t459;
t547 = pkin(8) * t489;
t301 = t309 + t547;
t465 = cos(qJ(5));
t531 = t465 * t489;
t369 = t391 * t459 + t393 * t460;
t461 = sin(qJ(5));
t536 = t369 * t461;
t330 = t531 + t536;
t452 = t467 * pkin(2);
t543 = pkin(1) + t452;
t418 = t543 * qJD(1);
t380 = pkin(3) * t391 + qJD(4) - t418;
t339 = -pkin(4) * t489 + t380;
t458 = qJ(2) + qJ(3);
t441 = pkin(9) + qJ(5) + t458;
t432 = sin(t441);
t433 = cos(t441);
t464 = sin(qJ(1));
t468 = cos(qJ(1));
t497 = g(1) * t468 + g(2) * t464;
t517 = qJD(5) * t461;
t560 = g(3) * t432 + t301 * t517 - t339 * t330 + t497 * t433;
t514 = qJDD(1) * t467;
t516 = qJD(1) * qJD(2);
t507 = t467 * t516;
t515 = qJDD(1) * t463;
t554 = t507 + t515;
t355 = qJD(3) * t509 - t455 * t510 + t462 * t514 + t554 * t466;
t453 = qJDD(2) + qJDD(3);
t379 = qJDD(2) * pkin(2) - t550 * t554;
t508 = t463 * t516;
t381 = t550 * (-t508 + t514);
t475 = qJD(3) * t488 + t466 * t379 - t462 * t381;
t296 = pkin(3) * t453 - qJ(4) * t355 + qJD(4) * t393 + t475;
t403 = t462 * t467 + t463 * t466;
t378 = t455 * t403;
t495 = t462 * t515 - t466 * t514;
t356 = qJD(1) * t378 + t495;
t519 = qJD(3) * t462;
t551 = (qJD(3) * t400 + t381) * t466 + t462 * t379 - t410 * t519;
t298 = -qJ(4) * t356 - qJD(4) * t391 + t551;
t282 = t460 * t296 - t298 * t459;
t316 = t355 * t460 - t356 * t459;
t278 = pkin(4) * t453 - pkin(8) * t316 + t282;
t283 = t459 * t296 + t460 * t298;
t315 = -t355 * t459 - t356 * t460;
t279 = pkin(8) * t315 + t283;
t559 = -t461 * t278 - t465 * t279 + t560;
t449 = qJD(5) + t455;
t538 = t330 * t449;
t448 = qJDD(5) + t453;
t553 = -t465 * t369 + t461 * t489;
t558 = t448 * MDP(24) + (-t330 ^ 2 + t553 ^ 2) * MDP(21) - t330 * MDP(20) * t553;
t539 = t553 * t449;
t450 = sin(t458);
t451 = cos(t458);
t556 = -g(3) * t451 + t497 * t450;
t476 = -g(3) * t433 + t465 * t278 - t461 * t279 - t339 * t553 + t497 * t432;
t366 = t369 * pkin(8);
t502 = t408 * t462 - t398;
t358 = t502 + t540;
t526 = -t466 * t408 - t394;
t359 = t388 + t526;
t533 = t460 * t462;
t542 = pkin(2) * qJD(3);
t528 = -t460 * t358 + t359 * t459 + (-t459 * t466 - t533) * t542;
t535 = t459 * t462;
t527 = -t459 * t358 - t460 * t359 + (t460 * t466 - t535) * t542;
t525 = -t462 * t419 + t466 * t420;
t505 = -t465 * t315 + t316 * t461;
t286 = qJD(5) * t553 + t505;
t549 = pkin(3) * t393;
t548 = pkin(3) * t459;
t345 = t459 * t354;
t308 = t460 * t344 - t345;
t299 = pkin(4) * t455 + t308 + t366;
t532 = t465 * t299;
t530 = t547 + t528;
t529 = t366 - t527;
t402 = t462 * t463 - t466 * t467;
t511 = qJD(2) * t550;
t409 = t463 * t511;
t411 = t467 * t511;
t518 = qJD(3) * t466;
t482 = -t466 * t409 - t462 * t411 - t419 * t518 - t420 * t519;
t321 = -qJ(4) * t378 - qJD(4) * t402 + t482;
t377 = t455 * t402;
t474 = -qJD(3) * t525 + t409 * t462 - t466 * t411;
t322 = qJ(4) * t377 - qJD(4) * t403 + t474;
t292 = t460 * t321 + t459 * t322;
t314 = t460 * t353 - t345;
t501 = -t466 * t419 - t420 * t462;
t367 = -qJ(4) * t403 + t501;
t368 = -qJ(4) * t402 + t525;
t326 = t459 * t367 + t460 * t368;
t524 = pkin(3) * t451 + t452;
t456 = t463 ^ 2;
t523 = -t467 ^ 2 + t456;
t447 = t463 * t541;
t513 = qJD(5) * t531 + t461 * t315 + t465 * t316;
t506 = pkin(3) * t378 + t447;
t291 = -t321 * t459 + t460 * t322;
t313 = -t353 * t459 - t534;
t325 = t460 * t367 - t368 * t459;
t498 = pkin(3) * t402 - t543;
t444 = t466 * pkin(2) + pkin(3);
t386 = -pkin(2) * t535 + t460 * t444;
t343 = -pkin(4) * t369 - t549;
t496 = g(1) * t464 - g(2) * t468;
t494 = -t461 * t299 - t465 * t301;
t375 = -t402 * t459 + t403 * t460;
t306 = -pkin(8) * t375 + t325;
t374 = -t402 * t460 - t403 * t459;
t307 = pkin(8) * t374 + t326;
t493 = t306 * t465 - t307 * t461;
t492 = t306 * t461 + t307 * t465;
t491 = t308 * t489 - t309 * t369;
t490 = t465 * t374 - t375 * t461;
t333 = t374 * t461 + t375 * t465;
t438 = pkin(3) * t460 + pkin(4);
t487 = t438 * t461 + t465 * t548;
t486 = t438 * t465 - t461 * t548;
t485 = -0.2e1 * pkin(1) * t516 - pkin(6) * qJDD(2);
t389 = pkin(2) * t508 - qJDD(1) * t543;
t285 = t369 * t517 + t513;
t469 = qJD(2) ^ 2;
t479 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t469 + t496;
t470 = qJD(1) ^ 2;
t478 = pkin(1) * t470 - pkin(6) * qJDD(1) + t497;
t477 = t356 * pkin(3) + qJDD(4) + t389;
t473 = g(3) * t450 - t418 * t391 + t497 * t451 - t551;
t472 = -t393 * t391 * MDP(11) + (t285 - t538) * MDP(22) + (-t286 + t539) * MDP(23) + (t391 * t455 + t355) * MDP(13) + (-t495 + (-qJD(1) * t403 - t393) * t455) * MDP(14) + (-t391 ^ 2 + t393 ^ 2) * MDP(12) + t453 * MDP(15) + t558;
t471 = -t418 * t393 + t475 + t556;
t454 = -qJ(4) - t550;
t446 = pkin(2) * t521;
t407 = pkin(1) + t524;
t387 = pkin(2) * t533 + t444 * t459;
t382 = pkin(4) + t386;
t349 = -pkin(4) * t374 + t498;
t340 = t343 + t446;
t338 = -t377 * t460 - t378 * t459;
t337 = t377 * t459 - t378 * t460;
t320 = -pkin(4) * t337 + t506;
t303 = t366 + t314;
t302 = t313 - t547;
t293 = -pkin(4) * t315 + t477;
t290 = qJD(5) * t333 - t465 * t337 + t338 * t461;
t289 = qJD(5) * t490 + t337 * t461 + t338 * t465;
t288 = pkin(8) * t337 + t292;
t287 = -pkin(8) * t338 + t291;
t1 = [qJDD(1) * MDP(1) + t496 * MDP(2) + t497 * MDP(3) + (qJDD(1) * t456 + 0.2e1 * t463 * t507) * MDP(4) + 0.2e1 * (t463 * t514 - t523 * t516) * MDP(5) + (qJDD(2) * t463 + t467 * t469) * MDP(6) + (qJDD(2) * t467 - t463 * t469) * MDP(7) + (t463 * t485 + t467 * t479) * MDP(9) + (-t463 * t479 + t467 * t485) * MDP(10) + (t355 * t403 + t377 * t393) * MDP(11) + (-t355 * t402 - t356 * t403 + t377 * t391 + t378 * t393) * MDP(12) + (-t377 * t455 + t403 * t453) * MDP(13) + (-t378 * t455 - t402 * t453) * MDP(14) + (-t356 * t543 - t418 * t378 + t389 * t402 + t391 * t447 + t496 * t451 + t501 * t453 + t474 * t455) * MDP(16) + (-t355 * t543 + t418 * t377 + t389 * t403 - t393 * t447 - t496 * t450 - t525 * t453 - t482 * t455) * MDP(17) + (-t282 * t375 + t283 * t374 + t291 * t369 + t292 * t489 - t308 * t338 + t309 * t337 + t315 * t326 - t316 * t325 - t497) * MDP(18) + (t283 * t326 + t309 * t292 + t282 * t325 + t308 * t291 + t477 * t498 + t380 * t506 - g(1) * (-t407 * t464 - t454 * t468) - g(2) * (t407 * t468 - t454 * t464)) * MDP(19) + (t285 * t333 + t289 * t553) * MDP(20) + (t285 * t490 - t286 * t333 + t289 * t330 - t290 * t553) * MDP(21) + (t289 * t449 + t333 * t448) * MDP(22) + (-t290 * t449 + t448 * t490) * MDP(23) + (-t320 * t330 + t349 * t286 - t293 * t490 + t339 * t290 + (-qJD(5) * t492 + t287 * t465 - t288 * t461) * t449 + t493 * t448 + t496 * t433) * MDP(25) + (t320 * t553 + t349 * t285 + t293 * t333 + t339 * t289 - (qJD(5) * t493 + t287 * t461 + t288 * t465) * t449 - t492 * t448 - t496 * t432) * MDP(26); (t315 * t387 - t316 * t386 + t369 * t528 + t489 * t527 + t491) * MDP(18) + ((t382 * t465 - t387 * t461) * t448 + t340 * t330 + (t529 * t461 + t530 * t465) * t449 + ((-t382 * t461 - t387 * t465) * t449 + t494) * qJD(5) + t476) * MDP(25) + qJDD(2) * MDP(8) + (-g(3) * t467 + t463 * t478) * MDP(9) + (t283 * t387 + t282 * t386 - t380 * (t446 - t549) - g(3) * t524 - t497 * (-pkin(2) * t463 - pkin(3) * t450) + t527 * t309 + t528 * t308) * MDP(19) + (t526 * t455 + (t393 * t521 - t462 * t453 - t455 * t518) * pkin(2) + t473) * MDP(17) + (-t502 * t455 + (-t391 * t521 + t466 * t453 - t455 * t519) * pkin(2) + t471) * MDP(16) + (-t340 * t553 + (-t382 * t448 - t278 + (qJD(5) * t387 - t530) * t449) * t461 + (-qJD(5) * t299 - t387 * t448 - t279 + (-qJD(5) * t382 + t529) * t449) * t465 + t560) * MDP(26) + (g(3) * t463 + t467 * t478) * MDP(10) + t472 + MDP(6) * t515 + MDP(7) * t514 + (-t463 * t467 * MDP(4) + t523 * MDP(5)) * t470; (t486 * t448 + t343 * t330 - (t302 * t465 - t303 * t461) * t449 + (-t487 * t449 + t494) * qJD(5) + t476) * MDP(25) + (-t313 * t369 - t314 * t489 + (t315 * t459 - t316 * t460) * pkin(3) + t491) * MDP(18) + (-t455 * t488 + t471) * MDP(16) + (t503 * t455 + t473) * MDP(17) + (-t487 * t448 - t343 * t553 + (t302 * t461 + t303 * t465) * t449 + (-t449 * t486 - t532) * qJD(5) + t559) * MDP(26) + t472 + (-t308 * t313 - t309 * t314 + (t282 * t460 + t283 * t459 + t380 * t393 + t556) * pkin(3)) * MDP(19); (-t369 ^ 2 - t489 ^ 2) * MDP(18) + (-t308 * t369 - t309 * t489 + t477 - t496) * MDP(19) + (t286 + t539) * MDP(25) + (t285 + t538) * MDP(26); (t513 - t538) * MDP(22) + (-t505 + t539) * MDP(23) + (-t494 * t449 + t476) * MDP(25) + ((-t301 * t461 + t532) * t449 + t559) * MDP(26) + (MDP(22) * t536 - MDP(23) * t553 + t494 * MDP(25) - MDP(26) * t532) * qJD(5) + t558;];
tau = t1;
