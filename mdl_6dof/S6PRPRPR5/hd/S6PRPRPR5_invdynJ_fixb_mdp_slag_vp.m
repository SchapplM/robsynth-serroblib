% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:15
% EndTime: 2019-03-08 19:45:21
% DurationCPUTime: 5.70s
% Computational Cost: add. (2793->453), mult. (6516->572), div. (0->0), fcn. (5375->14), ass. (0->202)
t425 = cos(pkin(11));
t549 = cos(qJ(4));
t492 = t549 * t425;
t474 = qJD(2) * t492;
t428 = sin(qJ(4));
t422 = sin(pkin(11));
t514 = qJD(2) * t422;
t489 = t428 * t514;
t379 = -t474 + t489;
t430 = cos(qJ(6));
t427 = sin(qJ(6));
t512 = qJD(4) * t427;
t358 = -t430 * t379 + t512;
t388 = t422 * t549 + t428 * t425;
t559 = t388 * qJD(2);
t564 = qJD(6) + t559;
t568 = t358 * t564;
t360 = qJD(4) * t430 + t379 * t427;
t477 = t564 * t360;
t497 = MDP(14) - MDP(17);
t511 = qJD(4) * t428;
t488 = t422 * t511;
t484 = qJDD(2) * t549;
t501 = qJDD(2) * t428;
t493 = qJD(4) * t474 + t422 * t484 + t425 * t501;
t340 = qJD(2) * t488 - t493;
t334 = -qJDD(6) + t340;
t328 = t430 * t334;
t479 = t427 * t564;
t567 = -t479 * t564 - t328;
t530 = t422 * MDP(6);
t566 = t425 * MDP(5) - t530;
t543 = pkin(8) + qJ(3);
t395 = t543 * t422;
t396 = t543 * t425;
t355 = -t428 * t395 + t396 * t549;
t421 = pkin(11) + qJ(4);
t416 = sin(t421);
t423 = sin(pkin(10));
t429 = sin(qJ(2));
t431 = cos(qJ(2));
t541 = cos(pkin(10));
t542 = cos(pkin(6));
t469 = t542 * t541;
t375 = t423 * t429 - t431 * t469;
t483 = t423 * t542;
t377 = t429 * t541 + t431 * t483;
t473 = g(1) * t377 + g(2) * t375;
t424 = sin(pkin(6));
t526 = t424 * t431;
t560 = -g(3) * t526 + t473;
t444 = t560 * t416;
t466 = t492 * t526;
t487 = qJD(4) * t549;
t490 = qJD(1) * t526;
t529 = t422 * t428;
t520 = qJD(1) * t466 + t428 * (qJD(3) * t422 + qJD(4) * t396) - t490 * t529 - qJD(3) * t492 + t395 * t487;
t565 = -qJD(4) * t520 + qJDD(4) * t355 + t444;
t515 = qJD(1) * t429;
t491 = t424 * t515;
t392 = qJD(2) * qJ(3) + t491;
t481 = qJD(1) * t542;
t406 = t425 * t481;
t350 = t406 + (-pkin(8) * qJD(2) - t392) * t422;
t363 = t425 * t392 + t422 * t481;
t513 = qJD(2) * t425;
t351 = pkin(8) * t513 + t363;
t518 = t549 * t350 - t428 * t351;
t558 = qJD(5) - t518;
t516 = t422 ^ 2 + t425 ^ 2;
t562 = MDP(7) * t516;
t527 = t424 * t429;
t445 = t388 * t526;
t519 = -qJD(1) * t445 + qJD(3) * t388 + qJD(4) * t355;
t468 = qJD(3) - t490;
t496 = MDP(15) - MDP(18);
t486 = qJD(2) * t515;
t557 = t424 * t486 + qJDD(3);
t310 = t428 * t350 + t549 * t351;
t303 = -qJD(4) * qJ(5) - t310;
t547 = pkin(5) * t379;
t298 = -t303 - t547;
t550 = pkin(4) + pkin(9);
t556 = t550 * t334 + (t298 - t310 + t547) * t564;
t354 = t395 * t549 + t428 * t396;
t417 = cos(t421);
t553 = qJD(4) * t519 + qJDD(4) * t354 - t417 * t560;
t552 = t379 ^ 2;
t551 = t559 ^ 2;
t384 = t388 * qJD(4);
t467 = t422 * t501 - t425 * t484;
t341 = qJD(2) * t384 + t467;
t548 = pkin(4) * t341;
t544 = g(3) * t431;
t540 = qJ(5) * t379;
t539 = qJDD(2) * pkin(2);
t538 = qJDD(4) * pkin(4);
t509 = qJD(6) * t430;
t494 = t430 * qJDD(4) + t427 * t341 + t379 * t509;
t510 = qJD(6) * t427;
t304 = -qJD(4) * t510 + t494;
t537 = t304 * t430;
t387 = -t492 + t529;
t415 = pkin(3) * t425 + pkin(2);
t459 = -qJ(5) * t388 - t415;
t316 = t387 * t550 + t459;
t536 = t316 * t334;
t535 = t334 * t427;
t534 = t358 * t379;
t533 = t360 * t379;
t532 = t379 * t559;
t531 = t387 * t427;
t528 = t423 * t424;
t524 = t427 * t431;
t523 = qJDD(1) - g(3);
t522 = -pkin(5) * t384 - t520;
t383 = -t425 * t487 + t488;
t521 = -t383 * pkin(5) + t519;
t503 = qJDD(2) * qJ(3);
t504 = qJDD(1) * t424;
t361 = t429 * t504 + t503 + (qJD(3) + t490) * qJD(2);
t480 = qJDD(1) * t542;
t332 = t425 * t361 + t422 * t480;
t508 = t310 * qJD(4);
t506 = pkin(5) * t559 + t558;
t502 = qJDD(2) * t425;
t500 = qJDD(4) * qJ(5);
t495 = g(3) * t527;
t485 = t431 * t504;
t482 = t424 * t541;
t478 = t430 * t564;
t476 = qJDD(4) * t427 - t430 * t341;
t404 = t425 * t480;
t322 = t404 + (-pkin(8) * qJDD(2) - t361) * t422;
t323 = pkin(8) * t502 + t332;
t475 = -t428 * t322 - t549 * t323 - t350 * t487 + t351 * t511;
t376 = t423 * t431 + t429 * t469;
t378 = -t429 * t483 + t431 * t541;
t472 = g(1) * t378 + g(2) * t376;
t463 = qJ(5) * t383 - qJD(5) * t388;
t471 = -t384 * t550 - t463 + t491;
t315 = pkin(4) * t384 + t463;
t470 = -t315 + t491;
t297 = -qJD(4) * t550 + t506;
t370 = -qJD(2) * t415 + t468;
t439 = -qJ(5) * t559 + t370;
t301 = t379 * t550 + t439;
t291 = t297 * t430 - t301 * t427;
t292 = t297 * t427 + t301 * t430;
t362 = -t392 * t422 + t406;
t464 = t362 * t422 - t363 * t425;
t454 = -t485 + t557;
t364 = t454 - t539;
t458 = -t364 + t473;
t455 = MDP(3) + t566;
t373 = -t422 * t527 + t425 * t542;
t374 = t422 * t542 + t425 * t527;
t327 = t428 * t373 + t374 * t549;
t453 = t384 * t427 + t387 * t509;
t452 = -t322 * t549 + t428 * t323 + t350 * t511 + t351 * t487;
t343 = t376 * t417 - t416 * t482;
t345 = t378 * t417 + t416 * t528;
t368 = t416 * t542 + t417 * t527;
t451 = -g(1) * t345 - g(2) * t343 - g(3) * t368;
t450 = -t478 * t564 + t535;
t447 = qJDD(5) + t452;
t443 = t485 + t560;
t331 = -t361 * t422 + t404;
t442 = -t331 * t422 + t332 * t425 - t472;
t293 = -qJD(4) * qJD(5) + t475 - t500;
t441 = t451 - t475;
t352 = -qJDD(2) * t415 + t454;
t290 = -pkin(5) * t341 - t293;
t324 = t388 * pkin(5) + t354;
t440 = t290 * t387 + t298 * t384 + t324 * t334 + t472;
t342 = t376 * t416 + t417 * t482;
t344 = t378 * t416 - t417 * t528;
t367 = t416 * t527 - t417 * t542;
t438 = g(1) * t344 + g(2) * t342 + g(3) * t367 - t452;
t437 = t290 + (t550 * t564 + t540) * t564 + t451;
t436 = qJ(5) * t340 + t352;
t435 = -qJD(5) * t559 + t436;
t313 = pkin(4) * t379 + t439;
t434 = t313 * t559 + qJDD(5) - t438;
t433 = qJD(2) ^ 2;
t401 = t424 * t524;
t391 = -qJD(2) * pkin(2) + t468;
t371 = qJD(4) * t379;
t335 = pkin(4) * t387 + t459;
t333 = pkin(4) * t559 + t540;
t326 = -t373 * t549 + t428 * t374;
t325 = -t387 * pkin(5) + t355;
t307 = qJD(2) * t445 + qJD(4) * t327;
t306 = -t373 * t487 - qJD(2) * t466 + (qJD(4) * t374 + t514 * t526) * t428;
t305 = t360 * qJD(6) + t476;
t302 = -qJD(4) * pkin(4) + t558;
t296 = t435 + t548;
t295 = t341 * t550 + t435;
t294 = t447 - t538;
t289 = -t340 * pkin(5) - qJDD(4) * t550 + t447;
t288 = t430 * t289;
t1 = [t523 * MDP(1) + (t331 * t373 + t332 * t374 - g(3)) * MDP(8) + (t306 * t379 + t307 * t559 - t326 * t340 - t327 * t341) * MDP(16) + (-t293 * t327 + t294 * t326 + t302 * t307 + t303 * t306 - g(3)) * MDP(19) + ((t307 * t430 - t326 * t510) * t564 - (t326 * t430 + t401) * t334 - t306 * t358 + t327 * t305) * MDP(25) + (-(t307 * t427 + t326 * t509) * t564 + t326 * t535 - t306 * t360 + t327 * t304) * MDP(26) + (-t373 * t422 + t374 * t425) * MDP(7) * qJDD(2) + ((-qJDD(2) * MDP(4) - t455 * t433 + (MDP(19) * t313 + MDP(8) * t391 + t496 * t559 + t497 * t379 + (-t427 * MDP(25) - t430 * MDP(26)) * t564) * qJD(2)) * t429 + ((-t362 * t514 + t363 * t513 - t364) * MDP(8) - t296 * MDP(19) + t564 * MDP(25) * t509 + (-t510 * t564 - t328) * MDP(26) - t497 * t341 + t496 * t340 + (-MDP(4) + t562) * t433 + t455 * qJDD(2)) * t431) * t424 + t496 * (qJD(4) * t306 - qJDD(4) * t327) - t497 * (qJD(4) * t307 + qJDD(4) * t326); qJDD(2) * MDP(2) + t443 * MDP(3) + (-t523 * t527 + t472) * MDP(4) + (-t495 + t442 + (qJD(2) * t468 + t503) * t516) * MDP(7) + (-t464 * qJD(3) + t458 * pkin(2) + t442 * qJ(3) + (-g(3) * (pkin(2) * t431 + qJ(3) * t429) + (-t391 * t429 + t431 * t464) * qJD(1)) * t424) * MDP(8) + (-t340 * t388 - t383 * t559) * MDP(9) + (t340 * t387 - t341 * t388 + t379 * t383 - t384 * t559) * MDP(10) + (-qJD(4) * t383 + qJDD(4) * t388) * MDP(11) + (-qJD(4) * t384 - qJDD(4) * t387) * MDP(12) + (-t341 * t415 + t352 * t387 + t370 * t384 - t379 * t491 - t553) * MDP(14) + (t340 * t415 + t352 * t388 - t370 * t383 - t491 * t559 - t565) * MDP(15) + (t293 * t387 + t294 * t388 - t302 * t383 + t303 * t384 - t340 * t354 - t341 * t355 + t379 * t520 + t519 * t559 - t472 - t495) * MDP(16) + (-t296 * t387 - t313 * t384 - t335 * t341 + t379 * t470 + t553) * MDP(17) + (-t296 * t388 + t313 * t383 + t335 * t340 + t470 * t559 + t565) * MDP(18) + (-t293 * t355 + t294 * t354 + t296 * t335 + t313 * t315 - t472 * t543 + t520 * t303 + t519 * t302 + (-g(3) * t543 - t313 * qJD(1)) * t527 + (-t544 * t424 + t473) * (pkin(4) * t417 + qJ(5) * t416 + t415)) * MDP(19) + (t304 * t531 + t360 * t453) * MDP(20) + ((-t358 * t427 + t360 * t430) * t384 + (t537 - t305 * t427 + (-t358 * t430 - t360 * t427) * qJD(6)) * t387) * MDP(21) + (t304 * t388 - t334 * t531 - t360 * t383 + t453 * t564) * MDP(22) + (-t387 * t328 - t305 * t388 + t358 * t383 + (t384 * t430 - t387 * t510) * t564) * MDP(23) + (-t334 * t388 - t383 * t564) * MDP(24) + (t288 * t388 - t291 * t383 + t325 * t305 + (-t295 * t388 + t416 * t473 + t536) * t427 - t440 * t430 - g(3) * (t416 * t524 + t429 * t430) * t424 + (t471 * t427 + t521 * t430) * t564 + t522 * t358 + ((-t316 * t430 - t324 * t427) * t564 - t292 * t388 + t298 * t531) * qJD(6)) * MDP(25) + (t292 * t383 + t325 * t304 + t522 * t360 + (t536 - (qJD(6) * t297 + t295) * t388 + t298 * qJD(6) * t387 + (-qJD(6) * t324 + t471) * t564 + t444) * t430 + (-(-qJD(6) * t301 + t289) * t388 + t495 + (qJD(6) * t316 - t521) * t564 + t440) * t427) * MDP(26) + t566 * (t424 * (t486 - t544) + t458 + t539); -MDP(5) * t502 + qJDD(2) * t530 - t433 * t562 + (qJD(2) * t464 - t443 - t539 + t557) * MDP(8) + ((-t379 - t489) * qJD(4) + t493) * MDP(15) + (-t551 - t552) * MDP(16) + (t340 + t371) * MDP(18) + (t548 - t303 * t379 + (-qJD(5) - t302) * t559 + t436 - t560) * MDP(19) + (t450 + t534) * MDP(25) + (t533 - t567) * MDP(26) + t497 * (0.2e1 * t559 * qJD(4) + t467); MDP(9) * t532 + (t551 - t552) * MDP(10) + ((t379 - t489) * qJD(4) + t493) * MDP(11) - t467 * MDP(12) + qJDD(4) * MDP(13) + (-t370 * t559 + t438 + t508) * MDP(14) + (qJD(4) * t518 + t370 * t379 - t441) * MDP(15) + (pkin(4) * t340 - qJ(5) * t341 + (-t303 - t310) * t559 + (t302 - t558) * t379) * MDP(16) + (t333 * t379 + t434 - t508 - 0.2e1 * t538) * MDP(17) + (0.2e1 * t500 - t313 * t379 + t333 * t559 + (0.2e1 * qJD(5) - t518) * qJD(4) + t441) * MDP(18) + (-t293 * qJ(5) - t294 * pkin(4) - t313 * t333 - t302 * t310 - g(1) * (-pkin(4) * t344 + qJ(5) * t345) - g(2) * (-pkin(4) * t342 + qJ(5) * t343) - g(3) * (-pkin(4) * t367 + qJ(5) * t368) - t558 * t303) * MDP(19) + (-t427 * t477 + t537) * MDP(20) + ((-t305 - t477) * t430 + (-t304 + t568) * t427) * MDP(21) + (t533 + t567) * MDP(22) + (t450 - t534) * MDP(23) + t564 * t379 * MDP(24) + (qJ(5) * t305 + t291 * t379 + t506 * t358 + t437 * t427 + t430 * t556) * MDP(25) + (qJ(5) * t304 - t292 * t379 + t506 * t360 - t427 * t556 + t437 * t430) * MDP(26); (-t340 + t371) * MDP(16) + (qJDD(4) - t532) * MDP(17) + (-qJD(4) ^ 2 - t551) * MDP(18) + (t303 * qJD(4) + t434 - t538) * MDP(19) + (-qJD(4) * t358 - t328) * MDP(25) + (-qJD(4) * t360 + t535) * MDP(26) + (-MDP(25) * t479 - MDP(26) * t478) * t564; t360 * t358 * MDP(20) + (-t358 ^ 2 + t360 ^ 2) * MDP(21) + (t494 + t568) * MDP(22) + (-t476 + t477) * MDP(23) - t334 * MDP(24) + (-t427 * t295 + t288 + t292 * t564 - t298 * t360 - g(1) * (t344 * t430 - t377 * t427) - g(2) * (t342 * t430 - t375 * t427) - g(3) * (t367 * t430 + t401)) * MDP(25) + (-t430 * t295 - t427 * t289 + t291 * t564 + t298 * t358 - g(1) * (-t344 * t427 - t377 * t430) - g(2) * (-t342 * t427 - t375 * t430) - g(3) * (-t367 * t427 + t430 * t526)) * MDP(26) + (-MDP(22) * t512 - MDP(23) * t360 - MDP(25) * t292 - MDP(26) * t291) * qJD(6);];
tau  = t1;
