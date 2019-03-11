% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:50
% EndTime: 2019-03-09 02:18:57
% DurationCPUTime: 5.28s
% Computational Cost: add. (4833->404), mult. (10766->516), div. (0->0), fcn. (8582->18), ass. (0->191)
t442 = sin(pkin(11));
t448 = sin(qJ(4));
t444 = cos(pkin(11));
t452 = cos(qJ(4));
t524 = t452 * t444;
t399 = t442 * t448 - t524;
t389 = t399 * qJD(1);
t400 = t442 * t452 + t444 * t448;
t390 = t400 * qJD(1);
t447 = sin(qJ(5));
t451 = cos(qJ(5));
t348 = t451 * t389 + t390 * t447;
t450 = cos(qJ(6));
t510 = qJD(6) * t450;
t570 = t348 * t450 + t510;
t439 = pkin(11) + qJ(4);
t435 = qJ(5) + t439;
t421 = sin(t435);
t441 = qJ(1) + pkin(10);
t432 = sin(t441);
t434 = cos(t441);
t479 = g(1) * t434 + g(2) * t432;
t569 = t479 * t421;
t468 = -t389 * t447 + t451 * t390;
t515 = qJD(1) * t448;
t495 = t442 * t515;
t503 = qJDD(1) * t452;
t504 = qJDD(1) * t448;
t514 = qJD(4) * t452;
t516 = qJD(1) * t444;
t496 = t442 * t503 + t444 * t504 + t514 * t516;
t356 = -qJD(4) * t495 + t496;
t392 = t400 * qJD(4);
t413 = t444 * t503;
t476 = -t442 * t504 + t413;
t357 = qJD(1) * t392 - t476;
t512 = qJD(5) * t451;
t513 = qJD(5) * t447;
t320 = t451 * t356 - t447 * t357 - t389 * t512 - t390 * t513;
t436 = qJDD(4) + qJDD(5);
t440 = qJD(4) + qJD(5);
t446 = sin(qJ(6));
t497 = t450 * t320 + t446 * t436 + t440 * t510;
t511 = qJD(6) * t446;
t300 = -t468 * t511 + t497;
t299 = t300 * t450;
t344 = t440 * t446 + t450 * t468;
t417 = t450 * t436;
t301 = qJD(6) * t344 + t320 * t446 - t417;
t342 = -t450 * t440 + t446 * t468;
t568 = -t446 * t301 - t570 * t342 + t299;
t298 = t300 * t446;
t321 = qJD(5) * t468 + t356 * t447 + t451 * t357;
t318 = qJDD(6) + t321;
t315 = t446 * t318;
t532 = t348 * t440;
t534 = t468 * t440;
t536 = t344 * t468;
t561 = qJD(6) + t348;
t567 = t436 * MDP(20) + (-t321 + t534) * MDP(19) - t348 ^ 2 * MDP(17) + (t348 * MDP(16) + MDP(17) * t468 - MDP(27) * t561) * t468 + (t320 + t532) * MDP(18) + (t570 * t344 + t298) * MDP(23) + (t570 * t561 + t315 - t536) * MDP(25);
t443 = sin(pkin(10));
t419 = pkin(1) * t443 + qJ(3);
t408 = t419 * qJD(1);
t430 = t444 * qJD(2);
t367 = t430 + (-pkin(7) * qJD(1) - t408) * t442;
t378 = t442 * qJD(2) + t444 * t408;
t368 = pkin(7) * t516 + t378;
t555 = t452 * t367 - t368 * t448;
t331 = -pkin(8) * t390 + t555;
t330 = qJD(4) * pkin(4) + t331;
t470 = -t367 * t448 - t368 * t452;
t332 = -pkin(8) * t389 - t470;
t539 = t332 * t447;
t310 = t330 * t451 - t539;
t306 = -pkin(5) * t440 - t310;
t566 = t306 * t348;
t397 = qJD(1) * qJD(3) + qJDD(1) * t419;
t428 = t444 * qJDD(2);
t542 = pkin(7) * qJDD(1);
t363 = t428 + (-t397 - t542) * t442;
t370 = t442 * qJDD(2) + t444 * t397;
t364 = t444 * t542 + t370;
t489 = t452 * t363 - t448 * t364;
t308 = qJDD(4) * pkin(4) - pkin(8) * t356 + qJD(4) * t470 + t489;
t471 = t448 * t363 + t452 * t364;
t309 = -pkin(8) * t357 + qJD(4) * t555 + t471;
t538 = t332 * t451;
t311 = t330 * t447 + t538;
t549 = qJD(5) * t311 - t451 * t308 + t309 * t447;
t289 = -pkin(5) * t436 + t549;
t422 = cos(t435);
t544 = g(3) * t422;
t564 = t289 + t544;
t445 = cos(pkin(10));
t423 = -pkin(1) * t445 - pkin(2);
t406 = -pkin(3) * t444 + t423;
t387 = qJD(1) * t406 + qJD(3);
t355 = pkin(4) * t389 + t387;
t415 = g(3) * t421;
t550 = t451 * (qJD(5) * t330 + t309) + t308 * t447 - t332 * t513;
t560 = t348 * t355 + t422 * t479 + t415 - t550;
t557 = pkin(5) * t468;
t537 = t342 * t468;
t556 = t561 * (pkin(9) * t561 + t557);
t543 = pkin(7) + t419;
t393 = t543 * t442;
t394 = t543 * t444;
t520 = -t448 * t393 + t452 * t394;
t307 = pkin(9) * t440 + t311;
t322 = pkin(5) * t348 - pkin(9) * t468 + t355;
t292 = -t307 * t446 + t322 * t450;
t554 = -t292 * t468 + t306 * t511 + t450 * t569;
t293 = t307 * t450 + t322 * t446;
t553 = t293 * t468 + t306 * t510 + t446 * t564;
t552 = -t468 * t355 - t544 - t549 + t569;
t460 = -t393 * t514 + qJD(3) * t524 + (-qJD(3) * t442 - qJD(4) * t394) * t448;
t335 = -pkin(8) * t392 + t460;
t391 = t399 * qJD(4);
t457 = -t400 * qJD(3) - qJD(4) * t520;
t336 = pkin(8) * t391 + t457;
t488 = -t452 * t393 - t394 * t448;
t340 = -pkin(8) * t400 + t488;
t341 = -pkin(8) * t399 + t520;
t472 = t340 * t451 - t341 * t447;
t294 = qJD(5) * t472 + t335 * t451 + t336 * t447;
t324 = t340 * t447 + t341 * t451;
t360 = t451 * t399 + t400 * t447;
t361 = -t399 * t447 + t400 * t451;
t366 = pkin(4) * t399 + t406;
t326 = pkin(5) * t360 - pkin(9) * t361 + t366;
t333 = -qJD(5) * t360 - t391 * t451 - t392 * t447;
t288 = pkin(9) * t436 + t550;
t485 = qJD(6) * t322 + t288;
t548 = t289 * t361 + t306 * t333 - t324 * t318 - (qJD(6) * t326 + t294) * t561 - t360 * t485;
t547 = pkin(4) * t392;
t541 = t306 * t361;
t540 = t326 * t318;
t535 = t344 * t446;
t529 = t432 * t446;
t528 = t432 * t450;
t527 = t434 * t446;
t526 = t434 * t450;
t525 = t444 * MDP(5);
t316 = t450 * t318;
t334 = qJD(5) * t361 - t391 * t447 + t451 * t392;
t523 = t300 * t360 + t344 * t334;
t519 = t442 ^ 2 + t444 ^ 2;
t517 = qJD(1) * t442;
t505 = qJDD(1) * t423;
t502 = t361 * t315;
t501 = t361 * t316;
t487 = t561 * t446;
t385 = qJDD(1) * t406 + qJDD(3);
t339 = pkin(4) * t357 + t385;
t291 = pkin(5) * t321 - pkin(9) * t320 + t339;
t484 = qJD(6) * t307 - t291;
t425 = pkin(4) * t447 + pkin(9);
t482 = pkin(4) * t390 + pkin(9) * t348 + qJD(6) * t425 + t557;
t312 = t331 * t447 + t538;
t481 = pkin(4) * t513 - t312;
t313 = t331 * t451 - t539;
t480 = -pkin(4) * t512 + t313;
t478 = g(1) * t432 - g(2) * t434;
t449 = sin(qJ(1));
t453 = cos(qJ(1));
t477 = g(1) * t449 - g(2) * t453;
t475 = -t301 * t360 - t334 * t342;
t474 = -t318 * t425 + t566;
t473 = t333 * t440 + t361 * t436;
t369 = -t397 * t442 + t428;
t469 = -t369 * t442 + t370 * t444;
t467 = t316 - (t348 * t446 + t511) * t561;
t465 = -t333 * t446 - t361 * t510;
t464 = -t333 * t450 + t361 * t511;
t463 = -pkin(9) * t318 + t310 * t561 + t566;
t433 = cos(t439);
t431 = sin(t439);
t426 = -pkin(4) * t451 - pkin(5);
t405 = qJDD(3) + t505;
t377 = -t408 * t442 + t430;
t374 = t422 * t526 + t529;
t373 = -t422 * t527 + t528;
t372 = -t422 * t528 + t527;
t371 = t422 * t529 + t526;
t359 = -qJD(4) * t392 - qJDD(4) * t399;
t358 = -qJD(4) * t391 + qJDD(4) * t400;
t319 = -t334 * t440 - t360 * t436;
t302 = pkin(5) * t334 - pkin(9) * t333 + t547;
t295 = qJD(5) * t324 + t335 * t447 - t336 * t451;
t290 = t450 * t291;
t1 = [(-qJD(4) * t460 - qJDD(4) * t520 + t406 * t356 + t385 * t400 - t387 * t391 - t431 * t478) * MDP(15) + (t356 * t400 - t390 * t391) * MDP(9) + (-t294 * t440 + t320 * t366 - t324 * t436 + t333 * t355 + t339 * t361 - t421 * t478 + t468 * t547) * MDP(22) + (-t320 * t360 - t321 * t361 - t333 * t348 - t334 * t468) * MDP(17) + (t320 * t361 + t333 * t468) * MDP(16) + (t397 * t519 + t469 - t479) * MDP(7) + (-t442 * MDP(6) + t525) * (-t405 + t478 - t505) + (-t295 * t440 + t321 * t366 + t334 * t355 + t339 * t360 + t348 * t547 + t422 * t478 + t436 * t472) * MDP(21) + (-t464 * t561 + t501 + t523) * MDP(25) + (t465 * t561 + t475 - t502) * MDP(26) + (t318 * t360 + t334 * t561) * MDP(27) + (-g(1) * t371 - g(2) * t373 - t293 * t334 + t295 * t344 - t472 * t300 + (-(-qJD(6) * t324 + t302) * t561 - t540 + t484 * t360 - qJD(6) * t541) * t446 + t548 * t450) * MDP(29) + (-g(1) * t372 - g(2) * t374 + t290 * t360 + t292 * t334 + t295 * t342 - t472 * t301 + (t302 * t561 + t540 + (-t307 * t360 - t324 * t561 + t541) * qJD(6)) * t450 + t548 * t446) * MDP(28) + qJDD(1) * MDP(1) + (qJD(4) * t457 + qJDD(4) * t488 + t406 * t357 + t385 * t399 + t387 * t392 + t433 * t478) * MDP(14) + t477 * MDP(2) + t473 * MDP(18) + (t405 * t423 - g(1) * (-pkin(1) * t449 - pkin(2) * t432 + qJ(3) * t434) - g(2) * (pkin(1) * t453 + pkin(2) * t434 + qJ(3) * t432) + t469 * t419 + (-t377 * t442 + t378 * t444) * qJD(3)) * MDP(8) + (t477 + (t443 ^ 2 + t445 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t299 * t361 - t344 * t464) * MDP(23) + ((-t342 * t450 - t535) * t333 + (-t298 - t301 * t450 + (t342 * t446 - t344 * t450) * qJD(6)) * t361) * MDP(24) + (g(1) * t453 + g(2) * t449) * MDP(3) + (-t356 * t399 - t357 * t400 + t389 * t391 - t390 * t392) * MDP(10) + t358 * MDP(11) + t359 * MDP(12) + t319 * MDP(19); (qJDD(2) - g(3)) * MDP(4) + (t369 * t444 + t370 * t442 - g(3)) * MDP(8) + t359 * MDP(14) - t358 * MDP(15) + t319 * MDP(21) - t473 * MDP(22) + (-t475 - t502) * MDP(28) + (-t501 + t523) * MDP(29) + (MDP(28) * t465 + MDP(29) * t464) * t561; (t377 * t517 - t378 * t516 + qJDD(3) - t478) * MDP(8) - t413 * MDP(14) + t496 * MDP(15) + (t321 + t534) * MDP(21) + (t320 - t532) * MDP(22) + (t467 - t537) * MDP(28) + (-t450 * t561 ^ 2 - t315 - t536) * MDP(29) - t519 * MDP(7) * qJD(1) ^ 2 + (-t525 + t423 * MDP(8) + (MDP(14) * t448 + MDP(6)) * t442) * qJDD(1) + ((t444 * t515 + t452 * t517 + t390) * MDP(14) + (-t389 - t495) * MDP(15)) * qJD(4); (t496 + (t389 - t495) * qJD(4)) * MDP(11) + (g(3) * t431 + t387 * t389 + t433 * t479 - t471) * MDP(15) + (-t535 * t561 + t568) * MDP(24) + (t312 * t440 + (-t348 * t390 + t436 * t451 - t440 * t513) * pkin(4) + t552) * MDP(21) + (t313 * t440 + (-t390 * t468 - t436 * t447 - t440 * t512) * pkin(4) + t560) * MDP(22) + (t426 * t301 - t564 * t450 + t474 * t446 + t481 * t342 + (t446 * t480 - t450 * t482) * t561 + t554) * MDP(28) + (t426 * t300 + t474 * t450 - t446 * t569 + t481 * t344 + (t446 * t482 + t450 * t480) * t561 + t553) * MDP(29) + (t467 + t537) * MDP(26) + qJDD(4) * MDP(13) + t390 * t389 * MDP(9) + (-g(3) * t433 - t387 * t390 + t431 * t479 + t489) * MDP(14) + t476 * MDP(12) + (-t389 ^ 2 + t390 ^ 2) * MDP(10) + t567; (t311 * t440 + t552) * MDP(21) + (t310 * t440 + t560) * MDP(22) + (-t344 * t487 + t568) * MDP(24) + (-t487 * t561 + t316 + t537) * MDP(26) + (-pkin(5) * t301 - t311 * t342 + t463 * t446 + (-t564 - t556) * t450 + t554) * MDP(28) + (-pkin(5) * t300 - t311 * t344 + t463 * t450 + (-t569 + t556) * t446 + t553) * MDP(29) + t567; t344 * t342 * MDP(23) + (-t342 ^ 2 + t344 ^ 2) * MDP(24) + (t342 * t561 + t497) * MDP(25) + (t344 * t561 + t417) * MDP(26) + t318 * MDP(27) + (-g(1) * t373 + g(2) * t371 + t293 * t561 - t306 * t344 + t290) * MDP(28) + (g(1) * t374 - g(2) * t372 + t292 * t561 + t306 * t342) * MDP(29) + ((-t288 + t415) * MDP(29) + (-MDP(26) * t468 - MDP(28) * t307 - MDP(29) * t322) * qJD(6)) * t450 + (-qJD(6) * t468 * MDP(25) + (-qJD(6) * t440 - t320) * MDP(26) + (-t485 + t415) * MDP(28) + t484 * MDP(29)) * t446;];
tau  = t1;
