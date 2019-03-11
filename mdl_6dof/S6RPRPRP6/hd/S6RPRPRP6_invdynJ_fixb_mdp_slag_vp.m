% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:55
% EndTime: 2019-03-09 03:20:03
% DurationCPUTime: 6.50s
% Computational Cost: add. (4035->455), mult. (9197->546), div. (0->0), fcn. (6745->10), ass. (0->206)
t436 = cos(pkin(9));
t556 = cos(qJ(3));
t496 = t556 * t436;
t478 = qJD(1) * t496;
t435 = sin(pkin(9));
t440 = sin(qJ(3));
t531 = t435 * t440;
t495 = qJD(1) * t531;
t386 = -t478 + t495;
t439 = sin(qJ(5));
t442 = cos(qJ(5));
t394 = t435 * t556 + t440 * t436;
t391 = t394 * qJD(3);
t490 = qJDD(1) * t556;
t503 = qJDD(1) * t440;
t474 = t435 * t503 - t436 * t490;
t447 = qJD(1) * t391 + t474;
t510 = qJD(5) * t442;
t511 = qJD(5) * t439;
t316 = qJD(3) * t511 - t442 * qJDD(3) - t386 * t510 - t439 * t447;
t359 = qJD(3) * t439 - t442 * t386;
t569 = t394 * qJD(1);
t576 = qJD(5) + t569;
t482 = t576 * t359;
t580 = t316 - t482;
t431 = pkin(9) + qJ(3);
t426 = cos(t431);
t441 = sin(qJ(1));
t443 = cos(qJ(1));
t552 = g(1) * t443;
t476 = g(2) * t441 + t552;
t425 = sin(t431);
t549 = g(3) * t425;
t584 = -t476 * t426 - t549;
t361 = qJD(3) * t442 + t386 * t439;
t572 = t442 * t576;
t583 = t361 * t572;
t513 = qJD(3) * t440;
t494 = t435 * t513;
t498 = qJD(3) * t478 + t435 * t490 + t436 * t503;
t353 = qJD(1) * t494 - t498;
t351 = -qJDD(5) + t353;
t422 = pkin(2) * t436 + pkin(1);
t398 = -qJD(1) * t422 + qJD(2);
t456 = -qJ(4) * t569 + t398;
t557 = pkin(3) + pkin(8);
t321 = t386 * t557 + t456;
t547 = pkin(7) + qJ(2);
t401 = t547 * t435;
t395 = qJD(1) * t401;
t402 = t547 * t436;
t396 = qJD(1) * t402;
t518 = -t556 * t395 - t440 * t396;
t567 = qJD(4) - t518;
t470 = pkin(4) * t569 + t567;
t326 = -qJD(3) * t557 + t470;
t304 = t321 * t442 + t326 * t439;
t544 = qJDD(1) * pkin(1);
t424 = qJDD(2) - t544;
t570 = qJD(3) * t569;
t448 = -t474 - t570;
t504 = qJDD(1) * t436;
t310 = -pkin(2) * t504 - pkin(3) * t448 + t353 * qJ(4) - qJD(4) * t569 + t424;
t300 = pkin(8) * t447 + t310;
t505 = qJD(1) * qJD(2);
t561 = qJDD(1) * t547 + t505;
t367 = t561 * t435;
t368 = t561 * t436;
t493 = qJD(3) * t556;
t463 = t367 * t556 + t440 * t368 - t395 * t513 + t396 * t493;
t458 = qJDD(4) + t463;
t307 = -t353 * pkin(4) - qJDD(3) * t557 + t458;
t487 = -t300 * t439 + t442 * t307;
t451 = -t304 * qJD(5) + t487;
t290 = -pkin(5) * t351 + qJ(6) * t316 - qJD(6) * t361 + t451;
t297 = -qJ(6) * t359 + t304;
t582 = t297 * t576 + t290;
t512 = qJD(5) * t361;
t317 = qJDD(3) * t439 - t442 * t447 + t512;
t500 = -t442 * t300 - t439 * t307 - t326 * t510;
t462 = -t321 * t511 - t500;
t291 = -qJ(6) * t317 - qJD(6) * t359 + t462;
t303 = -t321 * t439 + t442 * t326;
t296 = -qJ(6) * t361 + t303;
t295 = pkin(5) * t576 + t296;
t581 = -t295 * t576 + t291;
t345 = t442 * t351;
t485 = t439 * t576;
t579 = -t485 * t576 - t345;
t571 = g(1) * t441 - g(2) * t443;
t467 = -t424 + t571;
t480 = t440 * t367 - t556 * t368 + t395 * t493 + t396 * t513;
t578 = qJD(3) * t518 + t480 - t584;
t393 = -t496 + t531;
t468 = -qJ(4) * t394 - t422;
t334 = t393 * t557 + t468;
t356 = t401 * t556 + t440 * t402;
t343 = t394 * pkin(4) + t356;
t521 = t442 * t334 + t439 * t343;
t516 = t426 * pkin(3) + t425 * qJ(4);
t566 = qJ(2) * qJDD(1);
t526 = t442 * t443;
t529 = t439 * t441;
t378 = t425 * t526 - t529;
t527 = t441 * t442;
t528 = t439 * t443;
t380 = t425 * t527 + t528;
t420 = g(3) * t426;
t565 = -g(1) * t378 - g(2) * t380 + t442 * t420;
t355 = -t440 * t395 + t556 * t396;
t336 = -pkin(4) * t386 + t355;
t434 = qJD(3) * qJ(4);
t329 = t336 + t434;
t563 = t329 * t576 + t351 * t557;
t339 = t440 * (qJD(2) * t435 + qJD(3) * t402) - qJD(2) * t496 + t401 * t493;
t357 = -t440 * t401 + t402 * t556;
t562 = qJD(3) * t339 - qJDD(3) * t357 - t425 * t571;
t560 = t361 ^ 2;
t559 = t386 ^ 2;
t558 = t569 ^ 2;
t555 = pkin(5) * t439;
t437 = -qJ(6) - pkin(8);
t548 = pkin(3) - t437;
t546 = qJ(4) * t386;
t545 = qJ(6) * t442;
t543 = qJDD(3) * pkin(3);
t542 = t316 * t442;
t541 = t359 * t386;
t540 = t361 * t386;
t539 = t386 * t569;
t538 = t393 * t439;
t537 = t425 * t441;
t536 = t425 * t443;
t535 = t426 * t437;
t534 = t426 * t441;
t533 = t426 * t443;
t532 = t435 * MDP(5);
t530 = t439 * t351;
t524 = qJ(6) + t557;
t523 = -t296 + t295;
t332 = t557 * t569 + t546;
t522 = t442 * t332 + t439 * t336;
t331 = t442 * t336;
t508 = qJD(6) * t442;
t520 = t511 * t524 - t508 + pkin(5) * t386 - t331 - (-qJ(6) * t569 - t332) * t439;
t400 = t524 * t442;
t519 = -qJD(5) * t400 - qJD(6) * t439 - t545 * t569 - t522;
t423 = pkin(5) * t442 + pkin(4);
t515 = t423 + t547;
t514 = t435 ^ 2 + t436 ^ 2;
t509 = qJD(5) * t557;
t507 = t355 * qJD(3);
t502 = t426 * t555;
t499 = pkin(3) * t533 + qJ(4) * t536 + t443 * t422;
t497 = -g(1) * t536 - g(2) * t537 + t420;
t492 = qJ(4) + t555;
t489 = t514 * qJD(1) ^ 2;
t488 = -qJ(6) * t393 - t334;
t486 = t576 ^ 2;
t479 = 0.2e1 * t514;
t477 = -g(1) * t534 + g(2) * t533;
t390 = -t436 * t493 + t494;
t472 = qJ(4) * t390 - qJD(4) * t394;
t432 = qJDD(3) * qJ(4);
t433 = qJD(3) * qJD(4);
t311 = -t432 - t433 + t480;
t465 = t391 * t439 + t393 * t510;
t464 = -t391 * t442 + t393 * t511;
t318 = t391 * t557 + t472;
t340 = qJD(2) * t394 + qJD(3) * t357;
t323 = -t390 * pkin(4) + t340;
t461 = t442 * t318 + t439 * t323 - t334 * t511 + t343 * t510;
t460 = -t572 * t576 + t530;
t454 = -t463 - t497;
t453 = -t340 * qJD(3) - t356 * qJDD(3) - t477;
t308 = -pkin(4) * t447 - t311;
t452 = t308 + t584;
t450 = t479 * t505 - t476;
t322 = -pkin(4) * t391 - t339;
t338 = pkin(3) * t386 + t456;
t449 = t338 * t569 + qJDD(4) - t454;
t294 = t317 * pkin(5) + qJDD(6) + t308;
t405 = qJ(4) * t533;
t403 = qJ(4) * t534;
t399 = t524 * t439;
t397 = -qJDD(1) * t422 + qJDD(2);
t381 = -t425 * t529 + t526;
t379 = t425 * t528 + t527;
t372 = qJD(3) * t386;
t358 = t359 ^ 2;
t352 = pkin(3) * t393 + t468;
t350 = pkin(3) * t569 + t546;
t349 = -t434 - t355;
t348 = -qJD(3) * pkin(3) + t567;
t344 = -t393 * pkin(4) + t357;
t342 = t442 * t343;
t333 = pkin(3) * t391 + t472;
t320 = t442 * t323;
t314 = t442 * t317;
t313 = pkin(5) * t359 + qJD(6) + t329;
t312 = t458 - t543;
t309 = t393 * t545 + t521;
t301 = pkin(5) * t394 + t439 * t488 + t342;
t293 = -qJ(6) * t464 + t393 * t508 + t461;
t292 = -pkin(5) * t390 + t320 + t488 * t510 + (-qJ(6) * t391 - qJD(5) * t343 - qJD(6) * t393 - t318) * t439;
t1 = [((-t359 * t439 + t361 * t442) * t391 + (-t542 - t317 * t439 + (-t359 * t442 - t361 * t439) * qJD(5)) * t393) * MDP(20) + (-t292 * t361 - t293 * t359 + t301 * t316 - t309 * t317 + (-t295 * t439 + t297 * t442) * t391 + (-t290 * t439 + t291 * t442 + (-t295 * t442 - t297 * t439) * qJD(5)) * t393 - t477) * MDP(26) + t476 * MDP(3) + (-t316 * t538 + t361 * t465) * MDP(19) + (t436 * MDP(4) - t532) * (t467 + t544) + (t353 * t422 - t390 * t398 + t394 * t397 + t562) * MDP(14) + (-t310 * t393 - t333 * t386 - t338 * t391 + t352 * t448 - t453) * MDP(16) + (t398 * t391 + t397 * t393 - t422 * t447 + t453) * MDP(13) + (t291 * t309 + t297 * t293 + t290 * t301 + t295 * t292 + t294 * (-t393 * t423 + t357) + t313 * (pkin(5) * t464 + t322) - g(1) * (t515 * t443 + (-t425 * t492 - t426 * t548 - t422) * t441) - g(2) * ((t425 * t555 - t535) * t443 + t515 * t441 + t499)) * MDP(27) + (g(1) * t380 - g(2) * t378 + t304 * t390 + t308 * t538 - t344 * t316 + t322 * t361 + t329 * t465 + t351 * t521 - t394 * t462 - t461 * t576) * MDP(25) + (-t317 * t394 - t345 * t393 + t359 * t390 - t464 * t576) * MDP(22) + (-t316 * t394 - t361 * t390 - t393 * t530 + t465 * t576) * MDP(21) + ((-t318 * t439 + t320) * t576 - (-t334 * t439 + t342) * t351 + t487 * t394 - t303 * t390 + t322 * t359 + t344 * t317 - g(1) * t381 - g(2) * t379 + (-t308 * t393 - t329 * t391) * t442 + (-t304 * t394 + t329 * t538 - t521 * t576) * qJD(5)) * MDP(24) + (-t351 * t394 - t390 * t576) * MDP(23) + (t311 * t393 + t312 * t394 + t339 * t386 + t340 * t569 - t348 * t390 + t349 * t391 - t356 * t353 - t357 * t447 - t476) * MDP(15) + (t353 * t393 + t390 * t386 - t391 * t569 + t394 * t448) * MDP(9) + (-t310 * t394 - t333 * t569 + t338 * t390 + t352 * t353 - t562) * MDP(17) + (-t353 * t394 - t390 * t569) * MDP(8) + (t310 * t352 + t338 * t333 - t311 * t357 + t349 * t339 + t312 * t356 + t348 * t340 - t547 * t552 - g(2) * t499 + (-g(1) * (-t422 - t516) - g(2) * t547) * t441) * MDP(18) + t571 * MDP(2) + (pkin(1) * t467 + (t514 * t566 + t450) * qJ(2)) * MDP(7) + (t479 * t566 + t450) * MDP(6) + (-qJD(3) * t391 - qJDD(3) * t393) * MDP(11) + (-qJD(3) * t390 + qJDD(3) * t394) * MDP(10) + qJDD(1) * MDP(1); -MDP(4) * t504 + qJDD(1) * t532 - MDP(6) * t489 + (-qJ(2) * t489 - t467) * MDP(7) + ((-t386 - t495) * qJD(3) + t498) * MDP(14) + (-t558 - t559) * MDP(15) + (t353 + t372) * MDP(17) + (-t348 * t569 - t349 * t386 + t310 - t571) * MDP(18) + (t460 + t541) * MDP(24) + (t540 - t579) * MDP(25) + (-t439 * t580 - t314 + t583) * MDP(26) + (t313 * t386 - t439 * t582 + t442 * t581 - t571) * MDP(27) + (-MDP(13) + MDP(16)) * (-t474 - 0.2e1 * t570); MDP(8) * t539 + (t558 - t559) * MDP(9) + (-qJD(3) * t495 + t498) * MDP(10) - t474 * MDP(11) + qJDD(3) * MDP(12) + (-t398 * t569 + t454 + t507) * MDP(13) + t578 * MDP(14) + (pkin(3) * t353 - qJ(4) * t447 + (-t349 - t355) * t569) * MDP(15) + (t449 - t507 - 0.2e1 * t543) * MDP(16) + (t350 * t569 + 0.2e1 * t432 + 0.2e1 * t433 - t578) * MDP(17) + (-t311 * qJ(4) - t312 * pkin(3) - t338 * t350 - t348 * t355 - g(1) * (-pkin(3) * t536 + t405) - g(2) * (-pkin(3) * t537 + t403) - g(3) * t516 - t567 * t349) * MDP(18) + (-t361 * t485 - t542) * MDP(19) + (-t314 - t583 + (t316 + t482) * t439) * MDP(20) + (t540 + t579) * MDP(21) + (t460 - t541) * MDP(22) + (qJ(4) * t317 - t331 * t576 + t470 * t359 + t563 * t442 + ((t332 + t509) * t576 + t452) * t439) * MDP(24) + (-qJ(4) * t316 + t522 * t576 + t470 * t361 - t563 * t439 + (t509 * t576 + t452) * t442) * MDP(25) + (-t316 * t400 + t317 * t399 - t519 * t359 - t520 * t361 - t439 * t581 - t442 * t582 - t497) * MDP(26) + (-t291 * t399 - t290 * t400 + t294 * t492 - g(1) * (t443 * t502 + t405) - g(2) * (t441 * t502 + t403) - g(3) * (t516 - t535) + (-g(3) * t555 + t476 * t548) * t425 + (pkin(5) * t572 + t470) * t313 + t519 * t297 + t520 * t295) * MDP(27) + (qJD(3) * MDP(10) + t398 * MDP(14) + (t348 - t567) * MDP(15) + t350 * MDP(16) - t338 * MDP(17) + t576 * MDP(23) + t303 * MDP(24) - t304 * MDP(25)) * t386; (-t353 + t372) * MDP(15) + (qJDD(3) - t539) * MDP(16) + (-qJD(3) ^ 2 - t558) * MDP(17) + (t349 * qJD(3) + t449 - t543) * MDP(18) + (-qJD(3) * t359 - t345) * MDP(24) - qJD(3) * t361 * MDP(25) + (-qJD(3) * t313 + t497) * MDP(27) + (-MDP(25) * t486 + MDP(26) * t580 + MDP(27) * t582) * t442 + (t351 * MDP(25) + (t361 * t569 - t317 + t512) * MDP(26) + t581 * MDP(27) - MDP(24) * t486) * t439; t361 * t359 * MDP(19) + (-t358 + t560) * MDP(20) - t580 * MDP(21) + (t361 * t576 - t317) * MDP(22) - t351 * MDP(23) + (t304 * t576 - t329 * t361 + t451 + t565) * MDP(24) + (g(1) * t379 - g(2) * t381 + t303 * t576 + t329 * t359 + (qJD(5) * t321 - t420) * t439 + t500) * MDP(25) + (pkin(5) * t316 - t359 * t523) * MDP(26) + (t523 * t297 + (-t313 * t361 + t290 + t565) * pkin(5)) * MDP(27); (-t358 - t560) * MDP(26) + (-g(1) * t533 - g(2) * t534 + t295 * t361 + t297 * t359 + t294 - t549) * MDP(27);];
tau  = t1;
