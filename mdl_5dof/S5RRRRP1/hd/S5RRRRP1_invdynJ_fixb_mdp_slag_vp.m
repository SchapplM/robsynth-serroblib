% Calculate vector of inverse dynamics joint torques for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:53
% EndTime: 2021-01-15 23:53:04
% DurationCPUTime: 4.72s
% Computational Cost: add. (4804->398), mult. (11358->484), div. (0->0), fcn. (8236->12), ass. (0->201)
t468 = sin(qJ(3));
t471 = cos(qJ(2));
t584 = cos(qJ(3));
t534 = qJD(1) * t584;
t469 = sin(qJ(2));
t549 = qJD(1) * t469;
t398 = -t468 * t549 + t471 * t534;
t463 = qJD(2) + qJD(3);
t598 = t398 * t463;
t537 = t584 * t469;
t564 = t468 * t471;
t409 = t537 + t564;
t596 = t463 * t409;
t481 = t596 * qJD(1);
t511 = t468 * t469 - t584 * t471;
t475 = -t511 * qJDD(1) - t481;
t399 = -qJD(1) * t564 - t469 * t534;
t467 = sin(qJ(4));
t583 = cos(qJ(4));
t366 = t583 * t398 + t399 * t467;
t362 = t366 ^ 2;
t462 = qJDD(2) + qJDD(3);
t455 = qJDD(4) + t462;
t541 = qJDD(1) * t471;
t349 = qJDD(1) * t537 + t468 * t541 + t598;
t508 = -t467 * t398 + t583 * t399;
t496 = t508 * qJD(4) - t467 * t349 + t583 * t475;
t533 = qJD(4) * t583;
t547 = qJD(4) * t467;
t503 = t583 * t349 + t398 * t533 + t399 * t547 + t467 * t475;
t456 = qJD(4) + t463;
t569 = t508 * t456;
t570 = t366 * t456;
t586 = t508 ^ 2;
t597 = MDP(18) * t366 * t508 + t455 * MDP(22) + (-t362 + t586) * MDP(19) + (t496 - t569) * MDP(21) + (t503 - t570) * MDP(20);
t461 = t471 * pkin(2);
t594 = -pkin(1) - t461;
t574 = qJ(5) * t366;
t425 = t594 * qJD(1);
t381 = -pkin(3) * t398 + t425;
t529 = -pkin(4) * t366 + qJD(5);
t331 = t381 + t529;
t571 = t331 * t508;
t359 = t508 * qJ(5);
t385 = t511 * pkin(3) + t594;
t543 = qJD(1) * qJD(2);
t531 = t469 * t543;
t443 = pkin(2) * t531;
t330 = pkin(3) * t481 + t385 * qJDD(1) + t443;
t302 = -pkin(4) * t496 + qJDD(5) + t330;
t470 = sin(qJ(1));
t472 = cos(qJ(1));
t515 = g(1) * t472 + g(2) * t470;
t466 = qJ(2) + qJ(3);
t460 = qJ(4) + t466;
t448 = sin(t460);
t530 = t471 * t543;
t542 = qJDD(1) * t469;
t585 = pkin(6) + pkin(7);
t379 = qJDD(2) * pkin(2) + t585 * (-t530 - t542);
t380 = t585 * (-t531 + t541);
t427 = t585 * t471;
t415 = qJD(1) * t427;
t404 = t584 * t415;
t426 = t585 * t469;
t413 = qJD(1) * t426;
t575 = qJD(2) * pkin(2);
t406 = -t413 + t575;
t510 = -t468 * t406 - t404;
t495 = t510 * qJD(3) + t584 * t379 - t468 * t380;
t307 = t462 * pkin(3) - t349 * pkin(8) + t495;
t532 = t584 * qJD(3);
t548 = qJD(3) * t468;
t491 = t468 * t379 + t584 * t380 + t406 * t532 - t415 * t548;
t311 = t475 * pkin(8) + t491;
t393 = t399 * pkin(8);
t400 = t468 * t415;
t523 = t584 * t406 - t400;
t347 = t393 + t523;
t335 = pkin(3) * t463 + t347;
t576 = t398 * pkin(8);
t348 = -t510 + t576;
t492 = t467 * t307 + t583 * t311 + t335 * t533 - t348 * t547;
t449 = cos(t460);
t561 = t472 * t449;
t567 = t449 * t470;
t489 = g(1) * t561 + g(2) * t567 + g(3) * t448 - t492;
t483 = -t381 * t366 + t489;
t339 = t583 * t348;
t509 = -t467 * t335 - t339;
t593 = t583 * t307 - t467 * t311;
t497 = t509 * qJD(4) + t593;
t563 = t470 * t448;
t568 = t448 * t472;
t578 = g(3) * t449;
t588 = g(1) * t568 + g(2) * t563 - t578;
t487 = t497 + t588;
t479 = t381 * t508 + t487;
t521 = -t584 * t426 - t468 * t427;
t360 = -pkin(8) * t409 + t521;
t552 = -t468 * t426 + t584 * t427;
t361 = -t511 * pkin(8) + t552;
t556 = t467 * t360 + t583 * t361;
t451 = t584 * pkin(2) + pkin(3);
t565 = t467 * t468;
t370 = t451 * t533 + (-t468 * t547 + (t583 * t584 - t565) * qJD(3)) * pkin(2);
t536 = t583 * t468;
t395 = pkin(2) * t536 + t467 * t451;
t592 = -t370 * t456 - t395 * t455;
t581 = pkin(3) * t456;
t591 = -t467 * pkin(3) * t455 - t533 * t581;
t444 = t455 * pkin(4);
t590 = qJD(5) * t508 + t444;
t573 = t503 * qJ(5);
t587 = -t573 + t590;
t582 = pkin(3) * t399;
t572 = t496 * qJ(5);
t337 = t467 * t348;
t458 = cos(t466);
t562 = t470 * t458;
t560 = t472 * t458;
t528 = t583 * t335 - t337;
t309 = t528 + t359;
t305 = pkin(4) * t456 + t309;
t559 = t305 - t309;
t558 = t583 * t347 - t337;
t522 = t413 * t468 - t404;
t351 = t522 - t576;
t553 = -t584 * t413 - t400;
t352 = t393 + t553;
t557 = t467 * t351 + t583 * t352;
t317 = t359 + t557;
t555 = t370 - t317;
t525 = t583 * t351 - t352 * t467;
t316 = t525 - t574;
t371 = -t451 * t547 + (-t468 * t533 + (-t584 * t467 - t536) * qJD(3)) * pkin(2);
t554 = t371 - t316;
t551 = pkin(3) * t458 + pkin(4) * t449;
t464 = t469 ^ 2;
t550 = -t471 ^ 2 + t464;
t546 = t366 * qJD(5);
t454 = t469 * t575;
t540 = t461 + t551;
t538 = qJD(2) * t585;
t527 = -t347 * t467 - t339;
t524 = t583 * t360 - t361 * t467;
t519 = -pkin(2) * t565 + t583 * t451;
t518 = -g(1) * t563 + g(2) * t568;
t517 = g(1) * t567 - g(2) * t561;
t334 = -pkin(4) * t508 - t582;
t457 = sin(t466);
t516 = -pkin(3) * t457 - pkin(4) * t448;
t514 = g(1) * t470 - g(2) * t472;
t310 = -t509 + t574;
t513 = t305 * t366 - t310 * t508;
t512 = -0.2e1 * pkin(1) * t543 - pkin(6) * qJDD(2);
t414 = t469 * t538;
t416 = t471 * t538;
t505 = -t584 * t414 - t468 * t416 - t426 * t532 - t427 * t548;
t326 = -pkin(8) * t596 + t505;
t378 = t463 * t511;
t493 = -t552 * qJD(3) + t468 * t414 - t584 * t416;
t327 = t378 * pkin(8) + t493;
t504 = t583 * t326 + t467 * t327 + t360 * t533 - t361 * t547;
t500 = t583 * t511;
t473 = qJD(2) ^ 2;
t499 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t473 + t514;
t474 = qJD(1) ^ 2;
t498 = pkin(1) * t474 - pkin(6) * qJDD(1) + t515;
t377 = t583 * t409 - t467 * t511;
t494 = -t556 * qJD(4) - t467 * t326 + t583 * t327;
t488 = t399 * t398 * MDP(11) + (t349 - t598) * MDP(13) + (-t399 * t463 + t475) * MDP(14) + (-t398 ^ 2 + t399 ^ 2) * MDP(12) + t462 * MDP(15) + t597;
t484 = t489 - t572;
t482 = g(1) * t560 + g(2) * t562 + g(3) * t457 - t425 * t398 - t491;
t480 = t487 - t573;
t478 = -g(3) * t458 + t425 * t399 + t457 * t515 + t495;
t369 = pkin(3) * t596 + t454;
t476 = -t331 * t366 + t484 - t546;
t459 = -qJ(5) - pkin(8) - t585;
t453 = pkin(2) * t549;
t450 = t583 * pkin(3) + pkin(4);
t394 = qJDD(1) * t594 + t443;
t391 = pkin(4) + t519;
t386 = pkin(1) + t540;
t382 = t453 - t582;
t376 = t409 * t467 + t500;
t358 = t371 * t456;
t340 = t376 * pkin(4) + t385;
t332 = t334 + t453;
t322 = t377 * qJD(4) - t467 * t378 + t583 * t596;
t321 = qJD(4) * t500 + t583 * t378 + t409 * t547 + t467 * t596;
t320 = -qJ(5) * t376 + t556;
t319 = -qJ(5) * t377 + t524;
t318 = t322 * pkin(4) + t369;
t313 = t359 + t558;
t312 = t527 - t574;
t301 = t321 * qJ(5) - t377 * qJD(5) + t494;
t300 = -qJ(5) * t322 - qJD(5) * t376 + t504;
t299 = t492 + t546 + t572;
t298 = t497 + t587;
t1 = [(t349 * t409 + t378 * t399) * MDP(11) + (-t381 * t321 + t330 * t377 - t369 * t508 + t385 * t503 - t556 * t455 - t504 * t456 + t518) * MDP(24) + (-t300 * t456 + t302 * t377 - t318 * t508 - t320 * t455 - t321 * t331 + t340 * t503 + t518) * MDP(26) + (t321 * t508 + t377 * t503) * MDP(18) + (t349 * t594 - t425 * t378 + t394 * t409 - t399 * t454 - t514 * t457 - t552 * t462 - t505 * t463) * MDP(17) + (-t378 * t463 + t409 * t462) * MDP(13) + (g(1) * t562 - g(2) * t560 + t394 * t511 - t398 * t454 + t425 * t596 + t521 * t462 + t493 * t463 - t475 * t594) * MDP(16) + 0.2e1 * (t469 * t541 - t550 * t543) * MDP(5) + (qJDD(1) * t464 + 0.2e1 * t469 * t530) * MDP(4) + (t301 * t456 + t302 * t376 - t318 * t366 + t319 * t455 + t322 * t331 - t340 * t496 + t517) * MDP(25) + (t381 * t322 + t330 * t376 - t366 * t369 - t385 * t496 + t524 * t455 + t494 * t456 + t517) * MDP(23) + (-t321 * t366 + t322 * t508 - t376 * t503 + t377 * t496) * MDP(19) + (-t298 * t377 - t299 * t376 + t300 * t366 + t301 * t508 + t305 * t321 - t310 * t322 - t319 * t503 + t320 * t496 - t515) * MDP(27) + (-t511 * t462 - t463 * t596) * MDP(14) + (-t349 * t511 - t378 * t398 + t399 * t596 + t475 * t409) * MDP(12) + (t512 * t469 + t499 * t471) * MDP(9) + (-t499 * t469 + t512 * t471) * MDP(10) + t515 * MDP(3) + qJDD(1) * MDP(1) + (t299 * t320 + t310 * t300 + t298 * t319 + t305 * t301 + t302 * t340 + t331 * t318 - g(1) * (-t386 * t470 - t459 * t472) - g(2) * (t386 * t472 - t459 * t470)) * MDP(28) + (qJDD(2) * t469 + t471 * t473) * MDP(6) + (qJDD(2) * t471 - t469 * t473) * MDP(7) + (-t321 * t456 + t377 * t455) * MDP(20) + (-t322 * t456 - t376 * t455) * MDP(21) + t514 * MDP(2); (t553 * t463 + (t399 * t549 - t468 * t462 - t463 * t532) * pkin(2) + t482) * MDP(17) + (t382 * t508 + t557 * t456 + t483 + t592) * MDP(24) + (t317 * t456 + t332 * t508 + t476 + t592) * MDP(26) + (-t522 * t463 + (t398 * t549 + t584 * t462 - t463 * t548) * pkin(2) + t478) * MDP(16) + (t366 * t382 + t519 * t455 - t525 * t456 + t358 + t479) * MDP(23) + (-t316 * t456 + t332 * t366 + t391 * t455 + t358 + t480 + t571 + t590) * MDP(25) + t488 + (g(3) * t469 + t498 * t471) * MDP(10) + (-g(3) * t471 + t498 * t469) * MDP(9) + qJDD(2) * MDP(8) + MDP(6) * t542 + MDP(7) * t541 + (t366 * t555 - t391 * t503 + t395 * t496 + t508 * t554 + t513) * MDP(27) + (t299 * t395 + t298 * t391 - t331 * t332 - g(3) * t540 - t515 * (-pkin(2) * t469 + t516) + t555 * t310 + t554 * t305) * MDP(28) + (-t469 * t471 * MDP(4) + t550 * MDP(5)) * t474; (t523 * t463 + t482) * MDP(17) + (t558 * t456 - t508 * t582 + t483 + t591) * MDP(24) + (t313 * t456 + t334 * t508 + t476 + t591) * MDP(26) + (-t510 * t463 + t478) * MDP(16) + (-t527 * t456 + (-t366 * t399 + t583 * t455 - t456 * t547) * pkin(3) + t479) * MDP(23) + t488 + (-t312 * t508 - t313 * t366 - t450 * t503 + (t496 * t467 + (t366 * t583 - t467 * t508) * qJD(4)) * pkin(3) + t513) * MDP(27) + (-t312 * t456 + t571 + t334 * t366 + t450 * t455 + (-t339 + (-t335 - t581) * t467) * qJD(4) + t587 + t588 + t593) * MDP(25) + (t298 * t450 - t310 * t313 - t305 * t312 - t331 * t334 - g(3) * t551 - t515 * t516 + (t299 * t467 + (-t305 * t467 + t583 * t310) * qJD(4)) * pkin(3)) * MDP(28); (-t509 * t456 + t479) * MDP(23) + (t528 * t456 + t483) * MDP(24) + (t310 * t456 + 0.2e1 * t444 - (-t331 - t529) * t508 + t480) * MDP(25) + (-t586 * pkin(4) + t309 * t456 - (qJD(5) + t331) * t366 + t484) * MDP(26) + (-pkin(4) * t503 + t366 * t559) * MDP(27) + (t559 * t310 + (t515 * t448 + t298 + t571 - t578) * pkin(4)) * MDP(28) + t597; (-t496 - t569) * MDP(25) + (t503 + t570) * MDP(26) + (-t362 - t586) * MDP(27) + (-t305 * t508 - t310 * t366 + t302 - t514) * MDP(28);];
tau = t1;
