% Calculate vector of inverse dynamics joint torques for
% S5RRRRP6
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
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:36
% EndTime: 2021-01-16 00:10:51
% DurationCPUTime: 6.03s
% Computational Cost: add. (5160->473), mult. (11407->583), div. (0->0), fcn. (7962->10), ass. (0->215)
t475 = cos(qJ(2));
t477 = pkin(7) + pkin(6);
t438 = t477 * t475;
t427 = qJD(1) * t438;
t471 = sin(qJ(3));
t418 = t471 * t427;
t472 = sin(qJ(2));
t437 = t477 * t472;
t425 = qJD(1) * t437;
t602 = cos(qJ(3));
t385 = -t602 * t425 - t418;
t541 = t602 * qJD(3);
t610 = -pkin(2) * t541 + t385;
t468 = qJ(2) + qJ(3);
t463 = sin(t468);
t476 = cos(qJ(1));
t578 = t463 * t476;
t473 = sin(qJ(1));
t579 = t463 * t473;
t608 = g(1) * t578 + g(2) * t579;
t470 = sin(qJ(4));
t474 = cos(qJ(4));
t505 = -t471 * t472 + t602 * t475;
t553 = qJD(2) + qJD(3);
t388 = t553 * t505;
t483 = t388 * qJD(1);
t552 = qJDD(2) + qJDD(3);
t612 = t470 * t483 - t474 * t552;
t542 = qJD(1) * t602;
t576 = t471 * t475;
t417 = -qJD(1) * t576 - t472 * t542;
t424 = t602 * t472 + t576;
t496 = t424 * qJDD(1);
t481 = t496 + t483;
t529 = qJD(4) * t553;
t558 = qJD(4) * t470;
t500 = t417 * t558 + t470 * t552 + (t481 + t529) * t474;
t395 = -t417 * t470 - t474 * t553;
t560 = qJD(1) * t472;
t415 = t471 * t560 - t475 * t542;
t411 = qJD(4) + t415;
t586 = t395 * t411;
t611 = t500 - t586;
t499 = t474 * t417 - t470 * t553;
t609 = t499 * qJD(4);
t419 = t602 * t427;
t384 = -t471 * t425 + t419;
t559 = qJD(3) * t471;
t525 = pkin(2) * t559 - t384;
t593 = t474 * pkin(4);
t459 = pkin(3) + t593;
t464 = cos(t468);
t469 = -qJ(5) - pkin(8);
t532 = t464 * t459 - t463 * t469;
t592 = t475 * pkin(2);
t461 = pkin(1) + t592;
t607 = -pkin(8) * t424 - t461;
t584 = t415 * t470;
t606 = -qJ(5) * t584 + t474 * qJD(5);
t378 = -pkin(3) * t417 + pkin(8) * t415;
t369 = pkin(2) * t560 + t378;
t605 = t470 * t369 + t610 * t474;
t521 = g(1) * t476 + g(2) * t473;
t575 = t473 * t464;
t404 = t470 * t575 + t474 * t476;
t573 = t476 * t464;
t574 = t473 * t474;
t406 = -t470 * t573 + t574;
t594 = g(3) * t470;
t604 = -g(1) * t406 + g(2) * t404 + t463 * t594;
t339 = t470 * t496 - t609 + t612;
t603 = t499 ^ 2;
t389 = t553 * t424;
t534 = t602 * qJDD(1);
t555 = qJDD(1) * t472;
t517 = t471 * t555 - t475 * t534;
t357 = t389 * qJD(1) + t517;
t356 = qJDD(4) + t357;
t601 = pkin(4) * t356;
t455 = g(3) * t463;
t595 = g(3) * t464;
t591 = qJD(2) * pkin(2);
t590 = qJ(5) * t500;
t589 = qJ(5) * t339;
t588 = t500 * t470;
t420 = -t425 + t591;
t381 = t602 * t420 - t418;
t370 = -t553 * pkin(3) - t381;
t587 = t370 * t415;
t585 = t499 * t411;
t583 = t415 * t474;
t582 = t424 * t470;
t581 = t424 * t474;
t577 = t470 * t356;
t465 = t474 * qJ(5);
t400 = -t471 * t437 + t602 * t438;
t392 = t474 * t400;
t458 = pkin(2) * t471 + pkin(8);
t572 = -qJ(5) - t458;
t436 = t461 * qJD(1);
t367 = pkin(3) * t415 + pkin(8) * t417 - t436;
t382 = t471 * t420 + t419;
t371 = t553 * pkin(8) + t382;
t340 = t474 * t367 - t371 * t470;
t326 = qJ(5) * t499 + t340;
t322 = pkin(4) * t411 + t326;
t571 = -t326 + t322;
t570 = t470 * t378 + t474 * t381;
t531 = qJD(4) * t572;
t568 = t470 * t531 - t605 + t606;
t366 = t474 * t369;
t519 = -t417 * pkin(4) + t415 * t465;
t567 = t474 * t531 - t366 - t519 + (-qJD(5) + t610) * t470;
t380 = -pkin(3) * t505 + t607;
t566 = t470 * t380 + t392;
t535 = qJD(4) * t469;
t565 = t470 * t535 - t570 + t606;
t374 = t474 * t378;
t564 = t474 * t535 - t374 - t519 + (-qJD(5) + t381) * t470;
t401 = pkin(4) * t584;
t549 = pkin(4) * t558;
t563 = t549 + t401 + t525;
t562 = t608 * t474;
t466 = t472 ^ 2;
t561 = -t475 ^ 2 + t466;
t557 = qJD(4) * t474;
t556 = qJD(1) * qJD(2);
t554 = qJDD(1) * t475;
t551 = t472 * t591;
t548 = qJD(4) * pkin(8) * t411;
t349 = pkin(3) * t389 - pkin(8) * t388 + t551;
t544 = qJD(2) * t477;
t426 = t472 * t544;
t428 = t475 * t544;
t506 = -t602 * t437 - t471 * t438;
t352 = t506 * qJD(3) - t602 * t426 - t471 * t428;
t546 = t470 * t349 + t474 * t352 + t380 * t557;
t545 = g(1) * t573 + g(2) * t575 + t455;
t543 = t424 * t557;
t363 = t370 * t557;
t540 = t472 * t556;
t539 = t475 * t556;
t538 = pkin(4) * t470 + t477;
t391 = qJDD(2) * pkin(2) + t477 * (-t539 - t555);
t394 = t477 * (-t540 + t554);
t526 = -t602 * t391 + t471 * t394 + t420 * t559 + t427 * t541;
t334 = -t552 * pkin(3) + t526;
t316 = t339 * pkin(4) + qJDD(5) + t334;
t537 = -t316 - t595;
t536 = -t334 - t595;
t533 = t395 * pkin(4) + qJD(5);
t530 = t411 * t474;
t453 = pkin(2) * t540;
t325 = t357 * pkin(3) - pkin(8) * t483 + t607 * qJDD(1) + t453;
t489 = t471 * t391 + t602 * t394 + t420 * t541 - t427 * t559;
t333 = t552 * pkin(8) + t489;
t527 = t470 * t325 + t474 * t333 + t367 * t557 - t371 * t558;
t460 = -t602 * pkin(2) - pkin(3);
t355 = -t401 + t382;
t524 = -t355 + t549;
t523 = -g(1) * t404 - g(2) * t406;
t405 = -t464 * t574 + t470 * t476;
t407 = t470 * t473 + t474 * t573;
t522 = -g(1) * t405 - g(2) * t407;
t520 = g(1) * t473 - g(2) * t476;
t518 = -pkin(8) * t356 + t587;
t341 = t367 * t470 + t371 * t474;
t327 = -qJ(5) * t395 + t341;
t515 = -t322 * t474 - t327 * t470;
t514 = -t356 * t458 + t587;
t513 = t459 * t463 + t464 * t469;
t512 = -qJ(5) * t388 - qJD(5) * t424;
t441 = t464 * t594;
t511 = t334 * t470 - t341 * t417 + t363 + t441;
t510 = t340 * t417 + t370 * t558 + t562;
t509 = t521 * t463;
t508 = -0.2e1 * pkin(1) * t556 - pkin(6) * qJDD(2);
t507 = t461 + t532;
t504 = t388 * t470 + t543;
t503 = t388 * t474 - t424 * t558;
t498 = t470 * t509;
t350 = t370 + t533;
t495 = t316 * t470 - t327 * t417 + t441 + (t557 + t583) * t350;
t478 = qJD(2) ^ 2;
t494 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t478 + t520;
t479 = qJD(1) ^ 2;
t493 = pkin(1) * t479 - pkin(6) * qJDD(1) + t521;
t492 = g(1) * t407 - g(2) * t405 + t474 * t455 - t527;
t491 = -t436 * t417 - t526 - t595 + t608;
t324 = t474 * t325;
t490 = -t341 * qJD(4) - t333 * t470 + t324;
t488 = t322 * t417 + t537 * t474 + t562 + (t558 + t584) * t350;
t353 = t400 * qJD(3) - t471 * t426 + t602 * t428;
t487 = (t611 * t474 + (-t339 + t585) * t470) * MDP(19) + (-t499 * t530 + t588) * MDP(18) + (-t411 ^ 2 * t470 + t356 * t474 - t395 * t417) * MDP(21) + (t411 * t530 - t417 * t499 + t577) * MDP(20) + (t415 * t553 + t481) * MDP(13) + (-t517 + (-qJD(1) * t424 - t417) * t553) * MDP(14) + (-t415 ^ 2 + t417 ^ 2) * MDP(12) + t552 * MDP(15) + (-MDP(11) * t415 + t411 * MDP(22)) * t417;
t309 = qJD(5) * t499 + t490 - t590 + t601;
t311 = -qJD(5) * t395 + t527 - t589;
t486 = t515 * qJD(4) - t309 * t470 + t311 * t474 - t322 * t583 - t327 * t584 - t545;
t485 = t490 + t604;
t484 = -t436 * t415 - t489 + t545;
t435 = pkin(8) * t474 + t465;
t434 = t469 * t470;
t433 = t460 - t593;
t422 = t458 * t474 + t465;
t421 = t572 * t470;
t412 = -t461 * qJDD(1) + t453;
t393 = t395 ^ 2;
t376 = t474 * t380;
t368 = pkin(4) * t582 - t506;
t346 = t474 * t349;
t342 = -qJ(5) * t582 + t566;
t336 = -pkin(4) * t505 - t400 * t470 - t424 * t465 + t376;
t330 = t504 * pkin(4) + t353;
t314 = -qJ(5) * t543 + (-qJD(4) * t400 + t512) * t470 + t546;
t313 = pkin(4) * t389 - t352 * t470 + t346 + t512 * t474 + (-t392 + (qJ(5) * t424 - t380) * t470) * qJD(4);
t1 = [(-t417 * t388 + t481 * t424) * MDP(11) + (t309 * t336 + t311 * t342 + t322 * t313 + t327 * t314 + t316 * t368 + t350 * t330 + (-g(1) * t538 - g(2) * t507) * t476 + (g(1) * t507 - g(2) * t538) * t473) * MDP(28) + (qJDD(1) * t466 + 0.2e1 * t472 * t539) * MDP(4) + (-g(1) * t579 + g(2) * t578 - t352 * t553 - t436 * t388 - t400 * t552 + t412 * t424 - t417 * t551 - t461 * t481) * MDP(17) + (t508 * t472 + t494 * t475) * MDP(9) + (-t494 * t472 + t508 * t475) * MDP(10) + t520 * MDP(2) + t521 * MDP(3) + (qJDD(2) * t472 + t475 * t478) * MDP(6) + (qJDD(2) * t475 - t472 * t478) * MDP(7) + (t388 * t553 + t424 * t552) * MDP(13) + ((-t395 * t474 + t470 * t499) * t388 + (-t588 - t339 * t474 + (t395 * t470 + t474 * t499) * qJD(4)) * t424) * MDP(19) + (-t499 * t503 + t500 * t581) * MDP(18) + (t313 * t499 - t314 * t395 - t336 * t500 - t339 * t342 + t520 * t463 + t515 * t388 + (-t309 * t474 - t311 * t470 + (t322 * t470 - t327 * t474) * qJD(4)) * t424) * MDP(27) + (-t424 * t357 - t388 * t415 + t417 * t389 + t481 * t505) * MDP(12) + (-t309 * t505 + t313 * t411 + t316 * t582 + t322 * t389 + t330 * t395 + t336 * t356 + t339 * t368 + t504 * t350 + t522) * MDP(25) + (-t356 * t505 + t389 * t411) * MDP(22) + (t339 * t505 - t389 * t395 - t504 * t411 - t424 * t577) * MDP(21) + (-t389 * t553 + t505 * t552) * MDP(14) + (t356 * t581 - t389 * t499 + t503 * t411 - t500 * t505) * MDP(20) + (t311 * t505 - t314 * t411 + t316 * t581 - t327 * t389 - t330 * t499 - t342 * t356 + t503 * t350 + t368 * t500 + t523) * MDP(26) + ((-t400 * t557 + t346) * t411 + t376 * t356 - (-t371 * t557 + t324) * t505 + t340 * t389 + t353 * t395 - t506 * t339 + t424 * t363 + ((-qJD(4) * t380 - t352) * t411 - t400 * t356 - (-qJD(4) * t367 - t333) * t505 + t334 * t424 + t370 * t388) * t470 + t522) * MDP(23) + (-t353 * t553 - t461 * t357 - t436 * t389 - t412 * t505 + t415 * t551 + t520 * t464 + t506 * t552) * MDP(16) + (-(-t400 * t558 + t546) * t411 - t566 * t356 + t527 * t505 - t341 * t389 - t353 * t499 - t506 * t500 + t334 * t581 + t503 * t370 + t523) * MDP(24) + 0.2e1 * (t472 * t554 - t561 * t556) * MDP(5) + qJDD(1) * MDP(1); (t460 * t339 + t536 * t474 + t514 * t470 + t525 * t395 + (-t458 * t557 + t610 * t470 - t366) * t411 + t510) * MDP(23) + (t460 * t500 + t514 * t474 - t498 - t525 * t499 + (t458 * t558 + t605) * t411 + t511) * MDP(24) + (t384 * t553 + (-t415 * t560 + t602 * t552 - t553 * t559) * pkin(2) + t491) * MDP(16) + (t311 * t422 + t309 * t421 + t316 * t433 - g(3) * (t532 + t592) + t563 * t350 + t568 * t327 + t567 * t322 + t521 * (pkin(2) * t472 + t513)) * MDP(28) + t487 + (-t356 * t422 - t568 * t411 + t433 * t500 - t499 * t563 + t495 - t498) * MDP(26) + qJDD(2) * MDP(8) + MDP(7) * t554 + (-g(3) * t475 + t493 * t472) * MDP(9) + (g(3) * t472 + t493 * t475) * MDP(10) + MDP(6) * t555 + (t339 * t433 + t356 * t421 + t563 * t395 + t567 * t411 + t488) * MDP(25) + (t385 * t553 + (t417 * t560 - t471 * t552 - t553 * t541) * pkin(2) + t484) * MDP(17) + (-t339 * t422 - t568 * t395 - t421 * t500 + t499 * t567 + t486) * MDP(27) + (-t472 * t475 * MDP(4) + t561 * MDP(5)) * t479; (-pkin(3) * t339 - t374 * t411 - t382 * t395 + (t381 * t411 + t518) * t470 + (t536 - t548) * t474 + t510) * MDP(23) + (-pkin(3) * t500 + t570 * t411 + t382 * t499 + t518 * t474 + (-t509 + t548) * t470 + t511) * MDP(24) + (-g(3) * t532 + t309 * t434 + t311 * t435 - t316 * t459 + t564 * t322 + t565 * t327 + t524 * t350 + t521 * t513) * MDP(28) + t487 + (-t500 * t459 + t355 * t499 - t356 * t435 - t565 * t411 + (-pkin(4) * t609 - t509) * t470 + t495) * MDP(26) + (t382 * t553 + t491) * MDP(16) + (-t339 * t459 + t356 * t434 + t524 * t395 + t564 * t411 + t488) * MDP(25) + (t381 * t553 + t484) * MDP(17) + (-t339 * t435 - t565 * t395 - t434 * t500 + t499 * t564 + t486) * MDP(27); -t499 * t395 * MDP(18) + (-t393 + t603) * MDP(19) + (t500 + t586) * MDP(20) + (-t339 - t585) * MDP(21) + t356 * MDP(22) + (t341 * t411 + t370 * t499 + t485) * MDP(23) + (t340 * t411 + t370 * t395 + t492) * MDP(24) + (0.2e1 * t601 - t590 + t327 * t411 - (-t350 - t533) * t499 + t485) * MDP(25) + (-pkin(4) * t603 + t589 + t326 * t411 + (qJD(5) + t350) * t395 + t492) * MDP(26) + (-pkin(4) * t500 - t571 * t395) * MDP(27) + (t571 * t327 + (t350 * t499 + t309 + t604) * pkin(4)) * MDP(28); t611 * MDP(26) + (-t393 - t603) * MDP(27) + (-t322 * t499 + t327 * t395 - t537 - t608) * MDP(28) + (-t417 * t557 - t585 + (t471 * t554 + t472 * t534 + t529) * t470 + t612) * MDP(25);];
tau = t1;
