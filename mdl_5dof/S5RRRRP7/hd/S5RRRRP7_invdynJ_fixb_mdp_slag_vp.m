% Calculate vector of inverse dynamics joint torques for
% S5RRRRP7
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
%   see S5RRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:57
% EndTime: 2019-12-31 21:58:08
% DurationCPUTime: 6.71s
% Computational Cost: add. (5245->472), mult. (11621->590), div. (0->0), fcn. (8081->10), ass. (0->206)
t467 = sin(qJ(2));
t466 = sin(qJ(3));
t470 = cos(qJ(2));
t563 = t466 * t470;
t599 = cos(qJ(3));
t417 = t599 * t467 + t563;
t543 = qJD(2) + qJD(3);
t384 = t543 * t417;
t533 = t599 * t470;
t545 = qJDD(1) * t467;
t359 = qJD(1) * t384 - qJDD(1) * t533 + t466 * t545;
t358 = qJDD(4) + t359;
t468 = sin(qJ(1));
t471 = cos(qJ(1));
t518 = g(1) * t471 + g(2) * t468;
t469 = cos(qJ(4));
t548 = qJD(4) * t469;
t532 = qJD(1) * t599;
t551 = qJD(1) * t467;
t409 = t466 * t551 - t470 * t532;
t572 = t409 * t469;
t616 = t548 + t572;
t465 = sin(qJ(4));
t549 = qJD(4) * t465;
t615 = t409 * t465 + t549;
t464 = qJ(2) + qJ(3);
t460 = sin(t464);
t569 = t460 * t471;
t570 = t460 * t468;
t612 = g(1) * t569 + g(2) * t570;
t515 = t469 * pkin(4) + t465 * qJ(5);
t461 = cos(t464);
t611 = t518 * t461;
t600 = pkin(7) + pkin(6);
t405 = qJD(4) + t409;
t561 = t469 * t358;
t608 = (t405 * t549 - t561) * pkin(8);
t546 = qJD(1) * qJD(2);
t530 = t467 * t546;
t447 = pkin(2) * t530;
t504 = -t466 * t467 + t533;
t383 = t543 * t504;
t478 = t383 * qJD(1);
t458 = t470 * pkin(2) + pkin(1);
t604 = -pkin(8) * t417 - t458;
t324 = t359 * pkin(3) - pkin(8) * t478 + t604 * qJDD(1) + t447;
t529 = t470 * t546;
t386 = qJDD(2) * pkin(2) + t600 * (-t529 - t545);
t544 = qJDD(1) * t470;
t388 = t600 * (-t530 + t544);
t428 = t600 * t467;
t418 = qJD(1) * t428;
t587 = qJD(2) * pkin(2);
t414 = -t418 + t587;
t429 = t600 * t470;
t420 = qJD(1) * t429;
t531 = qJD(3) * t599;
t550 = qJD(3) * t466;
t487 = t466 * t386 + t599 * t388 + t414 * t531 - t420 * t550;
t542 = qJDD(2) + qJDD(3);
t332 = t542 * pkin(8) + t487;
t411 = -qJD(1) * t563 - t467 * t532;
t427 = t458 * qJD(1);
t367 = pkin(3) * t409 + pkin(8) * t411 - t427;
t413 = t599 * t420;
t379 = t466 * t414 + t413;
t370 = t543 * pkin(8) + t379;
t498 = t465 * t324 + t469 * t332 + t367 * t548 - t370 * t549;
t586 = qJ(5) * t358;
t309 = qJD(5) * t405 + t498 + t586;
t523 = -t469 * t324 + t465 * t332 + t367 * t549 + t370 * t548;
t598 = pkin(4) * t358;
t311 = qJDD(5) + t523 - t598;
t607 = t309 * t469 + t311 * t465;
t381 = -t466 * t418 + t413;
t521 = pkin(2) * t550 - t381;
t377 = -pkin(3) * t504 + t604;
t394 = -t466 * t428 + t599 * t429;
t555 = t465 * t377 + t469 * t394;
t606 = t615 * pkin(4) - t616 * qJ(5) - qJD(5) * t465;
t605 = t461 * pkin(3) + t460 * pkin(8);
t477 = t417 * qJDD(1) + t478;
t496 = t469 * t411 - t465 * t543;
t341 = -qJD(4) * t496 + t477 * t465 - t469 * t542;
t602 = t496 ^ 2;
t601 = t405 ^ 2;
t597 = pkin(4) * t411;
t596 = pkin(8) * t358;
t451 = g(3) * t460;
t592 = g(3) * t461;
t591 = g(3) * t465;
t590 = g(3) * t470;
t588 = pkin(8) * qJD(4);
t522 = -t599 * t386 + t466 * t388 + t414 * t550 + t420 * t531;
t333 = -t542 * pkin(3) + t522;
t525 = t469 * t543;
t340 = -qJD(4) * t525 - t411 * t549 - t465 * t542 - t469 * t477;
t313 = t341 * pkin(4) + t340 * qJ(5) + qJD(5) * t496 + t333;
t585 = t313 * t465;
t412 = t466 * t420;
t378 = t599 * t414 - t412;
t369 = -t543 * pkin(3) - t378;
t389 = -t411 * t465 - t525;
t337 = t389 * pkin(4) + qJ(5) * t496 + t369;
t584 = t337 * t409;
t583 = t340 * t465;
t582 = t341 * t469;
t343 = t367 * t465 + t370 * t469;
t581 = t343 * t405;
t456 = pkin(2) * t466 + pkin(8);
t580 = t358 * t456;
t579 = t369 * t409;
t578 = t389 * t405;
t577 = t389 * t465;
t576 = t496 * t389;
t575 = t496 * t405;
t574 = t496 * t469;
t571 = t417 * t469;
t565 = t465 * t358;
t564 = t465 * t468;
t562 = t468 * t469;
t560 = t469 * t471;
t559 = t471 * t465;
t558 = -t606 - t521;
t375 = -pkin(3) * t411 + pkin(8) * t409;
t557 = t465 * t375 + t469 * t378;
t368 = pkin(2) * t551 + t375;
t382 = -t599 * t418 - t412;
t556 = t465 * t368 + t469 * t382;
t554 = -t379 + t606;
t553 = t612 * t469;
t462 = t467 ^ 2;
t552 = -t470 ^ 2 + t462;
t342 = t367 * t469 - t370 * t465;
t547 = qJD(5) - t342;
t541 = t599 * pkin(2);
t540 = t467 * t587;
t538 = t461 * t591 - t612 * t465;
t537 = t451 + t611;
t536 = qJD(2) * t600;
t535 = t465 * t599;
t534 = t469 * t599;
t528 = -t313 - t592;
t527 = -t333 - t592;
t526 = t405 * t469;
t524 = pkin(2) * t531;
t400 = t461 * t564 + t560;
t402 = t461 * t559 - t562;
t520 = -g(1) * t400 + g(2) * t402;
t401 = t461 * t562 - t559;
t403 = t461 * t560 + t564;
t519 = g(1) * t401 - g(2) * t403;
t517 = g(1) * t468 - g(2) * t471;
t514 = pkin(4) * t465 - qJ(5) * t469;
t327 = -pkin(4) * t405 + t547;
t328 = qJ(5) * t405 + t343;
t512 = t327 * t469 - t328 * t465;
t511 = t579 - t580;
t510 = pkin(3) + t515;
t509 = t458 + t605;
t508 = -t327 * t411 + t337 * t549 + t553;
t507 = t342 * t411 + t369 * t549 + t553;
t506 = -0.2e1 * pkin(1) * t546 - pkin(6) * qJDD(2);
t505 = -t599 * t428 - t466 * t429;
t502 = t383 * t465 + t417 * t548;
t501 = -t383 * t469 + t417 * t549;
t351 = pkin(3) * t384 - pkin(8) * t383 + t540;
t419 = t467 * t536;
t421 = t470 * t536;
t354 = t505 * qJD(3) - t599 * t419 - t466 * t421;
t497 = t465 * t351 + t469 * t354 + t377 * t548 - t394 * t549;
t495 = t333 * t465 - t343 * t411 + t369 * t548 + t538;
t493 = t328 * t411 - t337 * t572 - t538 - t585;
t492 = -t456 * t549 + t469 * t524;
t473 = qJD(2) ^ 2;
t491 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t473 + t517;
t474 = qJD(1) ^ 2;
t490 = pkin(1) * t474 - pkin(6) * qJDD(1) + t518;
t489 = g(1) * t402 + g(2) * t400 + t460 * t591 - t523;
t488 = -t427 * t411 - t522 - t592 + t612;
t486 = t512 * qJD(4) + t607;
t485 = -t583 - t582 + (-t574 + t577) * qJD(4);
t484 = -t337 * t496 + qJDD(5) - t489;
t483 = t616 * t327 - t615 * t328 - t537 + t607;
t482 = -g(1) * t403 - g(2) * t401 - t469 * t451 + t498;
t355 = t394 * qJD(3) - t466 * t419 + t599 * t421;
t481 = ((-t340 - t578) * t469 + (-t341 + t575) * t465) * MDP(19) + (-t496 * t526 - t583) * MDP(18) + (-t389 * t411 - t601 * t465 + t561) * MDP(21) + (t405 * t526 - t411 * t496 + t565) * MDP(20) + (t409 * t543 + t477) * MDP(13) + (-t411 * t543 - t359) * MDP(14) + (-t409 ^ 2 + t411 ^ 2) * MDP(12) + t542 * MDP(15) + (-MDP(11) * t409 + t405 * MDP(22)) * t411;
t480 = -t427 * t409 - t487 + t537;
t476 = -g(3) * (t515 * t461 + t605) + t518 * t460 * t510 - t611 * pkin(8);
t457 = -t541 - pkin(3);
t415 = -t541 - t510;
t406 = -t458 * qJDD(1) + t447;
t399 = t411 * qJ(5);
t356 = -pkin(4) * t496 + qJ(5) * t389;
t353 = t514 * t417 - t505;
t345 = pkin(4) * t504 - t377 * t469 + t394 * t465;
t344 = -qJ(5) * t504 + t555;
t339 = -t375 * t469 + t378 * t465 + t597;
t338 = -t399 + t557;
t336 = -t368 * t469 + t382 * t465 + t597;
t335 = -t399 + t556;
t320 = -t340 + t578;
t316 = t514 * t383 + (t515 * qJD(4) - qJD(5) * t469) * t417 + t355;
t315 = -pkin(4) * t384 + t555 * qJD(4) - t351 * t469 + t354 * t465;
t314 = qJ(5) * t384 - qJD(5) * t504 + t497;
t1 = [(qJDD(2) * t467 + t470 * t473) * MDP(6) + (qJDD(2) * t470 - t467 * t473) * MDP(7) + (-g(1) * t570 + g(2) * t569 - t354 * t543 - t427 * t383 - t394 * t542 + t406 * t417 - t411 * t540 - t458 * t477) * MDP(17) + ((-t389 * t469 + t465 * t496) * t383 + (t583 - t582 + (t574 + t577) * qJD(4)) * t417) * MDP(19) + (-t314 * t389 - t315 * t496 - t340 * t345 - t341 * t344 + t517 * t460 + t512 * t383 + (-t309 * t465 + t311 * t469 + (-t327 * t465 - t328 * t469) * qJD(4)) * t417) * MDP(26) + (-t340 * t571 + t496 * t501) * MDP(18) + (-t358 * t504 + t384 * t405) * MDP(22) + (-t384 * t543 + t504 * t542) * MDP(14) + (t341 * t504 - t384 * t389 - t405 * t502 - t417 * t565) * MDP(21) + (t311 * t504 - t315 * t405 + t316 * t389 - t327 * t384 + t337 * t502 + t341 * t353 - t345 * t358 + t417 * t585 + t519) * MDP(25) + (t340 * t504 - t384 * t496 - t405 * t501 + t417 * t561) * MDP(20) + (-t309 * t504 - t313 * t571 + t314 * t405 + t316 * t496 + t328 * t384 + t337 * t501 + t340 * t353 + t344 * t358 - t520) * MDP(27) + (-t355 * t543 - t458 * t359 - t427 * t384 - t406 * t504 + t409 * t540 + t517 * t461 + t505 * t542) * MDP(16) + (t333 * t571 + t340 * t505 - t343 * t384 - t355 * t496 - t555 * t358 - t501 * t369 - t497 * t405 + t498 * t504 + t520) * MDP(24) + (t523 * t504 + t342 * t384 + t355 * t389 - t505 * t341 + ((-qJD(4) * t394 + t351) * t405 + t377 * t358 + t369 * qJD(4) * t417) * t469 + ((-qJD(4) * t377 - t354) * t405 - t394 * t358 + t333 * t417 + t369 * t383) * t465 + t519) * MDP(23) + (-t417 * t359 - t383 * t409 + t411 * t384 + t477 * t504) * MDP(12) + (t309 * t344 + t328 * t314 + t313 * t353 + t337 * t316 + t311 * t345 + t327 * t315 - g(1) * (-pkin(4) * t401 - qJ(5) * t400) - g(2) * (pkin(4) * t403 + qJ(5) * t402) + (-g(1) * t600 - g(2) * t509) * t471 + (g(1) * t509 - g(2) * t600) * t468) * MDP(28) + (t467 * t506 + t470 * t491) * MDP(9) + (-t467 * t491 + t470 * t506) * MDP(10) + t517 * MDP(2) + t518 * MDP(3) + (qJDD(1) * t462 + 0.2e1 * t467 * t529) * MDP(4) + (t383 * t543 + t417 * t542) * MDP(13) + 0.2e1 * (t467 * t544 - t552 * t546) * MDP(5) + qJDD(1) * MDP(1) + (-t411 * t383 + t477 * t417) * MDP(11); t481 + (t382 * t543 + (t411 * t551 - t466 * t542 - t543 * t531) * pkin(2) + t480) * MDP(17) + (t313 * t415 - t328 * t335 - t327 * t336 - t558 * t337 + (-t590 + t518 * t467 + (t327 * t535 + t328 * t534) * qJD(3)) * pkin(2) + t486 * t456 + t476) * MDP(28) + (-t457 * t340 + t511 * t469 - t521 * t496 + (-t492 + t556) * t405 + t495) * MDP(24) + (t467 * t490 - t590) * MDP(9) + (g(3) * t467 + t470 * t490) * MDP(10) + (t457 * t341 + t527 * t469 + t511 * t465 + t521 * t389 + ((-qJD(4) * t456 - t368) * t469 + (-t524 + t382) * t465) * t405 + t507) * MDP(23) + (t335 * t389 + t336 * t496 + (-t389 * t534 - t496 * t535) * qJD(3) * pkin(2) + t485 * t456 + t483) * MDP(26) + (t415 * t340 + (-qJD(4) * t337 + t580) * t469 - t558 * t496 + (-t335 + t492) * t405 + t493) * MDP(27) + (t415 * t341 + t528 * t469 + (-t580 + t584) * t465 - t558 * t389 + (-t456 * t548 - t465 * t524 + t336) * t405 + t508) * MDP(25) + MDP(7) * t544 + MDP(6) * t545 + (t381 * t543 + (-t409 * t551 + t599 * t542 - t543 * t550) * pkin(2) + t488) * MDP(16) + qJDD(2) * MDP(8) + (-t467 * t470 * MDP(4) + t552 * MDP(5)) * t474; t481 + (t378 * t543 + t480) * MDP(17) + (t486 * pkin(8) - t313 * t510 - t327 * t339 - t328 * t338 + t554 * t337 + t476) * MDP(28) + (t379 * t543 + t488) * MDP(16) + (pkin(8) * t485 + t338 * t389 + t339 * t496 + t483) * MDP(26) + (-t337 * t548 - t338 * t405 - t340 * t510 + t496 * t554 + t493 - t608) * MDP(27) + (t339 * t405 - t341 * t510 + (t584 - t596) * t465 + t554 * t389 + (-t405 * t588 + t528) * t469 + t508) * MDP(25) + (pkin(3) * t340 + t369 * t572 + t379 * t496 + t557 * t405 + t495 + t608) * MDP(24) + (-pkin(3) * t341 - t379 * t389 + (t378 * t405 + t579 - t596) * t465 + ((-t375 - t588) * t405 + t527) * t469 + t507) * MDP(23); -MDP(18) * t576 + (-t389 ^ 2 + t602) * MDP(19) + t320 * MDP(20) + (-t341 - t575) * MDP(21) + t358 * MDP(22) + (t369 * t496 + t489 + t581) * MDP(23) + (t342 * t405 + t369 * t389 - t482) * MDP(24) + (-t356 * t389 - t484 + t581 + 0.2e1 * t598) * MDP(25) + (pkin(4) * t340 - qJ(5) * t341 - (t328 - t343) * t496 + (t327 - t547) * t389) * MDP(26) + (0.2e1 * t586 - t337 * t389 - t356 * t496 + (0.2e1 * qJD(5) - t342) * t405 + t482) * MDP(27) + (t309 * qJ(5) - t311 * pkin(4) - t337 * t356 - t327 * t343 - g(1) * (-pkin(4) * t402 + qJ(5) * t403) - g(2) * (-pkin(4) * t400 + qJ(5) * t401) + t514 * t451 + t547 * t328) * MDP(28); t320 * MDP(26) + (-t601 - t602) * MDP(27) + (-t328 * t405 + t484 - t598) * MDP(28) + (-t576 - t358) * MDP(25);];
tau = t1;
