% Calculate vector of inverse dynamics joint torques for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRPPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:42
% EndTime: 2019-12-31 19:24:52
% DurationCPUTime: 6.91s
% Computational Cost: add. (3919->513), mult. (10902->645), div. (0->0), fcn. (7982->8), ass. (0->222)
t491 = cos(pkin(5));
t488 = sin(pkin(8));
t492 = sin(qJ(2));
t576 = qJD(1) * t492;
t554 = t488 * t576;
t542 = t491 * t554;
t450 = qJD(2) * t542;
t490 = cos(pkin(8));
t494 = cos(qJ(2));
t565 = qJD(1) * qJD(2);
t549 = t494 * t565;
t564 = qJDD(1) * t492;
t514 = t549 + t564;
t489 = sin(pkin(5));
t562 = qJDD(2) * t489;
t563 = qJDD(1) * t494;
t630 = t491 * t563 + t562;
t633 = t488 * t630 + t514 * t490 - t450;
t546 = qJ(3) * t491 + pkin(7);
t530 = qJD(2) * t546;
t569 = qJD(3) * t492;
t401 = -t491 * t569 - t494 * t530;
t448 = t546 * t492;
t366 = qJDD(2) * pkin(2) + qJD(1) * t401 - qJDD(1) * t448;
t596 = t489 * t494;
t558 = qJ(3) * t596;
t612 = pkin(2) * t492;
t521 = -t558 + t612;
t399 = qJD(2) * t521 - t489 * t569;
t598 = t489 * t492;
t473 = qJ(3) * t598;
t578 = t494 * pkin(2) + t473;
t626 = -pkin(1) - t578;
t367 = qJD(1) * t399 + qJDD(1) * t626;
t623 = t490 * (t366 * t491 + t367 * t489);
t593 = t491 * t492;
t631 = t490 * t593;
t629 = t492 * t565 - t563;
t627 = MDP(17) + MDP(20);
t624 = t490 * (-t448 * t491 + t489 * t626);
t429 = t488 * t494 + t631;
t416 = t429 * qJD(2);
t580 = t630 * t490;
t520 = t488 * t564 - t580;
t362 = qJD(1) * t416 + t520;
t575 = qJD(1) * t494;
t553 = t491 * t575;
t573 = qJD(2) * t490;
t385 = -t489 * t573 - t490 * t553 + t554;
t622 = t362 * qJ(5) + t385 * qJD(5);
t438 = t521 * qJD(1);
t531 = qJD(1) * t546;
t439 = t492 * t531;
t440 = t494 * t531;
t600 = t488 * t491;
t601 = t488 * t489;
t351 = t438 * t601 - t490 * t439 - t440 * t600;
t570 = qJD(3) * t490;
t582 = -qJD(4) * t491 + t351 + (qJ(4) * t576 - t570) * t489;
t591 = t491 * t494;
t428 = t488 * t591 + t490 * t492;
t363 = qJDD(1) * t428 + t488 * t562 + t490 * t549 - t450;
t444 = -t491 * qJD(2) + t489 * t575;
t620 = pkin(4) * t363 + qJD(5) * t444;
t619 = t489 ^ 2 + t491 ^ 2;
t618 = MDP(15) + MDP(19);
t617 = MDP(16) - MDP(21);
t495 = cos(qJ(1));
t586 = t495 * t491;
t493 = sin(qJ(1));
t590 = t492 * t493;
t435 = t489 * t590 - t586;
t589 = t492 * t495;
t592 = t491 * t493;
t436 = t489 * t589 + t592;
t616 = -g(1) * t436 - g(2) * t435 + g(3) * t596;
t400 = qJD(3) * t591 - t492 * t530;
t341 = t399 * t601 + t490 * t400 + t401 * t600;
t572 = qJD(2) * t492;
t330 = -t489 * (qJ(4) * t572 - qJD(4) * t494) - t341;
t574 = qJD(2) * t489;
t446 = t553 + t574;
t405 = pkin(7) * t575 + qJ(3) * t446;
t422 = qJD(2) * pkin(2) - t439;
t423 = t626 * qJD(1);
t344 = -t488 * t405 + (t422 * t491 + t423 * t489) * t490;
t441 = t444 ^ 2;
t398 = qJDD(2) * t491 + t629 * t489;
t611 = pkin(3) * t398;
t608 = g(1) * t493;
t606 = g(3) * t492;
t605 = pkin(3) + qJ(5);
t386 = qJD(1) * t428 + t488 * t574;
t604 = t385 * t386;
t603 = t399 * t490;
t602 = t438 * t490;
t599 = t489 * t490;
t597 = t489 * t493;
t595 = t489 * t495;
t594 = t490 * t491;
t588 = t493 * t494;
t587 = t494 * t495;
t414 = t429 * qJD(1);
t370 = t491 * t438 + t440 * t489;
t415 = t490 * t575 - t542;
t518 = -qJ(4) * t415 + t370;
t567 = qJD(4) * t488;
t585 = -t414 * t605 + (-qJD(5) * t490 - t567) * t489 - t518;
t449 = t546 * t494;
t433 = t488 * t449;
t584 = pkin(3) * t596 + t433;
t419 = t488 * t439;
t535 = t440 * t594 - t419;
t551 = t605 * t492;
t571 = qJD(3) * t488;
t552 = t489 * t571;
t583 = -qJD(5) * t491 + t552 - pkin(4) * t415 - (-qJD(1) * t551 - t602) * t489 - t535;
t581 = pkin(4) * t414 - t582;
t432 = pkin(2) * t600 + qJ(3) * t599;
t579 = t495 * pkin(7) + qJ(3) * t586;
t486 = t492 ^ 2;
t577 = -t494 ^ 2 + t486;
t568 = qJD(4) * t386;
t357 = t446 * qJD(3) - t629 * pkin(7) + (-t491 * t629 + t562) * qJ(3);
t353 = t488 * t357;
t561 = qJDD(4) + t353;
t560 = pkin(2) * t590;
t559 = pkin(2) * t589;
t557 = t488 * t593;
t556 = t490 * t591;
t321 = t490 * t357 + t366 * t600 + t367 * t601;
t345 = t490 * t405 + t422 * t600 + t423 * t601;
t356 = -t448 * t600 + t490 * t449 + t601 * t626;
t555 = -pkin(2) * t490 - pkin(3);
t547 = -qJ(4) * t488 - pkin(2);
t358 = t491 * t399 - t401 * t489;
t371 = t448 * t489 + t491 * t626;
t375 = t488 * t588 + (t491 * t590 + t595) * t490;
t377 = t429 * t495 - t490 * t597;
t541 = g(1) * t375 - g(2) * t377;
t376 = -t488 * t595 + t490 * t588 - t493 * t557;
t378 = t488 * t597 + t490 * t587 - t495 * t557;
t540 = g(1) * t376 - g(2) * t378;
t539 = g(1) * t435 - g(2) * t436;
t538 = g(1) * t495 + g(2) * t493;
t537 = -g(2) * t495 + t608;
t388 = t488 * t400;
t536 = -t401 * t594 + t388;
t406 = t488 * t590 - t493 * t556;
t407 = t428 * t493;
t461 = t493 * t558;
t534 = -t407 * pkin(3) - qJ(4) * t406 + t461;
t408 = t488 * t589 - t495 * t556;
t409 = t428 * t495;
t463 = t495 * t558;
t533 = -t409 * pkin(3) - qJ(4) * t408 + t463;
t410 = -qJ(4) * t491 - t432;
t365 = -t422 * t489 + t491 * t423 + qJD(3);
t333 = -t366 * t489 + t491 * t367 + qJDD(3);
t318 = -t398 * qJ(4) + t444 * qJD(4) - t321;
t529 = pkin(2) * t587 + t493 * pkin(7) + qJ(3) * t592 + (pkin(1) + t473) * t495;
t524 = pkin(4) * t596 - t612;
t523 = -t376 * pkin(3) - qJ(4) * t375 + t579;
t337 = qJ(4) * t444 - t345;
t430 = t490 * t494 - t557;
t522 = t430 * pkin(3) + qJ(4) * t429 + t578;
t519 = -0.2e1 * pkin(1) * t565 - pkin(7) * qJDD(2);
t517 = -qJ(4) * t428 + t371;
t349 = qJ(4) * t596 - t356;
t511 = g(1) * t408 + g(2) * t406 - g(3) * t429;
t510 = g(1) * t409 + g(2) * t407 - g(3) * t430;
t509 = -qJ(4) * t386 + t365;
t316 = -pkin(4) * t362 + qJDD(5) - t318;
t508 = t378 * pkin(3) + qJ(4) * t377 + t529;
t507 = t626 * t608;
t417 = -qJD(2) * t557 + t494 * t573;
t506 = -qJ(4) * t417 - qJD(4) * t428 + t358;
t505 = -t494 * t538 - t606;
t504 = qJD(4) - t344;
t496 = qJD(2) ^ 2;
t503 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t496 + t537;
t497 = qJD(1) ^ 2;
t502 = pkin(1) * t497 - pkin(7) * qJDD(1) + t538;
t501 = t561 - t623;
t427 = t488 * t492 - t556;
t499 = -g(1) * t377 - g(2) * t375 - g(3) * t427 + t561 - t611;
t317 = pkin(3) * t362 - qJ(4) * t363 + t333 - t568;
t498 = t317 + t616;
t470 = qJ(3) * t601;
t431 = pkin(2) * t594 - t470;
t412 = (-pkin(3) * t490 + t547) * t489;
t411 = t491 * t555 + t470;
t383 = (-t490 * t605 + t547) * t489;
t382 = pkin(4) * t599 - t410;
t372 = pkin(4) * t601 + t470 + (-qJ(5) + t555) * t491;
t369 = t446 * t490 - t619 * t554;
t368 = t619 * t490 * t576 + t446 * t488;
t355 = -t433 + t624;
t352 = t584 - t624;
t350 = t419 + (t438 * t489 - t440 * t491) * t490;
t348 = pkin(3) * t427 + t517;
t347 = (-pkin(3) * t576 - t602) * t489 + t535;
t343 = pkin(3) * t414 + t518;
t342 = -pkin(4) * t427 - t349;
t340 = -t388 + (t399 * t489 + t401 * t491) * t490;
t339 = t427 * t605 + t517;
t338 = t448 * t594 + pkin(4) * t428 + (qJ(5) * t494 - t490 * t626) * t489 + t584;
t336 = pkin(3) * t444 + t504;
t334 = (-pkin(3) * t572 - t603) * t489 + t536;
t329 = pkin(3) * t385 + t509;
t328 = pkin(3) * t416 + t506;
t327 = -pkin(4) * t385 + qJD(5) - t337;
t326 = -pkin(4) * t416 - t330;
t325 = pkin(4) * t417 + (-qJD(2) * t551 + qJD(5) * t494 - t603) * t489 + t536;
t324 = t385 * t605 + t509;
t323 = pkin(4) * t386 + t444 * t605 + t504;
t322 = qJD(5) * t427 + t416 * t605 + t506;
t320 = -t353 + t623;
t319 = t501 - t611;
t315 = -t398 * t605 + t501 + t620;
t314 = t317 + t622;
t1 = [qJDD(1) * MDP(1) + t537 * MDP(2) + t538 * MDP(3) + (qJDD(1) * t486 + 0.2e1 * t492 * t549) * MDP(4) + 0.2e1 * (t492 * t563 - t565 * t577) * MDP(5) + (qJDD(2) * t492 + t494 * t496) * MDP(6) + (qJDD(2) * t494 - t492 * t496) * MDP(7) + (t492 * t519 + t494 * t503) * MDP(9) + (-t492 * t503 + t494 * t519) * MDP(10) + (t333 * t427 - t340 * t444 + t355 * t398 + t358 * t385 + t362 * t371 + t365 * t416 + (-t320 * t494 + t344 * t572) * t489 + t540) * MDP(11) + (t333 * t428 + t341 * t444 - t356 * t398 + t358 * t386 + t363 * t371 + t365 * t417 + (t321 * t494 - t345 * t572) * t489 - t541) * MDP(12) + (-t320 * t428 - t321 * t427 - t340 * t386 - t341 * t385 - t344 * t417 - t345 * t416 - t355 * t363 - t356 * t362 + t539) * MDP(13) + (-g(1) * t579 - g(2) * t529 + t320 * t355 + t321 * t356 + t333 * t371 + t344 * t340 + t345 * t341 + t365 * t358 - t507) * MDP(14) + (t318 * t427 + t319 * t428 + t330 * t385 + t334 * t386 + t336 * t417 + t337 * t416 + t349 * t362 + t352 * t363 + t539) * MDP(15) + (-t317 * t427 - t328 * t385 - t329 * t416 - t334 * t444 - t348 * t362 + t352 * t398 + (-t319 * t494 + t336 * t572) * t489 - t540) * MDP(16) + (-t317 * t428 - t328 * t386 - t329 * t417 + t330 * t444 - t348 * t363 - t349 * t398 + (t318 * t494 - t337 * t572) * t489 + t541) * MDP(17) + (-g(1) * t523 - g(2) * t508 + t317 * t348 + t318 * t349 + t319 * t352 + t329 * t328 + t337 * t330 + t336 * t334 - t507) * MDP(18) + (t315 * t428 - t316 * t427 + t323 * t417 + t325 * t386 - t326 * t385 - t327 * t416 + t338 * t363 - t342 * t362 + t539) * MDP(19) + (-t314 * t428 - t322 * t386 - t324 * t417 - t326 * t444 - t339 * t363 + t342 * t398 + (-t316 * t494 + t327 * t572) * t489 + t541) * MDP(20) + (t314 * t427 + t322 * t385 + t324 * t416 + t325 * t444 - t338 * t398 + t339 * t362 + (t315 * t494 - t323 * t572) * t489 + t540) * MDP(21) + (t314 * t339 + t324 * t322 + t315 * t338 + t323 * t325 + t316 * t342 + t327 * t326 - g(1) * (-pkin(4) * t435 - qJ(5) * t376 + t523) - g(2) * (pkin(4) * t436 + qJ(5) * t378 + t508) - t507) * MDP(22); MDP(6) * t564 + MDP(7) * t563 + qJDD(2) * MDP(8) + (-g(3) * t494 + t492 * t502) * MDP(9) + (t494 * t502 + t606) * MDP(10) + (t320 * t491 + t350 * t444 - t365 * t414 - t370 * t385 + t398 * t431 + (-pkin(2) * t362 - t333 * t490 - t344 * t576 + t444 * t571) * t489 + t510) * MDP(11) + (-t321 * t491 - t351 * t444 - t365 * t415 - t370 * t386 - t398 * t432 + (-pkin(2) * t363 + t333 * t488 + t345 * t576 + t444 * t570) * t489 - t511) * MDP(12) + (t344 * t415 + t345 * t414 + t350 * t386 + t351 * t385 - t362 * t432 - t363 * t431 + (-t320 * t488 + t321 * t490 + (-t385 * t490 + t386 * t488) * qJD(3) + t505) * t489) * MDP(13) + (t321 * t432 + t320 * t431 - t345 * t351 - t344 * t350 - t365 * t370 - g(1) * (t463 - t559) - g(2) * (t461 - t560) - g(3) * t578 + (-t333 * pkin(2) + (-t344 * t488 + t345 * t490) * qJD(3)) * t489) * MDP(14) + (-t336 * t415 - t337 * t414 - t347 * t386 + t362 * t410 + t363 * t411 + t582 * t385 + (-t318 * t490 + (qJD(3) * t386 + t319) * t488 + t505) * t489) * MDP(15) + (t319 * t491 + t329 * t414 + t343 * t385 + t347 * t444 - t362 * t412 + t398 * t411 + (-t336 * t576 + t317 * t490 + (-qJD(3) * t444 + qJD(4) * t385) * t488) * t489 - t510) * MDP(16) + (-t318 * t491 + t329 * t415 + t343 * t386 - t363 * t412 - t398 * t410 + t582 * t444 + (t337 * t576 + (-t317 + t568) * t488) * t489 + t511) * MDP(17) + (t317 * t412 + t318 * t410 + t319 * t411 - g(1) * (t533 - t559) - g(2) * (t534 - t560) - g(3) * t522 + t582 * t337 + (-t347 + t552) * t336 + (-t489 * t567 - t343) * t329) * MDP(18) + (-t323 * t415 + t327 * t414 - t362 * t382 + t363 * t372 + t583 * t386 - t581 * t385 + (t315 * t488 + t316 * t490 + t505) * t489) * MDP(19) + (t316 * t491 + t324 * t415 - t363 * t383 + t382 * t398 + (-t314 * t488 - t327 * t576) * t489 - t581 * t444 - t585 * t386 + t511) * MDP(20) + (-t315 * t491 - t324 * t414 + t362 * t383 - t372 * t398 + (-t314 * t490 + t323 * t576) * t489 + t583 * t444 + t585 * t385 + t510) * MDP(21) + (t314 * t383 + t315 * t372 + t316 * t382 - g(1) * (-qJ(5) * t409 + t495 * t524 + t533) - g(2) * (-qJ(5) * t407 + t493 * t524 + t534) - g(3) * (pkin(4) * t598 + qJ(5) * t430 + t522) + t581 * t327 + t585 * t324 + t583 * t323) * MDP(22) + (-MDP(4) * t492 * t494 + MDP(5) * t577) * t497; (t344 * t368 - t345 * t369 + t333 + t616) * MDP(14) + (-t336 * t368 + t337 * t369 + t498) * MDP(18) + (-t323 * t368 - t327 * t369 + t498 + t622) * MDP(22) + (MDP(13) + t618) * (-t368 * t386 + t369 * t385) + (-t368 * t444 + t488 * t514 + t565 * t631 - t580) * (MDP(11) - t617) + (-t369 * t444 + t633) * (MDP(12) - t627); (t329 * t386 - t337 * t444 + t499) * MDP(18) + (-qJ(5) * t398 + t324 * t386 + t327 * t444 + t499 + t620) * MDP(22) + (-MDP(18) - MDP(22)) * t623 + t617 * (t398 - t604) + t627 * (-t386 ^ 2 - t441) + (-t385 * t444 + t633) * t618; (-t386 * t444 - t429 * t565 - t520) * MDP(19) + (t398 + t604) * MDP(20) + (-t385 ^ 2 - t441) * MDP(21) + (-g(1) * t378 - g(2) * t376 - g(3) * t428 - t323 * t444 - t324 * t385 + t316) * MDP(22);];
tau = t1;
