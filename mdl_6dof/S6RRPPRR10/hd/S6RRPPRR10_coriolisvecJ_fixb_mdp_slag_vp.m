% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:12
% EndTime: 2019-03-09 09:37:25
% DurationCPUTime: 7.66s
% Computational Cost: add. (4432->481), mult. (10383->645), div. (0->0), fcn. (7289->8), ass. (0->214)
t518 = sin(qJ(2));
t588 = qJD(1) * t518;
t499 = qJD(5) + t588;
t491 = qJD(6) + t499;
t514 = cos(pkin(10));
t521 = cos(qJ(2));
t587 = qJD(1) * t521;
t570 = t514 * t587;
t513 = sin(pkin(10));
t586 = qJD(2) * t513;
t465 = t570 + t586;
t568 = t513 * t587;
t585 = qJD(2) * t514;
t467 = -t568 + t585;
t517 = sin(qJ(5));
t520 = cos(qJ(5));
t405 = t520 * t465 + t467 * t517;
t519 = cos(qJ(6));
t404 = t465 * t517 - t467 * t520;
t516 = sin(qJ(6));
t610 = t404 * t516;
t634 = -t519 * t405 + t610;
t633 = t491 * t634;
t469 = t513 * t520 + t514 * t517;
t454 = t469 * qJD(5);
t532 = t518 * t469;
t593 = -qJD(1) * t532 - t454;
t617 = pkin(3) + pkin(7);
t632 = t404 * t499;
t631 = t405 * t499;
t539 = t404 * t519 + t405 * t516;
t630 = t491 * t539;
t571 = t514 * t588;
t580 = qJD(5) * t520;
t581 = qJD(5) * t517;
t606 = t513 * t517;
t592 = -t513 * t581 + t514 * t580 + t520 * t571 - t588 * t606;
t537 = -t514 * t520 + t606;
t501 = pkin(7) * t588;
t628 = qJD(3) + t501;
t515 = -pkin(2) - qJ(4);
t565 = -qJ(3) * t518 - pkin(1);
t464 = t515 * t521 + t565;
t430 = t464 * qJD(1);
t576 = pkin(3) * t588 + t628;
t436 = qJD(2) * t515 + t576;
t377 = -t430 * t513 + t514 * t436;
t353 = pkin(4) * t588 - pkin(8) * t467 + t377;
t378 = t514 * t430 + t513 * t436;
t360 = -pkin(8) * t465 + t378;
t331 = t353 * t517 + t360 * t520;
t329 = -pkin(9) * t405 + t331;
t579 = qJD(6) * t516;
t327 = t329 * t579;
t502 = pkin(7) * t587;
t503 = pkin(3) * t587;
t476 = t502 + t503;
t510 = qJD(2) * qJ(3);
t620 = qJD(4) + t510;
t449 = t476 + t620;
t413 = pkin(4) * t465 + t449;
t354 = pkin(5) * t405 + t413;
t627 = -t354 * t634 + t327;
t575 = qJD(1) * qJD(2);
t567 = t518 * t575;
t551 = t520 * t567;
t552 = t517 * t567;
t367 = -t465 * t580 - t467 * t581 + t513 * t551 + t514 * t552;
t498 = pkin(2) * t567;
t544 = -qJ(3) * t521 + qJ(4) * t518;
t582 = qJD(3) * t518;
t525 = qJD(2) * t544 - qJD(4) * t521 - t582;
t403 = qJD(1) * t525 + t498;
t566 = t521 * t575;
t497 = pkin(7) * t566;
t444 = t497 + (-qJD(4) + t503) * qJD(2);
t365 = -t403 * t513 + t514 * t444;
t605 = t513 * t518;
t535 = pkin(4) * t521 - pkin(8) * t605;
t528 = t535 * qJD(2);
t341 = qJD(1) * t528 + t365;
t366 = t514 * t403 + t513 * t444;
t584 = qJD(2) * t518;
t553 = pkin(8) * t514 * t584;
t352 = qJD(1) * t553 + t366;
t560 = t520 * t341 - t352 * t517;
t527 = -qJD(5) * t331 + t560;
t320 = pkin(5) * t566 - pkin(9) * t367 + t527;
t368 = -qJD(5) * t404 + t513 * t552 - t514 * t551;
t531 = t517 * t341 + t520 * t352 + t353 * t580 - t360 * t581;
t321 = -pkin(9) * t368 + t531;
t561 = t519 * t320 - t516 * t321;
t626 = t354 * t539 + t561;
t625 = MDP(30) * t566 + (t539 ^ 2 - t634 ^ 2) * MDP(27) + t634 * MDP(26) * t539;
t624 = -0.2e1 * t575;
t511 = t518 ^ 2;
t623 = MDP(5) * (-t521 ^ 2 + t511);
t488 = t617 * t518;
t472 = t514 * t488;
t391 = pkin(4) * t518 + t472 + (pkin(8) * t521 - t464) * t513;
t412 = t514 * t464 + t513 * t488;
t602 = t514 * t521;
t395 = -pkin(8) * t602 + t412;
t595 = t517 * t391 + t520 * t395;
t622 = t593 * t519;
t410 = -t469 * t516 - t519 * t537;
t616 = -pkin(8) + t515;
t478 = t616 * t513;
t479 = t616 * t514;
t591 = t520 * t478 + t517 * t479;
t505 = pkin(2) * t588;
t445 = qJD(1) * t544 + t505;
t398 = -t445 * t513 + t514 * t476;
t373 = qJD(1) * t535 + t398;
t399 = t514 * t445 + t513 * t476;
t382 = pkin(8) * t571 + t399;
t621 = qJD(4) * t537 - qJD(5) * t591 - t520 * t373 + t517 * t382;
t572 = -pkin(4) * t514 - pkin(3);
t577 = -t572 * t588 + t628;
t583 = qJD(2) * t521;
t619 = t518 * (-t449 + t620) - t515 * t583;
t558 = t516 * t367 + t519 * t368;
t325 = -qJD(6) * t539 + t558;
t618 = qJD(4) * t469 + t517 * t373 + t520 * t382 + t478 * t581 - t479 * t580;
t615 = qJD(2) * pkin(2);
t330 = t520 * t353 - t360 * t517;
t328 = pkin(9) * t404 + t330;
t326 = pkin(5) * t499 + t328;
t614 = t326 * t519;
t613 = t329 * t519;
t612 = t377 * t521;
t611 = t378 * t521;
t475 = t617 * t584;
t509 = qJD(2) * qJD(3);
t443 = -qJD(1) * t475 + t509;
t609 = t443 * t513;
t523 = qJD(1) ^ 2;
t607 = t511 * t523;
t604 = t514 * t518;
t522 = qJD(2) ^ 2;
t601 = t518 * t522;
t600 = t521 * t522;
t599 = t521 * t523;
t500 = t513 * pkin(4) + qJ(3);
t409 = t519 * t469 - t516 * t537;
t598 = -qJD(6) * t409 - t516 * t592 + t622;
t597 = qJD(6) * t410 + t516 * t593 + t519 * t592;
t504 = pkin(2) * t584;
t416 = t504 + t525;
t477 = t617 * t583;
t381 = t514 * t416 + t513 * t477;
t594 = pkin(5) * t592 + t577;
t489 = t617 * t521;
t485 = -pkin(2) * t521 + t565;
t460 = qJD(1) * t485;
t578 = qJD(6) * t519;
t574 = t518 * t599;
t573 = t519 * t367 - t516 * t368 - t405 * t578;
t453 = pkin(4) * t602 + t489;
t564 = t592 * t499;
t563 = pkin(1) * t624;
t562 = qJD(3) - t615;
t380 = -t416 * t513 + t514 * t477;
t364 = t380 + t528;
t370 = t553 + t381;
t559 = t520 * t364 - t370 * t517;
t556 = t520 * t391 - t395 * t517;
t555 = -t478 * t517 + t520 * t479;
t554 = qJD(6) * t326 + t321;
t490 = t518 * t566;
t550 = t499 * t593 - t537 * t566;
t386 = -pkin(9) * t469 + t591;
t548 = pkin(5) * t587 + pkin(9) * t593 + qJD(6) * t386 - t621;
t385 = pkin(9) * t537 + t555;
t547 = pkin(9) * t592 - qJD(6) * t385 + t618;
t546 = qJD(6) * t537 - t592;
t319 = t326 * t516 + t613;
t441 = t469 * t521;
t336 = pkin(5) * t518 + pkin(9) * t441 + t556;
t440 = t537 * t521;
t337 = pkin(9) * t440 + t595;
t543 = t336 * t516 + t337 * t519;
t542 = t365 * t514 + t366 * t513;
t541 = -t377 * t513 + t378 * t514;
t538 = t519 * t440 + t441 * t516;
t389 = t440 * t516 - t441 * t519;
t536 = -0.2e1 * qJD(2) * t460;
t533 = -qJ(3) * t583 - t582;
t431 = qJD(1) * t533 + t498;
t447 = t504 + t533;
t534 = pkin(7) * t522 + qJD(1) * t447 + t431;
t433 = (-pkin(7) + t572) * t584;
t530 = t517 * t364 + t520 * t370 + t391 * t580 - t395 * t581;
t324 = t404 * t579 + t573;
t415 = qJD(1) * t433 + t509;
t483 = pkin(7) * t567 - t509;
t484 = t501 + t562;
t487 = -t502 - t510;
t524 = -t483 * t521 + (t484 * t521 + (t487 + t502) * t518) * qJD(2);
t473 = -qJ(3) * t587 + t505;
t434 = t460 * t588;
t427 = pkin(5) * t469 + t500;
t411 = -t464 * t513 + t472;
t402 = -pkin(5) * t440 + t453;
t394 = t454 * t521 - t537 * t584;
t393 = qJD(2) * t532 + qJD(5) * t440;
t359 = -pkin(5) * t394 + t433;
t338 = pkin(5) * t368 + t415;
t333 = qJD(6) * t389 + t393 * t516 - t519 * t394;
t332 = qJD(6) * t538 + t393 * t519 + t394 * t516;
t323 = pkin(9) * t394 + t530;
t322 = pkin(5) * t583 - pkin(9) * t393 - qJD(5) * t595 + t559;
t318 = -t329 * t516 + t614;
t1 = [t623 * t624 + MDP(6) * t600 + ((t322 * t519 - t323 * t516) * t491 + t561 * t518 - t359 * t634 + t402 * t325 - t338 * t538 + t354 * t333 + (-t319 * t518 - t491 * t543) * qJD(6) + ((t336 * t519 - t337 * t516) * qJD(1) + t318) * t583) * MDP(31) + (-t325 * t518 - t333 * t491 + (qJD(1) * t538 + t634) * t583) * MDP(29) + (t324 * t538 - t325 * t389 + t332 * t634 + t333 * t539) * MDP(27) + (-t521 * t609 - t467 * t475 + (-qJD(1) * t381 - t366) * t518 + (t449 * t605 - t611 + (-t412 * t521 + t489 * t605) * qJD(1)) * qJD(2)) * MDP(16) + (t443 * t602 - t465 * t475 + (qJD(1) * t380 + t365) * t518 + (-t449 * t604 + t612 + (t411 * t521 - t489 * t604) * qJD(1)) * qJD(2)) * MDP(15) + (t365 * t411 + t366 * t412 + t377 * t380 + t378 * t381 + t443 * t489 - t449 * t475) * MDP(18) + 0.2e1 * MDP(4) * t490 + t524 * MDP(11) + (pkin(7) * t524 + t431 * t485 + t447 * t460) * MDP(14) + (-t368 * t518 + t394 * t499 + (qJD(1) * t440 - t405) * t583) * MDP(22) + (t499 * t583 + t490) * MDP(23) + (t491 * t583 + t490) * MDP(30) + (-t380 * t467 - t381 * t465 + (t365 * t513 - t366 * t514) * t521 + ((-t411 * t513 + t412 * t514) * qJD(1) + t541) * t584) * MDP(17) - MDP(7) * t601 + (pkin(7) * t601 + t521 * t563) * MDP(10) + (-pkin(7) * t600 + t518 * t563) * MDP(9) + (-t518 * t534 + t521 * t536) * MDP(13) + (t518 * t536 + t521 * t534) * MDP(12) + (t559 * t499 + t560 * t518 + t433 * t405 + t453 * t368 - t415 * t440 - t413 * t394 + (-t331 * t518 - t499 * t595) * qJD(5) + (qJD(1) * t556 + t330) * t583) * MDP(24) + (t324 * t389 - t332 * t539) * MDP(26) + (t324 * t518 + t332 * t491 + (qJD(1) * t389 - t539) * t583) * MDP(28) + (t402 * t324 + t327 * t518 + t354 * t332 + t338 * t389 - t359 * t539 + (-(-qJD(6) * t337 + t322) * t491 - t320 * t518) * t516 + (-(qJD(6) * t336 + t323) * t491 - t554 * t518) * t519 + (-qJD(1) * t543 - t319) * t583) * MDP(32) + (t367 * t440 + t368 * t441 - t393 * t405 - t394 * t404) * MDP(20) + (-t367 * t441 - t393 * t404) * MDP(19) + (t367 * t518 + t393 * t499 + (-qJD(1) * t441 - t404) * t583) * MDP(21) + (-t530 * t499 - t531 * t518 - t433 * t404 + t453 * t367 - t415 * t441 + t413 * t393 + (-qJD(1) * t595 - t331) * t583) * MDP(25); -MDP(4) * t574 + (-qJ(3) * t483 - qJD(3) * t487 - t460 * t473 + (-t487 * t518 + (-t484 - t615) * t521) * qJD(1) * pkin(7)) * MDP(14) + (t443 * t514 + t576 * t467 + (t399 * t518 + t513 * t619 + t611) * qJD(1)) * MDP(16) + (t609 + t576 * t465 + (-t398 * t518 - t514 * t619 - t612) * qJD(1)) * MDP(15) + (t500 * t368 + t577 * t405 + t592 * t413 + t415 * t469 + t621 * t499) * MDP(24) + (-t367 * t537 - t404 * t593) * MDP(19) + (-t367 * t469 + t368 * t537 + t404 * t592 - t405 * t593) * MDP(20) - t564 * MDP(22) + t434 * MDP(12) + t550 * MDP(21) + (t398 * t467 + t399 * t465 + (qJD(4) * t467 - t378 * t588 - t365) * t514 + (qJD(4) * t465 + t377 * t588 - t366) * t513) * MDP(17) + (qJ(3) * t443 - t377 * t398 - t378 * t399 + t542 * t515 + t576 * t449 + (-t377 * t514 - t378 * t513) * qJD(4)) * MDP(18) + (t500 * t367 - t404 * t577 + t593 * t413 - t415 * t537 + t618 * t499) * MDP(25) + (t427 * t325 + t338 * t409 + t597 * t354 - t594 * t634) * MDP(31) + (t324 * t410 - t539 * t598) * MDP(26) + (t427 * t324 + t338 * t410 + t598 * t354 - t539 * t594) * MDP(32) + (-t324 * t409 - t325 * t410 + t539 * t597 + t598 * t634) * MDP(27) + t523 * t623 + ((-t487 - t510) * t518 + (-t484 + t562) * t521) * qJD(1) * MDP(11) + (0.2e1 * t509 + (t460 * t521 + t473 * t518) * qJD(1)) * MDP(13) + ((t516 * t547 - t519 * t548) * MDP(31) - t597 * MDP(29) + (t516 * t548 + t519 * t547) * MDP(32) + t598 * MDP(28)) * t491 + ((qJD(2) * t555 - t330) * MDP(24) + (-qJD(2) * t469 + t405) * MDP(22) - t473 * MDP(12) + t404 * MDP(21) + (-qJD(2) * t591 + t331) * MDP(25) + ((t385 * t519 - t386 * t516) * qJD(2) - t318) * MDP(31) + (-qJD(2) * t409 - t634) * MDP(29) + (-(t385 * t516 + t386 * t519) * qJD(2) + t319) * MDP(32) + (qJD(2) * t410 + t539) * MDP(28) - t499 * MDP(23) - t491 * MDP(30)) * t587 + (MDP(9) * t518 * t523 + MDP(10) * t599) * pkin(1); MDP(12) * t574 + (-t522 - t607) * MDP(13) + (qJD(2) * t487 + t434 + t497) * MDP(14) + (-t513 * t607 + (-t465 + t570) * qJD(2)) * MDP(15) + (-t514 * t607 + (-t467 - t568) * qJD(2)) * MDP(16) + (-t465 * t514 + t467 * t513) * MDP(17) * t588 + (-qJD(2) * t449 + t541 * t588 + t542) * MDP(18) + (-qJD(2) * t405 + t550) * MDP(24) + (-t564 + (-t469 * t587 + t404) * qJD(2)) * MDP(25) + ((-t469 * t578 + t516 * t546 + t622) * t491 + (t410 * t587 + t634) * qJD(2)) * MDP(31) + ((t546 * t519 + (qJD(6) * t469 - t593) * t516) * t491 + (-t409 * t587 + t539) * qJD(2)) * MDP(32); (-t465 ^ 2 - t467 ^ 2) * MDP(17) + (t377 * t467 + t378 * t465 + t443) * MDP(18) + (t368 - t632) * MDP(24) + (t367 - t631) * MDP(25) + (t325 - t630) * MDP(31) + (t324 + t633) * MDP(32) + ((t467 - t585) * MDP(15) + (-t465 + t586) * MDP(16)) * t588; -t404 * t405 * MDP(19) + (t404 ^ 2 - t405 ^ 2) * MDP(20) + (t367 + t631) * MDP(21) + (-t368 - t632) * MDP(22) + MDP(23) * t566 + (t331 * t499 + t404 * t413 + t527) * MDP(24) + (t330 * t499 + t405 * t413 - t531) * MDP(25) + (t324 - t633) * MDP(28) + (-t325 - t630) * MDP(29) + (-(-t328 * t516 - t613) * t491 - t319 * qJD(6) + (-t404 * t634 - t491 * t579 + t519 * t566) * pkin(5) + t626) * MDP(31) + ((-t329 * t491 - t320) * t516 + (t328 * t491 - t554) * t519 + (-t404 * t539 - t491 * t578 - t516 * t566) * pkin(5) + t627) * MDP(32) + t625; (t573 - t633) * MDP(28) + (-t558 - t630) * MDP(29) + (t319 * t491 + t626) * MDP(31) + (t318 * t491 - t516 * t320 - t519 * t321 + t627) * MDP(32) + (MDP(28) * t610 + MDP(29) * t539 - MDP(31) * t319 - MDP(32) * t614) * qJD(6) + t625;];
tauc  = t1;
