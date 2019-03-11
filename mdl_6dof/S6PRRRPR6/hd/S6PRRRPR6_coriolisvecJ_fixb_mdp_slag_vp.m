% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:57
% EndTime: 2019-03-08 23:37:09
% DurationCPUTime: 6.98s
% Computational Cost: add. (3200->509), mult. (7838->685), div. (0->0), fcn. (5564->10), ass. (0->223)
t472 = sin(qJ(4));
t476 = cos(qJ(4));
t549 = t476 * qJD(3);
t473 = sin(qJ(3));
t562 = qJD(2) * t473;
t424 = t472 * t562 - t549;
t560 = qJD(3) * t472;
t426 = t476 * t562 + t560;
t471 = sin(qJ(6));
t475 = cos(qJ(6));
t368 = t424 * t471 + t426 * t475;
t496 = -t475 * t424 + t426 * t471;
t550 = t473 * MDP(27);
t520 = qJD(3) * t550;
t636 = MDP(23) * t368 * t496 + (t368 ^ 2 - t496 ^ 2) * MDP(24) - qJD(2) * t520;
t477 = cos(qJ(3));
t561 = qJD(2) * t477;
t633 = qJD(4) - t561;
t546 = -qJD(6) + t633;
t635 = t368 * t546;
t508 = pkin(3) * t473 - pkin(9) * t477;
t433 = t508 * qJD(3);
t438 = -pkin(3) * t477 - pkin(9) * t473 - pkin(2);
t474 = sin(qJ(2));
t556 = qJD(4) * t476;
t469 = sin(pkin(6));
t565 = qJD(1) * t469;
t478 = cos(qJ(2));
t582 = t477 * t478;
t634 = -(t472 * t474 + t476 * t582) * t565 + t472 * t433 + t438 * t556;
t632 = qJD(6) - qJD(4);
t534 = t474 * t565;
t435 = qJD(2) * pkin(8) + t534;
t470 = cos(pkin(6));
t564 = qJD(1) * t477;
t401 = -t473 * t435 + t470 * t564;
t509 = qJD(3) * pkin(3) + t401;
t487 = qJ(5) * t426 + t509;
t612 = pkin(4) + pkin(5);
t337 = -t424 * t612 + t487;
t446 = t633 * qJD(5);
t545 = qJD(2) * qJD(3);
t524 = t473 * t545;
t454 = qJ(5) * t524;
t563 = qJD(2) * t469;
t531 = t478 * t563;
t559 = qJD(3) * t473;
t361 = -t435 * t559 + (qJD(3) * t470 + t531) * t564;
t590 = t470 * t473;
t452 = qJD(1) * t590;
t402 = t477 * t435 + t452;
t390 = qJD(3) * pkin(9) + t402;
t400 = (t433 + t534) * qJD(2);
t533 = t478 * t565;
t403 = qJD(2) * t438 - t533;
t557 = qJD(4) * t472;
t489 = -t476 * t361 + t390 * t557 - t472 * t400 - t403 * t556;
t322 = t446 + t454 - t489;
t525 = t473 * t556;
t558 = qJD(3) * t477;
t486 = t472 * t558 + t525;
t544 = qJD(3) * qJD(4);
t397 = qJD(2) * t486 + t472 * t544;
t320 = pkin(10) * t397 + t322;
t523 = t477 * t545;
t526 = t473 * t557;
t396 = qJD(2) * t526 + (-t523 - t544) * t476;
t514 = -t472 * t361 - t390 * t556 + t476 * t400 - t403 * t557;
t321 = pkin(10) * t396 - t524 * t612 - t514;
t519 = t471 * t320 - t475 * t321;
t630 = t337 * t368 + t519;
t628 = MDP(5) * t473;
t467 = t473 ^ 2;
t627 = MDP(6) * (-t477 ^ 2 + t467);
t626 = t546 * t496;
t515 = pkin(4) * t524;
t323 = -t514 - t515;
t347 = t476 * t390 + t472 * t403;
t448 = t633 * qJ(5);
t339 = t448 + t347;
t625 = -t339 * t633 + t323;
t583 = t476 * t477;
t460 = pkin(8) * t583;
t624 = -(t472 * t582 - t474 * t476) * t565 + qJD(4) * t460 - t433 * t476 + t438 * t557;
t623 = qJD(5) * t472 + t402;
t622 = t477 * t549 - t526;
t621 = qJ(5) * t559 + t634;
t620 = MDP(28) * t471 + MDP(29) * t475;
t346 = -t472 * t390 + t476 * t403;
t547 = qJD(5) - t346;
t619 = MDP(17) + MDP(19);
t618 = MDP(18) - MDP(21);
t548 = pkin(10) * t426 - t547;
t329 = -t612 * t633 - t548;
t552 = qJD(6) * t475;
t538 = t475 * t320 + t471 * t321 + t329 * t552;
t616 = -t337 * t496 + t538;
t607 = qJ(5) * t472;
t615 = -t476 * t612 - t607;
t613 = t426 ^ 2;
t611 = pkin(9) - pkin(10);
t610 = pkin(10) * t473;
t609 = qJD(2) * pkin(2);
t608 = qJ(5) * t424;
t606 = qJ(5) * t476;
t513 = t473 * t531;
t362 = qJD(1) * t513 + qJD(3) * t452 + t435 * t558;
t484 = -qJ(5) * t396 + qJD(5) * t426 - t362;
t325 = pkin(4) * t397 - t484;
t605 = t325 * t472;
t604 = t325 * t476;
t348 = pkin(4) * t424 - t487;
t602 = t348 * t426;
t601 = t362 * t472;
t600 = t362 * t476;
t599 = t509 * t472;
t598 = t396 * t472;
t597 = t424 * t426;
t596 = t424 * t633;
t595 = t426 * t633;
t594 = t633 * t472;
t593 = t633 * t476;
t592 = t469 * t474;
t591 = t469 * t478;
t336 = pkin(10) * t424 + t347;
t330 = t336 + t448;
t589 = t471 * t330;
t588 = t471 * t476;
t587 = t472 * t475;
t586 = t472 * t477;
t585 = t473 * t476;
t480 = qJD(3) ^ 2;
t584 = t473 * t480;
t581 = t477 * t480;
t535 = -pkin(8) * t472 - pkin(4);
t541 = pkin(10) * t583;
t580 = pkin(10) * t526 + (-t541 + (-pkin(5) + t535) * t473) * qJD(3) + t624;
t579 = -(-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t585 - (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t472) * t477 - t621;
t578 = -qJD(5) * t477 + (-t473 * t549 - t477 * t557) * pkin(8) + t621;
t577 = t535 * t559 + t624;
t501 = pkin(4) * t472 - t606;
t576 = -t501 * t633 + t623;
t428 = t471 * t472 + t475 * t476;
t490 = t428 * t477;
t575 = -qJD(2) * t490 - t428 * t632;
t530 = t472 * t561;
t553 = qJD(6) * t471;
t574 = t471 * t556 + t472 * t552 - t476 * t553 - t561 * t588 + (t530 - t557) * t475;
t430 = t508 * qJD(2);
t573 = t476 * t401 + t472 * t430;
t492 = -t472 * t612 + t606;
t572 = t492 * t633 + t623;
t570 = t472 * t438 + t460;
t554 = qJD(5) * t476;
t551 = t473 * MDP(16);
t543 = pkin(9) * t594;
t542 = pkin(9) * t593;
t540 = pkin(9) * t559;
t539 = pkin(9) * t549;
t444 = t611 * t476;
t537 = -t475 * t396 + t471 * t397 + t424 * t552;
t350 = qJ(5) * t562 + t573;
t532 = t474 * t563;
t527 = t633 * t557;
t521 = qJD(3) * t551;
t518 = -t396 * t471 - t475 * t397;
t517 = -t472 * t401 + t430 * t476;
t459 = pkin(8) * t586;
t516 = t438 * t476 - t459;
t512 = t424 * t533;
t511 = t426 * t533;
t394 = -qJ(5) * t477 + t570;
t436 = -t533 - t609;
t505 = -t436 - t533;
t504 = (-t473 * t612 - t541) * qJD(2) - t517 + t632 * t444;
t443 = t611 * t472;
t503 = pkin(10) * t530 - qJD(6) * t443 + t611 * t557 + t350;
t502 = pkin(4) * t476 + t607;
t500 = qJ(5) * t475 - t471 * t612;
t499 = qJ(5) * t471 + t475 * t612;
t316 = t471 * t329 + t475 * t330;
t338 = -pkin(4) * t633 + t547;
t498 = t338 * t476 - t339 * t472;
t466 = t477 * pkin(4);
t363 = pkin(5) * t477 + t459 + t466 + (-t438 - t610) * t476;
t369 = t472 * t610 + t394;
t497 = t363 * t471 + t369 * t475;
t495 = -t587 + t588;
t494 = qJD(2) * t467 + t477 * t633;
t493 = pkin(8) + t501;
t415 = t477 * t592 + t590;
t378 = t415 * t472 + t476 * t591;
t379 = t415 * t476 - t472 * t591;
t414 = -t470 * t477 + t473 * t592;
t491 = t347 * t633 + t514;
t409 = t428 * t473;
t327 = -t426 * t553 + t537;
t488 = -pkin(8) + t492;
t483 = qJD(3) * (-t505 - t609);
t482 = t346 * t633 + t489;
t328 = qJD(6) * t368 + t518;
t481 = qJD(2) ^ 2;
t437 = -pkin(3) - t502;
t422 = pkin(3) - t615;
t408 = t471 * t585 - t473 * t587;
t407 = t493 * t473;
t395 = t466 - t516;
t393 = t488 * t473;
t377 = qJD(3) * t415 + t513;
t376 = -qJD(3) * t414 + t477 * t531;
t373 = pkin(4) * t426 + t608;
t355 = -t426 * t612 - t608;
t354 = -t396 + t596;
t353 = (qJD(4) * t502 - t554) * t473 + t493 * t558;
t352 = -pkin(4) * t562 - t517;
t343 = qJD(6) * t409 + t471 * t622 - t475 * t486;
t342 = -t473 * t495 * t632 + qJD(3) * t490;
t340 = (qJD(4) * t615 + t554) * t473 + t488 * t558;
t333 = -qJD(4) * t378 + t376 * t476 + t472 * t532;
t332 = qJD(4) * t379 + t376 * t472 - t476 * t532;
t324 = -t397 * t612 + t484;
t315 = t329 * t475 - t589;
t1 = [(t332 * t426 - t333 * t424 - t378 * t396 - t379 * t397) * MDP(20) + (t322 * t379 + t323 * t378 + t325 * t414 + t332 * t338 + t333 * t339 + t348 * t377) * MDP(22) + (-(t332 * t475 - t333 * t471 + (-t378 * t471 - t379 * t475) * qJD(6)) * t546 - t377 * t496 - t414 * t328) * MDP(28) + ((t332 * t471 + t333 * t475 + (t378 * t475 - t379 * t471) * qJD(6)) * t546 - t377 * t368 - t414 * t327) * MDP(29) + (-t377 * MDP(10) - t376 * MDP(11) + ((-t618 + t620) * t379 + (-MDP(28) * t475 + MDP(29) * t471 - t619) * t378) * t562) * qJD(3) + ((-MDP(10) * t473 - MDP(11) * t477) * t478 * t545 + (-t478 * MDP(4) + (-MDP(10) * t477 + MDP(11) * t473 - MDP(3)) * t474) * t481) * t469 + t619 * (-t332 * t633 + t377 * t424 + t414 * t397) + t618 * (-t333 * t633 + t377 * t426 - t396 * t414); -0.2e1 * t545 * t627 + (-pkin(8) * t581 + t473 * t483) * MDP(10) + (pkin(8) * t584 + t477 * t483) * MDP(11) + (-t396 * t585 + t426 * t622) * MDP(12) + ((-t424 * t476 - t426 * t472) * t558 + (t598 - t397 * t476 + (t424 * t472 - t426 * t476) * qJD(4)) * t473) * MDP(13) + (-t633 * t526 + t396 * t477 + (t426 * t473 + t476 * t494) * qJD(3)) * MDP(14) + (-t633 * t525 + t397 * t477 + (-t424 * t473 - t472 * t494) * qJD(3)) * MDP(15) + (t633 - t561) * t521 + (-t624 * t633 + ((pkin(8) * t424 - t599) * qJD(3) - t514) * t477 + (-t512 - t509 * t556 + pkin(8) * t397 + t601 + (pkin(8) * t594 + qJD(2) * t516 + t346) * qJD(3)) * t473) * MDP(17) + (-t634 * t633 + (-t509 * t549 + (qJD(3) * t426 + t527) * pkin(8) - t489) * t477 + (-t511 + t509 * t557 - pkin(8) * t396 + t600 + (pkin(8) * t593 - qJD(2) * t570 - t347) * qJD(3)) * t473) * MDP(18) + (t353 * t424 + t397 * t407 + (t348 * t560 + t323) * t477 - t577 * t633 + (-t512 + t348 * t556 + t605 + (-qJD(2) * t395 - t338) * qJD(3)) * t473) * MDP(19) + (-t394 * t397 - t395 * t396 + t577 * t426 - t578 * t424 + t498 * t558 + (-t322 * t472 + t323 * t476 + (-t338 * t472 - t339 * t476) * qJD(4)) * t473) * MDP(20) + (-t353 * t426 + t396 * t407 + (-t348 * t549 - t322) * t477 + t578 * t633 + (t511 + t348 * t557 - t604 + (qJD(2) * t394 + t339) * qJD(3)) * t473) * MDP(21) + (t322 * t394 + t323 * t395 + t325 * t407 + (-t473 * t533 + t353) * t348 + t578 * t339 + t577 * t338) * MDP(22) + (t327 * t409 + t342 * t368) * MDP(23) + (-t327 * t408 - t328 * t409 - t342 * t496 - t343 * t368) * MDP(24) + (t327 * t477 - t342 * t546 + (-qJD(2) * t409 - t368) * t559) * MDP(25) + (-t328 * t477 + t343 * t546 + (qJD(2) * t408 + t496) * t559) * MDP(26) + (t546 - t561) * t520 + (-t519 * t477 + t340 * t496 + t393 * t328 + t324 * t408 + t337 * t343 - (t471 * t579 + t475 * t580) * t546 + (-t316 * t477 + t497 * t546) * qJD(6) + (t496 * t533 + (-(t363 * t475 - t369 * t471) * qJD(2) - t315) * qJD(3)) * t473) * MDP(28) + (-(-t330 * t553 + t538) * t477 + t340 * t368 + t393 * t327 + t324 * t409 + t337 * t342 - ((-qJD(6) * t363 + t579) * t475 + (qJD(6) * t369 - t580) * t471) * t546 + (t368 * t533 + (qJD(2) * t497 + t316) * qJD(3)) * t473) * MDP(29) + MDP(7) * t581 - MDP(8) * t584 + 0.2e1 * t523 * t628; (qJD(3) * t402 - t436 * t562 - t362) * MDP(10) + t505 * t561 * MDP(11) + (t426 * t593 - t598) * MDP(12) + ((-t396 - t596) * t476 + (-t397 - t595) * t472) * MDP(13) + t633 * t556 * MDP(14) - t527 * MDP(15) + (-pkin(3) * t397 - t600 - t517 * t633 - t402 * t424 + (-t542 - t599) * qJD(4)) * MDP(17) + (pkin(3) * t396 + t601 + t573 * t633 - t402 * t426 + (-t476 * t509 + t543) * qJD(4)) * MDP(18) + (-t604 + t352 * t633 + t397 * t437 - t576 * t424 + (t348 * t472 - t542) * qJD(4)) * MDP(19) + (t350 * t424 - t352 * t426 + (t322 + t633 * t338 + (qJD(4) * t426 - t397) * pkin(9)) * t476 + ((qJD(4) * t424 - t396) * pkin(9) + t625) * t472) * MDP(20) + (-t605 - t350 * t633 + t396 * t437 + t576 * t426 + (-t348 * t476 - t543) * qJD(4)) * MDP(21) + (t325 * t437 - t338 * t352 - t339 * t350 - t576 * t348 + (qJD(4) * t498 + t322 * t476 + t323 * t472) * pkin(9)) * MDP(22) + (-t327 * t495 + t368 * t575) * MDP(23) + (-t327 * t428 + t328 * t495 - t368 * t574 - t496 * t575) * MDP(24) + (-t575 * t546 + (qJD(3) * t495 + t368) * t562) * MDP(25) + (t574 * t546 + (qJD(3) * t428 - t496) * t562) * MDP(26) + (t324 * t428 + t422 * t328 - (t471 * t503 - t475 * t504) * t546 + t572 * t496 + t574 * t337 + (-(t443 * t475 - t444 * t471) * qJD(3) + t315) * t562) * MDP(28) + (-t324 * t495 + t422 * t327 - (t471 * t504 + t475 * t503) * t546 + t572 * t368 + t575 * t337 + ((t443 * t471 + t444 * t475) * qJD(3) - t316) * t562) * MDP(29) + (-t477 * t628 + t627) * t481 + ((-t633 * t583 + (-t426 + t560) * t473) * MDP(14) + (t633 * t586 + (t424 + t549) * t473) * MDP(15) - t633 * t551 + (-t346 * t473 + (t477 * t509 - t540) * t472) * MDP(17) + (t509 * t583 + (t347 - t539) * t473) * MDP(18) + (t338 * t473 + (-t348 * t477 - t540) * t472) * MDP(19) + (t348 * t583 + (-t339 + t539) * t473) * MDP(21) - t546 * t550) * qJD(2); MDP(12) * t597 + (-t424 ^ 2 + t613) * MDP(13) + t354 * MDP(14) + (-t397 + t595) * MDP(15) + qJD(2) * t521 + (t426 * t509 + t491) * MDP(17) + (-t424 * t509 + t482) * MDP(18) + (-t373 * t424 + t491 + 0.2e1 * t515 - t602) * MDP(19) + (pkin(4) * t396 - qJ(5) * t397 + (t339 - t347) * t426 + (t338 - t547) * t424) * MDP(20) + (-t348 * t424 + t373 * t426 + 0.2e1 * t446 + 0.2e1 * t454 - t482) * MDP(21) + (-pkin(4) * t323 + qJ(5) * t322 - t338 * t347 + t339 * t547 - t348 * t373) * MDP(22) + (-t327 + t626) * MDP(25) + (t328 + t635) * MDP(26) + (t499 * t524 - t355 * t496 - (-t475 * t336 + t471 * t548) * t546 + (t500 * t546 + t316) * qJD(6) + t630) * MDP(28) + (t500 * t524 - t355 * t368 - (t336 * t471 + t475 * t548) * t546 + (-t499 * t546 - t589) * qJD(6) + t616) * MDP(29) - t636; (-t524 + t597) * MDP(19) + t354 * MDP(20) + (-t633 ^ 2 - t613) * MDP(21) + (t602 + t625) * MDP(22) + (-t426 * t496 - t475 * t524) * MDP(28) + (-t368 * t426 + t471 * t524) * MDP(29) - t620 * t546 ^ 2; (t537 - t626) * MDP(25) + (-t518 - t635) * MDP(26) + (-t316 * t546 - t630) * MDP(28) + (-t315 * t546 - t616) * MDP(29) + ((-MDP(26) * t426 - MDP(28) * t330) * t475 + (-MDP(25) * t426 - MDP(26) * t424 - MDP(28) * t329 + MDP(29) * t330) * t471) * qJD(6) + t636;];
tauc  = t1;
