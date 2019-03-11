% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:25
% EndTime: 2019-03-09 11:45:37
% DurationCPUTime: 7.68s
% Computational Cost: add. (10130->479), mult. (25928->599), div. (0->0), fcn. (19779->8), ass. (0->216)
t509 = cos(qJ(5));
t507 = sin(qJ(4));
t504 = sin(pkin(10));
t505 = cos(pkin(10));
t508 = sin(qJ(2));
t510 = cos(qJ(2));
t483 = t504 * t510 + t505 * t508;
t473 = t483 * qJD(2);
t518 = qJD(1) * t473;
t559 = qJD(1) * qJD(2);
t548 = t510 * t559;
t549 = t508 * t559;
t463 = -t504 * t549 + t505 * t548;
t611 = cos(qJ(4));
t554 = t611 * t463;
t515 = -t507 * t518 + t554;
t482 = -t504 * t508 + t505 * t510;
t472 = t482 * qJD(1);
t474 = t483 * qJD(1);
t617 = t611 * t472 - t507 * t474;
t632 = t617 * qJD(4);
t514 = t515 + t632;
t527 = -t507 * t472 - t474 * t611;
t558 = qJD(2) + qJD(4);
t539 = qJD(5) * t558;
t506 = sin(qJ(5));
t567 = qJD(5) * t506;
t364 = -t527 * t567 + (-t514 - t539) * t509;
t362 = t364 * t506;
t418 = t506 * t558 - t509 * t527;
t428 = qJD(5) - t617;
t542 = t507 * t463 + t611 * t518;
t633 = qJD(4) * t527;
t392 = t542 - t633;
t389 = t506 * t392;
t566 = qJD(5) * t509;
t574 = t428 * t566 + t389;
t636 = t617 * t509;
t530 = -t428 * t636 + t574;
t513 = t506 * t514;
t568 = qJD(5) * t418;
t365 = t513 + t568;
t415 = -t506 * t527 - t509 * t558;
t644 = t566 - t636;
t531 = -t506 * t365 - t644 * t415;
t563 = t617 * qJD(2);
t564 = t527 * qJD(2);
t595 = t418 * t527;
t642 = t428 * t506;
t626 = -t364 * t509 - t418 * t642;
t646 = (t531 + t626) * MDP(21) + (t530 + t595) * MDP(22) - t617 ^ 2 * MDP(14) + (MDP(13) * t617 + MDP(14) * t527 + MDP(24) * t428) * t527 + (-t564 - t542) * MDP(16) + (-t563 + t515) * MDP(15) + (t644 * t418 - t362) * MDP(20);
t391 = t509 * t392;
t643 = -t428 * t642 + t391;
t638 = pkin(5) * t527;
t637 = qJ(6) * t527;
t598 = t415 * t527;
t606 = -qJ(3) - pkin(7);
t547 = qJD(2) * t606;
t468 = qJD(3) * t510 + t508 * t547;
t452 = t468 * qJD(1);
t469 = -qJD(3) * t508 + t510 * t547;
t453 = t469 * qJD(1);
t414 = -t452 * t504 + t505 * t453;
t400 = -pkin(8) * t463 + t414;
t417 = t505 * t452 + t504 * t453;
t401 = -pkin(8) * t518 + t417;
t492 = t606 * t510;
t488 = qJD(1) * t492;
t477 = t504 * t488;
t491 = t606 * t508;
t487 = qJD(1) * t491;
t605 = qJD(2) * pkin(2);
t481 = t487 + t605;
t435 = t505 * t481 + t477;
t607 = pkin(8) * t474;
t410 = qJD(2) * pkin(3) + t435 - t607;
t587 = t505 * t488;
t436 = t504 * t481 - t587;
t608 = pkin(8) * t472;
t413 = t436 + t608;
t550 = qJD(4) * t611;
t569 = qJD(4) * t507;
t338 = -t611 * t400 + t507 * t401 + t410 * t569 + t413 * t550;
t329 = pkin(5) * t365 + qJ(6) * t364 - qJD(6) * t418 + t338;
t373 = t507 * t410 + t611 * t413;
t370 = pkin(9) * t558 + t373;
t556 = -pkin(2) * t510 - pkin(1);
t537 = t556 * qJD(1);
t489 = qJD(3) + t537;
t441 = -pkin(3) * t472 + t489;
t377 = -pkin(4) * t617 + pkin(9) * t527 + t441;
t343 = -t370 * t506 + t377 * t509;
t560 = qJD(6) - t343;
t332 = -pkin(5) * t428 + t560;
t372 = t410 * t611 - t507 * t413;
t369 = -pkin(4) * t558 - t372;
t347 = t415 * pkin(5) - t418 * qJ(6) + t369;
t631 = -t329 * t509 - t332 * t527 + t347 * t567;
t344 = t370 * t509 + t377 * t506;
t333 = qJ(6) * t428 + t344;
t603 = t329 * t506;
t630 = t333 * t527 - t603;
t629 = t441 * t527 - t338;
t628 = -t338 * t509 + t343 * t527 + t369 * t567;
t627 = t338 * t506 - t344 * t527 + t369 * t566;
t517 = -t507 * t400 - t401 * t611 - t410 * t550 + t413 * t569;
t625 = -t441 * t617 + t517;
t535 = pkin(5) * t506 - qJ(6) * t509;
t624 = -pkin(5) * t567 + qJ(6) * t566 + t506 * qJD(6) + t535 * t617;
t395 = -pkin(4) * t527 - pkin(9) * t617;
t622 = -0.2e1 * t559;
t621 = MDP(4) * t508;
t620 = MDP(5) * (t508 ^ 2 - t510 ^ 2);
t496 = pkin(2) * t549;
t585 = t507 * t483;
t355 = t496 - (t472 * t550 - t474 * t569 + t554) * pkin(9) + t392 * pkin(4) + (pkin(3) * t483 + pkin(9) * t585) * t559;
t522 = t506 * t355 - t370 * t567 + t377 * t566 - t509 * t517;
t604 = qJ(6) * t392;
t324 = qJD(6) * t428 + t522 + t604;
t538 = t509 * t355 - t370 * t566 - t377 * t567 + t506 * t517;
t609 = pkin(5) * t392;
t326 = -t538 - t609;
t619 = t324 * t509 + t326 * t506;
t438 = t507 * t482 + t483 * t611;
t454 = -pkin(3) * t482 + t556;
t526 = t482 * t611 - t585;
t386 = -pkin(4) * t526 - pkin(9) * t438 + t454;
t444 = t505 * t491 + t492 * t504;
t426 = -pkin(8) * t483 + t444;
t445 = t504 * t491 - t505 * t492;
t427 = pkin(8) * t482 + t445;
t388 = t507 * t426 + t427 * t611;
t575 = t506 * t386 + t509 * t388;
t439 = -t487 * t504 + t587;
t420 = t439 - t608;
t440 = t505 * t487 + t477;
t421 = t440 - t607;
t380 = t507 * t420 + t421 * t611;
t498 = pkin(2) * t505 + pkin(3);
t610 = pkin(2) * t504;
t615 = t611 * t498 - t507 * t610;
t458 = t615 * qJD(4);
t618 = -t458 + t380;
t520 = t507 * t498 + t610 * t611;
t572 = t520 * qJD(4) + t420 * t611 - t507 * t421;
t616 = t428 * t567 - t391;
t614 = t611 * t426 - t507 * t427;
t613 = t418 ^ 2;
t612 = t428 ^ 2;
t602 = t333 * t428;
t601 = t344 * t428;
t600 = t365 * t509;
t467 = pkin(9) + t520;
t599 = t392 * t467;
t597 = t415 * t506;
t596 = t418 * t415;
t594 = t418 * t509;
t592 = t617 * t506;
t589 = t438 * t509;
t588 = t458 * t509;
t511 = qJD(2) ^ 2;
t584 = t508 * t511;
t583 = t510 * t511;
t512 = qJD(1) ^ 2;
t582 = t510 * t512;
t580 = t373 + t624;
t579 = -t572 + t624;
t577 = t509 * t372 + t506 * t395;
t570 = qJD(1) * t508;
t449 = pkin(2) * t570 + pkin(3) * t474;
t381 = t395 + t449;
t576 = t509 * t380 + t506 * t381;
t425 = t505 * t468 + t504 * t469;
t501 = t508 * t605;
t557 = t332 * t566 + t619;
t551 = t467 * t567;
t450 = pkin(3) * t473 + t501;
t546 = pkin(1) * t622;
t340 = t380 * t506 - t381 * t509 + t638;
t543 = t458 * t506 - t340;
t424 = -t468 * t504 + t505 * t469;
t536 = t509 * pkin(5) + t506 * qJ(6);
t534 = t332 * t509 - t333 * t506;
t533 = -t369 * t617 - t599;
t532 = -t372 * t506 + t395 * t509;
t490 = -pkin(4) - t536;
t529 = t617 * t642 - t616;
t476 = t482 * qJD(2);
t396 = qJD(4) * t526 - t507 * t473 + t476 * t611;
t525 = t396 * t506 + t438 * t566;
t524 = -t396 * t509 + t438 * t567;
t523 = t347 * t418 - t538;
t403 = -pkin(8) * t476 + t424;
t404 = -pkin(8) * t473 + t425;
t351 = qJD(4) * t614 + t507 * t403 + t611 * t404;
t397 = qJD(4) * t438 + t473 * t611 + t507 * t476;
t359 = pkin(4) * t397 - pkin(9) * t396 + t450;
t521 = t509 * t351 + t506 * t359 + t386 * t566 - t388 * t567;
t519 = t574 * pkin(9);
t516 = qJD(5) * t534 + t619;
t352 = qJD(4) * t388 - t403 * t611 + t507 * t404;
t466 = -pkin(4) - t615;
t443 = t490 - t615;
t442 = pkin(3) * t518 + t496;
t382 = pkin(5) * t418 + qJ(6) * t415;
t360 = t438 * t535 - t614;
t349 = pkin(5) * t526 - t386 * t509 + t388 * t506;
t348 = -qJ(6) * t526 + t575;
t345 = t415 * t428 - t364;
t342 = -t532 + t638;
t341 = t577 - t637;
t339 = t576 - t637;
t330 = t535 * t396 + (qJD(5) * t536 - qJD(6) * t509) * t438 + t352;
t328 = -pkin(5) * t397 + qJD(5) * t575 + t351 * t506 - t359 * t509;
t327 = qJ(6) * t397 - qJD(6) * t526 + t521;
t1 = [(-t438 * t392 + t396 * t617 + t397 * t527 + t514 * t526) * MDP(14) + (-t396 * t527 + t438 * t514) * MDP(13) + (-t538 * t526 + t343 * t397 + t352 * t415 - t614 * t365 + ((-qJD(5) * t388 + t359) * t428 + t386 * t392 + t369 * qJD(5) * t438) * t509 + ((-qJD(5) * t386 - t351) * t428 - t388 * t392 + t338 * t438 + t369 * t396) * t506) * MDP(25) + (-t327 * t415 + t328 * t418 - t348 * t365 - t349 * t364 + t534 * t396 + (-t324 * t506 + t326 * t509 + (-t332 * t506 - t333 * t509) * qJD(5)) * t438) * MDP(28) + MDP(6) * t583 + (-t414 * t483 + t417 * t482 - t424 * t474 + t425 * t472 - t435 * t476 - t436 * t473 - t444 * t463 - t445 * t518) * MDP(11) + ((-t415 * t509 - t418 * t506) * t396 + (t362 - t600 + (-t594 + t597) * qJD(5)) * t438) * MDP(21) + (t326 * t526 - t328 * t428 + t330 * t415 - t332 * t397 + t347 * t525 - t349 * t392 + t360 * t365 + t438 * t603) * MDP(27) + t620 * t622 + 0.2e1 * t548 * t621 + (-pkin(7) * t583 + t508 * t546) * MDP(9) + (t364 * t526 + t391 * t438 + t397 * t418 - t428 * t524) * MDP(22) - MDP(7) * t584 + (pkin(7) * t584 + t510 * t546) * MDP(10) + (-t324 * t526 + t327 * t428 - t329 * t589 - t330 * t418 + t333 * t397 + t347 * t524 + t348 * t392 + t360 * t364) * MDP(29) + (t338 * t589 - t344 * t397 + t352 * t418 + t364 * t614 - t369 * t524 - t392 * t575 - t428 * t521 + t522 * t526) * MDP(26) + (t365 * t526 - t389 * t438 - t397 * t415 - t428 * t525) * MDP(23) + (-t364 * t589 - t418 * t524) * MDP(20) + (t414 * t444 + t417 * t445 + t435 * t424 + t436 * t425 + (t489 + t537) * t501) * MDP(12) + (t454 * t392 + t441 * t397 - t442 * t526 - t450 * t617) * MDP(18) + (t441 * t396 + t442 * t438 - t450 * t527 + t454 * t514) * MDP(19) + (t324 * t348 + t326 * t349 + t327 * t333 + t328 * t332 + t329 * t360 + t330 * t347) * MDP(30) + (-t392 * t526 + t397 * t428) * MDP(24) + (MDP(15) * t396 - t397 * MDP(16) - t352 * MDP(18) - MDP(19) * t351) * t558; -t582 * t621 + (t449 * t527 + t558 * t618 + t625) * MDP(19) + (-t435 * t439 - t436 * t440 + (t414 * t505 + t417 * t504 - t489 * t570) * pkin(2)) * MDP(12) + (t339 * t415 - t340 * t418 + (-t332 * t617 - t415 * t458 + (-t365 + t568) * t467) * t509 + (t333 * t617 - t364 * t467 + t418 * t458 + (t415 * t467 - t333) * qJD(5)) * t506 + t557) * MDP(28) + (-t466 * t364 + t533 * t509 + t572 * t418 + (t551 + t576 - t588) * t428 + t627) * MDP(26) + (t466 * t365 + t533 * t506 + t572 * t415 + ((-qJD(5) * t467 - t381) * t509 + t618 * t506) * t428 + t628) * MDP(25) + (t449 * t617 - t558 * t572 + t629) * MDP(18) + (t364 * t443 + t579 * t418 + t599 * t509 + (-t339 - t551 + (-t347 + t458) * t509) * t428 + t630) * MDP(29) + (t365 * t443 + (-t347 * t617 - t599) * t506 - t579 * t415 + (-t467 * t566 - t543) * t428 + t631) * MDP(27) + (t329 * t443 - t579 * t347 + (-t339 + t588) * t333 + t543 * t332 + t516 * t467) * MDP(30) + (t529 - t598) * MDP(23) + ((t436 + t439) * t474 + (-t440 + t435) * t472 + (-t505 * t463 - t504 * t518) * pkin(2)) * MDP(11) + t512 * t620 + (MDP(9) * t508 * t512 + MDP(10) * t582) * pkin(1) + t646; (-t472 ^ 2 - t474 ^ 2) * MDP(11) + (t435 * t474 - t436 * t472 + t496) * MDP(12) + (-t564 + t542 - 0.2e1 * t633) * MDP(18) + (t563 + t515 + 0.2e1 * t632) * MDP(19) + (t529 + t598) * MDP(25) + (-t509 * t612 - t389 + t595) * MDP(26) + (t598 + t643) * MDP(27) + (t531 - t626) * MDP(28) + (t530 - t595) * MDP(29) + (t347 * t527 + (-t326 + t602) * t509 + (t332 * t428 + t324) * t506) * MDP(30); (t373 * t558 + t629) * MDP(18) + (t372 * t558 + t625) * MDP(19) + (-t598 + t643) * MDP(23) + (-pkin(4) * t365 - t369 * t592 - t373 * t415 - t428 * t532 - t519 + t628) * MDP(25) + (pkin(4) * t364 + pkin(9) * t616 - t369 * t636 - t373 * t418 + t577 * t428 + t627) * MDP(26) + (t342 * t428 - t347 * t592 + t365 * t490 - t415 * t580 - t519 + t631) * MDP(27) + (-t332 * t636 + t341 * t415 - t342 * t418 - t333 * t642 + (-t362 - t600 + (t594 + t597) * qJD(5)) * pkin(9) + t557) * MDP(28) + (t364 * t490 + (-pkin(9) * t567 - t341) * t428 + t580 * t418 + (pkin(9) * t392 - t347 * t428) * t509 + t630) * MDP(29) + (pkin(9) * t516 + t329 * t490 - t332 * t342 - t333 * t341 - t347 * t580) * MDP(30) + t646; MDP(20) * t596 + (-t415 ^ 2 + t613) * MDP(21) + t345 * MDP(22) + (t418 * t428 - t506 * t539 + t527 * t566 - t513) * MDP(23) + t392 * MDP(24) + (-t369 * t418 + t538 + t601) * MDP(25) + (t343 * t428 + t369 * t415 - t522) * MDP(26) + (-t382 * t415 - t523 + t601 + 0.2e1 * t609) * MDP(27) + (pkin(5) * t364 - qJ(6) * t365 + (t333 - t344) * t418 + (t332 - t560) * t415) * MDP(28) + (0.2e1 * t604 - t347 * t415 + t382 * t418 + (0.2e1 * qJD(6) - t343) * t428 + t522) * MDP(29) + (-pkin(5) * t326 + qJ(6) * t324 - t332 * t344 + t333 * t560 - t347 * t382) * MDP(30); (-t392 + t596) * MDP(27) + t345 * MDP(28) + (-t612 - t613) * MDP(29) + (t523 - t602 - t609) * MDP(30);];
tauc  = t1;
