% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:21
% EndTime: 2019-03-09 10:19:34
% DurationCPUTime: 7.84s
% Computational Cost: add. (8109->470), mult. (20734->635), div. (0->0), fcn. (15929->10), ass. (0->220)
t555 = sin(pkin(10));
t560 = sin(qJ(2));
t608 = qJD(1) * t560;
t557 = cos(pkin(10));
t563 = cos(qJ(2));
t627 = t557 * t563;
t521 = qJD(1) * t627 - t555 * t608;
t515 = qJD(4) - t521;
t561 = cos(qJ(6));
t536 = t555 * t563 + t557 * t560;
t524 = t536 * qJD(1);
t559 = sin(qJ(4));
t562 = cos(qJ(4));
t603 = t562 * qJD(2);
t490 = t524 * t559 - t603;
t492 = qJD(2) * t559 + t524 * t562;
t554 = sin(pkin(11));
t556 = cos(pkin(11));
t579 = -t490 * t556 - t492 * t554;
t623 = t561 * t579;
t431 = t490 * t554 - t492 * t556;
t558 = sin(qJ(6));
t643 = t431 * t558;
t388 = t623 + t643;
t507 = qJD(6) + t515;
t645 = t388 * t507;
t602 = qJD(1) * qJD(2);
t594 = t563 * t602;
t595 = t560 * t602;
t509 = -t555 * t595 + t557 * t594;
t606 = qJD(4) * t559;
t446 = qJD(4) * t603 + t562 * t509 - t524 * t606;
t522 = t536 * qJD(2);
t508 = qJD(1) * t522;
t597 = -pkin(2) * t563 - pkin(1);
t584 = t597 * qJD(1);
t540 = qJD(3) + t584;
t451 = -pkin(3) * t521 - pkin(8) * t524 + t540;
t649 = -qJ(3) - pkin(7);
t541 = t649 * t560;
t538 = qJD(1) * t541;
t648 = qJD(2) * pkin(2);
t530 = t538 + t648;
t542 = t649 * t563;
t539 = qJD(1) * t542;
t628 = t557 * t539;
t477 = t555 * t530 - t628;
t470 = qJD(2) * pkin(8) + t477;
t413 = t451 * t559 + t470 * t562;
t545 = pkin(2) * t595;
t450 = pkin(3) * t508 - pkin(8) * t509 + t545;
t442 = t562 * t450;
t593 = qJD(2) * t649;
t517 = qJD(3) * t563 + t560 * t593;
t500 = t517 * qJD(1);
t518 = -qJD(3) * t560 + t563 * t593;
t501 = t518 * qJD(1);
t445 = t500 * t557 + t501 * t555;
t566 = -qJD(4) * t413 - t445 * t559 + t442;
t355 = pkin(4) * t508 - qJ(5) * t446 - qJD(5) * t492 + t566;
t447 = t492 * qJD(4) + t509 * t559;
t605 = qJD(4) * t562;
t568 = t562 * t445 + t559 * t450 + t451 * t605 - t470 * t606;
t358 = -qJ(5) * t447 - qJD(5) * t490 + t568;
t340 = t556 * t355 - t358 * t554;
t401 = t446 * t556 - t447 * t554;
t336 = pkin(5) * t508 - pkin(9) * t401 + t340;
t341 = t554 * t355 + t556 * t358;
t400 = -t446 * t554 - t447 * t556;
t337 = pkin(9) * t400 + t341;
t412 = t562 * t451 - t470 * t559;
t395 = -qJ(5) * t492 + t412;
t381 = pkin(4) * t515 + t395;
t396 = -qJ(5) * t490 + t413;
t630 = t556 * t396;
t363 = t554 * t381 + t630;
t656 = pkin(9) * t579;
t350 = t363 + t656;
t604 = qJD(6) * t558;
t349 = t350 * t604;
t527 = t555 * t539;
t476 = t530 * t557 + t527;
t469 = -qJD(2) * pkin(3) - t476;
t427 = pkin(4) * t490 + qJD(5) + t469;
t382 = -pkin(5) * t579 + t427;
t669 = -t558 * t336 - t561 * t337 - t382 * t388 + t349;
t658 = -t561 * t431 + t558 * t579;
t668 = t508 * MDP(26) + (-t388 ^ 2 + t658 ^ 2) * MDP(23) - t388 * MDP(22) * t658;
t535 = t554 * t562 + t556 * t559;
t661 = t515 * t535;
t577 = t554 * t559 - t556 * t562;
t660 = t515 * t577;
t647 = t658 * t507;
t462 = pkin(2) * t608 + pkin(3) * t524 - pkin(8) * t521;
t457 = t562 * t462;
t482 = t538 * t557 + t527;
t547 = pkin(2) * t555 + pkin(8);
t620 = qJ(5) + t547;
t589 = qJD(4) * t620;
t666 = -pkin(4) * t524 - t457 + (qJ(5) * t521 - t589) * t562 + (-qJD(5) + t482) * t559;
t613 = t559 * t462 + t562 * t482;
t635 = t521 * t559;
t665 = -qJ(5) * t635 - qJD(5) * t562 + t559 * t589 + t613;
t664 = t606 - t635;
t591 = t561 * t336 - t558 * t337;
t663 = -t382 * t658 + t591;
t662 = pkin(9) * t431;
t479 = t535 * t561 - t558 * t577;
t618 = qJD(6) * t479 - t660 * t558 + t561 * t661;
t534 = t555 * t560 - t627;
t526 = t534 * qJD(2);
t596 = t536 * t605;
t659 = -t526 * t559 + t596;
t657 = -0.2e1 * t602;
t655 = MDP(4) * t560;
t654 = MDP(5) * (t560 ^ 2 - t563 ^ 2);
t617 = t554 * t665 + t556 * t666;
t616 = t554 * t666 - t556 * t665;
t481 = t538 * t555 - t628;
t652 = pkin(4) * t664 - t481;
t578 = -t535 * t558 - t561 * t577;
t619 = qJD(6) * t578 - t558 * t661 - t561 * t660;
t651 = -t479 * t508 - t507 * t619;
t590 = -t561 * t400 + t401 * t558;
t347 = qJD(6) * t658 + t590;
t650 = pkin(4) * t554;
t646 = t388 * t524;
t644 = t658 * t524;
t642 = t446 * t559;
t640 = t490 * t515;
t639 = t490 * t524;
t638 = t492 * t515;
t637 = t492 * t524;
t633 = t536 * t559;
t632 = t536 * t562;
t390 = t554 * t396;
t626 = t559 * t508;
t564 = qJD(2) ^ 2;
t625 = t560 * t564;
t362 = t556 * t381 - t390;
t348 = pkin(5) * t515 + t362 + t662;
t624 = t561 * t348;
t485 = t541 * t555 - t542 * t557;
t483 = t562 * t485;
t498 = t562 * t508;
t622 = t563 * t564;
t565 = qJD(1) ^ 2;
t621 = t563 * t565;
t601 = t560 * t648;
t463 = pkin(3) * t522 + pkin(8) * t526 + t601;
t458 = t562 * t463;
t461 = t517 * t557 + t518 * t555;
t475 = pkin(3) * t534 - pkin(8) * t536 + t597;
t576 = qJ(5) * t526 - qJD(5) * t536;
t370 = pkin(4) * t522 - t461 * t559 + t458 + t576 * t562 + (-t483 + (qJ(5) * t536 - t475) * t559) * qJD(4);
t598 = t562 * t461 + t559 * t463 + t475 * t605;
t374 = -qJ(5) * t596 + (-qJD(4) * t485 + t576) * t559 + t598;
t345 = t554 * t370 + t556 * t374;
t366 = t556 * t395 - t390;
t468 = t562 * t475;
t405 = pkin(4) * t534 - qJ(5) * t632 - t485 * t559 + t468;
t612 = t559 * t475 + t483;
t416 = -qJ(5) * t633 + t612;
t377 = t554 * t405 + t556 * t416;
t611 = pkin(5) * t661 + t652;
t531 = t620 * t559;
t532 = t620 * t562;
t472 = -t554 * t531 + t556 * t532;
t607 = qJD(4) * t536;
t599 = qJD(6) * t623 + t558 * t400 + t561 * t401;
t549 = -pkin(2) * t557 - pkin(3);
t592 = pkin(1) * t657;
t344 = t556 * t370 - t374 * t554;
t365 = -t395 * t554 - t630;
t376 = t556 * t405 - t416 * t554;
t444 = t500 * t555 - t557 * t501;
t460 = t517 * t555 - t557 * t518;
t471 = -t556 * t531 - t532 * t554;
t484 = -t557 * t541 - t542 * t555;
t588 = t515 * t562;
t587 = -t507 * t618 + t578 * t508;
t585 = pkin(4) * t633 + t484;
t440 = -pkin(9) * t577 + t472;
t583 = pkin(5) * t524 - pkin(9) * t660 + qJD(6) * t440 - t617;
t439 = -pkin(9) * t535 + t471;
t582 = pkin(9) * t661 - qJD(6) * t439 - t616;
t339 = t558 * t348 + t561 * t350;
t581 = t444 * t536 - t485 * t508;
t465 = t535 * t536;
t466 = t577 * t536;
t580 = -t561 * t465 + t466 * t558;
t420 = -t465 * t558 - t466 * t561;
t575 = -pkin(4) * t562 + t549;
t574 = pkin(4) * t659 + t460;
t573 = -t515 * t664 + t498;
t548 = pkin(4) * t556 + pkin(5);
t572 = t548 * t558 + t561 * t650;
t571 = t548 * t561 - t558 * t650;
t407 = pkin(4) * t447 + t444;
t569 = -t526 * t562 - t536 * t606;
t346 = t431 * t604 + t599;
t567 = t469 * t515 - t547 * t508;
t489 = pkin(5) * t577 + t575;
t480 = t508 * t534;
t422 = pkin(5) * t465 + t585;
t421 = pkin(4) * t492 - pkin(5) * t431;
t418 = -t526 * t577 + t535 * t607;
t417 = t526 * t535 + t577 * t607;
t378 = -pkin(5) * t417 + t574;
t375 = -pkin(5) * t400 + t407;
t369 = -pkin(9) * t465 + t377;
t364 = pkin(5) * t534 + pkin(9) * t466 + t376;
t361 = qJD(6) * t420 - t561 * t417 - t418 * t558;
t360 = qJD(6) * t580 + t417 * t558 - t418 * t561;
t352 = t366 + t662;
t351 = t365 - t656;
t343 = pkin(9) * t417 + t345;
t342 = pkin(5) * t522 + pkin(9) * t418 + t344;
t338 = -t350 * t558 + t624;
t1 = [t654 * t657 + (-pkin(7) * t622 + t560 * t592) * MDP(9) + (pkin(7) * t625 + t563 * t592) * MDP(10) + (-t445 * t534 + t460 * t524 + t461 * t521 + t476 * t526 - t477 * t522 + t484 * t509 + t581) * MDP(11) + (t444 * t484 + t445 * t485 - t476 * t460 + t477 * t461 + (t540 + t584) * t601) * MDP(12) + (t446 * t632 + t492 * t569) * MDP(13) + (-(-t490 * t562 - t492 * t559) * t526 + (-t642 - t447 * t562 + (t490 * t559 - t492 * t562) * qJD(4)) * t536) * MDP(14) + (t446 * t534 + t492 * t522 + t498 * t536 + t515 * t569) * MDP(15) + (-t447 * t534 - t490 * t522 - t515 * t659 - t536 * t626) * MDP(16) + (t515 * t522 + t480) * MDP(17) + ((-t485 * t605 + t458) * t515 + t468 * t508 + (-t470 * t605 + t442) * t534 + t412 * t522 + t460 * t490 + t484 * t447 + t469 * t596 + ((-qJD(4) * t475 - t461) * t515 + (-qJD(4) * t451 - t445) * t534 - t469 * t526 + t581) * t559) * MDP(18) + (-(-t485 * t606 + t598) * t515 - t612 * t508 - t568 * t534 - t413 * t522 + t460 * t492 + t484 * t446 + t444 * t632 + t569 * t469) * MDP(19) + (t340 * t466 - t341 * t465 + t344 * t431 + t345 * t579 + t362 * t418 + t363 * t417 - t376 * t401 + t377 * t400) * MDP(20) + (t340 * t376 + t341 * t377 + t362 * t344 + t363 * t345 + t407 * t585 + t427 * t574) * MDP(21) + (t346 * t420 + t360 * t658) * MDP(22) + (t346 * t580 - t347 * t420 + t360 * t388 - t361 * t658) * MDP(23) + (t346 * t534 + t360 * t507 + t420 * t508 + t522 * t658) * MDP(24) + (-t347 * t534 - t361 * t507 + t388 * t522 + t508 * t580) * MDP(25) + (t507 * t522 + t480) * MDP(26) + ((t342 * t561 - t343 * t558) * t507 + (t364 * t561 - t369 * t558) * t508 + t591 * t534 + t338 * t522 - t378 * t388 + t422 * t347 - t375 * t580 + t382 * t361 + ((-t364 * t558 - t369 * t561) * t507 - t339 * t534) * qJD(6)) * MDP(27) + (-t339 * t522 + t422 * t346 + t349 * t534 + t382 * t360 + t375 * t420 + t378 * t658 + (-(-qJD(6) * t369 + t342) * t507 - t364 * t508 - t336 * t534) * t558 + (-(qJD(6) * t364 + t343) * t507 - t369 * t508 - (qJD(6) * t348 + t337) * t534) * t561) * MDP(28) + MDP(6) * t622 - MDP(7) * t625 + 0.2e1 * t594 * t655; -t621 * t655 + t565 * t654 + ((t476 - t482) * t521 + (-t508 * t555 - t509 * t557) * pkin(2)) * MDP(11) + (t476 * t481 - t477 * t482 + (-t444 * t557 + t445 * t555 - t540 * t608) * pkin(2)) * MDP(12) + (t492 * t588 + t642) * MDP(13) + ((t446 - t640) * t562 + (-t447 - t638) * t559) * MDP(14) + (t515 * t588 + t626 - t637) * MDP(15) + (t573 + t639) * MDP(16) + (-t444 * t562 + t549 * t447 - t481 * t490 + (-t547 * t605 - t457) * t515 + (t482 * t515 + t567) * t559) * MDP(18) + (t444 * t559 + t549 * t446 - t481 * t492 + (t547 * t606 + t613) * t515 + t567 * t562) * MDP(19) + (-t340 * t535 - t341 * t577 + t362 * t660 - t363 * t661 + t400 * t472 - t401 * t471 + t431 * t617 + t579 * t616) * MDP(20) + (t340 * t471 + t341 * t472 + t617 * t362 + t616 * t363 + t407 * t575 + t427 * t652) * MDP(21) + (t346 * t479 + t619 * t658) * MDP(22) + (t346 * t578 - t347 * t479 + t388 * t619 - t618 * t658) * MDP(23) + (-t644 - t651) * MDP(24) + (t587 - t646) * MDP(25) + ((t439 * t561 - t440 * t558) * t508 + t489 * t347 - t375 * t578 + (t558 * t582 - t561 * t583) * t507 - t611 * t388 + t618 * t382) * MDP(27) + (-(t439 * t558 + t440 * t561) * t508 + t489 * t346 + t375 * t479 + (t558 * t583 + t561 * t582) * t507 + t611 * t658 + t619 * t382) * MDP(28) - ((-t477 + t481) * MDP(11) + t515 * MDP(17) + t412 * MDP(18) - t413 * MDP(19) + t507 * MDP(26) + t338 * MDP(27) - t339 * MDP(28)) * t524 + (MDP(9) * t560 * t565 + MDP(10) * t621) * pkin(1); (-t521 ^ 2 - t524 ^ 2) * MDP(11) + (t476 * t524 - t477 * t521 + t545) * MDP(12) + (t573 - t639) * MDP(18) + (-t515 ^ 2 * t562 - t626 - t637) * MDP(19) + (t400 * t535 + t401 * t577 - t431 * t661 - t579 * t660) * MDP(20) + (-t340 * t577 + t341 * t535 - t362 * t661 - t363 * t660 - t427 * t524) * MDP(21) + (t587 + t646) * MDP(27) + (-t644 + t651) * MDP(28); t492 * t490 * MDP(13) + (-t490 ^ 2 + t492 ^ 2) * MDP(14) + (t446 + t640) * MDP(15) + (-t447 + t638) * MDP(16) + t508 * MDP(17) + (t413 * t515 - t469 * t492 + t566) * MDP(18) + (t412 * t515 + t469 * t490 - t568) * MDP(19) + ((t400 * t554 - t401 * t556) * pkin(4) + (t362 - t366) * t579 + (-t363 - t365) * t431) * MDP(20) + (-t362 * t365 - t363 * t366 + (t340 * t556 + t341 * t554 - t427 * t492) * pkin(4)) * MDP(21) + (t346 - t645) * MDP(24) + (-t347 + t647) * MDP(25) + (t571 * t508 - (t351 * t561 - t352 * t558) * t507 + t421 * t388 + (-t507 * t572 - t339) * qJD(6) + t663) * MDP(27) + (-t572 * t508 + (t351 * t558 + t352 * t561) * t507 - t421 * t658 + (-t507 * t571 - t624) * qJD(6) + t669) * MDP(28) + t668; (-t431 ^ 2 - t579 ^ 2) * MDP(20) + (-t362 * t431 - t363 * t579 + t407) * MDP(21) + (t347 + t647) * MDP(27) + (t346 + t645) * MDP(28); (t599 - t645) * MDP(24) + (-t590 + t647) * MDP(25) + (t339 * t507 + t663) * MDP(27) + (t338 * t507 + t669) * MDP(28) + (MDP(24) * t643 - MDP(25) * t658 - MDP(27) * t339 - MDP(28) * t624) * qJD(6) + t668;];
tauc  = t1;
