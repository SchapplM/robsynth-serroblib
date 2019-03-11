% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:41
% EndTime: 2019-03-09 13:14:52
% DurationCPUTime: 6.96s
% Computational Cost: add. (10762->422), mult. (28393->556), div. (0->0), fcn. (23059->10), ass. (0->216)
t538 = cos(qJ(6));
t599 = qJD(6) * t538;
t535 = sin(qJ(5));
t539 = cos(qJ(5));
t532 = sin(pkin(11));
t533 = cos(pkin(11));
t537 = sin(qJ(2));
t541 = cos(qJ(2));
t509 = -t532 * t537 + t533 * t541;
t499 = t509 * qJD(1);
t510 = t532 * t541 + t533 * t537;
t501 = t510 * qJD(1);
t536 = sin(qJ(4));
t540 = cos(qJ(4));
t557 = t499 * t536 + t540 * t501;
t558 = t499 * t540 - t536 * t501;
t412 = t535 * t557 - t539 * t558;
t676 = t412 * t538;
t680 = t599 + t676;
t500 = t510 * qJD(2);
t490 = qJD(1) * t500;
t595 = qJD(1) * qJD(2);
t590 = t541 * t595;
t591 = t537 * t595;
t491 = -t532 * t591 + t533 * t590;
t409 = qJD(4) * t557 + t540 * t490 + t491 * t536;
t637 = -qJ(3) - pkin(7);
t588 = qJD(2) * t637;
t496 = qJD(3) * t541 + t537 * t588;
t478 = t496 * qJD(1);
t497 = -qJD(3) * t537 + t541 * t588;
t479 = t497 * qJD(1);
t440 = -t478 * t532 + t533 * t479;
t423 = -pkin(8) * t491 + t440;
t441 = t533 * t478 + t532 * t479;
t424 = -pkin(8) * t490 + t441;
t519 = t637 * t541;
t515 = qJD(1) * t519;
t504 = t532 * t515;
t518 = t637 * t537;
t514 = qJD(1) * t518;
t636 = qJD(2) * pkin(2);
t508 = t514 + t636;
t460 = t533 * t508 + t504;
t638 = pkin(8) * t501;
t434 = qJD(2) * pkin(3) + t460 - t638;
t619 = t533 * t515;
t461 = t532 * t508 - t619;
t639 = pkin(8) * t499;
t439 = t461 + t639;
t604 = qJD(4) * t536;
t549 = -(qJD(4) * t434 + t424) * t540 - t536 * t423 + t439 * t604;
t357 = -pkin(9) * t409 - t549;
t581 = t540 * t434 - t439 * t536;
t657 = pkin(9) * t557;
t388 = t581 - t657;
t529 = qJD(2) + qJD(4);
t386 = pkin(4) * t529 + t388;
t603 = qJD(4) * t540;
t408 = -t536 * t490 + t540 * t491 + t499 * t603 - t501 * t604;
t563 = -t434 * t536 - t439 * t540;
t547 = qJD(4) * t563 + t540 * t423 - t536 * t424;
t358 = -pkin(9) * t408 + t547;
t656 = pkin(9) * t558;
t389 = -t563 + t656;
t602 = qJD(5) * t535;
t584 = t358 * t535 - t389 * t602;
t342 = t539 * (qJD(5) * t386 + t357) + t584;
t593 = -pkin(2) * t541 - pkin(1);
t568 = t593 * qJD(1);
t516 = qJD(3) + t568;
t466 = -pkin(3) * t499 + t516;
t425 = -pkin(4) * t558 + t466;
t668 = t412 * t425;
t679 = -t342 + t668;
t534 = sin(qJ(6));
t600 = qJD(6) * t534;
t601 = qJD(5) * t539;
t374 = t539 * t408 - t535 * t409 - t557 * t602 + t558 * t601;
t528 = qJD(5) + t529;
t611 = t538 * t374 + t528 * t599;
t645 = t535 * t558 + t539 * t557;
t354 = -t600 * t645 + t611;
t353 = t354 * t538;
t400 = t528 * t534 + t538 * t645;
t634 = t374 * t534;
t355 = t400 * qJD(6) + t634;
t622 = t645 * t534;
t398 = -t538 * t528 + t622;
t678 = -t534 * t355 - t680 * t398 + t353;
t352 = t354 * t534;
t375 = qJD(5) * t645 + t408 * t535 + t539 * t409;
t371 = t534 * t375;
t672 = -qJD(6) - t412;
t612 = -t599 * t672 + t371;
t624 = t412 * t528;
t627 = t645 * t528;
t629 = t400 * t645;
t677 = (-t375 + t627) * MDP(23) - t412 ^ 2 * MDP(21) + (t412 * MDP(20) + MDP(21) * t645 + MDP(31) * t672) * t645 + (t374 + t624) * MDP(22) + (t680 * t400 + t352) * MDP(27) + (-t672 * t676 + t612 - t629) * MDP(29);
t632 = t389 * t535;
t363 = t386 * t539 - t632;
t361 = -pkin(5) * t528 - t363;
t669 = t361 * t412;
t631 = t389 * t539;
t364 = t386 * t535 + t631;
t585 = t357 * t535 - t539 * t358;
t343 = qJD(5) * t364 + t585;
t653 = t645 * t425;
t675 = -t343 - t653;
t464 = -t514 * t532 + t619;
t442 = t464 - t639;
t465 = t533 * t514 + t504;
t443 = t465 - t638;
t523 = pkin(2) * t533 + pkin(3);
t641 = pkin(2) * t532;
t554 = t523 * t540 - t536 * t641;
t674 = -t554 * qJD(4) + t536 * t442 + t540 * t443;
t495 = t523 * t536 + t540 * t641;
t673 = -t495 * qJD(4) - t540 * t442 + t443 * t536;
t385 = pkin(5) * t645 + pkin(10) * t412;
t373 = t538 * t375;
t630 = t398 * t645;
t666 = t672 * t534;
t670 = (-t666 * t672 + t373 + t630) * MDP(30) + (t400 * t666 + t678) * MDP(28) + t677;
t664 = -t656 - t673;
t663 = -t657 + t674;
t362 = pkin(10) * t528 + t364;
t378 = pkin(5) * t412 - pkin(10) * t645 + t425;
t349 = t362 * t538 + t378 * t534;
t650 = t343 * t534 + t349 * t645 + t361 * t599;
t565 = t362 * t534 - t378 * t538;
t651 = -t343 * t538 + t361 * t600 + t565 * t645;
t620 = t558 * t529;
t621 = t557 * t529;
t662 = (t557 ^ 2 - t558 ^ 2) * MDP(14) - t558 * MDP(13) * t557 + (t408 - t620) * MDP(15) + (-t409 + t621) * MDP(16);
t659 = -0.2e1 * t595;
t640 = pkin(4) * t557;
t655 = MDP(4) * t537;
t654 = MDP(5) * (t537 ^ 2 - t541 ^ 2);
t652 = -t600 * t672 - t373;
t647 = -t466 * t557 + t547;
t646 = -t466 * t558 + t549;
t463 = t509 * t536 + t510 * t540;
t503 = t509 * qJD(2);
t419 = qJD(4) * t463 + t540 * t500 + t503 * t536;
t451 = -t496 * t532 + t533 * t497;
t430 = -pkin(8) * t503 + t451;
t452 = t533 * t496 + t532 * t497;
t431 = -pkin(8) * t500 + t452;
t468 = t533 * t518 + t519 * t532;
t453 = -pkin(8) * t510 + t468;
t469 = t532 * t518 - t533 * t519;
t454 = pkin(8) * t509 + t469;
t552 = t536 * t430 + t540 * t431 + t453 * t603 - t454 * t604;
t369 = -pkin(9) * t419 + t552;
t556 = t540 * t509 - t510 * t536;
t418 = qJD(4) * t556 - t500 * t536 + t503 * t540;
t562 = -t453 * t536 - t454 * t540;
t545 = qJD(4) * t562 + t540 * t430 - t431 * t536;
t370 = -pkin(9) * t418 + t545;
t392 = -pkin(9) * t463 + t453 * t540 - t454 * t536;
t393 = pkin(9) * t556 - t562;
t564 = t392 * t539 - t393 * t535;
t344 = qJD(5) * t564 + t369 * t539 + t370 * t535;
t377 = t392 * t535 + t393 * t539;
t420 = t463 * t535 - t539 * t556;
t380 = -qJD(5) * t420 + t418 * t539 - t419 * t535;
t421 = t463 * t539 + t535 * t556;
t481 = -pkin(3) * t509 + t593;
t435 = -pkin(4) * t556 + t481;
t383 = pkin(5) * t420 - pkin(10) * t421 + t435;
t644 = t343 * t421 + t361 * t380 - t377 * t375 + (qJD(6) * t383 + t344) * t672 - (qJD(6) * t378 + t342) * t420;
t635 = t361 * t421;
t633 = t383 * t375;
t628 = t400 * t534;
t542 = qJD(2) ^ 2;
t618 = t537 * t542;
t617 = t541 * t542;
t543 = qJD(1) ^ 2;
t616 = t541 * t543;
t494 = pkin(4) + t554;
t560 = t494 * t539 - t495 * t535;
t613 = -qJD(5) * t560 + t535 * t664 + t539 * t663;
t559 = t494 * t535 + t495 * t539;
t610 = qJD(5) * t559 - t535 * t663 + t539 * t664;
t605 = qJD(1) * t537;
t527 = t537 * t636;
t522 = pkin(2) * t591;
t467 = pkin(3) * t490 + t522;
t474 = pkin(3) * t500 + t527;
t473 = pkin(2) * t605 + pkin(3) * t501;
t589 = -pkin(4) * t528 - t386;
t587 = pkin(1) * t659;
t426 = t473 + t640;
t447 = pkin(10) + t559;
t572 = qJD(6) * t447 + t385 + t426;
t524 = pkin(4) * t535 + pkin(10);
t571 = qJD(6) * t524 + t385 + t640;
t365 = t388 * t535 + t631;
t570 = pkin(4) * t602 - t365;
t366 = t388 * t539 - t632;
t569 = -pkin(4) * t601 + t366;
t567 = -t375 * t524 + t669;
t566 = -t375 * t447 + t669;
t394 = pkin(4) * t409 + t467;
t397 = pkin(4) * t419 + t474;
t555 = t412 * t666 - t652;
t553 = t380 * t538 - t421 * t600;
t525 = -pkin(4) * t539 - pkin(5);
t446 = -pkin(5) - t560;
t381 = qJD(5) * t421 + t418 * t535 + t539 * t419;
t350 = pkin(5) * t381 - pkin(10) * t380 + t397;
t347 = pkin(5) * t375 - pkin(10) * t374 + t394;
t346 = t538 * t347;
t345 = qJD(5) * t377 + t369 * t535 - t370 * t539;
t1 = [(t418 * MDP(15) - t419 * MDP(16) + t545 * MDP(18) - t552 * MDP(19)) * t529 + (t481 * t408 + t466 * t418 + t467 * t463 + t474 * t557) * MDP(19) + (t408 * t463 + t418 * t557) * MDP(13) + (t481 * t409 + t466 * t419 - t467 * t556 - t474 * t558) * MDP(18) + (t408 * t556 - t409 * t463 + t418 * t558 - t419 * t557) * MDP(14) + (t374 * t435 + t380 * t425 + t394 * t421 + t397 * t645) * MDP(26) + (t374 * t421 + t380 * t645) * MDP(20) + (-t374 * t420 - t375 * t421 - t380 * t412 - t381 * t645) * MDP(21) + (-t440 * t510 + t441 * t509 - t451 * t501 + t452 * t499 - t460 * t503 - t461 * t500 - t468 * t491 - t469 * t490) * MDP(11) + MDP(6) * t617 + (t354 * t420 + t373 * t421 + t381 * t400 - t553 * t672) * MDP(29) + (-t421 * t371 - t355 * t420 - t381 * t398 - (-t380 * t534 - t421 * t599) * t672) * MDP(30) + (t375 * t420 - t381 * t672) * MDP(31) + (t345 * t400 - t349 * t381 - t564 * t354 + ((-qJD(6) * t377 + t350) * t672 - t633 - (-qJD(6) * t362 + t347) * t420 - qJD(6) * t635) * t534 + t644 * t538) * MDP(33) + (t345 * t398 + t346 * t420 - t565 * t381 - t564 * t355 + (-t350 * t672 + t633 + (-t362 * t420 + t377 * t672 + t635) * qJD(6)) * t538 + t644 * t534) * MDP(32) + (t375 * t435 + t381 * t425 + t394 * t420 + t397 * t412) * MDP(25) + (MDP(22) * t380 - MDP(23) * t381 - MDP(25) * t345 - MDP(26) * t344) * t528 + 0.2e1 * t590 * t655 + t654 * t659 + (t440 * t468 + t441 * t469 + t460 * t451 + t461 * t452 + (t516 + t568) * t527) * MDP(12) + (-pkin(7) * t617 + t537 * t587) * MDP(9) - MDP(7) * t618 + (pkin(7) * t618 + t541 * t587) * MDP(10) + (t353 * t421 + t400 * t553) * MDP(27) + ((-t398 * t538 - t628) * t380 + (-t352 - t355 * t538 + (t398 * t534 - t400 * t538) * qJD(6)) * t421) * MDP(28); t662 + ((t461 + t464) * t501 + (t460 - t465) * t499 + (-t490 * t532 - t491 * t533) * pkin(2)) * MDP(11) + (t446 * t355 + t566 * t534 + t610 * t398 - (t534 * t613 - t538 * t572) * t672 + t651) * MDP(32) + (t446 * t354 + t566 * t538 + t610 * t400 - (t534 * t572 + t538 * t613) * t672 + t650) * MDP(33) + (MDP(9) * t537 * t543 + MDP(10) * t616) * pkin(1) + (t628 * t672 + t678) * MDP(28) - t616 * t655 + (-t473 * t557 + t674 * t529 + t646) * MDP(19) + (t473 * t558 + t673 * t529 + t647) * MDP(18) + (-t412 * t426 - t528 * t610 + t675) * MDP(25) + (-t426 * t645 + t528 * t613 + t679) * MDP(26) + t543 * t654 + (-t460 * t464 - t461 * t465 + (t440 * t533 + t441 * t532 - t516 * t605) * pkin(2)) * MDP(12) + (t555 + t630) * MDP(30) + t677; (-t499 ^ 2 - t501 ^ 2) * MDP(11) + (t460 * t501 - t461 * t499 + t522) * MDP(12) + (t409 + t621) * MDP(18) + (t408 + t620) * MDP(19) + (t375 + t627) * MDP(25) + (t374 - t624) * MDP(26) + (t555 - t630) * MDP(32) + (-t538 * t672 ^ 2 - t371 - t629) * MDP(33); (-t529 * t563 + t647) * MDP(18) + (t529 * t581 + t646) * MDP(19) + (-t412 * t640 + t365 * t528 - t653 + (t535 * t589 - t631) * qJD(5) - t585) * MDP(25) + (-t645 * t640 + t366 * t528 + t668 + (qJD(5) * t589 - t357) * t539 - t584) * MDP(26) + (t525 * t355 + t567 * t534 + t570 * t398 - (t534 * t569 - t538 * t571) * t672 + t651) * MDP(32) + (t525 * t354 + t567 * t538 + t570 * t400 - (t534 * t571 + t538 * t569) * t672 + t650) * MDP(33) + t662 + t670; (t364 * t528 + t675) * MDP(25) + (t363 * t528 + t679) * MDP(26) + (-pkin(5) * t355 + (-t363 * t534 + t385 * t538) * t672 - t364 * t398 + t534 * t669 - t612 * pkin(10) + t651) * MDP(32) + (-pkin(5) * t354 - (t363 * t538 + t385 * t534) * t672 - t364 * t400 + t361 * t676 + t652 * pkin(10) + t650) * MDP(33) + t670; t400 * t398 * MDP(27) + (-t398 ^ 2 + t400 ^ 2) * MDP(28) + (-t398 * t672 + t611) * MDP(29) + (-t400 * t672 - t634) * MDP(30) + t375 * MDP(31) + (-t342 * t534 - t349 * t672 - t361 * t400 + t346) * MDP(32) + (-t342 * t538 - t347 * t534 + t361 * t398 + t565 * t672) * MDP(33) + (-MDP(29) * t622 - MDP(30) * t400 - MDP(32) * t349 + MDP(33) * t565) * qJD(6);];
tauc  = t1;
