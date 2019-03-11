% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:30
% EndTime: 2019-03-09 10:26:44
% DurationCPUTime: 7.54s
% Computational Cost: add. (10740->496), mult. (32285->681), div. (0->0), fcn. (26382->12), ass. (0->236)
t549 = cos(pkin(11));
t556 = cos(qJ(2));
t547 = sin(pkin(6));
t627 = qJD(1) * t547;
t609 = t556 * t627;
t528 = t549 * t609;
t546 = sin(pkin(11));
t553 = sin(qJ(2));
t610 = t553 * t627;
t512 = -t546 * t610 + t528;
t508 = qJD(4) - t512;
t550 = cos(pkin(6));
t659 = pkin(1) * t553;
t618 = t550 * t659;
t640 = t547 * t556;
t657 = pkin(8) + qJ(3);
t509 = t657 * t640 + t618;
t498 = t509 * qJD(1);
t486 = t546 * t498;
t658 = pkin(1) * t556;
t617 = t550 * t658;
t533 = qJD(1) * t617;
t607 = t657 * t553;
t591 = t547 * t607;
t497 = -qJD(1) * t591 + t533;
t446 = t497 * t549 - t486;
t576 = t546 * t556 + t549 * t553;
t515 = t576 * t627;
t461 = pkin(2) * t610 + pkin(3) * t515 - pkin(9) * t512;
t555 = cos(qJ(4));
t450 = t555 * t461;
t552 = sin(qJ(4));
t539 = pkin(2) * t546 + pkin(9);
t636 = qJ(5) + t539;
t596 = qJD(4) * t636;
t675 = -pkin(4) * t515 - t450 + (qJ(5) * t512 - t596) * t555 + (-qJD(5) + t446) * t552;
t634 = t555 * t446 + t552 * t461;
t644 = t512 * t552;
t674 = -qJ(5) * t644 - qJD(5) * t555 + t552 * t596 + t634;
t624 = qJD(4) * t552;
t673 = t624 - t644;
t554 = cos(qJ(6));
t551 = sin(qJ(6));
t626 = qJD(1) * t550;
t535 = qJD(2) + t626;
t478 = t515 * t552 - t555 * t535;
t545 = sin(pkin(12));
t548 = cos(pkin(12));
t577 = -t515 * t555 - t535 * t552;
t578 = -t478 * t545 - t548 * t577;
t651 = t578 * t551;
t408 = -t554 * t508 + t651;
t597 = -t548 * t478 + t545 * t577;
t661 = qJD(6) - t597;
t672 = t408 * t661;
t410 = t508 * t551 + t554 * t578;
t671 = t410 * t661;
t525 = t545 * t555 + t548 * t552;
t633 = t508 * t525;
t524 = t545 * t552 - t548 * t555;
t670 = t508 * t524;
t594 = t661 * t554;
t625 = qJD(2) * t547;
t608 = t553 * t625;
t589 = qJD(1) * t608;
t507 = qJD(2) * t528 - t546 * t589;
t623 = qJD(4) * t555;
t432 = t555 * t507 - t515 * t624 + t535 * t623;
t433 = -qJD(4) * t577 + t507 * t552;
t402 = t432 * t545 + t548 * t433;
t637 = t551 * t402;
t669 = -t661 * t594 - t637;
t403 = t432 * t548 - t433 * t545;
t518 = t576 * t547;
t513 = qJD(2) * t518;
t506 = qJD(1) * t513;
t603 = t403 * t551 - t554 * t506;
t375 = t410 * qJD(6) + t603;
t668 = t375 * t524 + t408 * t633;
t667 = pkin(2) * t608;
t666 = MDP(4) * t553;
t665 = MDP(5) * (t553 ^ 2 - t556 ^ 2);
t496 = (pkin(2) + t658) * t550 - t591;
t456 = t546 * t496 + t549 * t509;
t444 = pkin(9) * t550 + t456;
t638 = t549 * t556;
t641 = t547 * t553;
t517 = t546 * t641 - t547 * t638;
t588 = (-pkin(2) * t556 - pkin(1)) * t547;
t463 = pkin(3) * t517 - pkin(9) * t518 + t588;
t635 = t555 * t444 + t552 * t463;
t631 = t674 * t545 + t548 * t675;
t629 = t545 * t675 - t674 * t548;
t639 = t549 * t498;
t445 = t497 * t546 + t639;
t662 = pkin(4) * t673 - t445;
t481 = pkin(2) * t535 + t497;
t430 = t546 * t481 + t639;
t425 = pkin(9) * t535 + t430;
t574 = qJD(1) * t588;
t520 = qJD(3) + t574;
t441 = -pkin(3) * t512 - pkin(9) * t515 + t520;
t399 = t425 * t555 + t441 * t552;
t531 = qJD(2) * t533;
t560 = (-qJD(2) * t607 + qJD(3) * t556) * t547;
t472 = qJD(1) * t560 + t531;
t483 = -qJD(2) * t509 - qJD(3) * t641;
t473 = t483 * qJD(1);
t415 = t472 * t549 + t473 * t546;
t530 = pkin(2) * t589;
t442 = pkin(3) * t506 - pkin(9) * t507 + t530;
t602 = -t415 * t552 + t555 * t442;
t558 = -qJD(4) * t399 + t602;
t351 = pkin(4) * t506 - qJ(5) * t432 + qJD(5) * t577 + t558;
t568 = t555 * t415 - t425 * t624 + t441 * t623 + t552 * t442;
t355 = -qJ(5) * t433 - qJD(5) * t478 + t568;
t340 = t351 * t548 - t355 * t545;
t338 = -pkin(5) * t506 - t340;
t538 = pkin(4) * t545 + pkin(10);
t660 = t661 * (-pkin(4) * t577 + pkin(5) * t578 - pkin(10) * t597 + qJD(6) * t538) + t338;
t621 = qJD(6) * t554;
t611 = t554 * t403 + t551 * t506 + t508 * t621;
t622 = qJD(6) * t551;
t374 = -t578 * t622 + t611;
t656 = t374 * t551;
t387 = -qJ(5) * t478 + t399;
t654 = t387 * t545;
t653 = t408 * t578;
t652 = t410 * t578;
t650 = t478 * t508;
t649 = t478 * t515;
t648 = t577 * t508;
t647 = t577 * t515;
t646 = t506 * t552;
t643 = t525 * t554;
t542 = t547 ^ 2;
t557 = qJD(1) ^ 2;
t642 = t542 * t557;
t384 = t548 * t387;
t400 = t554 * t402;
t341 = t545 * t351 + t548 * t355;
t489 = t518 * t552 - t550 * t555;
t514 = (-t546 * t553 + t638) * t625;
t454 = -qJD(4) * t489 + t514 * t555;
t490 = t518 * t555 + t550 * t552;
t534 = qJD(2) * t617;
t482 = t534 + t560;
t427 = t482 * t549 + t483 * t546;
t462 = pkin(3) * t513 - pkin(9) * t514 + t667;
t599 = -t427 * t552 + t555 * t462;
t362 = pkin(4) * t513 - qJ(5) * t454 - qJD(4) * t635 - qJD(5) * t490 + t599;
t453 = qJD(4) * t490 + t514 * t552;
t567 = t555 * t427 - t444 * t624 + t552 * t462 + t463 * t623;
t368 = -qJ(5) * t453 - qJD(5) * t489 + t567;
t345 = t545 * t362 + t548 * t368;
t398 = -t425 * t552 + t555 * t441;
t386 = qJ(5) * t577 + t398;
t381 = pkin(4) * t508 + t386;
t357 = t545 * t381 + t384;
t598 = -t444 * t552 + t555 * t463;
t390 = pkin(4) * t517 - qJ(5) * t490 + t598;
t396 = -qJ(5) * t489 + t635;
t370 = t545 * t390 + t548 * t396;
t630 = pkin(5) * t515 - t631;
t616 = t542 * t659;
t614 = t525 * t637;
t613 = t525 * t400;
t612 = t556 * t642;
t541 = -pkin(2) * t549 - pkin(3);
t606 = qJD(1) * qJD(2) * t542;
t605 = t636 * t552;
t339 = pkin(10) * t506 + t341;
t414 = t472 * t546 - t549 * t473;
t394 = pkin(4) * t433 + t414;
t350 = pkin(5) * t402 - pkin(10) * t403 + t394;
t604 = -t339 * t551 + t554 * t350;
t601 = -t554 * t515 + t551 * t670;
t600 = t515 * t551 + t554 * t670;
t426 = t482 * t546 - t549 * t483;
t429 = t481 * t549 - t486;
t455 = t496 * t549 - t546 * t509;
t595 = t508 * t555;
t590 = t556 * t606;
t587 = t374 * t524 + t410 * t633;
t575 = -pkin(4) * t555 + t541;
t471 = pkin(5) * t524 - pkin(10) * t525 + t575;
t585 = pkin(10) * t515 - qJD(6) * t471 - t629;
t523 = t636 * t555;
t477 = t548 * t523 - t545 * t605;
t584 = -pkin(5) * t633 - pkin(10) * t670 + qJD(6) * t477 - t662;
t583 = t339 * t554 + t350 * t551;
t353 = pkin(10) * t508 + t357;
t424 = -pkin(3) * t535 - t429;
t407 = pkin(4) * t478 + qJD(5) + t424;
t376 = -pkin(5) * t597 - pkin(10) * t578 + t407;
t347 = t353 * t554 + t376 * t551;
t582 = t353 * t551 - t376 * t554;
t344 = t362 * t548 - t368 * t545;
t365 = pkin(10) * t517 + t370;
t434 = t548 * t489 + t490 * t545;
t435 = -t489 * t545 + t490 * t548;
t443 = -pkin(3) * t550 - t455;
t561 = pkin(4) * t489 + t443;
t379 = pkin(5) * t434 - pkin(10) * t435 + t561;
t581 = t365 * t554 + t379 * t551;
t580 = -t365 * t551 + t379 * t554;
t356 = t381 * t548 - t654;
t369 = t390 * t548 - t396 * t545;
t412 = t435 * t554 + t517 * t551;
t411 = t435 * t551 - t517 * t554;
t573 = t555 * t506 - t508 * t673;
t572 = t400 + (t551 * t597 - t622) * t661;
t571 = pkin(4) * t453 + t426;
t570 = -pkin(8) * t640 - t618;
t569 = -pkin(8) * t589 + t531;
t566 = t424 * t508 - t539 * t506;
t565 = t525 * t621 - t601;
t564 = t525 * t622 + t600;
t562 = t570 * t535;
t352 = -pkin(5) * t508 - t356;
t360 = t386 * t548 - t654;
t559 = -t538 * t402 + (t352 + t360) * t661;
t540 = -pkin(4) * t548 - pkin(5);
t476 = t523 * t545 + t548 * t605;
t406 = -t453 * t545 + t454 * t548;
t405 = t548 * t453 + t454 * t545;
t378 = qJD(6) * t412 + t406 * t551 - t513 * t554;
t377 = -qJD(6) * t411 + t406 * t554 + t513 * t551;
t364 = -pkin(5) * t517 - t369;
t359 = t386 * t545 + t384;
t358 = pkin(5) * t405 - pkin(10) * t406 + t571;
t343 = pkin(10) * t513 + t345;
t342 = -pkin(5) * t513 - t344;
t337 = -qJD(6) * t347 + t604;
t336 = -qJD(6) * t582 + t583;
t1 = [(-t414 * t455 + t415 * t456 - t429 * t426 + t430 * t427 + (t520 + t574) * t667) * MDP(12) + (t432 * t490 - t454 * t577) * MDP(13) + (-t432 * t489 - t433 * t490 + t453 * t577 - t454 * t478) * MDP(14) + (t432 * t517 + t454 * t508 + t490 * t506 - t513 * t577) * MDP(15) + (-t433 * t517 - t453 * t508 - t478 * t513 - t489 * t506) * MDP(16) + (t506 * t517 + t508 * t513) * MDP(17) + (t599 * t508 + t598 * t506 + t602 * t517 + t398 * t513 + t426 * t478 + t443 * t433 + t414 * t489 + t424 * t453 + (-t399 * t517 - t508 * t635) * qJD(4)) * MDP(18) + (-t399 * t513 + t414 * t490 + t424 * t454 - t426 * t577 + t443 * t432 - t506 * t635 - t508 * t567 - t517 * t568) * MDP(19) + (-t340 * t435 - t341 * t434 - t344 * t578 + t345 * t597 - t356 * t406 - t357 * t405 - t369 * t403 - t370 * t402) * MDP(20) + (t340 * t369 + t341 * t370 + t356 * t344 + t357 * t345 + t394 * t561 + t407 * t571) * MDP(21) + (t374 * t412 + t377 * t410) * MDP(22) + (-t374 * t411 - t375 * t412 - t377 * t408 - t378 * t410) * MDP(23) + (t374 * t434 + t377 * t661 + t402 * t412 + t405 * t410) * MDP(24) + (-t375 * t434 - t378 * t661 - t402 * t411 - t405 * t408) * MDP(25) + (t402 * t434 + t405 * t661) * MDP(26) + ((-qJD(6) * t581 - t343 * t551 + t358 * t554) * t661 + t580 * t402 + t337 * t434 - t582 * t405 + t342 * t408 + t364 * t375 + t338 * t411 + t352 * t378) * MDP(27) + (-(qJD(6) * t580 + t343 * t554 + t358 * t551) * t661 - t581 * t402 - t336 * t434 - t347 * t405 + t342 * t410 + t364 * t374 + t338 * t412 + t352 * t377) * MDP(28) - 0.2e1 * t606 * t665 + (t562 + (t550 * t570 - 0.2e1 * t616) * qJD(1)) * qJD(2) * MDP(9) + (-0.2e1 * pkin(1) * t590 - (-pkin(8) * t608 + t534) * t535 - t569 * t550) * MDP(10) + (t414 * t518 - t415 * t517 + t426 * t515 + t427 * t512 - t429 * t514 - t430 * t513 - t455 * t507 - t456 * t506) * MDP(11) + 0.2e1 * t590 * t666 + (MDP(6) * t556 * t625 - MDP(7) * t608) * (t535 + t626); -t508 * t515 * MDP(17) + t642 * t665 + (t557 * t616 + (qJD(2) * t570 - t562) * qJD(1)) * MDP(9) + (pkin(1) * t612 + (-pkin(8) * t610 + t533) * t535 - t569) * MDP(10) + ((t430 - t445) * t515 + (t429 - t446) * t512 + (-t506 * t546 - t507 * t549) * pkin(2)) * MDP(11) + (t429 * t445 - t430 * t446 + (-t414 * t549 + t415 * t546 - t520 * t610) * pkin(2)) * MDP(12) + (t432 * t552 - t577 * t595) * MDP(13) + ((t432 - t650) * t555 + (-t433 + t648) * t552) * MDP(14) + (t508 * t595 + t646 + t647) * MDP(15) + (t573 + t649) * MDP(16) + (-t398 * t515 - t414 * t555 + t541 * t433 - t445 * t478 + (-t539 * t623 - t450) * t508 + (t446 * t508 + t566) * t552) * MDP(18) + (t399 * t515 + t414 * t552 + t541 * t432 + t445 * t577 + (t539 * t624 + t634) * t508 + t566 * t555) * MDP(19) + (-t340 * t525 - t341 * t524 + t356 * t670 - t357 * t633 - t477 * t402 + t403 * t476 - t578 * t631 + t597 * t629) * MDP(20) + (-t340 * t476 + t341 * t477 + t631 * t356 + t629 * t357 + t394 * t575 + t407 * t662) * MDP(21) + (t374 * t643 - t410 * t564) * MDP(22) + (t601 * t410 + t600 * t408 + (-t656 - t375 * t554 + (t408 * t551 - t410 * t554) * qJD(6)) * t525) * MDP(23) + (-t564 * t661 + t587 + t613) * MDP(24) + (-t565 * t661 - t614 - t668) * MDP(25) + (t402 * t524 + t633 * t661) * MDP(26) + ((t471 * t554 - t477 * t551) * t402 + t337 * t524 + t476 * t375 + t338 * t551 * t525 + (t551 * t585 - t554 * t584) * t661 + t630 * t408 - t633 * t582 + t565 * t352) * MDP(27) + (-(t471 * t551 + t477 * t554) * t402 - t336 * t524 + t476 * t374 + t338 * t643 + (t551 * t584 + t554 * t585) * t661 + t630 * t410 - t633 * t347 - t564 * t352) * MDP(28) - t612 * t666 + (MDP(6) * t609 - MDP(7) * t610) * (qJD(2) - t535); (-t512 ^ 2 - t515 ^ 2) * MDP(11) + (t429 * t515 - t430 * t512 + t530) * MDP(12) + (t573 - t649) * MDP(18) + (-t508 ^ 2 * t555 - t646 + t647) * MDP(19) + (-t402 * t525 + t403 * t524 + t578 * t633 - t597 * t670) * MDP(20) + (-t340 * t524 + t341 * t525 - t356 * t633 - t357 * t670 - t407 * t515) * MDP(21) + (-t614 + t668) * MDP(27) + (t587 - t613) * MDP(28) + (-t565 * MDP(27) + t564 * MDP(28)) * t661; -t577 * t478 * MDP(13) + (-t478 ^ 2 + t577 ^ 2) * MDP(14) + (t432 + t650) * MDP(15) + (-t433 - t648) * MDP(16) + t506 * MDP(17) + (t399 * t508 + t424 * t577 + t558) * MDP(18) + (t398 * t508 + t424 * t478 - t568) * MDP(19) + ((-t402 * t545 - t403 * t548) * pkin(4) + (t356 - t360) * t597 + (t357 - t359) * t578) * MDP(20) + (t356 * t359 - t357 * t360 + (t340 * t548 + t341 * t545 + t407 * t577) * pkin(4)) * MDP(21) + (t410 * t594 + t656) * MDP(22) + ((t374 - t672) * t554 + (-t375 - t671) * t551) * MDP(23) + (-t652 - t669) * MDP(24) + (t572 + t653) * MDP(25) - t661 * t578 * MDP(26) + (-t359 * t408 + t540 * t375 + t559 * t551 - t554 * t660 + t578 * t582) * MDP(27) + (t347 * t578 - t359 * t410 + t540 * t374 + t551 * t660 + t559 * t554) * MDP(28); (-t578 ^ 2 - t597 ^ 2) * MDP(20) + (t356 * t578 - t357 * t597 + t394) * MDP(21) + (t572 - t653) * MDP(27) + (-t652 + t669) * MDP(28); t410 * t408 * MDP(22) + (-t408 ^ 2 + t410 ^ 2) * MDP(23) + (t611 + t672) * MDP(24) + (-t603 + t671) * MDP(25) + t402 * MDP(26) + (t347 * t661 - t352 * t410 + t604) * MDP(27) + (t352 * t408 - t582 * t661 - t583) * MDP(28) + (-MDP(24) * t651 - MDP(25) * t410 - MDP(27) * t347 + MDP(28) * t582) * qJD(6);];
tauc  = t1;
