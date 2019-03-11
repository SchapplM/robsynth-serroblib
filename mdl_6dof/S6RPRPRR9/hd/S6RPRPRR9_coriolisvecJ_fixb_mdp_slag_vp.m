% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:23
% EndTime: 2019-03-09 04:05:38
% DurationCPUTime: 8.85s
% Computational Cost: add. (14273->532), mult. (49792->759), div. (0->0), fcn. (43794->14), ass. (0->232)
t537 = sin(pkin(6));
t546 = cos(qJ(3));
t536 = sin(pkin(7));
t668 = cos(pkin(6));
t610 = t668 * t536;
t591 = t546 * t610;
t539 = cos(pkin(12));
t540 = cos(pkin(7));
t648 = t540 * t546;
t621 = t539 * t648;
t682 = t537 * t621 + t591;
t545 = cos(qJ(5));
t632 = qJD(1) * t537;
t616 = t539 * t632;
t517 = t536 * t616;
t607 = qJD(1) * t668;
t570 = t540 * t607 - t517;
t562 = -qJD(3) - t570;
t498 = t545 * t562;
t542 = sin(qJ(5));
t535 = sin(pkin(12));
t543 = sin(qJ(3));
t649 = t540 * t543;
t487 = t537 * (t535 * t546 + t539 * t649) + t543 * t610;
t481 = qJD(1) * t487;
t534 = sin(pkin(13));
t538 = cos(pkin(13));
t617 = t535 * t632;
t677 = qJD(1) * t682 - t543 * t617;
t577 = t481 * t538 + t534 * t677;
t431 = t542 * t577 + t498;
t603 = -t481 * t534 + t538 * t677;
t678 = t603 - qJD(5);
t681 = t431 * t678;
t680 = t542 * t678;
t503 = (-t534 * t543 + t538 * t546) * t536;
t561 = (-t535 * t648 - t539 * t543) * t537;
t496 = qJD(1) * t561;
t650 = t539 * t546;
t497 = (-t535 * t649 + t650) * t632;
t679 = -qJD(3) * t503 + t496 * t534 + t497 * t538;
t430 = qJD(6) + t431;
t676 = (t535 ^ 2 + t539 ^ 2) * MDP(6) * t537 ^ 2;
t653 = t537 * t539;
t559 = (t540 * t653 + t610) * pkin(9);
t618 = pkin(1) * t668;
t634 = qJ(2) * t653 + t535 * t618;
t484 = t559 + t634;
t609 = t668 * t540;
t507 = t536 * t653 - t609;
t527 = t539 * t618;
t655 = t535 * t537;
t553 = t668 * pkin(2) + (-pkin(9) * t540 - qJ(2)) * t655;
t488 = t527 + t553;
t499 = (-pkin(9) * t535 * t536 - pkin(2) * t539 - pkin(1)) * t537;
t576 = t488 * t540 + t499 * t536;
t425 = -pkin(3) * t507 - qJ(4) * t487 - t484 * t543 + t546 * t576;
t654 = t535 * t543;
t486 = t537 * t654 - t682;
t671 = -t484 * t546 - t543 * t576;
t428 = -qJ(4) * t486 - t671;
t384 = t425 * t534 + t428 * t538;
t382 = -pkin(10) * t507 + t384;
t457 = t486 * t538 + t487 * t534;
t458 = -t486 * t534 + t487 * t538;
t459 = -t488 * t536 + t499 * t540;
t571 = pkin(3) * t486 + t459;
t389 = pkin(4) * t457 - pkin(10) * t458 + t571;
t675 = t382 * t545 + t389 * t542;
t480 = t487 * qJD(3);
t468 = t677 * qJD(3);
t469 = qJD(1) * t480;
t440 = t468 * t534 + t469 * t538;
t596 = pkin(1) * t607;
t506 = qJ(2) * t616 + t535 * t596;
t470 = qJD(1) * t559 + t506;
t524 = t539 * t596;
t475 = qJD(1) * t553 + t524;
t492 = qJD(1) * t499 + qJD(2);
t657 = t492 * t536;
t635 = t475 * t648 + t546 * t657;
t589 = -t470 * t543 + t635;
t423 = -qJ(4) * t481 + t589;
t411 = -pkin(3) * t562 + t423;
t555 = -t543 * (t475 * t540 + t657) - t546 * t470;
t424 = qJ(4) * t677 - t555;
t652 = t538 * t424;
t376 = t411 * t534 + t652;
t374 = -pkin(10) * t562 + t376;
t456 = -t475 * t536 + t492 * t540;
t429 = -pkin(3) * t677 + qJD(4) + t456;
t386 = -pkin(4) * t603 - pkin(10) * t577 + t429;
t355 = t374 * t545 + t386 * t542;
t622 = qJD(1) * qJD(2);
t612 = t537 * t622;
t595 = t535 * t612;
t574 = t540 * t595;
t594 = t539 * t612;
t630 = qJD(3) * t546;
t613 = t540 * t630;
t614 = t536 * t630;
t619 = t475 * t613 + t492 * t614 + t546 * t594;
t550 = (-qJD(3) * t470 - t574) * t543 + t619;
t394 = -qJ(4) * t469 + qJD(4) * t677 + t550;
t552 = t555 * qJD(3);
t573 = -qJ(4) * t468 - qJD(4) * t481;
t367 = t538 * t394 + (-t543 * t594 - t546 * t574 + t552 + t573) * t534;
t441 = t468 * t538 - t469 * t534;
t511 = t536 * t595;
t464 = pkin(3) * t469 + t511;
t393 = pkin(4) * t440 - pkin(10) * t441 + t464;
t605 = t367 * t542 - t393 * t545;
t670 = -qJD(5) * t355 - t605;
t344 = -pkin(5) * t440 - t670;
t433 = -t542 * t562 + t545 * t577;
t672 = t430 * (pkin(5) * t433 + pkin(11) * t430) + t344;
t669 = MDP(4) * t535 + MDP(5) * t539;
t628 = qJD(5) * t542;
t403 = -qJD(5) * t498 + t441 * t545 - t577 * t628;
t541 = sin(qJ(6));
t544 = cos(qJ(6));
t624 = qJD(6) * t544;
t620 = t403 * t544 + t440 * t541 - t624 * t678;
t625 = qJD(6) * t541;
t362 = -t433 * t625 + t620;
t666 = t362 * t541;
t660 = t433 * t541;
t407 = t544 * t678 + t660;
t665 = t407 * t430;
t409 = t433 * t544 - t541 * t678;
t664 = t409 * t430;
t663 = t409 * t603;
t601 = t430 * t544;
t662 = t431 * t577;
t661 = t433 * t577;
t659 = t603 * t545;
t419 = t534 * t424;
t646 = t542 * t441;
t404 = qJD(5) * t433 + t646;
t647 = t541 * t404;
t645 = t544 * t404;
t378 = t423 * t534 + t652;
t643 = t378 + t678 * (pkin(5) * t542 - pkin(11) * t545);
t379 = t423 * t538 - t419;
t406 = pkin(3) * t481 + pkin(4) * t577 - pkin(10) * t603;
t642 = t379 * t545 + t406 * t542;
t640 = t440 * t545 - t603 * t680;
t504 = (t534 * t546 + t538 * t543) * t536;
t490 = t504 * t542 - t540 * t545;
t598 = t536 * t617;
t639 = qJD(5) * t490 + t542 * t598 + t545 * t679;
t638 = -qJD(3) * t504 - t496 * t538 + t497 * t534;
t491 = t504 * t545 + t540 * t542;
t636 = qJD(5) * t491 - t542 * t679 + t545 * t598;
t631 = qJD(2) * t537;
t629 = qJD(5) * t541;
t627 = qJD(5) * t544;
t626 = qJD(5) * t545;
t530 = -pkin(3) * t538 - pkin(4);
t615 = t535 * t631;
t516 = t536 * t615;
t611 = pkin(3) * t480 + t516;
t565 = t367 * t545 - t374 * t628 + t386 * t626 + t393 * t542;
t343 = pkin(11) * t440 + t565;
t558 = qJD(2) * t561;
t554 = qJD(1) * t558;
t548 = t554 + t552;
t366 = t394 * t534 - t538 * (t548 + t573);
t350 = pkin(5) * t404 - pkin(11) * t403 + t366;
t606 = -t343 * t541 + t350 * t544;
t551 = t488 * t613 + t499 * t614 + t631 * t650 + (-qJD(3) * t484 - t540 * t615) * t543;
t400 = -qJ(4) * t480 - qJD(4) * t486 + t551;
t479 = (t591 + (t621 - t654) * t537) * qJD(3);
t549 = qJD(3) * t671 + t558;
t401 = -qJ(4) * t479 - qJD(4) * t487 + t549;
t371 = t400 * t534 - t401 * t538;
t604 = t403 * t541 - t440 * t544;
t375 = t411 * t538 - t419;
t383 = t425 * t538 - t428 * t534;
t602 = t678 * t545;
t509 = -pkin(5) * t545 - pkin(11) * t542 + t530;
t600 = pkin(11) * t577 - qJD(6) * t509 + t642;
t588 = qJD(6) * t503 + t639;
t587 = qJD(6) * t491 + t638;
t585 = t343 * t544 + t350 * t541;
t353 = -pkin(11) * t678 + t355;
t373 = pkin(4) * t562 - t375;
t360 = pkin(5) * t431 - pkin(11) * t433 + t373;
t346 = t353 * t544 + t360 * t541;
t584 = t353 * t541 - t360 * t544;
t357 = pkin(11) * t457 + t675;
t381 = pkin(4) * t507 - t383;
t434 = t458 * t542 + t507 * t545;
t435 = t458 * t545 - t507 * t542;
t361 = pkin(5) * t434 - pkin(11) * t435 + t381;
t583 = t357 * t544 + t361 * t541;
t582 = -t357 * t541 + t361 * t544;
t372 = t400 * t538 + t401 * t534;
t448 = t479 * t534 + t480 * t538;
t451 = t479 * t538 - t480 * t534;
t397 = pkin(4) * t448 - pkin(10) * t451 + t611;
t581 = -t372 * t542 + t397 * t545;
t354 = -t374 * t542 + t386 * t545;
t579 = -t382 * t542 + t389 * t545;
t415 = t435 * t544 + t457 * t541;
t414 = t435 * t541 - t457 * t544;
t575 = (-qJ(2) * t617 + t524) * t535 - t506 * t539;
t567 = -t430 * t624 - t647;
t566 = t430 * t625 - t645;
t564 = t372 * t545 - t382 * t628 + t389 * t626 + t397 * t542;
t529 = pkin(3) * t534 + pkin(10);
t563 = -t373 * t678 - t440 * t529;
t417 = t541 * t577 + t544 * t659;
t560 = -t542 * t625 + t544 * t626 - t417;
t557 = qJD(5) * t407 + t567;
t352 = pkin(5) * t678 - t354;
t556 = -pkin(11) * t404 + (t352 + t354) * t430;
t416 = t541 * t659 - t544 * t577;
t413 = qJD(5) * t435 + t451 * t542;
t412 = -qJD(5) * t434 + t451 * t545;
t405 = t409 * t628;
t369 = qJD(6) * t415 + t412 * t541 - t448 * t544;
t368 = -qJD(6) * t414 + t412 * t544 + t448 * t541;
t363 = qJD(6) * t409 + t604;
t358 = -pkin(5) * t577 + t379 * t542 - t406 * t545;
t356 = -pkin(5) * t457 - t579;
t351 = pkin(5) * t413 - pkin(11) * t412 + t371;
t348 = -pkin(5) * t448 + qJD(5) * t675 - t581;
t347 = pkin(11) * t448 + t564;
t342 = -qJD(6) * t346 + t606;
t341 = -qJD(6) * t584 + t585;
t1 = [0.2e1 * t622 * t676 + ((-qJD(6) * t583 - t347 * t541 + t351 * t544) * t430 + t582 * t404 + t342 * t434 - t584 * t413 + t348 * t407 + t356 * t363 + t344 * t414 + t352 * t369) * MDP(29) + (t366 * t458 - t367 * t457 + t371 * t577 + t372 * t603 - t375 * t451 - t376 * t448 - t383 * t441 - t384 * t440) * MDP(15) + (-t355 * t448 + t366 * t435 + t371 * t433 + t373 * t412 + t381 * t403 - t440 * t675 - t457 * t565 + t564 * t678) * MDP(23) + (t456 * t479 + t459 * t468 + 0.2e1 * t481 * t516 + t507 * t550 + t551 * t562) * MDP(14) + (-t366 * t383 + t367 * t384 - t375 * t371 + t376 * t372 + t429 * t611 + t464 * t571) * MDP(16) + (-t581 * t678 + t579 * t440 - t605 * t457 + t354 * t448 + t371 * t431 + t381 * t404 + t366 * t434 + t373 * t413 + (-t355 * t457 + t675 * t678) * qJD(5)) * MDP(22) + (t456 * t480 + t459 * t469 + t486 * t511 - t507 * t548 - t516 * t677 - t549 * t562) * MDP(13) + (-t468 * t507 - t479 * t562) * MDP(10) + (t469 * t507 + t480 * t562) * MDP(11) + (-(qJD(6) * t582 + t347 * t544 + t351 * t541) * t430 - t583 * t404 - t341 * t434 - t346 * t413 + t348 * t409 + t356 * t362 + t344 * t415 + t352 * t368) * MDP(30) + (-t468 * t486 - t469 * t487 + t479 * t677 - t480 * t481) * MDP(9) + (t468 * t487 + t479 * t481) * MDP(8) + (t440 * t457 - t448 * t678) * MDP(21) + (-t404 * t457 + t413 * t678 - t431 * t448 - t434 * t440) * MDP(20) + (t403 * t457 - t412 * t678 + t433 * t448 + t435 * t440) * MDP(19) + (t362 * t415 + t368 * t409) * MDP(24) + (-t362 * t414 - t363 * t415 - t368 * t407 - t369 * t409) * MDP(25) + (t404 * t434 + t413 * t430) * MDP(28) + (-t363 * t434 - t369 * t430 - t404 * t414 - t407 * t413) * MDP(27) + (t362 * t434 + t368 * t430 + t404 * t415 + t409 * t413) * MDP(26) + (t403 * t435 + t412 * t433) * MDP(17) + (-t403 * t434 - t404 * t435 - t412 * t431 - t413 * t433) * MDP(18) + (((t539 * t634 + (qJ(2) * t655 - t527) * t535) * qJD(1) - t575) * MDP(7) - 0.2e1 * t669 * t607) * t631; t575 * MDP(7) * t632 + (t540 * t469 + t496 * t562 + (qJD(3) * t543 * t562 + t617 * t677) * t536) * MDP(13) + (t540 * t468 - t497 * t562 + (-t481 * t617 + t562 * t630) * t536) * MDP(14) + (-t440 * t504 - t441 * t503 - t577 * t638 - t603 * t679) * MDP(15) + (-t366 * t503 + t367 * t504 + t375 * t638 - t376 * t679 - t429 * t598 + t464 * t540) * MDP(16) + (-t404 * t503 - t431 * t638 - t440 * t490 + t636 * t678) * MDP(22) + (-t403 * t503 - t433 * t638 - t440 * t491 - t639 * t678) * MDP(23) + ((-t491 * t541 - t503 * t544) * t404 + t490 * t363 + (t541 * t588 - t544 * t587) * t430 + t636 * t407) * MDP(29) + (-(t491 * t544 - t503 * t541) * t404 + t490 * t362 + (t541 * t587 + t544 * t588) * t430 + t636 * t409) * MDP(30) + (t537 * t668 * t669 - t676) * qJD(1) ^ 2; -t481 * t677 * MDP(8) + (t481 ^ 2 - t677 ^ 2) * MDP(9) + (t562 * t677 + t468) * MDP(10) + (t481 * (qJD(3) - t517) + (t481 * t609 - t480) * qJD(1)) * MDP(11) + (-t456 * t481 - t555 * t570 + t554) * MDP(13) + (qJD(3) * t635 - t456 * t677 + t543 * t574 + t570 * t589 - t619) * MDP(14) + ((-t440 * t534 - t441 * t538) * pkin(3) + (t375 - t379) * t603 + (t376 - t378) * t577) * MDP(15) + (t375 * t378 - t376 * t379 + (-t366 * t538 + t367 * t534 - t429 * t481) * pkin(3)) * MDP(16) + (t403 * t542 - t433 * t602) * MDP(17) + ((t403 + t681) * t545 + (t433 * t678 - t404) * t542) * MDP(18) + (t440 * t542 + t602 * t678 - t661) * MDP(19) + (t628 * t678 + t640 + t662) * MDP(20) + t678 * t577 * MDP(21) + (-t354 * t577 - t378 * t431 + t530 * t404 + (-t366 - (-qJD(5) * t529 - t406) * t678) * t545 + (-t379 * t678 + t563) * t542) * MDP(22) + (t355 * t577 + t366 * t542 - t378 * t433 + t530 * t403 - (t529 * t628 + t642) * t678 + t563 * t545) * MDP(23) + (t362 * t542 * t544 + t409 * t560) * MDP(24) + (t407 * t417 + t409 * t416 + (-t407 * t544 - t409 * t541) * t626 + (-t666 - t363 * t544 + (t407 * t541 - t409 * t544) * qJD(6)) * t542) * MDP(25) + (-t362 * t545 + t405 + (t645 - t663) * t542 + t560 * t430) * MDP(26) + (t363 * t545 + (-t541 * t626 + t416) * t430 + (t407 * t678 + t567) * t542) * MDP(27) + (-t404 * t545 - t430 * t680) * MDP(28) + (t509 * t645 - t352 * t416 - t358 * t407 + (t541 * t600 - t544 * t643) * t430 + (t352 * t629 + t529 * t557 - t342) * t545 + (t352 * t624 + t344 * t541 + t584 * t603 + t529 * t363 + (t430 * t529 * t541 - t584) * qJD(5)) * t542) * MDP(29) + (-t509 * t647 - t352 * t417 - t358 * t409 + (t541 * t643 + t544 * t600) * t430 + (t352 * t627 + t341 + (qJD(5) * t409 + t566) * t529) * t545 + (-t352 * t625 + t344 * t544 + t346 * t603 + t529 * t362 + (t529 * t601 - t346) * qJD(5)) * t542) * MDP(30); (-t577 ^ 2 - t603 ^ 2) * MDP(15) + (t375 * t577 - t376 * t603 + t464) * MDP(16) + (t640 - t662) * MDP(22) - MDP(23) * t661 + t416 * t430 * MDP(29) + (t417 * t430 + t405) * MDP(30) + ((-t430 * t629 - t363) * MDP(29) + (-t430 * t627 - t362) * MDP(30) - t678 ^ 2 * MDP(23)) * t545 + (qJD(5) * t678 * MDP(22) - t440 * MDP(23) + (-t407 * t603 + t557) * MDP(29) + (t566 - t663) * MDP(30)) * t542; -t431 ^ 2 * MDP(18) + (t403 - t681) * MDP(19) - t646 * MDP(20) + t440 * MDP(21) + (-t355 * t678 + t670) * MDP(22) + (-t354 * t678 + t373 * t431 - t565) * MDP(23) + (t409 * t601 + t666) * MDP(24) + ((t362 - t665) * t544 + (-t363 - t664) * t541) * MDP(25) + (t430 * t601 + t647) * MDP(26) + (-t430 ^ 2 * t541 + t645) * MDP(27) + (-pkin(5) * t363 - t355 * t407 + t556 * t541 - t544 * t672) * MDP(29) + (-pkin(5) * t362 - t355 * t409 + t541 * t672 + t556 * t544) * MDP(30) + (t431 * MDP(17) + (-qJD(5) - t678) * MDP(20) - t373 * MDP(22) - t409 * MDP(26) + t407 * MDP(27) - t430 * MDP(28) + t584 * MDP(29) + t346 * MDP(30) + t433 * MDP(18)) * t433; t409 * t407 * MDP(24) + (-t407 ^ 2 + t409 ^ 2) * MDP(25) + (t620 + t665) * MDP(26) + (-t604 + t664) * MDP(27) + t404 * MDP(28) + (t346 * t430 - t352 * t409 + t606) * MDP(29) + (t352 * t407 - t430 * t584 - t585) * MDP(30) + (-MDP(26) * t660 - MDP(27) * t409 - MDP(29) * t346 + MDP(30) * t584) * qJD(6);];
tauc  = t1;
