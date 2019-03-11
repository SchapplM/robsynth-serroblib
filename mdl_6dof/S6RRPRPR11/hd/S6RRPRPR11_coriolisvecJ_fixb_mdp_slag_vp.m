% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:16
% EndTime: 2019-03-09 11:15:29
% DurationCPUTime: 7.60s
% Computational Cost: add. (5079->478), mult. (11682->649), div. (0->0), fcn. (7769->8), ass. (0->229)
t560 = cos(qJ(6));
t561 = cos(qJ(4));
t558 = sin(qJ(4));
t629 = qJD(2) * t558;
t562 = cos(qJ(2));
t630 = qJD(1) * t562;
t504 = t561 * t630 + t629;
t602 = t558 * t630;
t627 = qJD(2) * t561;
t506 = -t602 + t627;
t555 = sin(pkin(10));
t556 = cos(pkin(10));
t583 = -t504 * t556 - t506 * t555;
t444 = t504 * t555 - t506 * t556;
t557 = sin(qJ(6));
t662 = t444 * t557;
t392 = t560 * t583 + t662;
t559 = sin(qJ(2));
t631 = qJD(1) * t559;
t537 = qJD(4) + t631;
t531 = qJD(6) + t537;
t663 = t392 * t531;
t616 = qJD(1) * qJD(2);
t601 = t559 * t616;
t457 = -qJD(4) * t504 + t558 * t601;
t563 = -pkin(2) - pkin(8);
t599 = -qJ(3) * t559 - pkin(1);
t500 = t562 * t563 + t599;
t471 = t500 * qJD(1);
t542 = pkin(7) * t631;
t679 = qJD(3) + t542;
t618 = pkin(3) * t631 + t679;
t475 = qJD(2) * t563 + t618;
t420 = t471 * t561 + t475 * t558;
t536 = pkin(2) * t601;
t587 = pkin(8) * t559 - qJ(3) * t562;
t625 = qJD(3) * t559;
t568 = qJD(2) * t587 - t625;
t452 = qJD(1) * t568 + t536;
t600 = t562 * t616;
t535 = pkin(7) * t600;
t496 = pkin(3) * t600 + t535;
t593 = -t452 * t558 + t561 * t496;
t567 = -qJD(4) * t420 + t593;
t364 = pkin(4) * t600 - qJ(5) * t457 - qJD(5) * t506 + t567;
t623 = qJD(4) * t561;
t458 = qJD(2) * t623 - qJD(4) * t602 - t561 * t601;
t612 = -t561 * t452 - t475 * t623 - t558 * t496;
t624 = qJD(4) * t558;
t570 = -t471 * t624 - t612;
t366 = -qJ(5) * t458 - qJD(5) * t504 + t570;
t349 = t556 * t364 - t366 * t555;
t411 = t457 * t556 - t458 * t555;
t347 = pkin(5) * t600 - pkin(9) * t411 + t349;
t350 = t555 * t364 + t556 * t366;
t410 = -t457 * t555 - t458 * t556;
t348 = pkin(9) * t410 + t350;
t419 = -t471 * t558 + t561 * t475;
t405 = -qJ(5) * t506 + t419;
t397 = pkin(4) * t537 + t405;
t406 = -qJ(5) * t504 + t420;
t651 = t556 * t406;
t368 = t555 * t397 + t651;
t676 = pkin(9) * t583;
t359 = t368 + t676;
t620 = qJD(6) * t557;
t358 = t359 * t620;
t543 = pkin(7) * t630;
t512 = pkin(3) * t630 + t543;
t552 = qJD(2) * qJ(3);
t490 = t552 + t512;
t449 = pkin(4) * t504 + qJD(5) + t490;
t398 = -pkin(5) * t583 + t449;
t686 = -t557 * t347 - t560 * t348 - t398 * t392 + t358;
t678 = -t560 * t444 + t557 * t583;
t685 = MDP(28) * t600 + (-t392 ^ 2 + t678 ^ 2) * MDP(25) - t392 * MDP(24) * t678;
t664 = t678 * t531;
t547 = pkin(2) * t631;
t483 = qJD(1) * t587 + t547;
t592 = -t483 * t558 + t561 * t512;
t643 = qJ(5) - t563;
t650 = t558 * t559;
t683 = -(pkin(4) * t562 - qJ(5) * t650) * qJD(1) - t592 - qJD(5) * t561 + t624 * t643;
t516 = t643 * t561;
t609 = t561 * t631;
t635 = t561 * t483 + t558 * t512;
t682 = qJ(5) * t609 + qJD(4) * t516 + qJD(5) * t558 + t635;
t582 = t555 * t561 + t556 * t558;
t636 = t537 * t582;
t595 = t560 * t347 - t557 * t348;
t681 = -t398 * t678 + t595;
t667 = pkin(3) + pkin(7);
t680 = pkin(9) * t444;
t603 = t556 * t623;
t652 = t555 * t558;
t637 = -t555 * t624 + t556 * t609 - t631 * t652 + t603;
t677 = -0.2e1 * t616;
t553 = t559 ^ 2;
t554 = t562 ^ 2;
t675 = (t553 - t554) * MDP(5);
t640 = t555 * t682 + t556 * t683;
t639 = t555 * t683 - t556 * t682;
t672 = t636 * t560;
t524 = t667 * t559;
t634 = t561 * t500 + t558 * t524;
t581 = -t556 * t561 + t652;
t443 = -t557 * t582 - t560 * t581;
t622 = qJD(4) * t562;
t605 = t558 * t622;
t610 = pkin(4) * t561 + pkin(3);
t628 = qJD(2) * t559;
t671 = (-pkin(7) - t610) * t628 - pkin(4) * t605;
t670 = t559 * t627 + t605;
t594 = -t560 * t410 + t411 * t557;
t356 = qJD(6) * t678 + t594;
t669 = pkin(4) * t623 + t610 * t631 + t679;
t399 = t555 * t406;
t367 = t556 * t397 - t399;
t668 = -t349 * t581 + t350 * t582 - t367 * t636 + t368 * t637;
t666 = pkin(4) * t555;
t665 = qJD(2) * pkin(2);
t661 = t457 * t561;
t511 = t667 * t628;
t551 = qJD(2) * qJD(3);
t480 = -qJD(1) * t511 + t551;
t660 = t480 * t558;
t659 = t480 * t561;
t657 = t504 * t537;
t656 = t506 * t537;
t655 = t506 * t562;
t654 = t537 * t559;
t653 = t537 * t563;
t564 = qJD(2) ^ 2;
t649 = t559 * t564;
t357 = pkin(5) * t537 + t367 + t680;
t648 = t560 * t357;
t647 = t561 * t562;
t646 = t562 * t564;
t565 = qJD(1) ^ 2;
t645 = t562 * t565;
t644 = t558 * pkin(4) + qJ(3);
t546 = pkin(2) * t628;
t460 = t546 + t568;
t626 = qJD(2) * t562;
t513 = t667 * t626;
t495 = t561 * t513;
t596 = qJ(5) * t562 - t500;
t621 = qJD(5) * t562;
t379 = pkin(4) * t626 + t495 + t596 * t623 + (-qJ(5) * t628 - qJD(4) * t524 - t460 + t621) * t558;
t569 = t561 * t460 - t500 * t624 + t558 * t513 + t524 * t623;
t384 = qJ(5) * t670 - t561 * t621 + t569;
t354 = t555 * t379 + t556 * t384;
t442 = -t557 * t581 + t560 * t582;
t642 = -qJD(6) * t442 - t557 * t637 - t672;
t641 = t443 * qJD(6) - t557 * t636 + t560 * t637;
t373 = t556 * t405 - t399;
t508 = t561 * t524;
t433 = pkin(4) * t559 + t558 * t596 + t508;
t438 = -qJ(5) * t647 + t634;
t386 = t555 * t433 + t556 * t438;
t638 = pkin(5) * t637 + t669;
t515 = t643 * t558;
t451 = -t556 * t515 - t555 * t516;
t525 = t667 * t562;
t518 = -pkin(2) * t562 + t599;
t491 = qJD(1) * t518;
t619 = qJD(6) * t560;
t615 = t561 * t654;
t614 = t559 * t645;
t613 = t557 * t410 + t560 * t411 + t583 * t619;
t611 = pkin(4) * t647 + t525;
t608 = t558 * t628;
t606 = t537 * t623;
t604 = t561 * t622;
t598 = pkin(1) * t677;
t597 = qJD(3) - t665;
t353 = t556 * t379 - t384 * t555;
t372 = -t405 * t555 - t651;
t385 = t556 * t433 - t438 * t555;
t450 = t515 * t555 - t556 * t516;
t527 = t559 * t600;
t428 = -pkin(9) * t582 + t451;
t590 = pkin(5) * t630 - pkin(9) * t636 + qJD(6) * t428 - t640;
t427 = pkin(9) * t581 + t450;
t589 = pkin(9) * t637 - qJD(6) * t427 - t639;
t588 = qJD(6) * t581 - t637;
t345 = t557 * t357 + t560 * t359;
t479 = t582 * t562;
t375 = pkin(5) * t559 + pkin(9) * t479 + t385;
t478 = t581 * t562;
t376 = pkin(9) * t478 + t386;
t585 = t375 * t557 + t376 * t560;
t584 = t560 * t478 + t479 * t557;
t430 = t478 * t557 - t479 * t560;
t580 = -qJD(1) * t554 + t654;
t579 = -0.2e1 * qJD(2) * t491;
t578 = t537 * t558;
t541 = pkin(4) * t556 + pkin(5);
t577 = t541 * t557 + t560 * t666;
t576 = t541 * t560 - t557 * t666;
t571 = -qJ(3) * t626 - t625;
t468 = qJD(1) * t571 + t536;
t485 = t546 + t571;
t575 = pkin(7) * t564 + qJD(1) * t485 + t468;
t573 = t490 * t559 + t563 * t626;
t355 = t444 * t620 + t613;
t431 = pkin(4) * t458 + t480;
t514 = pkin(7) * t601 - t551;
t517 = t542 + t597;
t523 = -t543 - t552;
t566 = -t514 * t562 + (t517 * t562 + (t523 + t543) * t559) * qJD(2);
t529 = t561 * t600;
t509 = -qJ(3) * t630 + t547;
t472 = t491 * t631;
t470 = pkin(5) * t582 + t644;
t441 = -pkin(5) * t478 + t611;
t435 = -t555 * t670 - t556 * t608 + t562 * t603;
t434 = qJD(4) * t479 - t581 * t628;
t415 = pkin(4) * t506 - pkin(5) * t444;
t396 = -pkin(5) * t434 + t671;
t383 = -pkin(5) * t410 + t431;
t370 = qJD(6) * t430 - t560 * t434 - t435 * t557;
t369 = qJD(6) * t584 + t434 * t557 - t435 * t560;
t361 = t373 + t680;
t360 = t372 - t676;
t352 = pkin(9) * t434 + t354;
t351 = pkin(5) * t626 + pkin(9) * t435 + t353;
t344 = -t359 * t557 + t648;
t1 = [((-t504 * t558 + t506 * t561) * t628 + (-t661 + t458 * t558 + (t504 * t561 + t506 * t558) * qJD(4)) * t562) * MDP(16) + (-t537 * t604 + t457 * t559 + (t558 * t580 + t655) * qJD(2)) * MDP(17) + 0.2e1 * MDP(4) * t527 + (t349 * t385 + t350 * t386 + t367 * t353 + t368 * t354 + t431 * t611 + t449 * t671) * MDP(23) + ((-t460 * t558 + t495) * t537 - t511 * t504 + t525 * t458 + (-t490 * t627 + t593) * t559 + (-t420 * t559 - t537 * t634) * qJD(4) + (-t490 * t624 + t659 + ((-t500 * t558 + t508) * qJD(1) + t419) * qJD(2)) * t562) * MDP(20) - MDP(7) * t649 + (pkin(7) * t649 + t562 * t598) * MDP(10) + (-pkin(7) * t646 + t559 * t598) * MDP(9) + (t537 * t626 + t527) * MDP(19) + (t531 * t626 + t527) * MDP(28) + (-t457 * t558 * t562 + (-t604 + t608) * t506) * MDP(15) + (-t559 * t575 + t562 * t579) * MDP(13) + (t559 * t579 + t562 * t575) * MDP(12) + (t355 * t584 - t356 * t430 + t369 * t392 - t370 * t678) * MDP(25) + ((t351 * t560 - t352 * t557) * t531 + t595 * t559 - t396 * t392 + t441 * t356 - t383 * t584 + t398 * t370 + (-t345 * t559 - t531 * t585) * qJD(6) + ((t375 * t560 - t376 * t557) * qJD(1) + t344) * t626) * MDP(29) + (-t356 * t559 - t370 * t531 + (qJD(1) * t584 + t392) * t626) * MDP(27) + (t355 * t559 + t369 * t531 + (qJD(1) * t430 + t678) * t626) * MDP(26) + (t441 * t355 + t358 * t559 + t398 * t369 + t383 * t430 + t396 * t678 + (-(-qJD(6) * t376 + t351) * t531 - t347 * t559) * t557 + (-(qJD(6) * t375 + t352) * t531 - (qJD(6) * t357 + t348) * t559) * t560 + (-qJD(1) * t585 - t345) * t626) * MDP(30) + (t355 * t430 + t369 * t678) * MDP(24) + t675 * t677 + (-t569 * t537 - t511 * t506 + t525 * t457 + ((qJD(2) * t490 + qJD(4) * t471) * t558 + t612) * t559 + (-t490 * t623 - t660 + (-qJD(1) * t634 - t420) * qJD(2)) * t562) * MDP(21) + (t349 * t479 + t350 * t478 + t353 * t444 + t354 * t583 + t367 * t435 + t368 * t434 - t385 * t411 + t386 * t410) * MDP(22) + MDP(6) * t646 + t566 * MDP(11) + (pkin(7) * t566 + t468 * t518 + t485 * t491) * MDP(14) + (t537 * t605 - t458 * t559 + (-t504 * t562 + t561 * t580) * qJD(2)) * MDP(18); t565 * t675 + ((-t523 - t552) * t559 + (-t517 + t597) * t562) * qJD(1) * MDP(11) + t472 * MDP(12) + (0.2e1 * t551 + (t491 * t562 + t509 * t559) * qJD(1)) * MDP(13) + (-qJ(3) * t514 - qJD(3) * t523 - t491 * t509 + (-t523 * t559 + (-t517 - t665) * t562) * qJD(1) * pkin(7)) * MDP(14) + (-t506 * t578 + t661) * MDP(15) + ((-t458 - t656) * t561 + (-t457 + t657) * t558) * MDP(16) + (-t537 * t624 + t529 + (-t537 * t650 - t655) * qJD(1)) * MDP(17) + (-t606 + (-t615 + (t504 - t629) * t562) * qJD(1)) * MDP(18) + (qJ(3) * t458 + t660 - t592 * t537 + t618 * t504 + (t490 * t561 - t558 * t653) * qJD(4) + (-t419 * t562 + t561 * t573) * qJD(1)) * MDP(20) + (qJ(3) * t457 + t659 + t635 * t537 + t618 * t506 + (-t490 * t558 - t561 * t653) * qJD(4) + (t420 * t562 - t558 * t573) * qJD(1)) * MDP(21) + (t410 * t451 - t411 * t450 + t444 * t640 + t583 * t639 - t668) * MDP(22) + (t349 * t450 + t350 * t451 + t640 * t367 + t639 * t368 + t431 * t644 + t449 * t669) * MDP(23) + (t355 * t443 + t642 * t678) * MDP(24) + (-t355 * t442 - t356 * t443 + t392 * t642 - t641 * t678) * MDP(25) + (t470 * t356 + t383 * t442 - t392 * t638 + t641 * t398) * MDP(29) + (t470 * t355 + t383 * t443 + t642 * t398 + t638 * t678) * MDP(30) - MDP(4) * t614 + (t642 * MDP(26) - t641 * MDP(27) + (t557 * t589 - t560 * t590) * MDP(29) + (t557 * t590 + t560 * t589) * MDP(30)) * t531 + (-t509 * MDP(12) + (qJD(2) * t443 - t678) * MDP(26) + (-qJD(2) * t442 - t392) * MDP(27) + ((t427 * t560 - t428 * t557) * qJD(2) - t344) * MDP(29) + (-(t427 * t557 + t428 * t560) * qJD(2) + t345) * MDP(30) - t537 * MDP(19) - t531 * MDP(28)) * t630 + (MDP(9) * t559 * t565 + MDP(10) * t645) * pkin(1); MDP(12) * t614 + (-t553 * t565 - t564) * MDP(13) + (qJD(2) * t523 + t472 + t535) * MDP(14) + (-qJD(2) * t504 - t537 * t578 + t529) * MDP(20) + (-t606 - qJD(2) * t506 + (-t558 * t626 - t615) * qJD(1)) * MDP(21) + (t410 * t582 + t411 * t581 - t444 * t636 + t583 * t637) * MDP(22) + (-t449 * qJD(2) + t668) * MDP(23) + ((t557 * t588 - t582 * t619 - t672) * t531 + (t443 * t630 + t392) * qJD(2)) * MDP(29) + ((t588 * t560 + (qJD(6) * t582 + t636) * t557) * t531 + (-t442 * t630 - t678) * qJD(2)) * MDP(30); t506 * t504 * MDP(15) + (-t504 ^ 2 + t506 ^ 2) * MDP(16) + (t457 + t657) * MDP(17) + (-t458 + t656) * MDP(18) + MDP(19) * t600 + (t420 * t537 - t490 * t506 + t567) * MDP(20) + (t419 * t537 + t490 * t504 - t570) * MDP(21) + ((t410 * t555 - t411 * t556) * pkin(4) + (t367 - t373) * t583 + (-t368 - t372) * t444) * MDP(22) + (-t367 * t372 - t368 * t373 + (t349 * t556 + t350 * t555 - t449 * t506) * pkin(4)) * MDP(23) + (t355 - t663) * MDP(26) + (-t356 + t664) * MDP(27) + (t576 * t600 - (t360 * t560 - t361 * t557) * t531 + t415 * t392 + (-t531 * t577 - t345) * qJD(6) + t681) * MDP(29) + (-t577 * t600 + (t360 * t557 + t361 * t560) * t531 - t415 * t678 + (-t531 * t576 - t648) * qJD(6) + t686) * MDP(30) + t685; (-t444 ^ 2 - t583 ^ 2) * MDP(22) + (-t367 * t444 - t368 * t583 + t431) * MDP(23) + (t356 + t664) * MDP(29) + (t355 + t663) * MDP(30); (t613 - t663) * MDP(26) + (-t594 + t664) * MDP(27) + (t345 * t531 + t681) * MDP(29) + (t344 * t531 + t686) * MDP(30) + (MDP(26) * t662 - MDP(27) * t678 - MDP(29) * t345 - MDP(30) * t648) * qJD(6) + t685;];
tauc  = t1;
