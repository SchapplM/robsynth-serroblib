% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:31
% EndTime: 2019-03-10 00:58:44
% DurationCPUTime: 7.92s
% Computational Cost: add. (10536->483), mult. (26443->621), div. (0->0), fcn. (19685->8), ass. (0->236)
t571 = cos(qJ(5));
t646 = qJD(5) * t571;
t569 = sin(qJ(3));
t570 = sin(qJ(2));
t654 = qJD(1) * t570;
t632 = t569 * t654;
t573 = cos(qJ(3));
t574 = cos(qJ(2));
t653 = qJD(1) * t574;
t634 = t573 * t653;
t520 = -t632 + t634;
t521 = -t569 * t653 - t573 * t654;
t568 = sin(qJ(4));
t572 = cos(qJ(4));
t712 = t572 * t520 + t521 * t568;
t725 = t712 * t571;
t732 = t646 - t725;
t602 = t520 * t568 - t572 * t521;
t567 = sin(qJ(5));
t647 = qJD(5) * t567;
t564 = qJD(2) + qJD(3);
t679 = t569 * t574;
t535 = t570 * t573 + t679;
t707 = qJD(1) * t535;
t581 = t564 * t707;
t580 = t568 * t581;
t641 = qJD(1) * qJD(2);
t631 = t574 * t641;
t496 = qJD(3) * t634 - t564 * t632 + t573 * t631;
t677 = t572 * t496;
t579 = -t580 + t677;
t577 = t712 * qJD(4) + t579;
t628 = qJD(4) + t564;
t723 = -qJD(5) * t628 - t577;
t406 = t571 * t723 + t602 * t647;
t404 = t406 * t567;
t474 = t567 * t628 + t571 * t602;
t407 = t474 * qJD(5) + t567 * t577;
t623 = t568 * t496 + t572 * t581;
t419 = t602 * qJD(4) + t623;
t416 = t571 * t419;
t472 = t567 * t602 - t571 * t628;
t642 = qJD(5) - t712;
t414 = t567 * t419;
t670 = t642 * t646 + t414;
t730 = t642 * t567;
t731 = -t623 * MDP(21) - t712 ^ 2 * MDP(19) + (-t564 * t712 + t579) * MDP(20) + (-MDP(18) * t712 + MDP(19) * t602 + MDP(21) * t564 - MDP(29) * t642) * t602 + (t474 * t732 - t404) * MDP(25) + (-t474 * t602 - t642 * t725 + t670) * MDP(27) + (-t406 * t571 - t567 * t407 - t472 * t732 - t474 * t730) * MDP(26) + (t472 * t602 - t642 * t730 + t416) * MDP(28);
t702 = pkin(7) + pkin(8);
t635 = qJD(2) * t702;
t611 = qJD(1) * t635;
t530 = t574 * t611;
t545 = t702 * t574;
t538 = qJD(1) * t545;
t651 = qJD(3) * t569;
t620 = -t569 * t530 - t538 * t651;
t544 = t702 * t570;
t536 = qJD(1) * t544;
t528 = qJD(2) * pkin(2) - t536;
t529 = t570 * t611;
t714 = t573 * (qJD(3) * t528 - t529);
t427 = -pkin(9) * t581 + t620 + t714;
t517 = t521 * pkin(9);
t522 = t569 * t538;
t622 = t573 * t528 - t522;
t469 = t517 + t622;
t462 = pkin(3) * t564 + t469;
t613 = qJD(4) * t462 + t427;
t526 = t573 * t538;
t601 = -t528 * t569 - t526;
t621 = t569 * t529 - t573 * t530;
t587 = qJD(3) * t601 + t621;
t428 = -pkin(9) * t496 + t587;
t698 = pkin(9) * t520;
t470 = -t601 + t698;
t649 = qJD(4) * t568;
t624 = t568 * t428 - t470 * t649;
t375 = t613 * t572 + t624;
t558 = -pkin(2) * t574 - pkin(1);
t543 = t558 * qJD(1);
t504 = -pkin(3) * t520 + t543;
t684 = t504 * t712;
t729 = -t375 - t684;
t467 = t568 * t470;
t430 = t469 * t572 - t467;
t648 = qJD(4) * t572;
t724 = -pkin(3) * t648 + t430;
t619 = t536 * t569 - t526;
t477 = t619 - t698;
t657 = -t573 * t536 - t522;
t478 = t517 + t657;
t557 = pkin(2) * t573 + pkin(3);
t681 = t568 * t569;
t713 = t477 * t568 + t478 * t572 - t557 * t648 - (-t569 * t649 + (t572 * t573 - t681) * qJD(3)) * pkin(2);
t726 = t712 * t567;
t711 = qJ(6) * t726 + t571 * qJD(6);
t680 = t569 * t572;
t722 = (t569 * t648 + (t568 * t573 + t680) * qJD(3)) * pkin(2) + t477 * t572;
t455 = pkin(4) * t602 - pkin(10) * t712;
t563 = t571 * qJ(6);
t609 = pkin(5) * t602 - t563 * t712;
t719 = -0.2e1 * t641;
t716 = MDP(4) * t570;
t715 = MDP(5) * (t570 ^ 2 - t574 ^ 2);
t483 = -pkin(9) * t535 - t544 * t573 - t545 * t569;
t534 = t569 * t570 - t573 * t574;
t600 = t544 * t569 - t545 * t573;
t484 = -pkin(9) * t534 - t600;
t451 = -t483 * t572 + t568 * t484;
t662 = -t478 * t568 + t557 * t649 + t722;
t588 = t520 * t648 + t521 * t649 + t677;
t710 = -t580 + t588;
t701 = pkin(3) * t521;
t443 = t455 - t701;
t709 = t567 * t443 + t571 * t724;
t559 = pkin(2) * t654;
t438 = t443 + t559;
t708 = t567 * t438 + t571 * t713;
t376 = t568 * t427 - t572 * t428 + t462 * t649 + t470 * t648;
t594 = -t504 * t602 - t376;
t468 = t572 * t470;
t422 = t568 * t462 + t468;
t418 = pkin(10) * t628 + t422;
t435 = -pkin(4) * t712 - pkin(10) * t602 + t504;
t400 = -t418 * t567 + t571 * t435;
t421 = t572 * t462 - t467;
t417 = -pkin(4) * t628 - t421;
t598 = -t376 * t571 - t400 * t602 + t417 * t647;
t401 = t418 * t571 + t435 * t567;
t413 = t417 * t646;
t608 = t376 * t567 + t401 * t602 + t413;
t703 = t474 ^ 2;
t700 = pkin(3) * t572;
t699 = pkin(5) * t571;
t697 = -qJ(6) - pkin(10);
t696 = pkin(3) * qJD(4);
t386 = -qJ(6) * t474 + t400;
t381 = pkin(5) * t642 + t386;
t695 = t381 * t571;
t694 = t417 * t712;
t499 = -t534 * t568 + t535 * t572;
t686 = t499 * t567;
t685 = t499 * t571;
t682 = t543 * t521;
t575 = qJD(2) ^ 2;
t678 = t570 * t575;
t452 = t483 * t568 + t484 * t572;
t449 = t571 * t452;
t676 = t574 * t575;
t576 = qJD(1) ^ 2;
t675 = t574 * t576;
t514 = pkin(2) * t680 + t557 * t568 + pkin(10);
t674 = -qJ(6) - t514;
t555 = pkin(3) * t568 + pkin(10);
t673 = -qJ(6) - t555;
t671 = t381 - t386;
t669 = t571 * t421 + t567 * t455;
t498 = t572 * t534 + t535 * t568;
t508 = pkin(3) * t534 + t558;
t450 = pkin(4) * t498 - pkin(10) * t499 + t508;
t666 = t567 * t450 + t449;
t618 = qJD(5) * t674;
t665 = t567 * t618 - t708 + t711;
t434 = t571 * t438;
t664 = t571 * t618 - t434 - t609 + (-qJD(6) + t713) * t567;
t617 = qJD(5) * t673;
t661 = t567 * t617 - t709 + t711;
t442 = t571 * t443;
t660 = t571 * t617 - t442 - t609 + (-qJD(6) + t724) * t567;
t629 = qJD(5) * t697;
t659 = t567 * t629 - t669 + t711;
t625 = -t421 * t567 + t571 * t455;
t658 = -qJD(6) * t567 + t571 * t629 - t609 - t625;
t656 = (t647 - t726) * pkin(5);
t652 = qJD(2) * t570;
t650 = qJD(3) * t573;
t561 = pkin(2) * t652;
t590 = t535 * qJD(3);
t503 = qJD(2) * t535 + t590;
t537 = t570 * t635;
t539 = t574 * t635;
t593 = -t573 * t537 - t569 * t539 - t544 * t650 - t545 * t651;
t447 = -pkin(9) * t503 + t593;
t502 = t564 * t534;
t585 = qJD(3) * t600 + t537 * t569 - t573 * t539;
t448 = pkin(9) * t502 + t585;
t391 = -qJD(4) * t451 + t447 * t572 + t448 * t568;
t439 = -qJD(4) * t498 - t502 * t572 - t503 * t568;
t440 = qJD(4) * t499 - t502 * t568 + t572 * t503;
t493 = pkin(3) * t503 + t561;
t398 = pkin(4) * t440 - pkin(10) * t439 + t493;
t638 = t571 * t391 + t567 * t398 + t450 * t646;
t636 = -pkin(4) - t699;
t633 = t499 * t646;
t630 = -pkin(2) * t564 - t528;
t627 = pkin(1) * t719;
t513 = pkin(2) * t681 - t557 * t572 - pkin(4);
t429 = t469 * t568 + t468;
t610 = pkin(3) * t649 - t429;
t387 = -qJ(6) * t472 + t401;
t606 = -t387 * t567 - t695;
t605 = -t419 * t514 - t694;
t604 = -t419 * t555 - t694;
t599 = -qJ(6) * t439 - qJD(6) * t499;
t597 = -t543 * t520 - t620;
t596 = t439 * t567 + t633;
t595 = t439 * t571 - t499 * t647;
t371 = pkin(5) * t407 + t376;
t552 = qJD(2) * t559;
t388 = t419 * pkin(4) - t588 * pkin(10) + t552 + (pkin(10) * t568 + pkin(3)) * qJD(1) * (qJD(2) * t679 + t570 * t650 + t573 * t652 + t574 * t651);
t592 = t571 * t375 + t567 * t388 - t418 * t647 + t435 * t646;
t385 = t571 * t388;
t586 = -qJD(5) * t401 - t375 * t567 + t385;
t392 = qJD(4) * t452 + t447 * t568 - t448 * t572;
t365 = pkin(5) * t419 + qJ(6) * t406 - qJD(6) * t474 + t586;
t367 = -qJ(6) * t407 - qJD(6) * t472 + t592;
t366 = t367 * t571;
t583 = qJD(5) * t606 - t365 * t567 + t381 * t725 + t387 * t726 + t366;
t578 = t521 * t520 * MDP(11) + (-t520 * t564 + t496) * MDP(13) + (-t521 * t564 - t581) * MDP(14) + (-t520 ^ 2 + t521 ^ 2) * MDP(12) + t731;
t556 = -pkin(4) - t700;
t542 = pkin(10) * t571 + t563;
t541 = t697 * t567;
t532 = t555 * t571 + t563;
t531 = t673 * t567;
t507 = t559 - t701;
t506 = t514 * t571 + t563;
t505 = t674 * t567;
t479 = pkin(3) * t581 + t552;
t471 = t472 ^ 2;
t446 = t571 * t450;
t408 = t472 * pkin(5) + qJD(6) + t417;
t402 = -qJ(6) * t686 + t666;
t397 = t571 * t398;
t395 = pkin(5) * t498 - t452 * t567 - t499 * t563 + t446;
t369 = -qJ(6) * t633 + (-qJD(5) * t452 + t599) * t567 + t638;
t368 = pkin(5) * t440 - t391 * t567 + t397 + t599 * t571 + (-t449 + (qJ(6) * t499 - t450) * t567) * qJD(5);
t1 = [0.2e1 * t631 * t716 + (t602 * t439 + t499 * t710) * MDP(18) + (t504 * t439 + t479 * t499 + t493 * t602 + t508 * t577) * MDP(24) + (-(-t452 * t647 + t638) * t642 - t666 * t419 - t592 * t498 - t401 * t440 + t392 * t474 - t451 * t406 + t376 * t685 + t595 * t417) * MDP(31) + (-t407 * t498 - t414 * t499 - t440 * t472 - t596 * t642) * MDP(28) + (-t406 * t498 + t416 * t499 + t440 * t474 + t595 * t642) * MDP(27) + ((-t452 * t646 + t397) * t642 + t446 * t419 + (-t418 * t646 + t385) * t498 + t400 * t440 + t392 * t472 + t451 * t407 + t499 * t413 + ((-qJD(5) * t450 - t391) * t642 - t452 * t419 + (-qJD(5) * t435 - t375) * t498 + t376 * t499 + t417 * t439) * t567) * MDP(30) + (t419 * t498 + t440 * t642) * MDP(29) + (MDP(20) * t439 - MDP(21) * t440 - MDP(23) * t392 - MDP(24) * t391) * t628 + MDP(6) * t676 + (-t496 * t534 - t502 * t520 + t521 * t503 - t535 * t581) * MDP(12) + (-t499 * t419 + t439 * t712 - t602 * t440 - t498 * t710) * MDP(19) + (t508 * t419 + t504 * t440 + t479 * t498 - t493 * t712) * MDP(23) + (t558 * t496 - t543 * t502 + (-t521 + t707) * t561) * MDP(17) + (t496 * t535 + t502 * t521) * MDP(11) + (-t502 * MDP(13) - t503 * MDP(14) + MDP(16) * t585 - MDP(17) * t593) * t564 + (-t520 * t561 + t543 * t503 + (t558 * t590 + (t570 * pkin(2) * t534 + t535 * t558) * qJD(2)) * qJD(1)) * MDP(16) + ((-t472 * t571 - t474 * t567) * t439 + (t404 - t407 * t571 + (t472 * t567 - t474 * t571) * qJD(5)) * t499) * MDP(26) + (-t406 * t685 + t474 * t595) * MDP(25) + (t367 * t402 + t387 * t369 + t365 * t395 + t381 * t368 + t371 * (pkin(5) * t686 + t451) + t408 * (pkin(5) * t596 + t392)) * MDP(33) - MDP(7) * t678 + (pkin(7) * t678 + t574 * t627) * MDP(10) + (-pkin(7) * t676 + t570 * t627) * MDP(9) + t715 * t719 + (-t368 * t474 - t369 * t472 + t395 * t406 - t402 * t407 + t606 * t439 + (-t365 * t571 - t367 * t567 + (t381 * t567 - t387 * t571) * qJD(5)) * t499) * MDP(32); (t406 * t505 - t407 * t506 - t472 * t665 - t474 * t664 + t583) * MDP(32) + (t521 * t559 + t657 * t564 + (qJD(3) * t630 + t529) * t573 + t597) * MDP(17) + (t520 * t559 + t682 - t619 * t564 + (t569 * t630 - t526) * qJD(3) + t621) * MDP(16) + (-t513 * t406 + t605 * t571 + t662 * t474 + (t514 * t647 + t708) * t642 + t608) * MDP(31) + (t513 * t407 + t605 * t567 + t662 * t472 + (-t514 * t646 + t567 * t713 - t434) * t642 + t598) * MDP(30) + (-t507 * t602 + t628 * t713 + t729) * MDP(24) + t576 * t715 + t578 + (t367 * t506 + t365 * t505 + t371 * (t513 - t699) + t665 * t387 + t664 * t381 + ((qJD(4) * t557 - t478) * t568 + t656 + t722) * t408) * MDP(33) - t675 * t716 + (t507 * t712 - t628 * t662 + t594) * MDP(23) + (MDP(9) * t570 * t576 + MDP(10) * t675) * pkin(1); (t406 * t531 - t407 * t532 - t472 * t661 - t474 * t660 + t583) * MDP(32) + (-t556 * t406 + t604 * t571 + t610 * t474 + (t555 * t647 + t709) * t642 + t608) * MDP(31) + (-t564 * t601 + t587 + t682) * MDP(16) + (t564 * t622 + t597 - t714) * MDP(17) + (t556 * t407 + t604 * t567 + t610 * t472 + (-t555 * t646 + t567 * t724 - t442) * t642 + t598) * MDP(30) + (t602 * t701 - t684 + t430 * t628 + (-t628 * t696 - t613) * t572 - t624) * MDP(24) + t578 + (t367 * t532 + t365 * t531 + t371 * (t636 - t700) + (-t468 + (-t469 + t696) * t568 + t656) * t408 + t661 * t387 + t660 * t381) * MDP(33) + (t429 * t628 + (-t521 * t712 - t628 * t649) * pkin(3) + t594) * MDP(23); (t422 * t628 + t594) * MDP(23) + (t421 * t628 + t729) * MDP(24) + (-pkin(4) * t407 - pkin(10) * t670 - t417 * t726 - t422 * t472 - t625 * t642 + t598) * MDP(30) + (pkin(4) * t406 + t669 * t642 - t422 * t474 - t417 * t725 + (t642 * t647 - t416) * pkin(10) + t608) * MDP(31) + (t406 * t541 - t407 * t542 + t366 - t658 * t474 - t659 * t472 - t642 * t695 + (-t387 * t642 - t365) * t567) * MDP(32) + (t367 * t542 + t365 * t541 + t371 * t636 + (pkin(5) * t730 - t422) * t408 + t659 * t387 + t658 * t381) * MDP(33) + t731; t474 * t472 * MDP(25) + (-t471 + t703) * MDP(26) + (t472 * t642 - t406) * MDP(27) + (t474 * t642 + t567 * t723 - t602 * t646) * MDP(28) + t419 * MDP(29) + (t401 * t642 - t417 * t474 + t586) * MDP(30) + (t400 * t642 + t417 * t472 - t592) * MDP(31) + (pkin(5) * t406 - t472 * t671) * MDP(32) + (t671 * t387 + (-t408 * t474 + t365) * pkin(5)) * MDP(33); (-t471 - t703) * MDP(32) + (t381 * t474 + t387 * t472 + t371) * MDP(33);];
tauc  = t1;
