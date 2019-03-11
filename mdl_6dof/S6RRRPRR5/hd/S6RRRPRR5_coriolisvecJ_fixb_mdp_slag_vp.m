% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:23
% EndTime: 2019-03-09 18:23:37
% DurationCPUTime: 7.84s
% Computational Cost: add. (6081->501), mult. (14301->635), div. (0->0), fcn. (10430->8), ass. (0->231)
t571 = sin(qJ(6));
t572 = sin(qJ(5));
t575 = cos(qJ(6));
t576 = cos(qJ(5));
t527 = t571 * t576 + t572 * t575;
t573 = sin(qJ(3));
t574 = sin(qJ(2));
t577 = cos(qJ(2));
t684 = cos(qJ(3));
t530 = t573 * t577 + t574 * t684;
t644 = qJD(1) * t530;
t688 = qJD(5) + qJD(6);
t652 = (t644 + t688) * t527;
t567 = qJD(2) + qJD(3);
t663 = t573 * t574;
t604 = t567 * t663;
t628 = t684 * t577;
t611 = qJD(1) * t628;
t649 = t567 * t611;
t469 = qJD(1) * t604 - t649;
t698 = qJD(5) + t644;
t498 = qJD(6) + t698;
t638 = qJD(6) * t571;
t641 = qJD(5) * t572;
t661 = t575 * t576;
t689 = -t571 * t572 + t661;
t699 = -t571 * t641 - t572 * t638 + t689 * t644 + t661 * t688;
t705 = t527 * t469 - t498 * t699;
t643 = qJD(1) * t574;
t627 = t573 * t643;
t510 = -t611 + t627;
t484 = -t576 * t510 + t567 * t572;
t704 = t484 * t698;
t528 = -t628 + t663;
t465 = t689 * t528;
t703 = -t469 * t689 - t652 * t498;
t486 = t510 * t572 + t567 * t576;
t670 = t486 * t571;
t424 = t575 * t484 + t670;
t702 = t424 * t498;
t603 = t484 * t571 - t575 * t486;
t701 = t498 * t603;
t471 = pkin(3) * t644 + qJ(4) * t510;
t506 = t644 * pkin(9);
t438 = t471 + t506;
t686 = pkin(3) + pkin(9);
t639 = qJD(5) * t686;
t700 = t438 + t639;
t685 = -pkin(8) - pkin(7);
t543 = t685 * t577;
t536 = qJD(1) * t543;
t514 = t573 * t536;
t542 = t685 * t574;
t534 = qJD(1) * t542;
t677 = qJD(2) * pkin(2);
t518 = t534 + t677;
t473 = -t684 * t518 - t514;
t635 = qJD(4) + t473;
t680 = t644 * pkin(4);
t636 = t680 + t635;
t419 = -t567 * t686 + t636;
t559 = -pkin(2) * t577 - pkin(1);
t541 = t559 * qJD(1);
t585 = -qJ(4) * t644 + t541;
t422 = t510 * t686 + t585;
t386 = t419 * t572 + t422 * t576;
t373 = -pkin(10) * t484 + t386;
t371 = t373 * t638;
t517 = t684 * t536;
t474 = t573 * t518 - t517;
t681 = t510 * pkin(4);
t443 = t474 - t681;
t565 = t567 * qJ(4);
t433 = t443 + t565;
t406 = pkin(5) * t484 + t433;
t697 = t406 * t424 + t371;
t481 = t567 * t530;
t470 = t481 * qJD(1);
t640 = qJD(5) * t576;
t414 = t572 * t470 + t510 * t640 - t567 * t641;
t633 = qJD(1) * qJD(2);
t625 = t574 * t633;
t596 = pkin(2) * t625 + qJ(4) * t469 - qJD(4) * t644;
t383 = t470 * t686 + t596;
t629 = qJD(2) * t685;
t612 = qJD(1) * t629;
t521 = t574 * t612;
t522 = t577 * t612;
t626 = qJD(3) * t684;
t642 = qJD(3) * t573;
t409 = t518 * t642 + t573 * t521 - t684 * t522 - t536 * t626;
t399 = -pkin(4) * t469 + t409;
t621 = -t383 * t572 + t576 * t399;
t583 = -qJD(5) * t386 + t621;
t355 = -pkin(5) * t469 - pkin(10) * t414 + t583;
t464 = t576 * t470;
t415 = qJD(5) * t486 - t464;
t590 = t576 * t383 + t572 * t399 + t419 * t640 - t422 * t641;
t356 = -pkin(10) * t415 + t590;
t622 = t575 * t355 - t571 * t356;
t696 = t406 * t603 + t622;
t695 = (-t424 ^ 2 + t603 ^ 2) * MDP(30) - t469 * MDP(33) - t424 * MDP(29) * t603;
t694 = -0.2e1 * t633;
t693 = MDP(5) * (t574 ^ 2 - t577 ^ 2);
t466 = t527 * t528;
t601 = -qJ(4) * t530 + t559;
t448 = t528 * t686 + t601;
t489 = -t542 * t684 - t573 * t543;
t459 = t530 * pkin(4) + t489;
t455 = t572 * t459;
t655 = t576 * t448 + t455;
t476 = t534 * t684 + t514;
t651 = -pkin(2) * t626 - qJD(4) + t476;
t692 = -t510 * pkin(5) - pkin(10) * t641;
t490 = -t573 * t542 + t684 * t543;
t630 = -pkin(5) * t576 - pkin(4);
t691 = pkin(5) * t640 - t630 * t644;
t475 = t573 * t534 - t517;
t449 = t475 - t681;
t632 = pkin(2) * t642;
t690 = (-t449 + t632) * t576;
t619 = t414 * t571 + t575 * t415;
t366 = -qJD(6) * t603 + t619;
t687 = t644 ^ 2;
t683 = pkin(10) * t644;
t682 = pkin(10) * t528;
t558 = -pkin(2) * t684 - pkin(3);
t554 = -pkin(9) + t558;
t679 = -pkin(10) + t554;
t678 = -pkin(10) - t686;
t385 = t576 * t419 - t422 * t572;
t372 = -pkin(10) * t486 + t385;
t370 = pkin(5) * t698 + t372;
t676 = t370 * t575;
t675 = t373 * t575;
t674 = t414 * t576;
t535 = t574 * t629;
t537 = t577 * t629;
t428 = -t684 * t535 - t573 * t537 - t542 * t626 - t543 * t642;
t673 = t428 * t567;
t429 = -qJD(3) * t490 + t573 * t535 - t684 * t537;
t672 = t429 * t567;
t671 = t474 * t567;
t669 = t510 * t644;
t668 = t644 * t576;
t666 = t528 * t572;
t664 = t572 * t469;
t579 = qJD(2) ^ 2;
t662 = t574 * t579;
t463 = t576 * t469;
t660 = t577 * t579;
t580 = qJD(1) ^ 2;
t659 = t577 * t580;
t658 = t686 * t469;
t458 = pkin(2) * t643 + t471;
t430 = t458 + t506;
t657 = t576 * t430 + t572 * t449;
t656 = t576 * t438 + t572 * t443;
t650 = t691 - t651;
t648 = t680 - t651;
t647 = t635 + t691;
t637 = qJD(6) * t575;
t563 = t574 * t677;
t631 = t575 * t414 - t571 * t415 - t484 * t637;
t520 = t679 * t576;
t539 = t678 * t576;
t555 = pkin(2) * t573 + qJ(4);
t624 = -t448 - t682;
t623 = pkin(1) * t694;
t613 = t518 * t626 + t684 * t521 + t573 * t522 + t536 * t642;
t407 = -t567 * qJD(4) - t613;
t393 = -pkin(4) * t470 - t407;
t620 = -t386 * t510 + t393 * t576;
t618 = t698 ^ 2;
t617 = t698 * t433;
t616 = t698 * t486;
t615 = qJD(6) * t370 + t356;
t610 = -t475 + t632;
t519 = t679 * t572;
t608 = qJD(6) * t519 + (-t430 - t683) * t572 + t554 * t641 - t690 + t692;
t437 = t576 * t443;
t538 = t678 * t572;
t607 = qJD(6) * t538 + t437 + t692 + (-t683 - t700) * t572;
t493 = pkin(10) * t668;
t606 = -t520 * t688 - t572 * t632 + t493 + t657;
t605 = -t539 * t688 + t493 + t656;
t361 = t370 * t571 + t675;
t600 = t385 * t510 + t393 * t572 + (t640 + t668) * t433;
t599 = -t572 * t618 - t463;
t598 = t481 * t572 + t528 * t640;
t597 = -t481 * t576 + t528 * t641;
t480 = -qJD(2) * t628 - t577 * t626 + t604;
t595 = qJ(4) * t480 - qJD(4) * t530 + t563;
t452 = pkin(3) * t510 + t585;
t594 = t452 * t644 + t409;
t593 = -t541 * t644 - t409;
t592 = t541 * t510 - t613;
t392 = t481 * t686 + t595;
t404 = -t480 * pkin(4) + t429;
t589 = t576 * t392 + t572 * t404 - t448 * t641 + t459 * t640;
t365 = -t486 * t638 + t631;
t360 = -t373 * t571 + t676;
t369 = pkin(5) * t415 + t393;
t588 = t360 * t510 + t369 * t527 + t406 * t699;
t587 = -t361 * t510 + t369 * t689 - t406 * t652;
t586 = -t452 * t510 - t407;
t584 = -t576 * t618 + t664;
t403 = -pkin(4) * t481 - t428;
t582 = t649 + (t510 - t627) * t567;
t581 = MDP(11) * t669 + (-t365 * t527 - t366 * t689 + t424 * t652 + t603 * t699) * MDP(30) + (t365 * t689 + t603 * t652) * MDP(29) + ((-t415 - t616) * t576 + (-t414 + t704) * t572) * MDP(23) + (-t510 * t603 + t703) * MDP(31) + (-t424 * t510 + t705) * MDP(32) + (-t572 * t616 + t674) * MDP(22) + (t486 * t510 + t599) * MDP(24) + (-t484 * t510 + t584) * MDP(25) + t582 * MDP(13) + (-t510 ^ 2 + t687) * MDP(12) + (MDP(26) * t698 + t498 * MDP(33)) * t510;
t566 = t572 * pkin(5);
t556 = qJ(4) + t566;
t540 = t555 + t566;
t472 = pkin(3) * t528 + t601;
t461 = -t565 - t474;
t460 = -pkin(4) * t528 - t490;
t457 = -pkin(3) * t567 + t635;
t456 = t576 * t459;
t445 = t469 * t530;
t434 = t528 * t630 - t490;
t405 = pkin(3) * t481 + t595;
t402 = t576 * t404;
t400 = pkin(3) * t470 + t596;
t398 = t576 * t682 + t655;
t388 = pkin(5) * t530 + t572 * t624 + t456;
t379 = pkin(5) * t597 + t403;
t375 = t466 * t688 - t481 * t689;
t374 = t465 * t688 + t527 * t481;
t359 = -pkin(10) * t597 + t589;
t358 = -pkin(5) * t480 + t402 + (-pkin(10) * t481 - t392) * t572 + (t576 * t624 - t455) * qJD(5);
t1 = [(t361 * t480 + t434 * t365 + t369 * t466 + t371 * t530 + t406 * t374 - t379 * t603 + (-(-qJD(6) * t398 + t358) * t498 + t388 * t469 - t355 * t530) * t571 + (-(qJD(6) * t388 + t359) * t498 + t398 * t469 - t615 * t530) * t575) * MDP(35) + (t365 * t466 - t374 * t603) * MDP(29) + (t365 * t530 + t374 * t498 - t466 * t469 + t480 * t603) * MDP(31) + (-t400 * t530 - t405 * t644 + t452 * t480 + t469 * t472 - t673) * MDP(20) + (-t469 * t559 - t480 * t541 + 0.2e1 * t644 * t563 + t673) * MDP(17) + (-t480 * t644 - t445) * MDP(11) + (t407 * t528 + t409 * t530 + t428 * t510 + t429 * t644 - t457 * t480 + t461 * t481 - t469 * t489 + t470 * t490) * MDP(18) + (t469 * t528 - t470 * t530 + t480 * t510 - t481 * t644) * MDP(12) + ((t358 * t575 - t359 * t571) * t498 - (t388 * t575 - t398 * t571) * t469 + t622 * t530 - t360 * t480 + t379 * t424 + t434 * t366 - t369 * t465 + t406 * t375 + ((-t388 * t571 - t398 * t575) * t498 - t361 * t530) * qJD(6)) * MDP(34) + (t365 * t465 - t366 * t466 - t374 * t424 + t375 * t603) * MDP(30) + (-t366 * t530 - t375 * t498 + t424 * t480 - t465 * t469) * MDP(32) + (-MDP(13) * t480 - MDP(14) * t481) * t567 + ((-t392 * t572 + t402) * t698 - (-t448 * t572 + t456) * t469 + t621 * t530 - t385 * t480 + t403 * t484 + t460 * t415 + (-t393 * t528 - t433 * t481) * t576 + (-t386 * t530 + t433 * t666 - t655 * t698) * qJD(5)) * MDP(27) + (-t415 * t530 - t463 * t528 + t480 * t484 - t597 * t698) * MDP(25) + (t386 * t480 + t393 * t666 + t403 * t486 + t460 * t414 + t433 * t598 + t469 * t655 - t530 * t590 - t589 * t698) * MDP(28) + (t414 * t530 - t480 * t486 - t528 * t664 + t598 * t698) * MDP(24) + (-t480 * t698 - t445) * MDP(26) + (-t672 + t470 * t559 + t481 * t541 + (qJD(1) * t528 + t510) * t563) * MDP(16) + (-t400 * t528 - t405 * t510 - t452 * t481 - t470 * t472 + t672) * MDP(19) + ((-t484 * t572 + t486 * t576) * t481 + (t674 - t415 * t572 + (-t484 * t576 - t486 * t572) * qJD(5)) * t528) * MDP(23) + (-pkin(7) * t660 + t574 * t623) * MDP(9) + (t414 * t666 + t486 * t598) * MDP(22) + 0.2e1 * t577 * MDP(4) * t625 - MDP(7) * t662 + (pkin(7) * t662 + t577 * t623) * MDP(10) + MDP(6) * t660 + t693 * t694 + (-t480 * t498 - t445) * MDP(33) + (t400 * t472 + t405 * t452 + t407 * t490 + t409 * t489 + t428 * t461 + t429 * t457) * MDP(21); ((t519 * t575 + t520 * t571) * t469 + t540 * t365 + (t571 * t608 + t575 * t606) * t498 - t650 * t603 + t587) * MDP(35) + (-t407 * t555 + t409 * t558 - t452 * t458 + t457 * t610 + t461 * t651) * MDP(21) + (-t554 * t463 + t555 * t415 + t648 * t484 + ((-qJD(5) * t554 + t430) * t572 + t690) * t698 + t600) * MDP(27) + t581 - t574 * MDP(4) * t659 + t580 * t693 + (t555 * t414 + (-t554 * t640 + t657) * t698 + t648 * t486 + (t554 * t469 - t632 * t698 - t617) * t572 + t620) * MDP(28) + (t458 * t510 + t567 * t610 + t594) * MDP(19) + (t476 * t567 + (-t567 * t626 - t643 * t644) * pkin(2) + t592) * MDP(17) + (t475 * t567 + (-t510 * t643 - t567 * t642) * pkin(2) + t593) * MDP(16) + (-t469 * t558 - t470 * t555 + (-t461 + t610) * t644 + (t457 + t651) * t510) * MDP(18) + (t458 * t644 - t567 * t651 + t586) * MDP(20) + (-(-t519 * t571 + t520 * t575) * t469 + t540 * t366 + (t571 * t606 - t575 * t608) * t498 + t650 * t424 + t588) * MDP(34) + (MDP(9) * t574 * t580 + MDP(10) * t659) * pkin(1); (-pkin(3) * t409 - qJ(4) * t407 - t452 * t471 - t457 * t474 - t461 * t635) * MDP(21) + (t576 * t658 + qJ(4) * t415 + (t572 * t700 - t437) * t698 + t636 * t484 + t600) * MDP(27) + ((t538 * t575 + t539 * t571) * t469 + t556 * t365 + (t571 * t607 + t575 * t605) * t498 - t647 * t603 + t587) * MDP(35) + (-(-t538 * t571 + t539 * t575) * t469 + t556 * t366 + (t571 * t605 - t575 * t607) * t498 + t647 * t424 + t588) * MDP(34) + t581 + (qJ(4) * t414 + (t576 * t639 + t656) * t698 + t636 * t486 + (-t617 - t658) * t572 + t620) * MDP(28) + (t593 + t671) * MDP(16) + (pkin(3) * t469 - qJ(4) * t470 + (-t461 - t474) * t644 + (t457 - t635) * t510) * MDP(18) + (t471 * t510 + t594 - t671) * MDP(19) + (t471 * t644 + t567 * t635 + t586) * MDP(20) + (-t473 * t567 + t592) * MDP(17); t582 * MDP(18) - MDP(19) * t669 + (-t567 ^ 2 - t687) * MDP(20) + (t461 * t567 + t594) * MDP(21) + (-t484 * t567 + t599) * MDP(27) + (-t486 * t567 + t584) * MDP(28) + (-t567 * t424 + t703) * MDP(34) + (t567 * t603 + t705) * MDP(35); t486 * t484 * MDP(22) + (-t484 ^ 2 + t486 ^ 2) * MDP(23) + (t414 + t704) * MDP(24) + (t464 + (-qJD(5) + t698) * t486) * MDP(25) - t469 * MDP(26) + (t386 * t698 - t433 * t486 + t583) * MDP(27) + (t385 * t698 + t433 * t484 - t590) * MDP(28) + (t365 + t702) * MDP(31) + (-t366 - t701) * MDP(32) + (-(-t372 * t571 - t675) * t498 - t361 * qJD(6) + (-t424 * t486 - t469 * t575 - t498 * t638) * pkin(5) + t696) * MDP(34) + ((-t373 * t498 - t355) * t571 + (t372 * t498 - t615) * t575 + (t469 * t571 + t486 * t603 - t498 * t637) * pkin(5) + t697) * MDP(35) + t695; (t631 + t702) * MDP(31) + (-t619 - t701) * MDP(32) + (t361 * t498 + t696) * MDP(34) + (-t571 * t355 - t575 * t356 + t360 * t498 + t697) * MDP(35) + (-MDP(31) * t670 + MDP(32) * t603 - MDP(34) * t361 - MDP(35) * t676) * qJD(6) + t695;];
tauc  = t1;
