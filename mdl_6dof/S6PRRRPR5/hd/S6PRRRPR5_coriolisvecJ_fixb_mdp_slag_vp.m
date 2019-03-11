% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:37
% EndTime: 2019-03-08 23:28:52
% DurationCPUTime: 9.91s
% Computational Cost: add. (7498->499), mult. (20456->723), div. (0->0), fcn. (16940->14), ass. (0->233)
t565 = sin(pkin(7));
t572 = sin(qJ(3));
t683 = t565 * t572;
t556 = pkin(9) * t683;
t568 = cos(pkin(7));
t576 = cos(qJ(3));
t577 = cos(qJ(2));
t675 = t576 * t577;
t573 = sin(qJ(2));
t679 = t572 * t573;
t592 = -t568 * t679 + t675;
t566 = sin(pkin(6));
t663 = qJD(1) * t566;
t681 = t568 * t576;
t718 = t592 * t663 - (pkin(2) * t681 - t556) * qJD(3);
t608 = pkin(3) * t572 - pkin(10) * t576;
t591 = t608 * qJD(3);
t642 = t573 * t663;
t717 = (t591 - t642) * t565;
t571 = sin(qJ(4));
t575 = cos(qJ(4));
t682 = t565 * t576;
t512 = pkin(9) * t682 + (pkin(2) * t572 + pkin(10)) * t568;
t609 = -pkin(3) * t576 - pkin(10) * t572;
t513 = (-pkin(2) + t609) * t565;
t703 = t575 * t512 + t571 * t513;
t716 = qJD(4) * t703 - t571 * t718 - t717 * t575;
t653 = qJD(4) * t575;
t654 = qJD(4) * t571;
t715 = -t512 * t654 + t513 * t653 + t717 * t571 - t575 * t718;
t658 = qJD(2) * t576;
t640 = t565 * t658;
t714 = qJD(4) - t640;
t528 = -t575 * t568 + t571 * t683;
t655 = qJD(3) * t576;
t638 = t565 * t655;
t487 = -qJD(4) * t528 + t575 * t638;
t529 = t568 * t571 + t575 * t683;
t657 = qJD(3) * t572;
t639 = t565 * t657;
t713 = pkin(4) * t639 - qJ(5) * t487 - qJD(5) * t529 - t716;
t636 = t571 * t655;
t488 = qJD(4) * t529 + t565 * t636;
t712 = -qJ(5) * t488 - qJD(5) * t528 + t715;
t661 = qJD(2) * t565;
t532 = pkin(9) * t661 + t642;
t522 = t572 * t532;
t569 = cos(pkin(6));
t662 = qJD(1) * t569;
t643 = t565 * t662;
t696 = qJD(2) * pkin(2);
t543 = t577 * t663 + t696;
t686 = t543 * t568;
t462 = t576 * (t643 + t686) - t522;
t517 = t608 * t661;
t628 = -t462 * t571 + t575 * t517;
t697 = -qJ(5) - pkin(10);
t631 = qJD(4) * t697;
t676 = t575 * t576;
t711 = -(pkin(4) * t572 - qJ(5) * t676) * t661 - t628 - qJD(5) * t571 + t575 * t631;
t614 = t571 * t640;
t671 = t575 * t462 + t571 * t517;
t710 = -qJ(5) * t614 - qJD(5) * t575 - t571 * t631 + t671;
t574 = cos(qJ(6));
t570 = sin(qJ(6));
t659 = qJD(2) * t568;
t555 = qJD(3) + t659;
t641 = t572 * t661;
t615 = t571 * t641;
t502 = -t575 * t555 + t615;
t504 = t555 * t571 + t575 * t641;
t564 = sin(pkin(13));
t567 = cos(pkin(13));
t599 = -t502 * t564 + t567 * t504;
t690 = t599 * t570;
t440 = -t574 * t714 + t690;
t625 = -t567 * t502 - t504 * t564;
t701 = qJD(6) - t625;
t709 = t440 * t701;
t442 = t570 * t714 + t574 * t599;
t708 = t442 * t701;
t533 = t564 * t571 - t567 * t575;
t707 = t714 * t533;
t677 = t573 * t576;
t678 = t572 * t577;
t594 = t568 * t677 + t678;
t637 = t568 * t657;
t666 = pkin(2) * t637 + pkin(9) * t638 - t594 * t663;
t534 = t564 * t575 + t567 * t571;
t665 = t714 * t534;
t622 = t701 * t574;
t648 = qJD(2) * qJD(3);
t632 = t565 * t648;
t611 = t576 * t632;
t478 = -qJD(4) * t615 + t555 * t653 + t575 * t611;
t479 = (t572 * t653 + t636) * t661 + t555 * t654;
t433 = t478 * t564 + t567 * t479;
t680 = t570 * t433;
t706 = -t701 * t622 - t680;
t673 = -t712 * t564 + t713 * t567;
t672 = t713 * t564 + t712 * t567;
t670 = t710 * t564 + t711 * t567;
t668 = t711 * t564 - t710 * t567;
t704 = pkin(4) * t488 + t666;
t545 = t572 * t643;
t463 = t576 * t532 + t572 * t686 + t545;
t702 = -t463 + (-t614 + t654) * pkin(4);
t451 = pkin(10) * t555 + t463;
t554 = t568 * t662;
t474 = t554 + (qJD(2) * t609 - t543) * t565;
t423 = t451 * t575 + t474 * t571;
t582 = t592 * qJD(2);
t613 = t569 * t638;
t436 = (t543 * t681 - t522) * qJD(3) + (t566 * t582 + t613) * qJD(1);
t483 = (t591 + t642) * t661;
t579 = -qJD(4) * t423 - t571 * t436 + t575 * t483;
t612 = t572 * t632;
t375 = pkin(4) * t612 - qJ(5) * t478 - qJD(5) * t504 + t579;
t589 = t575 * t436 - t451 * t654 + t474 * t653 + t571 * t483;
t379 = -qJ(5) * t479 - qJD(5) * t502 + t589;
t366 = t375 * t567 - t379 * t564;
t364 = -pkin(5) * t612 - t366;
t559 = pkin(4) * t564 + pkin(11);
t700 = t701 * (pkin(4) * t504 + pkin(5) * t599 - pkin(11) * t625 + qJD(6) * t559) + t364;
t561 = t565 ^ 2;
t699 = (-MDP(5) * t572 * t576 + (t572 ^ 2 - t576 ^ 2) * MDP(6)) * t561;
t434 = t478 * t567 - t479 * t564;
t651 = qJD(6) * t574;
t645 = t574 * t434 + t570 * t612 + t651 * t714;
t652 = qJD(6) * t570;
t394 = -t599 * t652 + t645;
t694 = t394 * t570;
t415 = -qJ(5) * t502 + t423;
t693 = t415 * t564;
t692 = t440 * t599;
t691 = t442 * t599;
t689 = t502 * t714;
t688 = t504 * t714;
t687 = t534 * t574;
t685 = t714 * t571;
t684 = t714 * t575;
t411 = t567 * t415;
t430 = t574 * t433;
t674 = -pkin(5) * t639 - t673;
t367 = t564 * t375 + t567 * t379;
t422 = -t451 * t571 + t575 * t474;
t414 = -qJ(5) * t504 + t422;
t405 = pkin(4) * t714 + t414;
t383 = t564 * t405 + t411;
t624 = -t512 * t571 + t575 * t513;
t439 = -pkin(4) * t682 - qJ(5) * t529 + t624;
t446 = -qJ(5) * t528 + t703;
t402 = t564 * t439 + t567 * t446;
t669 = pkin(5) * t641 - t670;
t660 = qJD(2) * t566;
t656 = qJD(3) * t575;
t650 = qJD(3) - t555;
t646 = t570 * t682;
t644 = -pkin(4) * t575 - pkin(3);
t634 = t697 * t571;
t633 = qJD(1) * t660;
t365 = pkin(11) * t612 + t367;
t597 = t568 * t573 * t633;
t437 = qJD(3) * t545 + t532 * t655 + t543 * t637 + t576 * t597 + t633 * t678;
t420 = pkin(4) * t479 + t437;
t381 = pkin(5) * t433 - pkin(11) * t434 + t420;
t630 = -t365 * t570 + t574 * t381;
t629 = t434 * t570 - t574 * t612;
t627 = t570 * t707 - t574 * t641;
t626 = t570 * t641 + t574 * t707;
t616 = t565 * t573 * t660;
t475 = t567 * t528 + t529 * t564;
t476 = -t528 * t564 + t529 * t567;
t511 = t556 + (-pkin(2) * t576 - pkin(3)) * t568;
t581 = pkin(4) * t528 + t511;
t416 = pkin(5) * t475 - pkin(11) * t476 + t581;
t607 = -pkin(11) * t639 - qJD(6) * t416 - t672;
t482 = pkin(5) * t533 - pkin(11) * t534 + t644;
t606 = pkin(11) * t641 - qJD(6) * t482 - t668;
t398 = -pkin(11) * t682 + t402;
t444 = t487 * t564 + t567 * t488;
t445 = t487 * t567 - t488 * t564;
t605 = -pkin(5) * t444 + pkin(11) * t445 + qJD(6) * t398 - t704;
t551 = t697 * t575;
t492 = -t567 * t551 + t564 * t634;
t604 = -pkin(5) * t665 - pkin(11) * t707 + qJD(6) * t492 - t702;
t603 = t365 * t574 + t381 * t570;
t378 = pkin(11) * t714 + t383;
t450 = -pkin(3) * t555 - t462;
t435 = pkin(4) * t502 + qJD(5) + t450;
t390 = -pkin(5) * t625 - pkin(11) * t599 + t435;
t369 = t378 * t574 + t390 * t570;
t602 = t378 * t570 - t390 * t574;
t382 = t405 * t567 - t693;
t593 = t568 * t678 + t677;
t486 = t566 * t593 + t569 * t683;
t525 = -t565 * t566 * t577 + t568 * t569;
t452 = -t486 * t571 + t525 * t575;
t453 = t486 * t575 + t525 * t571;
t418 = t452 * t564 + t453 * t567;
t595 = t568 * t675 - t679;
t485 = -t566 * t595 - t569 * t682;
t601 = t418 * t574 + t485 * t570;
t600 = -t418 * t570 + t485 * t574;
t401 = t439 * t567 - t446 * t564;
t596 = t430 + (t570 * t625 - t652) * t701;
t454 = t476 * t570 + t574 * t682;
t587 = t534 * t651 - t627;
t586 = -t534 * t652 - t626;
t498 = -t543 * t565 + t554;
t584 = qJD(3) * (t498 * t565 - t561 * t696);
t377 = -pkin(5) * t714 - t382;
t387 = t414 * t567 - t693;
t580 = -t559 * t433 + (t377 + t387) * t701;
t578 = qJD(2) ^ 2;
t560 = -pkin(4) * t567 - pkin(5);
t491 = -t551 * t564 - t567 * t634;
t455 = t476 * t574 - t646;
t448 = t613 + (qJD(3) * t595 + t582) * t566;
t447 = t569 * t639 + (qJD(2) * t594 + qJD(3) * t593) * t566;
t417 = -t567 * t452 + t453 * t564;
t410 = qJD(4) * t452 + t448 * t575 + t571 * t616;
t409 = -qJD(4) * t453 - t448 * t571 + t575 * t616;
t408 = -qJD(6) * t646 + t445 * t570 + t476 * t651 - t574 * t639;
t407 = -qJD(6) * t454 + t445 * t574 + t570 * t639;
t397 = pkin(5) * t682 - t401;
t395 = qJD(6) * t442 + t629;
t386 = t414 * t564 + t411;
t385 = t409 * t564 + t410 * t567;
t384 = -t567 * t409 + t410 * t564;
t363 = -qJD(6) * t369 + t630;
t362 = -qJD(6) * t602 + t603;
t1 = [(-t447 * t555 + t525 * t612) * MDP(10) + (-t448 * t555 + t525 * t611) * MDP(11) + (t409 * t714 + t447 * t502 + t452 * t612 + t479 * t485) * MDP(17) + (-t410 * t714 + t447 * t504 - t453 * t612 + t478 * t485) * MDP(18) + (t384 * t599 + t385 * t625 + t417 * t434 - t418 * t433) * MDP(19) + (-t366 * t417 + t367 * t418 - t382 * t384 + t383 * t385 + t420 * t485 + t435 * t447) * MDP(20) + ((-qJD(6) * t601 - t385 * t570 + t447 * t574) * t701 + t600 * t433 + t384 * t440 + t417 * t395) * MDP(26) + (-(qJD(6) * t600 + t385 * t574 + t447 * t570) * t701 - t601 * t433 + t384 * t442 + t417 * t394) * MDP(27) + (-MDP(4) * t577 + (-MDP(3) + (-MDP(10) * t576 + MDP(11) * t572) * t561) * t573) * t578 * t566; (-t437 * t568 - t555 * t666 + t572 * t584) * MDP(10) + (-t436 * t568 + t555 * t718 + t576 * t584) * MDP(11) + (t478 * t529 + t487 * t504) * MDP(12) + (-t478 * t528 - t479 * t529 - t487 * t502 - t488 * t504) * MDP(13) + (t487 * t714 + (-t478 * t576 + (qJD(2) * t529 + t504) * t657) * t565) * MDP(14) + (-t488 * t714 + (t479 * t576 + (-qJD(2) * t528 - t502) * t657) * t565) * MDP(15) + (-t561 * t658 + t565 * t714) * MDP(16) * t657 + (t437 * t528 + t450 * t488 + t511 * t479 - t716 * t714 + t666 * t502 + (-t579 * t576 + (qJD(2) * t624 + t422) * t657) * t565) * MDP(17) + (t437 * t529 + t450 * t487 + t511 * t478 - t715 * t714 + t666 * t504 + (t589 * t576 + (-qJD(2) * t703 - t423) * t657) * t565) * MDP(18) + (-t366 * t476 - t367 * t475 - t382 * t445 - t383 * t444 - t401 * t434 - t402 * t433 - t599 * t673 + t625 * t672) * MDP(19) + (t366 * t401 + t367 * t402 + t673 * t382 + t672 * t383 + t420 * t581 + t435 * t704) * MDP(20) + (t394 * t455 + t407 * t442) * MDP(21) + (-t394 * t454 - t395 * t455 - t407 * t440 - t408 * t442) * MDP(22) + (t394 * t475 + t407 * t701 + t433 * t455 + t442 * t444) * MDP(23) + (-t395 * t475 - t408 * t701 - t433 * t454 - t440 * t444) * MDP(24) + (t433 * t475 + t444 * t701) * MDP(25) + ((-t398 * t570 + t416 * t574) * t433 + t363 * t475 - t602 * t444 + t397 * t395 + t364 * t454 + t377 * t408 + (t570 * t607 - t574 * t605) * t701 + t674 * t440) * MDP(26) + (-(t398 * t574 + t416 * t570) * t433 - t362 * t475 - t369 * t444 + t397 * t394 + t364 * t455 + t377 * t407 + (t570 * t605 + t574 * t607) * t701 + t674 * t442) * MDP(27) - 0.2e1 * t699 * t648 + (MDP(7) * t638 - MDP(8) * t639) * (t555 + t659); t650 * MDP(7) * t640 + (t463 * t555 - t437) * MDP(10) + (t462 * t555 + (qJD(3) * t532 + t597) * t572 + (-t498 * t661 - qJD(3) * t686 + (-qJD(3) * t565 * t569 - t577 * t660) * qJD(1)) * t576) * MDP(11) + (t478 * t571 + t504 * t684) * MDP(12) + ((t478 - t689) * t575 + (-t479 - t688) * t571) * MDP(13) + (t714 * t653 + (-t714 * t676 + (qJD(3) * t571 - t504) * t572) * t661) * MDP(14) + (-t714 * t654 + (t576 * t685 + (t502 + t656) * t572) * t661) * MDP(15) + (-pkin(3) * t479 - t437 * t575 - t628 * t714 - t463 * t502 + (-pkin(10) * t684 + t450 * t571) * qJD(4) + (-t422 * t572 + (-pkin(10) * t657 - t450 * t576) * t571) * t661) * MDP(17) + (-pkin(3) * t478 + t437 * t571 - t463 * t504 + t671 * t714 + (pkin(10) * t685 + t450 * t575) * qJD(4) + (-t450 * t676 + (-pkin(10) * t656 + t423) * t572) * t661) * MDP(18) + (-t366 * t534 - t367 * t533 + t382 * t707 - t665 * t383 - t492 * t433 + t434 * t491 - t670 * t599 + t668 * t625) * MDP(19) + (-t366 * t491 + t367 * t492 + t670 * t382 + t668 * t383 + t420 * t644 + t435 * t702) * MDP(20) + (t394 * t687 + t442 * t586) * MDP(21) + (t627 * t442 + t626 * t440 + (-t694 - t395 * t574 + (t440 * t570 - t442 * t574) * qJD(6)) * t534) * MDP(22) + (t394 * t533 + t430 * t534 + t442 * t665 + t586 * t701) * MDP(23) + (-t395 * t533 - t440 * t665 - t534 * t680 - t587 * t701) * MDP(24) + (t433 * t533 + t665 * t701) * MDP(25) + ((t482 * t574 - t492 * t570) * t433 + t363 * t533 + t491 * t395 + t364 * t570 * t534 + (t570 * t606 - t574 * t604) * t701 + t669 * t440 - t665 * t602 + t587 * t377) * MDP(26) + (-(t482 * t570 + t492 * t574) * t433 - t362 * t533 + t491 * t394 + t364 * t687 + (t570 * t604 + t574 * t606) * t701 + t669 * t442 - t665 * t369 + t586 * t377) * MDP(27) + (-t498 * MDP(10) - MDP(16) * t714 - MDP(8) * t650) * t641 + t699 * t578; t504 * t502 * MDP(12) + (-t502 ^ 2 + t504 ^ 2) * MDP(13) + (t478 + t689) * MDP(14) + (-t479 + t688) * MDP(15) + MDP(16) * t612 + (t423 * t714 - t450 * t504 + t579) * MDP(17) + (t422 * t714 + t450 * t502 - t589) * MDP(18) + ((-t433 * t564 - t434 * t567) * pkin(4) + (t382 - t387) * t625 + (t383 - t386) * t599) * MDP(19) + (t382 * t386 - t383 * t387 + (t366 * t567 + t367 * t564 - t435 * t504) * pkin(4)) * MDP(20) + (t442 * t622 + t694) * MDP(21) + ((t394 - t709) * t574 + (-t395 - t708) * t570) * MDP(22) + (-t691 - t706) * MDP(23) + (t596 + t692) * MDP(24) - t701 * t599 * MDP(25) + (-t386 * t440 + t560 * t395 + t580 * t570 - t574 * t700 + t599 * t602) * MDP(26) + (t369 * t599 - t386 * t442 + t560 * t394 + t570 * t700 + t580 * t574) * MDP(27); (-t599 ^ 2 - t625 ^ 2) * MDP(19) + (t382 * t599 - t383 * t625 + t420) * MDP(20) + (t596 - t692) * MDP(26) + (-t691 + t706) * MDP(27); t442 * t440 * MDP(21) + (-t440 ^ 2 + t442 ^ 2) * MDP(22) + (t645 + t709) * MDP(23) + (-t629 + t708) * MDP(24) + t433 * MDP(25) + (t369 * t701 - t377 * t442 + t630) * MDP(26) + (t377 * t440 - t602 * t701 - t603) * MDP(27) + (-MDP(23) * t690 - MDP(24) * t442 - MDP(26) * t369 + MDP(27) * t602) * qJD(6);];
tauc  = t1;
