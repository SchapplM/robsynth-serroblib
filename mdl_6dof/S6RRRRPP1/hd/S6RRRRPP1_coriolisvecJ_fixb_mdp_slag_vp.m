% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:59
% EndTime: 2019-03-09 20:47:12
% DurationCPUTime: 7.24s
% Computational Cost: add. (11516->515), mult. (27691->652), div. (0->0), fcn. (19738->8), ass. (0->235)
t594 = qJD(2) + qJD(3);
t600 = sin(qJ(3));
t603 = cos(qJ(2));
t702 = cos(qJ(3));
t645 = qJD(1) * t702;
t601 = sin(qJ(2));
t663 = qJD(1) * t601;
t718 = -t600 * t663 + t603 * t645;
t505 = t718 * t594;
t682 = t600 * t603;
t565 = t601 * t702 + t682;
t597 = sin(pkin(10));
t598 = cos(pkin(10));
t599 = sin(qJ(4));
t602 = cos(qJ(4));
t706 = -t597 * t599 + t598 * t602;
t494 = t706 * t565;
t562 = t597 * t602 + t598 * t599;
t546 = t562 * qJD(4);
t716 = -t562 * t718 + t546;
t660 = qJD(4) * t602;
t661 = qJD(4) * t599;
t715 = -t597 * t661 + t598 * t660 - t706 * t718;
t550 = -qJD(1) * t682 - t601 * t645;
t524 = -t550 * t599 - t602 * t594;
t526 = -t550 * t602 + t594 * t599;
t466 = t598 * t524 + t526 * t597;
t542 = qJD(4) - t718;
t717 = t466 * t542;
t687 = t718 * t599;
t709 = (t661 - t687) * pkin(4);
t619 = -t600 * t601 + t603 * t702;
t517 = t594 * t619;
t646 = t565 * t660;
t714 = t517 * t599 + t646;
t624 = -t524 * t597 + t598 * t526;
t713 = t624 ^ 2;
t657 = qJD(1) * qJD(2);
t712 = -0.2e1 * t657;
t711 = MDP(5) * (t601 ^ 2 - t603 ^ 2);
t703 = -pkin(8) - pkin(7);
t574 = t703 * t603;
t569 = qJD(1) * t574;
t553 = t600 * t569;
t573 = t703 * t601;
t567 = qJD(1) * t573;
t698 = qJD(2) * pkin(2);
t555 = t567 + t698;
t511 = t555 * t702 + t553;
t489 = -t594 * pkin(3) - t511;
t458 = t524 * pkin(4) + qJD(5) + t489;
t408 = t466 * pkin(5) - qJ(6) * t624 + t458;
t710 = t408 * t624;
t507 = -pkin(3) * t550 - pkin(9) * t718;
t487 = pkin(2) * t663 + t507;
t483 = t602 * t487;
t515 = t567 * t702 + t553;
t593 = t602 * qJ(5);
t631 = -t550 * pkin(4) - t593 * t718;
t432 = -t515 * t599 + t483 + t631;
t654 = qJ(5) * t687;
t666 = t599 * t487 + t602 * t515;
t440 = -t654 + t666;
t592 = t602 * qJD(5);
t644 = qJD(3) * t702;
t634 = pkin(2) * t644;
t627 = t602 * t634;
t587 = pkin(2) * t600 + pkin(9);
t678 = -qJ(5) - t587;
t636 = qJD(4) * t678;
t516 = t599 * t636 + t592 + t627;
t609 = (-t634 - qJD(5)) * t599 + t602 * t636;
t670 = (-t432 + t609) * t598 + (t440 - t516) * t597;
t637 = t602 * t507 - t511 * t599;
t435 = t631 + t637;
t667 = t599 * t507 + t602 * t511;
t442 = -t654 + t667;
t699 = -qJ(5) - pkin(9);
t641 = qJD(4) * t699;
t540 = t599 * t641 + t592;
t613 = -t599 * qJD(5) + t602 * t641;
t669 = (-t435 + t613) * t598 + (t442 - t540) * t597;
t554 = t702 * t569;
t514 = t600 * t567 - t554;
t662 = qJD(3) * t600;
t632 = pkin(2) * t662 - t514;
t708 = -t716 * pkin(5) + t715 * qJ(6) + qJD(6) * t562 - t709;
t707 = t702 * t573 + t600 * t574;
t705 = qJD(1) * t565;
t704 = t542 ^ 2;
t701 = pkin(5) * t550;
t700 = t602 * pkin(4);
t589 = -pkin(2) * t603 - pkin(1);
t572 = t589 * qJD(1);
t484 = -pkin(3) * t718 + pkin(9) * t550 + t572;
t512 = t600 * t555 - t554;
t490 = pkin(9) * t594 + t512;
t444 = t484 * t599 + t490 * t602;
t430 = -qJ(5) * t524 + t444;
t426 = t598 * t430;
t443 = t602 * t484 - t490 * t599;
t429 = -qJ(5) * t526 + t443;
t394 = t429 * t597 + t426;
t697 = t394 * t624;
t696 = t430 * t597;
t668 = t602 * t505 + t594 * t660;
t461 = t550 * t661 + t668;
t695 = t461 * t599;
t694 = t489 * t718;
t693 = t505 * t599;
t518 = t594 * t565;
t506 = t518 * qJD(1);
t692 = t506 * t599;
t691 = t506 * t602;
t689 = t524 * t542;
t688 = t526 * t542;
t686 = t565 * t599;
t685 = t565 * t602;
t604 = qJD(2) ^ 2;
t681 = t601 * t604;
t528 = t600 * t573 - t574 * t702;
t521 = t602 * t528;
t680 = t603 * t604;
t605 = qJD(1) ^ 2;
t679 = t603 * t605;
t642 = t601 * t657;
t449 = pkin(2) * t642 + pkin(3) * t506 - pkin(9) * t505;
t447 = t602 * t449;
t649 = qJD(2) * t703;
t633 = qJD(1) * t649;
t556 = t601 * t633;
t557 = t603 * t633;
t452 = t555 * t644 + t556 * t702 + t600 * t557 + t569 * t662;
t638 = -t452 * t599 + t447;
t376 = pkin(4) * t506 - qJ(5) * t461 - qJD(4) * t444 - qJD(5) * t526 + t638;
t462 = qJD(4) * t526 + t693;
t652 = t599 * t449 + t602 * t452 + t484 * t660;
t616 = -t490 * t661 + t652;
t380 = -qJ(5) * t462 - qJD(5) * t524 + t616;
t677 = -t598 * t376 + t597 * t380;
t366 = t597 * t376 + t598 * t380;
t655 = t601 * t698;
t457 = pkin(3) * t518 - pkin(9) * t517 + t655;
t455 = t602 * t457;
t568 = t601 * t649;
t570 = t603 * t649;
t470 = qJD(3) * t707 + t702 * t568 + t600 * t570;
t510 = -pkin(3) * t619 - pkin(9) * t565 + t589;
t623 = -qJ(5) * t517 - qJD(5) * t565;
t382 = pkin(4) * t518 - t470 * t599 + t455 + t623 * t602 + (-t521 + (qJ(5) * t565 - t510) * t599) * qJD(4);
t651 = t599 * t457 + t602 * t470 + t510 * t660;
t391 = -qJ(5) * t646 + (-qJD(4) * t528 + t623) * t599 + t651;
t370 = t597 * t382 + t598 * t391;
t401 = t597 * t432 + t598 * t440;
t539 = t550 * qJ(6);
t396 = -t539 + t401;
t460 = t598 * t516 + t597 * t609;
t676 = -t396 + t460;
t675 = t701 + t670;
t403 = t597 * t435 + t598 * t442;
t398 = -t539 + t403;
t486 = t598 * t540 + t597 * t613;
t674 = -t398 + t486;
t673 = t701 + t669;
t672 = t512 + t708;
t671 = t708 - t632;
t423 = pkin(4) * t542 + t429;
t393 = t597 * t423 + t426;
t499 = t602 * t510;
t441 = -pkin(4) * t619 - t528 * t599 - t565 * t593 + t499;
t665 = t599 * t510 + t521;
t448 = -qJ(5) * t686 + t665;
t410 = t597 * t441 + t598 * t448;
t659 = qJD(6) * t542;
t395 = t429 * t598 - t696;
t658 = qJD(6) - t395;
t656 = t702 * pkin(2);
t653 = t506 * qJ(6) + t366;
t650 = -pkin(3) - t700;
t481 = t489 * t660;
t643 = t699 * t599;
t640 = t678 * t599;
t639 = pkin(1) * t712;
t416 = t461 * t597 + t598 * t462;
t635 = t542 * t602;
t453 = t555 * t662 + t600 * t556 - t702 * t557 - t569 * t644;
t588 = -t656 - pkin(3);
t362 = -pkin(5) * t506 + t677;
t417 = t461 * t598 - t462 * t597;
t558 = t587 * t602 + t593;
t501 = t558 * t597 - t598 * t640;
t502 = t598 * t558 + t597 * t640;
t630 = -t502 * t416 + t417 * t501 - t460 * t466;
t571 = pkin(9) * t602 + t593;
t522 = t571 * t597 - t598 * t643;
t523 = t598 * t571 + t597 * t643;
t629 = -t523 * t416 + t417 * t522 - t486 * t466;
t628 = -t444 * t550 + t453 * t599 + t481;
t369 = t382 * t598 - t391 * t597;
t392 = t423 * t598 - t696;
t409 = t441 * t598 - t448 * t597;
t625 = -t506 * t587 - t694;
t622 = pkin(4) * t686 - t707;
t620 = t443 * t550 - t453 * t602 + t489 * t661;
t618 = t517 * t602 - t565 * t661;
t418 = pkin(4) * t462 + t453;
t617 = t550 * t572 - t453;
t471 = t600 * t568 - t570 * t702 + t573 * t662 - t574 * t644;
t373 = pkin(5) * t416 - qJ(6) * t417 - qJD(6) * t624 + t418;
t386 = -pkin(5) * t542 + qJD(6) - t392;
t615 = -t373 * t706 - t386 * t550 + t716 * t408;
t387 = qJ(6) * t542 + t393;
t614 = -t373 * t562 + t387 * t550 - t715 * t408;
t504 = -pkin(5) * t706 - t562 * qJ(6) + t650;
t612 = t714 * pkin(4) + t471;
t360 = t653 + t659;
t611 = t360 * t706 + t362 * t562 + t715 * t386 - t716 * t387;
t610 = t366 * t706 - t715 * t392 - t716 * t393 + t677 * t562;
t606 = t594 * t705;
t608 = ((t461 - t689) * t602 + (-t462 - t688) * t599) * MDP(19) + (t526 * t635 + t695) * MDP(18) + (-t524 * t550 - t599 * t704 + t691) * MDP(21) + (t526 * t550 + t542 * t635 + t692) * MDP(20) + (-t550 * t594 - t606) * MDP(14) + (t550 ^ 2 - t718 ^ 2) * MDP(12) + (MDP(11) * t718 + MDP(22) * t542) * t550;
t607 = -t572 * t718 - t452;
t585 = -pkin(4) * t598 - pkin(5);
t583 = pkin(4) * t597 + qJ(6);
t493 = t562 * t565;
t488 = -t656 + t504;
t434 = -t517 * t706 + t546 * t565;
t433 = -qJD(4) * t494 - t517 * t562;
t425 = t493 * pkin(5) - t494 * qJ(6) + t622;
t413 = pkin(4) * t526 + pkin(5) * t624 + qJ(6) * t466;
t407 = pkin(5) * t619 - t409;
t406 = -qJ(6) * t619 + t410;
t377 = -t433 * pkin(5) + t434 * qJ(6) - t494 * qJD(6) + t612;
t368 = -pkin(5) * t518 - t369;
t367 = qJ(6) * t518 - qJD(6) * t619 + t370;
t1 = [MDP(6) * t680 + (t505 * t589 + t517 * t572 + (-t550 + t705) * t655) * MDP(17) + (t506 * t589 + t518 * t572 + (-qJD(1) * t619 - t718) * t655) * MDP(16) + ((-t528 * t660 + t455) * t542 + t499 * t506 - (-t490 * t660 + t447) * t619 + t443 * t518 + t471 * t524 - t707 * t462 + t565 * t481 + ((-qJD(4) * t510 - t470) * t542 - t528 * t506 - (-qJD(4) * t484 - t452) * t619 + t453 * t565 + t489 * t517) * t599) * MDP(23) + (-t360 * t619 + t367 * t542 - t373 * t494 - t377 * t624 + t387 * t518 + t406 * t506 + t408 * t434 - t417 * t425) * MDP(29) + (-t506 * t619 + t518 * t542) * MDP(22) + (t362 * t619 - t368 * t542 + t373 * t493 + t377 * t466 - t386 * t518 - t407 * t506 - t408 * t433 + t416 * t425) * MDP(27) + (t505 * t619 - t506 * t565 + t517 * t718 + t518 * t550) * MDP(12) + (t505 * t565 - t517 * t550) * MDP(11) + (-pkin(7) * t680 + t601 * t639) * MDP(9) - MDP(7) * t681 + (pkin(7) * t681 + t603 * t639) * MDP(10) + (-t461 * t619 + t506 * t685 + t518 * t526 + t542 * t618) * MDP(20) + (t461 * t685 + t526 * t618) * MDP(18) + (-(-t528 * t661 + t651) * t542 - t665 * t506 + t616 * t619 - t444 * t518 + t471 * t526 - t707 * t461 + t453 * t685 + t618 * t489) * MDP(24) + (-t360 * t493 + t362 * t494 - t367 * t466 + t368 * t624 - t386 * t434 + t387 * t433 - t406 * t416 + t407 * t417) * MDP(28) + (-t366 * t493 - t369 * t624 - t370 * t466 + t392 * t434 + t393 * t433 - t409 * t417 - t410 * t416 + t494 * t677) * MDP(25) + (t366 * t410 + t392 * t369 + t393 * t370 - t409 * t677 + t418 * t622 + t458 * t612) * MDP(26) + t711 * t712 + 0.2e1 * t603 * MDP(4) * t642 + (t360 * t406 + t362 * t407 + t367 * t387 + t368 * t386 + t373 * t425 + t377 * t408) * MDP(30) + ((-t524 * t602 - t526 * t599) * t517 + (-t695 - t462 * t602 + (t524 * t599 - t526 * t602) * qJD(4)) * t565) * MDP(19) + (t462 * t619 - t506 * t686 - t518 * t524 - t714 * t542) * MDP(21) + (MDP(13) * t517 - MDP(14) * t518 - MDP(16) * t471 - MDP(17) * t470) * t594; t605 * t711 + t608 + (t588 * t462 + t625 * t599 + t632 * t524 + (-t587 * t660 - t483 + (-t634 + t515) * t599) * t542 + t620) * MDP(23) - t601 * MDP(4) * t679 + (t366 * t502 + t677 * t501 + t418 * (t588 - t700) + (t709 + t632) * t458 + (t460 - t401) * t393 + t670 * t392) * MDP(26) + (-t417 * t488 + t502 * t506 + t542 * t676 + t624 * t671 + t614) * MDP(29) + (t401 * t466 - t624 * t670 + t610 + t630) * MDP(25) + (t396 * t466 - t624 * t675 + t611 + t630) * MDP(28) + (t416 * t488 - t466 * t671 - t501 * t506 + t542 * t675 + t615) * MDP(27) + (t515 * t594 + (t550 * t663 - t594 * t644) * pkin(2) + t607) * MDP(17) + (t588 * t461 + t625 * t602 + t632 * t526 + (t587 * t661 - t627 + t666) * t542 + t628) * MDP(24) + (t360 * t502 + t362 * t501 + t373 * t488 - t386 * t675 + t387 * t676 - t408 * t671) * MDP(30) + (t514 * t594 + (-t594 * t662 + t663 * t718) * pkin(2) + t617) * MDP(16) + (MDP(9) * t601 * t605 + MDP(10) * t679) * pkin(1); t608 + (t512 * t594 + t617) * MDP(16) + (-t417 * t504 + t506 * t523 + t542 * t674 + t624 * t672 + t614) * MDP(29) + (t403 * t466 - t624 * t669 + t610 + t629) * MDP(25) + (t398 * t466 - t624 * t673 + t611 + t629) * MDP(28) + (t366 * t523 + t677 * t522 + t418 * t650 + (-t512 + t709) * t458 + (t486 - t403) * t393 + t669 * t392) * MDP(26) + (t416 * t504 - t466 * t672 - t506 * t522 + t542 * t673 + t615) * MDP(27) + (t360 * t523 + t362 * t522 + t373 * t504 - t386 * t673 + t387 * t674 - t408 * t672) * MDP(30) + (t511 * t594 + t607) * MDP(17) + (-pkin(3) * t462 - t637 * t542 - t512 * t524 - t489 * t687 + (-t542 * t660 - t692) * pkin(9) + t620) * MDP(23) + (-pkin(3) * t461 + t667 * t542 - t512 * t526 - t602 * t694 + (t542 * t661 - t691) * pkin(9) + t628) * MDP(24); t526 * t524 * MDP(18) + (-t524 ^ 2 + t526 ^ 2) * MDP(19) + (t668 + t689) * MDP(20) + (t688 - t693) * MDP(21) + (t444 * t542 - t489 * t526 + t638) * MDP(23) + (t443 * t542 + t489 * t524 - t652) * MDP(24) + (t393 * t624 - t697) * MDP(25) + (t392 * t394 - t393 * t395) * MDP(26) + (t394 * t542 - t677 - t710) * MDP(27) + (t387 * t624 - t416 * t583 + t417 * t585 - t697) * MDP(28) + (-t395 * t542 + t413 * t624 + t653 + 0.2e1 * t659) * MDP(29) + (t360 * t583 + t362 * t585 - t386 * t394 + t387 * t658 - t408 * t413) * MDP(30) + (MDP(22) + (pkin(5) - t585) * MDP(27) + t583 * MDP(29)) * t506 + ((MDP(21) * t550 - MDP(23) * t490) * t602 + (MDP(20) * t550 - MDP(21) * t594 - MDP(23) * t484 + MDP(24) * t490) * t599) * qJD(4) + ((-t416 * t597 - t417 * t598) * MDP(25) + (t366 * t597 - t458 * t526 - t598 * t677) * MDP(26)) * pkin(4) + ((-t392 + t395) * MDP(25) - t413 * MDP(27) + (t386 - t658) * MDP(28) - t408 * MDP(29)) * t466; (t392 * t624 + t393 * t466 + t418) * MDP(26) + (t542 * t624 + t416) * MDP(27) + (-t417 + t717) * MDP(29) + (-t386 * t624 + t387 * t466 + t373) * MDP(30) + (MDP(25) + MDP(28)) * (-t466 ^ 2 - t713); (t417 + t717) * MDP(28) + (-t704 - t713) * MDP(29) + (-t387 * t542 + t362 + t710) * MDP(30) + (t466 * t624 - t606) * MDP(27);];
tauc  = t1;
