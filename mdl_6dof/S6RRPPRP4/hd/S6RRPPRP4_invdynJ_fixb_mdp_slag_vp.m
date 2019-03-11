% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:59
% EndTime: 2019-03-09 08:40:12
% DurationCPUTime: 11.74s
% Computational Cost: add. (5308->662), mult. (11787->798), div. (0->0), fcn. (8089->8), ass. (0->262)
t592 = cos(pkin(9));
t595 = sin(qJ(2));
t694 = qJD(1) * t595;
t665 = t592 * t694;
t591 = sin(pkin(9));
t692 = qJD(2) * t591;
t523 = t665 + t692;
t598 = cos(qJ(2));
t681 = qJD(1) * qJD(2);
t660 = t598 * t681;
t679 = qJDD(1) * t595;
t619 = t660 + t679;
t483 = qJDD(2) * t591 + t592 * t619;
t667 = t591 * t694;
t683 = t592 * qJD(2);
t522 = -t667 + t683;
t594 = sin(qJ(5));
t597 = cos(qJ(5));
t574 = t592 * qJDD(2);
t606 = t619 * t591 - t574;
t685 = qJD(5) * t597;
t673 = -t597 * t483 + t522 * t685 - t594 * t606;
t686 = qJD(5) * t594;
t397 = t523 * t686 + t673;
t693 = qJD(1) * t598;
t559 = qJD(5) + t693;
t631 = t597 * t522 + t523 * t594;
t746 = t559 * t631;
t756 = t397 - t746;
t757 = t756 * MDP(27);
t719 = t522 * t594;
t448 = -t523 * t597 + t719;
t735 = t448 ^ 2;
t748 = t448 * t559;
t656 = t483 * t594 - t597 * t606;
t398 = -qJD(5) * t448 + t656;
t755 = -t398 - t748;
t526 = t591 * t594 + t592 * t597;
t498 = t526 * t598;
t507 = t526 * qJD(5);
t701 = -qJD(1) * t498 - t507;
t664 = t592 * t693;
t666 = t591 * t693;
t700 = t591 * t685 - t592 * t686 - t594 * t664 + t597 * t666;
t577 = t598 * qJDD(1);
t661 = t595 * t681;
t738 = t577 - t661;
t659 = -t591 * qJ(4) - pkin(2);
t754 = t592 * pkin(3) - t659;
t753 = MDP(11) + MDP(15);
t678 = MDP(24) + MDP(26);
t677 = MDP(25) - MDP(28);
t752 = -2 * pkin(1);
t751 = t559 ^ 2;
t750 = t631 ^ 2;
t578 = t595 * qJ(3);
t583 = t598 * pkin(2);
t671 = -pkin(1) - t583;
t627 = t671 - t578;
t517 = t627 * qJD(1);
t572 = pkin(7) * t693;
t540 = qJD(2) * qJ(3) + t572;
t455 = t517 * t592 - t591 * t540;
t434 = pkin(3) * t693 + qJD(4) - t455;
t411 = pkin(4) * t693 - pkin(8) * t523 + t434;
t456 = t591 * t517 + t592 * t540;
t443 = -qJ(4) * t693 + t456;
t415 = -pkin(8) * t522 + t443;
t387 = t411 * t594 + t415 * t597;
t384 = qJ(6) * t559 + t387;
t749 = t384 * t559;
t747 = t526 * t595;
t716 = t591 * t597;
t527 = -t592 * t594 + t716;
t621 = t595 * t527;
t596 = sin(qJ(1));
t599 = cos(qJ(1));
t643 = g(1) * t599 + g(2) * t596;
t745 = t595 * t643;
t640 = pkin(2) * t595 - qJ(3) * t598;
t505 = qJD(2) * t640 - qJD(3) * t595;
t446 = qJD(1) * t505 + qJDD(1) * t627;
t489 = pkin(7) * t738 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t409 = t446 * t592 - t591 * t489;
t403 = t738 * pkin(3) + qJDD(4) - t409;
t744 = t443 * t693 + t403;
t670 = -pkin(7) * t591 - pkin(3);
t713 = t592 * t598;
t610 = -pkin(8) * t713 + (-pkin(4) + t670) * t595;
t530 = t640 * qJD(1);
t718 = t530 * t592;
t424 = qJD(1) * t610 - t718;
t513 = t591 * t530;
t567 = qJ(4) * t694;
t714 = t592 * t595;
t715 = t591 * t598;
t622 = -pkin(7) * t714 + pkin(8) * t715;
t436 = qJD(1) * t622 + t513 + t567;
t726 = pkin(8) - qJ(3);
t536 = t726 * t591;
t537 = t726 * t592;
t630 = -t536 * t597 + t537 * t594;
t743 = -qJD(3) * t526 - qJD(5) * t630 + t594 * t424 + t597 * t436;
t465 = -t536 * t594 - t537 * t597;
t742 = qJD(3) * t527 - qJD(5) * t465 - t424 * t597 + t436 * t594;
t696 = t583 + t578;
t533 = -pkin(1) - t696;
t552 = pkin(7) * t715;
t582 = t598 * pkin(3);
t444 = pkin(4) * t598 + t552 + t582 + (-pkin(8) * t595 - t533) * t592;
t485 = pkin(7) * t713 + t591 * t533;
t467 = -qJ(4) * t598 + t485;
t717 = t591 * t595;
t454 = pkin(8) * t717 + t467;
t741 = t594 * t444 + t597 * t454;
t688 = qJD(4) * t591;
t699 = qJ(4) * t664 - t572;
t734 = -pkin(3) - pkin(4);
t655 = -t666 * t734 + t688 - t699;
t569 = pkin(7) * t679;
t500 = -qJDD(2) * pkin(2) + pkin(7) * t660 + qJDD(3) + t569;
t587 = g(3) * t598;
t657 = -t500 - t587;
t691 = qJD(2) * t595;
t739 = qJ(4) * t691 - qJD(4) * t598;
t515 = t523 * qJD(4);
t400 = pkin(3) * t606 - t483 * qJ(4) + t500 - t515;
t611 = t587 + t400;
t737 = t398 * pkin(5) + t397 * qJ(6) + qJD(6) * t448;
t720 = t505 * t592;
t418 = qJD(2) * t610 - t720;
t488 = t591 * t505;
t419 = qJD(2) * t622 + t488 + t739;
t736 = -qJD(5) * t741 + t418 * t597 - t419 * t594;
t520 = t523 ^ 2;
t601 = qJD(1) ^ 2;
t733 = pkin(1) * t601;
t525 = -qJDD(5) - t738;
t732 = pkin(5) * t525;
t731 = pkin(7) * t523;
t730 = g(1) * t596;
t593 = qJD(2) * pkin(2);
t725 = qJ(3) * t592;
t724 = qJ(6) * t525;
t723 = t387 * t559;
t721 = t631 * t448;
t712 = t595 * t596;
t711 = t595 * t599;
t710 = t596 * t592;
t709 = t596 * t598;
t708 = t598 * t599;
t707 = t599 * t591;
t706 = t700 * pkin(5) - t701 * qJ(6) - qJD(6) * t527 + t655;
t705 = qJ(6) * t694 - t743;
t704 = -pkin(5) * t694 - t742;
t410 = t591 * t446 + t592 * t489;
t663 = t598 * t683;
t698 = -qJ(4) * t663 - qJD(4) * t714;
t697 = g(1) * t712 - g(2) * t711;
t588 = t595 ^ 2;
t589 = t598 ^ 2;
t695 = t588 - t589;
t690 = qJD(2) * t598;
t689 = qJD(3) * t592;
t386 = t411 * t597 - t415 * t594;
t682 = qJD(6) - t386;
t680 = qJDD(1) * qJ(4);
t676 = pkin(7) * t691;
t675 = t597 * t715;
t674 = qJ(4) * t661 + t410;
t672 = g(1) * t708 + g(2) * t709 + g(3) * t595;
t669 = pkin(3) * t591 + pkin(7);
t532 = pkin(7) * t694 + qJD(3) - t593;
t662 = qJ(3) * t577;
t484 = t533 * t592 - t552;
t514 = t592 * pkin(4) + t754;
t390 = pkin(4) * t738 - pkin(8) * t483 + t403;
t637 = t591 * t679 - t574;
t393 = t637 * pkin(8) + (-t680 + (pkin(8) * t692 - qJD(4)) * qJD(1)) * t598 + t674;
t654 = -t597 * t390 + t594 * t393 + t411 * t686 + t415 * t685;
t653 = g(1) * t592 * t711 + g(2) * t595 * t710 + qJD(3) * t666 + t591 * t662;
t652 = pkin(3) * t713 + qJ(4) * t715 + t696;
t651 = t599 * pkin(1) + pkin(2) * t708 + t596 * pkin(7) + qJ(3) * t711;
t650 = t592 * t662;
t649 = t591 * t734 - pkin(7);
t648 = t670 * t595;
t509 = t591 * t709 + t592 * t599;
t510 = t592 * t709 - t707;
t439 = t509 * t597 - t510 * t594;
t511 = t598 * t707 - t710;
t512 = t591 * t596 + t592 * t708;
t441 = -t511 * t597 + t512 * t594;
t647 = g(1) * t439 + g(2) * t441;
t438 = t509 * t594 + t510 * t597;
t442 = t511 * t594 + t512 * t597;
t646 = g(1) * t438 - g(2) * t442;
t645 = -g(1) * t509 + g(2) * t511;
t644 = g(1) * t510 - g(2) * t512;
t642 = -g(2) * t599 + t730;
t584 = t599 * pkin(7);
t641 = -t510 * pkin(3) - qJ(4) * t509 + t584;
t495 = t594 * t714 - t595 * t716;
t639 = -pkin(5) * t495 + qJ(6) * t747;
t636 = t522 * t689 - t606 * t725 - t672;
t429 = -t522 * pkin(3) - t523 * qJ(4) + t532;
t633 = t444 * t597 - t454 * t594;
t629 = qJ(3) * t483 + qJD(3) * t523;
t600 = qJD(2) ^ 2;
t628 = qJDD(2) * t598 - t595 * t600;
t473 = -pkin(7) * t665 + t513;
t462 = -t592 * t676 + t488;
t624 = t591 * t660 - t574;
t617 = t594 * t390 + t597 * t393 + t411 * t685 - t415 * t686;
t616 = t594 * t418 + t597 * t419 + t444 * t685 - t454 * t686;
t549 = qJ(4) * t714;
t466 = t595 * t649 + t549;
t414 = pkin(4) * t522 - t429;
t614 = t512 * pkin(3) + qJ(4) * t511 + t651;
t613 = t627 * t730;
t612 = -g(1) * t511 - g(2) * t509 - g(3) * t717;
t480 = t596 * t747;
t482 = t599 * t747;
t609 = g(1) * t482 + g(2) * t480 - g(3) * t498 - t525 * t630;
t479 = t596 * t621;
t481 = t599 * t621;
t497 = t594 * t713 - t675;
t608 = -g(1) * t481 - g(2) * t479 - g(3) * t497 - t465 * t525;
t433 = t649 * t690 - t698;
t604 = g(1) * t441 - g(2) * t439 + g(3) * t495 - t654;
t603 = -g(1) * t442 - g(2) * t438 - g(3) * t747 + t617;
t385 = pkin(5) * t631 + qJ(6) * t448 + t414;
t602 = -t385 * t448 + qJDD(6) - t604;
t394 = -pkin(4) * t606 - t400;
t557 = qJ(3) * t708;
t553 = qJ(3) * t709;
t501 = t522 * t693;
t490 = t595 * t669 - t549;
t478 = pkin(3) * t666 - t699;
t472 = t667 * pkin(7) + t718;
t471 = -t484 + t582;
t461 = t591 * t676 + t720;
t460 = qJD(1) * t648 - t718;
t459 = t473 + t567;
t458 = t669 * t690 + t698;
t447 = qJD(2) * t648 - t720;
t432 = t462 + t739;
t428 = qJD(2) * t498 + qJD(5) * t621;
t427 = -qJD(2) * t675 + t507 * t595 + t594 * t663;
t423 = pkin(5) * t526 - qJ(6) * t527 + t514;
t406 = t466 - t639;
t405 = -pkin(5) * t448 + qJ(6) * t631;
t402 = -pkin(5) * t598 - t633;
t401 = qJ(6) * t598 + t741;
t399 = (-qJD(1) * qJD(4) - t680) * t598 + t674;
t382 = -pkin(5) * t559 + t682;
t381 = pkin(5) * t427 - qJ(6) * t428 - qJD(6) * t747 + t433;
t380 = pkin(5) * t691 - t736;
t379 = -qJ(6) * t691 + qJD(6) * t598 + t616;
t378 = t394 + t737;
t377 = qJDD(6) + t654 + t732;
t376 = qJD(6) * t559 + t617 - t724;
t1 = [(-t398 * t598 - t427 * t559 + t495 * t525 + t631 * t691) * MDP(22) + (-t377 * t598 + t378 * t495 - t380 * t559 + t381 * t631 + t382 * t691 + t385 * t427 + t398 * t406 + t402 * t525 + t646) * MDP(26) + ((-pkin(7) * qJDD(2) + t681 * t752) * t595 + (0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t600 + t642) * t598) * MDP(9) + (-pkin(7) * t628 + t619 * t752 - t697) * MDP(10) + (-t386 * t691 + t394 * t495 + t466 * t398 + t414 * t427 + t433 * t631 - t633 * t525 + t559 * t736 - t654 * t598 + t646) * MDP(24) + (t376 * t401 + t384 * t379 + t378 * t406 + t385 * t381 + t377 * t402 + t382 * t380 - g(1) * (-pkin(4) * t510 - pkin(5) * t438 + qJ(6) * t439 + t641) - g(2) * (pkin(4) * t512 + pkin(5) * t442 - pkin(8) * t711 + qJ(6) * t441 + t614) - (t595 * t726 + t671) * t730) * MDP(29) + (-t376 * t495 + t377 * t747 - t379 * t631 - t380 * t448 + t382 * t428 - t384 * t427 - t397 * t402 - t398 * t401 - t697) * MDP(27) + (t397 * t495 - t398 * t747 + t427 * t448 - t428 * t631) * MDP(20) + (-t397 * t598 + t428 * t559 + t448 * t691 - t525 * t747) * MDP(21) + (t376 * t598 - t378 * t747 + t379 * t559 + t381 * t448 - t384 * t691 - t385 * t428 + t397 * t406 - t401 * t525 - t647) * MDP(28) + (-t397 * t747 - t428 * t448) * MDP(19) + (t387 * t691 + t394 * t747 - t466 * t397 + t414 * t428 - t433 * t448 + t525 * t741 - t559 * t616 - t598 * t617 + t647) * MDP(25) + (-t461 * t523 + t462 * t522 - t484 * t483 + t485 * t574 + (-t409 * t595 - t455 * t690) * t592 + (-t410 * t595 - t456 * t690 - t485 * t619) * t591 + t697) * MDP(13) + (t432 * t522 + t447 * t523 + t467 * t574 + t471 * t483 + (t403 * t595 + t434 * t690) * t592 + (-t399 * t595 - t443 * t690 - t467 * t619) * t591 + t697) * MDP(16) + 0.2e1 * (t577 * t595 - t681 * t695) * MDP(5) + (t410 * t485 + t456 * t462 + t409 * t484 + t455 * t461 - g(1) * t584 - g(2) * t651 - t613 + (t500 * t595 + t532 * t690) * pkin(7)) * MDP(14) + (-t525 * t598 - t559 * t691) * MDP(23) + (-t458 * t522 - t490 * t574 + ((qJDD(1) * t490 + t400) * t591 + (-qJD(1) * t471 - t434) * qJD(2)) * t595 + (t429 * t692 + t471 * qJDD(1) + t403 + (t490 * t692 + t447) * qJD(1)) * t598 + t644) * MDP(15) + (-t458 * t523 - t483 * t490 + (-t400 * t592 + (qJD(1) * t467 + t443) * qJD(2)) * t595 + (-qJD(1) * t432 - qJDD(1) * t467 - t429 * t683 - t399) * t598 - t645) * MDP(17) + ((t500 * t591 + (qJD(1) * t484 + t455) * qJD(2) + t637 * pkin(7)) * t595 + (-t461 * qJD(1) - t484 * qJDD(1) - t409 + (t532 * t591 + (-t522 + t667) * pkin(7)) * qJD(2)) * t598 + t644) * MDP(11) + (qJDD(1) * t588 + 0.2e1 * t595 * t660) * MDP(4) + (qJDD(2) * t595 + t598 * t600) * MDP(6) + qJDD(1) * MDP(1) + (-g(1) * t641 - g(2) * t614 + t399 * t467 + t400 * t490 + t403 * t471 + t429 * t458 + t443 * t432 + t434 * t447 - t613) * MDP(18) + t642 * MDP(2) + t643 * MDP(3) + t628 * MDP(7) + ((pkin(7) * t483 + t500 * t592 + (-qJD(1) * t485 - t456) * qJD(2)) * t595 + (qJD(1) * t462 + qJDD(1) * t485 + t410 + (t532 * t592 + t731) * qJD(2)) * t598 + t645) * MDP(12); (-MDP(4) * t595 * t598 + t695 * MDP(5)) * t601 + (t378 * t526 - t382 * t694 + t385 * t700 + t398 * t423 - t559 * t704 + t631 * t706 + t609) * MDP(26) + (t525 * t526 - t559 * t700 - t631 * t694) * MDP(22) + (t650 - pkin(2) * t483 + (-t657 - t745) * t591 + ((-qJ(3) * t683 + t456) * t595 + (-t731 - t473 + (qJD(3) - t532) * t592) * t598) * qJD(1)) * MDP(12) + (t386 * t694 + t394 * t526 + t514 * t398 + t700 * t414 + t559 * t742 + t655 * t631 + t609) * MDP(24) + (-t459 * t522 - t460 * t523 + (-t434 * t693 + t399) * t592 + (t629 + t744) * t591 + t636) * MDP(16) + t559 * MDP(23) * t694 + (t397 * t526 - t398 * t527 + t448 * t700 - t631 * t701) * MDP(20) + (-t378 * t527 + t384 * t694 - t385 * t701 + t397 * t423 + t448 * t706 + t559 * t705 + t608) * MDP(28) + (-t397 * t527 - t448 * t701) * MDP(19) + (-t448 * t694 - t525 * t527 + t559 * t701) * MDP(21) + (-t387 * t694 + t394 * t527 - t514 * t397 + t701 * t414 - t448 * t655 + t559 * t743 - t608) * MDP(25) + (-t376 * t526 + t377 * t527 + t382 * t701 - t384 * t700 + t397 * t630 - t398 * t465 - t448 * t704 - t631 * t705 + t672) * MDP(27) + (t376 * t465 + t378 * t423 - t377 * t630 - g(1) * (-pkin(5) * t482 - pkin(8) * t708 + qJ(6) * t481 + t557) - g(2) * (-pkin(5) * t480 - pkin(8) * t709 + qJ(6) * t479 + t553) - g(3) * (pkin(4) * t713 + pkin(5) * t498 + qJ(6) * t497 + t652) + t706 * t385 + t705 * t384 + t704 * t382 + (g(3) * pkin(8) + t643 * (-t592 * t734 - t659)) * t595) * MDP(29) + (-t500 * pkin(2) - t456 * t473 - t455 * t472 - t532 * t572 - g(1) * (-pkin(2) * t711 + t557) - g(2) * (-pkin(2) * t712 + t553) - g(3) * t696 + (-t455 * t591 + t456 * t592) * qJD(3) + (-t409 * t591 + t410 * t592) * qJ(3)) * MDP(14) + (t472 * t523 - t473 * t522 + (t455 * t693 + t410) * t592 + (t456 * t693 - t409 + t629) * t591 + t636) * MDP(13) + (-t650 + t478 * t523 + t483 * t754 + (t515 - t611 + t745) * t591 + (-t443 * t595 + t459 * t598 + (qJ(3) * t691 + (-qJD(3) + t429) * t598) * t592) * qJD(1)) * MDP(17) + (-t754 * t637 - t611 * t592 + (t478 + t688) * t522 + (t434 * t595 - t460 * t598 + (-t429 * t598 + (-t598 * t754 - t578) * qJD(2)) * t591) * qJD(1) + t653) * MDP(15) + qJDD(2) * MDP(8) + ((-pkin(7) * qJDD(1) + t733) * t598 + t672) * MDP(10) + (-t569 - t587 + (t643 + t733) * t595) * MDP(9) + MDP(7) * t577 + MDP(6) * t679 + (t399 * t725 - t429 * t478 - t434 * t460 - g(1) * t557 - g(2) * t553 - g(3) * t652 + (-t459 + t689) * t443 + (qJ(3) * t403 + qJD(3) * t434 - qJD(4) * t429) * t591 + (-t400 + t745) * t754) * MDP(18) + (-pkin(2) * t637 + t657 * t592 + ((-qJ(3) * t692 - t455) * t595 + (pkin(7) * t522 + t472 + (-t532 - t593) * t591) * t598) * qJD(1) + t653) * MDP(11); (t455 * t523 - t456 * t522 - t657) * MDP(14) + (-t434 * t523 - t443 * t522 + t611) * MDP(18) + (t750 + t735) * MDP(27) + (t624 * pkin(4) - t382 * t448 - t384 * t631 + t611 - t737) * MDP(29) + (MDP(12) - MDP(17)) * (-t501 + t483) + (t678 * t719 + (t594 * t677 - t597 * t678) * t523) * qJD(5) + ((MDP(29) * pkin(4) + t753) * t591 * qJDD(1) + (-MDP(14) - MDP(18) - MDP(29)) * t643) * t595 + (MDP(16) + MDP(13)) * (-t522 ^ 2 - t520) + t677 * (t673 + t746) + t678 * (-t656 + t748) + t753 * (-t523 * t693 + t624); t738 * MDP(15) + (t501 + t483) * MDP(16) + (-t589 * t601 - t520) * MDP(17) + (t612 + t744) * MDP(18) + t612 * MDP(29) + (-t522 * MDP(15) + t429 * MDP(18) - t385 * MDP(29) + t448 * t677 - t631 * t678) * t523 + (t757 + (-t377 + t749) * MDP(29) - t678 * t525 - t677 * t751) * t597 + (t755 * MDP(27) + (t382 * t559 + t376) * MDP(29) + t677 * t525 - t678 * t751) * t594; -MDP(19) * t721 + (t735 - t750) * MDP(20) - t756 * MDP(21) + t755 * MDP(22) - t525 * MDP(23) + (t414 * t448 + t604 + t723) * MDP(24) + (t386 * t559 + t414 * t631 - t603) * MDP(25) + (-t405 * t631 - t602 + t723 - 0.2e1 * t732) * MDP(26) + (pkin(5) * t397 - qJ(6) * t398 - (t384 - t387) * t448 + (t382 - t682) * t631) * MDP(27) + (-0.2e1 * t724 - t385 * t631 - t405 * t448 + (0.2e1 * qJD(6) - t386) * t559 + t603) * MDP(28) + (t376 * qJ(6) - t377 * pkin(5) - t385 * t405 - t382 * t387 - g(1) * (-pkin(5) * t441 + qJ(6) * t442) - g(2) * (pkin(5) * t439 + qJ(6) * t438) - g(3) * t639 + t682 * t384) * MDP(29); (t525 - t721) * MDP(26) - t757 + (-t735 - t751) * MDP(28) + (t602 + t732 - t749) * MDP(29);];
tau  = t1;
