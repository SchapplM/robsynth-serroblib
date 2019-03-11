% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:53
% EndTime: 2019-03-09 06:11:14
% DurationCPUTime: 13.68s
% Computational Cost: add. (31040->657), mult. (60819->878), div. (0->0), fcn. (71890->8), ass. (0->352)
t432 = sin(qJ(5));
t435 = cos(qJ(5));
t483 = t432 * pkin(5) - qJ(6) * t435;
t489 = t435 * mrSges(7,1) + t432 * mrSges(7,3);
t639 = Ifges(7,5) * t435;
t394 = Ifges(7,3) * t432 + t639;
t641 = Ifges(6,4) * t432;
t404 = Ifges(6,1) * t435 - t641;
t659 = -t435 / 0.2e1;
t661 = t432 / 0.2e1;
t662 = -t432 / 0.2e1;
t401 = t432 * Ifges(7,1) - t639;
t399 = t435 * Ifges(6,2) + t641;
t565 = t432 * t399;
t421 = Ifges(7,5) * t432;
t700 = t435 * Ifges(7,3) - t421;
t566 = t432 * t700;
t657 = t435 / 0.2e1;
t424 = Ifges(6,4) * t435;
t699 = -t432 * Ifges(6,1) - t424;
t696 = -t566 / 0.2e1 - t565 / 0.2e1 + t401 * t657 + t699 * t659;
t704 = -Ifges(6,2) * t432 + t424;
t707 = Ifges(7,1) * t435 + t421;
t693 = t404 * t661 - t662 * t707 + t696 + (t394 - t704) * t659;
t459 = -t483 * t489 + t693;
t428 = t432 ^ 2;
t429 = t435 ^ 2;
t776 = t428 + t429;
t430 = sin(pkin(10));
t431 = cos(pkin(10));
t434 = sin(qJ(3));
t436 = cos(qJ(3));
t381 = -t434 * t430 + t431 * t436;
t382 = t430 * t436 + t434 * t431;
t433 = sin(qJ(4));
t655 = cos(qJ(4));
t350 = -t655 * t381 + t382 * t433;
t742 = t350 * t435;
t541 = mrSges(7,2) * t742;
t468 = t433 * t381 + t382 * t655;
t718 = t468 * mrSges(7,1);
t768 = -t541 - t718;
t677 = t768 / 0.2e1;
t750 = mrSges(6,1) * t468 + mrSges(6,3) * t742;
t775 = t677 - t750 / 0.2e1;
t651 = pkin(3) * t433;
t416 = pkin(9) + t651;
t582 = t416 * t435;
t739 = t432 * t350;
t751 = -mrSges(6,2) * t468 + mrSges(6,3) * t739;
t752 = mrSges(7,2) * t739 + mrSges(7,3) * t468;
t767 = t751 + t752;
t774 = t767 * t582;
t618 = t435 * mrSges(7,3);
t621 = t432 * mrSges(7,1);
t391 = -t618 + t621;
t578 = t432 * qJ(6);
t484 = t435 * pkin(5) + t578;
t384 = -pkin(4) - t484;
t543 = t655 * pkin(3);
t370 = -t543 + t384;
t688 = m(7) / 0.2e1;
t491 = (-t370 - t384) * t688;
t619 = t435 * mrSges(6,2);
t622 = t432 * mrSges(6,1);
t392 = t619 + t622;
t417 = -t543 - pkin(4);
t581 = t417 * t392;
t585 = t384 * t391;
t649 = pkin(4) * t392;
t671 = t370 / 0.2e1;
t773 = t649 / 0.2e1 - t391 * t671 - t585 / 0.2e1 - t581 / 0.2e1 + t483 * t491 - t459;
t645 = pkin(7) + qJ(2);
t386 = t645 * t430;
t387 = t645 * t431;
t355 = -t386 * t436 - t434 * t387;
t300 = -t382 * pkin(8) + t355;
t356 = -t434 * t386 + t387 * t436;
t301 = t381 * pkin(8) + t356;
t182 = t655 * t300 - t433 * t301;
t113 = t468 * t483 - t182;
t245 = t391 * t468;
t246 = t392 * t468;
t527 = -pkin(2) * t431 - pkin(1);
t360 = -pkin(3) * t381 + t527;
t625 = t468 * mrSges(5,1);
t744 = t350 * mrSges(5,2);
t499 = -t744 + t625;
t590 = t350 * qJ(6);
t709 = t300 * t433 + t655 * t301;
t557 = t435 * t709;
t647 = pkin(9) * t468;
t184 = pkin(4) * t350 + t360 - t647;
t573 = t432 * t184;
t94 = t557 + t573;
t69 = t94 + t590;
t93 = t184 * t435 - t432 * t709;
t71 = -pkin(5) * t350 - t93;
t727 = -t350 * t483 + t709;
t740 = t392 * t350;
t741 = t391 * t350;
t772 = t113 * t741 - t182 * t740 - t727 * t245 - t709 * t246 - t360 * t499 - t69 * t752 - t71 * t768 - t93 * t750 - t94 * t751;
t770 = qJ(6) * t752;
t769 = t776 * t655;
t490 = t435 * mrSges(6,1) - t432 * mrSges(6,2);
t736 = t709 * t490;
t743 = t709 * mrSges(5,1);
t754 = t727 * t489;
t762 = t182 * mrSges(5,2);
t766 = -t736 - t743 - t754 - t762;
t765 = -t754 / 0.2e1 - t762 / 0.2e1;
t763 = t744 / 0.2e1;
t724 = mrSges(7,2) + mrSges(6,3);
t722 = Ifges(7,2) + Ifges(6,3);
t761 = Ifges(5,2) + t722 - Ifges(5,1);
t760 = t113 * t727;
t759 = t182 * t432;
t758 = t182 * t433;
t599 = t182 * t709;
t756 = t384 * t727;
t755 = t435 * t182;
t587 = t468 * t432;
t542 = mrSges(6,3) * t587;
t253 = -mrSges(6,2) * t350 - t542;
t626 = t350 * mrSges(7,3);
t260 = -mrSges(7,2) * t587 + t626;
t711 = t253 + t260;
t648 = pkin(9) * t350;
t650 = pkin(4) * t468;
t262 = t648 + t650;
t708 = -t699 + t401;
t749 = -t742 * t708 / 0.4e1 - t736 / 0.2e1 - t743 / 0.2e1;
t675 = -t350 / 0.2e1;
t746 = -t718 / 0.2e1;
t745 = pkin(4) * t709;
t643 = t71 + t93;
t738 = t435 * t643;
t457 = t724 * (t429 / 0.2e1 + t428 / 0.2e1);
t735 = -t404 - t707;
t734 = -t489 - t490;
t733 = t769 * pkin(3);
t586 = t468 * t435;
t330 = Ifges(7,5) * t586;
t157 = t350 * Ifges(7,6) + Ifges(7,3) * t587 + t330;
t159 = t350 * Ifges(6,6) + t468 * t704;
t681 = -Ifges(6,6) / 0.2e1;
t535 = Ifges(7,6) / 0.2e1 + t681;
t492 = t535 * t350;
t732 = t157 / 0.2e1 - t159 / 0.2e1 + t492;
t579 = t429 * t350;
t580 = t428 * t350;
t731 = -t579 - t580;
t652 = pkin(3) * t382;
t190 = t652 + t262;
t97 = t190 * t435 - t759;
t98 = t432 * t190 + t755;
t480 = -t97 * t432 + t98 * t435;
t717 = qJ(6) * t468;
t75 = t98 + t717;
t726 = pkin(5) * t468;
t76 = -t97 - t726;
t482 = t76 * t432 + t75 * t435;
t705 = Ifges(7,4) * t435 + Ifges(7,6) * t432;
t706 = Ifges(6,5) * t435 - Ifges(6,6) * t432;
t703 = t706 + t705;
t536 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t537 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t728 = (-t535 * t432 - t537 * t435) * t350 + t536 * t468;
t725 = mrSges(6,1) + mrSges(7,1);
t723 = Ifges(7,4) + Ifges(6,5);
t628 = t468 * mrSges(5,3);
t716 = t384 * t468;
t714 = t537 * t350;
t713 = t245 + t246;
t395 = Ifges(6,5) * t432 + Ifges(6,6) * t435;
t397 = Ifges(7,4) * t432 - Ifges(7,6) * t435;
t698 = t397 / 0.4e1 + t395 / 0.4e1;
t418 = m(7) * qJ(6) + mrSges(7,3);
t257 = mrSges(6,1) * t350 - mrSges(6,3) * t586;
t258 = -mrSges(7,1) * t350 + mrSges(7,2) * t586;
t695 = t257 * t659 + t258 * t657 + t711 * t662;
t242 = t484 * t468;
t694 = -t483 * t245 / 0.2e1 + t242 * t489 / 0.2e1 + t182 * t392 / 0.2e1 - t113 * t391 / 0.2e1;
t692 = t382 ^ 2;
t691 = 2 * qJD(4);
t690 = m(6) / 0.2e1;
t689 = -m(7) / 0.2e1;
t687 = -pkin(4) / 0.2e1;
t686 = m(5) * pkin(3);
t685 = m(7) * pkin(3);
t684 = mrSges(6,1) / 0.2e1;
t683 = -mrSges(7,1) / 0.2e1;
t682 = -mrSges(6,2) / 0.2e1;
t680 = t76 / 0.2e1;
t679 = -t97 / 0.2e1;
t678 = t98 / 0.2e1;
t676 = t468 / 0.2e1;
t673 = t350 / 0.4e1;
t672 = -t468 / 0.2e1;
t670 = -t384 / 0.2e1;
t669 = t384 / 0.2e1;
t668 = -t489 / 0.2e1;
t667 = -t490 / 0.2e1;
t664 = -t399 / 0.4e1;
t663 = t417 / 0.2e1;
t660 = t432 / 0.4e1;
t658 = -t435 / 0.4e1;
t656 = t435 / 0.4e1;
t654 = m(7) * t483;
t653 = m(7) * t432;
t646 = pkin(9) * t435;
t644 = t69 - t94;
t640 = Ifges(5,5) * t350;
t637 = Ifges(5,6) * t468;
t156 = Ifges(7,6) * t468 - t394 * t350;
t158 = Ifges(6,6) * t468 - t350 * t704;
t160 = Ifges(7,4) * t468 - t350 * t707;
t162 = Ifges(6,5) * t468 - t350 * t404;
t500 = t382 * mrSges(4,1) + t381 * mrSges(4,2);
t163 = t350 * Ifges(6,5) + t404 * t468;
t558 = t435 * t163;
t161 = t350 * Ifges(7,4) + t468 * t707;
t559 = t435 * t161;
t576 = t432 * t159;
t577 = t432 * t157;
t3 = m(7) * (t69 * t75 + t71 * t76 + t760) + (mrSges(5,1) * t652 - t559 / 0.2e1 - t558 / 0.2e1 - t577 / 0.2e1 + t576 / 0.2e1 + (-t703 / 0.2e1 + Ifges(5,4)) * t350 + t761 * t468) * t350 + m(6) * (t93 * t97 + t94 * t98 - t599) - t709 * t628 + (mrSges(5,2) * t652 + t709 * mrSges(5,3) + t156 * t661 + t158 * t662 + t703 * t676 + (t160 + t162) * t657 - Ifges(5,4) * t468) * t468 + m(5) * t360 * t652 + (Ifges(4,4) * t381 + (Ifges(4,1) - Ifges(4,2)) * t382) * t381 - t692 * Ifges(4,4) + t527 * t500 + t98 * t253 + t97 * t257 + t76 * t258 + t75 * t260 - t772;
t630 = t3 * qJD(1);
t107 = t262 * t435 - t759;
t108 = t432 * t262 + t755;
t461 = t163 / 0.2e1 + t161 / 0.2e1 + t714;
t82 = t108 + t717;
t83 = -t107 - t726;
t4 = -t107 * t257 - t108 * t253 - t82 * t260 - t83 * t258 - m(6) * (t107 * t93 + t108 * t94 - t599) - m(7) * (t69 * t82 + t71 * t83 + t760) + (-Ifges(5,4) * t350 + t732 * t432 + t461 * t435) * t350 + (((Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t429 + ((Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t432 + (-Ifges(6,4) + Ifges(7,5)) * t435) * t432 - t761) * t350 + (Ifges(5,4) - t703) * t468) * t468 + t772;
t623 = t4 * qJD(1);
t419 = t435 * mrSges(7,2);
t616 = t69 * t432;
t247 = t700 * t468;
t248 = t399 * t468;
t249 = -Ifges(7,1) * t587 + t330;
t250 = t699 * t468;
t465 = t489 * t468;
t467 = t490 * t468;
t9 = -t94 * t258 - m(7) * (t113 * t242 + t71 * t94) - t242 * t245 - t113 * t465 + t94 * t257 + t182 * t467 + ((-t248 / 0.2e1 - t247 / 0.2e1 + t71 * mrSges(7,2) + t461) * t432 + (-t249 / 0.2e1 - t250 / 0.2e1 + t69 * mrSges(7,2) + t94 * mrSges(6,3) - t732) * t435) * t468 + (-m(7) * t69 - t542 - t711) * t93;
t613 = t9 * qJD(1);
t612 = t94 * t432;
t597 = t182 * t468;
t10 = (t381 ^ 2 + t692) * mrSges(4,3) + (t713 + t628) * t468 - (-t350 * mrSges(5,3) + t711 * t435 + (-t257 + t258) * t432) * t350 + m(7) * (t113 * t468 - (t432 * t71 + t435 * t69) * t350) + m(6) * (-t597 - (-t432 * t93 + t435 * t94) * t350) + m(5) * (-t350 * t709 - t597) + m(4) * (-t355 * t382 + t356 * t381) + (m(3) * qJ(2) + mrSges(3,3)) * (t430 ^ 2 + t431 ^ 2);
t609 = qJD(1) * t10;
t30 = m(7) * (-t113 * t586 + t350 * t69) + t350 * t260 - t245 * t586;
t608 = qJD(1) * t30;
t546 = t686 / 0.2e1;
t439 = -t350 * t433 * t546 + (t370 * t688 + t417 * t690 - t546 * t655 + t667 + t668) * t468 + t724 * (-t580 / 0.2e1 - t579 / 0.2e1) + (t688 + t690) * t731 * t416;
t444 = (t432 * t98 + t435 * t97) * t690 + (t432 * t75 - t435 * t76) * t688 + t750 * t657 + t768 * t659 + t382 * t546 + t767 * t661;
t11 = t439 - t444 - t499 - t500;
t607 = t11 * qJD(1);
t605 = t113 * t483;
t603 = t113 * t432;
t569 = t432 * t258;
t448 = (t432 * t643 + t435 * t644) * t689 + t257 * t661 - t569 / 0.2e1 + t711 * t659;
t449 = (t618 / 0.2e1 - t619 / 0.2e1 - t622 / 0.2e1 - t621 / 0.2e1 + t483 * t689) * t350;
t12 = -t449 + t448;
t601 = t12 * qJD(1);
t440 = -m(6) * (t107 * t435 + t108 * t432) / 0.2e1 + (t432 * t82 - t435 * t83) * t689 + t763 - t625 / 0.2e1 + t768 * t657 + t750 * t659 + t767 * t662;
t529 = t667 - mrSges(5,1) / 0.2e1;
t550 = t731 * pkin(9);
t443 = -t457 * t350 + (t668 + t529) * t468 + t763 + (t550 - t650) * t690 + (t550 + t716) * t688;
t15 = t440 + t443;
t600 = t15 * qJD(1);
t592 = t468 * t395;
t591 = t468 * t397;
t575 = t432 * t160;
t574 = t432 * t162;
t572 = t432 * t245;
t571 = t432 * t750;
t570 = t432 * t768;
t561 = t435 * t156;
t560 = t435 * t158;
t556 = t435 * t489;
t549 = t733 * pkin(9);
t505 = 0.2e1 * t675;
t237 = t505 * t653;
t547 = t237 * qJD(1);
t544 = t83 * t688;
t534 = mrSges(7,2) * t661;
t533 = mrSges(6,3) * t662;
t532 = t419 / 0.2e1;
t531 = mrSges(6,3) * t657;
t530 = t612 / 0.2e1;
t528 = -t572 / 0.2e1 - t556 * t672 + t350 * t532;
t526 = t432 * t655;
t525 = t435 * t655;
t521 = -t742 / 0.2e1;
t514 = t572 / 0.2e1;
t511 = -t739 / 0.2e1;
t504 = t664 - t700 / 0.4e1;
t503 = -t699 / 0.4e1 + t401 / 0.4e1;
t501 = t391 + t654;
t494 = -t543 / 0.2e1;
t493 = t543 / 0.2e1;
t481 = t83 * t432 + t82 * t435;
t479 = t432 * t494;
t478 = t435 * t493;
t445 = -t694 + t586 * t664 - t248 * t656 + t247 * t658 - t576 / 0.4e1 + t577 / 0.4e1 + mrSges(6,3) * t530 + t94 * t533 + t558 / 0.4e1 + t559 / 0.4e1 + t394 * t587 / 0.4e1 + t703 * t673 + (t250 + t249) * t660 + t643 * t532 + (t530 - t616 / 0.2e1) * mrSges(7,2) - (t704 + t708) * t587 / 0.4e1 + (-t700 - t735) * t586 / 0.4e1;
t438 = t445 + ((t612 - t616 + t738) * t688 - t468 * t457 + t695) * t416 + (t242 * t370 + t605) * t688 + (t489 * t671 + t490 * t663) * t468;
t447 = (-pkin(5) * t76 + qJ(6) * t75) * t689 + pkin(5) * t677 - t770 / 0.2e1 - t75 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t680 + mrSges(6,1) * t679 + mrSges(6,2) * t678;
t5 = t438 + t447 - t728;
t90 = t370 * t501 + t459 + t581;
t477 = t5 * qJD(1) + t90 * qJD(3);
t453 = -mrSges(5,2) * t543 + (-mrSges(5,1) + t734) * t651 + t724 * t733;
t460 = t769 * t416;
t115 = -(t370 * t433 + t460) * t685 - m(6) * (t417 * t433 + t460) * pkin(3) - t453;
t475 = -t107 * t432 + t108 * t435;
t437 = t775 * t416 * t432 + t765 + (t417 * t709 + t475 * t416 + (t525 * t94 - t526 * t93 - t758) * pkin(3)) * t690 + t749 + t713 * t651 / 0.2e1 + t774 / 0.2e1 + t711 * t478 + t257 * t479 - t741 * t671 - t740 * t663 + (t370 * t727 + t481 * t416 + (t113 * t433 + t525 * t69 + t526 * t71) * pkin(3)) * t688 + Ifges(5,6) * t672 + t565 * t673 + Ifges(5,5) * t675 + t493 * t569 + (-t704 * t656 - t394 * t658 + t566 / 0.4e1 + t735 * t660) * t350 + (Ifges(6,6) * t656 + Ifges(7,6) * t658 + t660 * t723 + t698) * t468 + t108 * t531 + t82 * t532 + t107 * t533 + t83 * t534;
t442 = (pkin(9) * t480 - t745) * t690 + (pkin(9) * t482 + t756) * t688 - t740 * t687 - t741 * t669 + t765;
t2 = -t437 + t98 * t531 + t75 * t532 + t97 * t533 + t76 * t534 + t442 - t561 / 0.4e1 + t560 / 0.4e1 + t574 / 0.4e1 + t575 / 0.4e1 + t592 / 0.4e1 + t591 / 0.4e1 - t637 / 0.2e1 - t640 / 0.2e1 + t767 * t646 / 0.2e1 - (-t565 / 0.4e1 - t566 / 0.4e1) * t350 + (t570 / 0.2e1 - t571 / 0.2e1) * pkin(9) + t749;
t476 = -t2 * qJD(1) - t115 * qJD(3);
t474 = t242 * t384 + t605;
t455 = m(7) * (-t603 + (t350 * t416 - t370 * t468) * t435);
t471 = m(7) * t680 + t746;
t26 = t514 + (mrSges(7,2) * t505 + t468 * t668) * t435 - t455 / 0.2e1 + t471;
t353 = (m(7) * t370 - t489) * t432;
t473 = qJD(1) * t26 + qJD(3) * t353;
t34 = t626 + 0.2e1 * (t590 / 0.2e1 + t557 / 0.4e1 + t573 / 0.4e1 - t94 / 0.4e1) * m(7);
t472 = qJD(1) * t34 + qJD(5) * t418;
t441 = (-pkin(5) * t526 + qJ(6) * t525) * t685 / 0.2e1 + t494 * t619 + mrSges(7,3) * t478 + t725 * t479;
t45 = t441 + t773;
t446 = (-pkin(5) * t83 + qJ(6) * t82) * t688 - pkin(5) * t768 / 0.2e1 + t770 / 0.2e1 + t107 * t684 + t108 * t682 + t82 * mrSges(7,3) / 0.2e1 + t83 * t683;
t8 = (-t706 / 0.4e1 - t705 / 0.4e1) * t350 + t474 * t689 + (t248 / 0.4e1 + t247 / 0.4e1 - t163 / 0.4e1 - t161 / 0.4e1 - t714 + (-t93 / 0.2e1 - t71 / 0.2e1) * mrSges(7,2) + (t643 * t689 - t258 / 0.2e1 + t257 / 0.2e1) * pkin(9)) * t435 + (t159 / 0.4e1 - t157 / 0.4e1 - t250 / 0.4e1 - t249 / 0.4e1 - t492 + (-t94 / 0.2e1 + t69 / 0.2e1) * mrSges(7,2) + (t644 * t688 + t260 / 0.2e1 + t253 / 0.2e1) * pkin(9)) * t432 + ((-t404 / 0.4e1 - t707 / 0.4e1 + pkin(4) * t684 + mrSges(7,1) * t670 - t504) * t435 + (t704 / 0.4e1 - t394 / 0.4e1 + pkin(4) * t682 + mrSges(7,3) * t670 + t503) * t432 + t457 * pkin(9) + t536) * t468 + t446 + t694;
t99 = t384 * t501 + t459 - t649;
t463 = t8 * qJD(1) + t45 * qJD(3) - t99 * qJD(4);
t263 = (-t489 + (t493 + t671 + t669) * m(7)) * t432;
t456 = m(7) * (-t603 + (t648 - t716) * t435);
t28 = -t541 + t514 + (t683 - t556 / 0.2e1) * t468 + t544 - t456 / 0.2e1;
t354 = (m(7) * t384 - t489) * t432;
t462 = qJD(1) * t28 + qJD(3) * t263 + qJD(4) * t354;
t452 = (-mrSges(7,2) * t578 - pkin(5) * t419 + t703) * qJD(5);
t451 = qJD(5) * (-m(7) * t484 + t734);
t385 = m(7) * t646 + t419;
t366 = m(7) * t582 + t419;
t264 = (m(7) * t493 + t489 + t491) * t432;
t236 = m(7) * t511 + t688 * t739;
t46 = t441 - t773;
t32 = 0.2e1 * t688 * t69 + t260;
t29 = t456 / 0.2e1 - t541 / 0.2e1 + t746 + t544 + t528;
t27 = t455 / 0.2e1 + mrSges(7,2) * t521 + t471 + t528;
t16 = -t440 + t443;
t14 = t439 + t444;
t13 = -t449 - t448;
t7 = t465 * t669 + t467 * t687 + t474 * t688 + t445 + t446 + ((-t432 * t644 + t738) * t688 + t695) * pkin(9) + t728 - t724 * t776 * t647 / 0.2e1;
t6 = Ifges(7,6) * t511 + t521 * t723 + t722 * t676 - t681 * t739 + t438 - t447;
t1 = t437 + (-Ifges(5,6) / 0.2e1 + t698) * t468 - (Ifges(5,5) / 0.2e1 + t503 * t435 + t504 * t432) * t350 + (mrSges(6,3) * t679 + mrSges(7,2) * t680 + t162 / 0.4e1 + t160 / 0.4e1 + t775 * pkin(9)) * t432 + (mrSges(6,3) * t678 + t75 * mrSges(7,2) / 0.2e1 - t156 / 0.4e1 + t158 / 0.4e1 + (t752 / 0.2e1 + t751 / 0.2e1) * pkin(9)) * t435 + t529 * t709 + t442;
t17 = [qJD(2) * t10 + qJD(3) * t3 - qJD(4) * t4 - qJD(5) * t9 + qJD(6) * t30, qJD(3) * t14 + qJD(4) * t16 + qJD(5) * t13 + qJD(6) * t236 + t609, t14 * qJD(2) + t1 * qJD(4) + t6 * qJD(5) + t27 * qJD(6) + t630 + (t766 + (-t655 * t709 + t758) * t686 + t591 / 0.2e1 + t592 / 0.2e1 - t640 + (m(7) * t727 - t741) * t370 + (m(6) * t709 - t740) * t417 + t574 / 0.2e1 + t575 / 0.2e1 - (-mrSges(5,3) * t543 + t696) * t350 - t637 + t774 + t482 * mrSges(7,2) + t480 * mrSges(6,3) + (m(6) * t480 + m(7) * t482 + t570 - t571) * t416 + Ifges(4,5) * t381 - Ifges(4,6) * t382 - t355 * mrSges(4,2) - t356 * mrSges(4,1) - t628 * t651 + t560 / 0.2e1 - t561 / 0.2e1) * qJD(3), -t623 + t16 * qJD(2) + t1 * qJD(3) + t7 * qJD(5) + t29 * qJD(6) + ((pkin(9) * t481 + t756) * t688 + (pkin(9) * t475 - t745) * t690) * t691 + ((t82 * mrSges(7,2) + t108 * mrSges(6,3) + pkin(9) * t752) * t435 + (t83 * mrSges(7,2) - t107 * mrSges(6,3) + pkin(9) * t768) * t432 + (t395 / 0.2e1 + t397 / 0.2e1 - Ifges(5,6) + (-pkin(9) * mrSges(6,2) - t535) * t435 + (-pkin(9) * mrSges(6,1) + t537) * t432) * t468 + (-Ifges(5,5) - t585 + t649 - t693) * t350 + t766) * qJD(4), t13 * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + t32 * qJD(6) - t613 + ((-m(7) * pkin(5) - t725) * t94 + (-mrSges(6,2) + t418) * t93 + ((-mrSges(7,2) * qJ(6) - Ifges(6,6) + Ifges(7,6)) * t435 + (mrSges(7,2) * pkin(5) - t723) * t432) * t468) * qJD(5), qJD(2) * t236 + qJD(3) * t27 + qJD(4) * t29 + qJD(5) * t32 + t608; -qJD(3) * t11 - qJD(4) * t15 - qJD(5) * t12 - qJD(6) * t237 - t609, 0, -t607, -t600, -t601 + (t618 - t619 - t654) * qJD(5) + (m(7) * qJD(6) - t725 * qJD(5)) * t432, qJD(5) * t653 - t547; qJD(2) * t11 - qJD(4) * t2 + qJD(5) * t5 - qJD(6) * t26 - t630, t607, -qJD(4) * t115 + qJD(5) * t90 - qJD(6) * t353, t453 * qJD(4) + t46 * qJD(5) + t264 * qJD(6) + ((t384 * t651 + t549) * t688 + (-pkin(4) * t651 + t549) * t690) * t691 + t476, t46 * qJD(4) + t366 * qJD(6) + t416 * t451 + t452 + t477, qJD(4) * t264 + qJD(5) * t366 - t473; qJD(2) * t15 + qJD(3) * t2 - qJD(5) * t8 - qJD(6) * t28 + t623, t600, -qJD(5) * t45 - qJD(6) * t263 - t476, qJD(5) * t99 - qJD(6) * t354, pkin(9) * t451 + t385 * qJD(6) + t452 - t463, qJD(5) * t385 - t462; qJD(2) * t12 - qJD(3) * t5 + qJD(4) * t8 + qJD(6) * t34 + t613, t601, qJD(4) * t45 - t477, t463, t418 * qJD(6), t472; qJD(2) * t237 + qJD(3) * t26 + qJD(4) * t28 - qJD(5) * t34 - t608, t547, qJD(4) * t263 + t473, t462, -t472, 0;];
Cq  = t17;
