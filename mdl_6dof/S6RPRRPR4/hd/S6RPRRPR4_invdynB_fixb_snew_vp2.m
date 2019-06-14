% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:27:54
% EndTime: 2019-05-05 22:28:18
% DurationCPUTime: 23.99s
% Computational Cost: add. (369066->364), mult. (899262->457), div. (0->0), fcn. (714086->12), ass. (0->148)
t727 = qJD(1) ^ 2;
t755 = cos(qJ(4));
t719 = cos(pkin(10));
t754 = pkin(2) * t719;
t717 = sin(pkin(10));
t753 = mrSges(3,2) * t717;
t714 = t719 ^ 2;
t752 = t714 * t727;
t723 = sin(qJ(1));
t726 = cos(qJ(1));
t704 = -g(1) * t726 - g(2) * t723;
t700 = -pkin(1) * t727 + qJDD(1) * qJ(2) + t704;
t748 = qJD(1) * qJD(2);
t746 = -g(3) * t719 - 0.2e1 * t717 * t748;
t673 = (-pkin(7) * qJDD(1) + t727 * t754 - t700) * t717 + t746;
t690 = -g(3) * t717 + (t700 + 0.2e1 * t748) * t719;
t747 = qJDD(1) * t719;
t674 = -pkin(2) * t752 + pkin(7) * t747 + t690;
t722 = sin(qJ(3));
t725 = cos(qJ(3));
t651 = t725 * t673 - t674 * t722;
t735 = t717 * t725 + t719 * t722;
t734 = -t717 * t722 + t719 * t725;
t698 = t734 * qJD(1);
t749 = qJD(3) * t698;
t688 = qJDD(1) * t735 + t749;
t699 = t735 * qJD(1);
t627 = (-t688 + t749) * pkin(8) + (t698 * t699 + qJDD(3)) * pkin(3) + t651;
t652 = t722 * t673 + t725 * t674;
t687 = -qJD(3) * t699 + qJDD(1) * t734;
t693 = qJD(3) * pkin(3) - pkin(8) * t699;
t697 = t698 ^ 2;
t630 = -pkin(3) * t697 + pkin(8) * t687 - qJD(3) * t693 + t652;
t721 = sin(qJ(4));
t620 = t721 * t627 + t755 * t630;
t680 = t721 * t698 + t699 * t755;
t646 = qJD(4) * t680 - t755 * t687 + t688 * t721;
t679 = -t755 * t698 + t699 * t721;
t662 = mrSges(5,1) * t679 + mrSges(5,2) * t680;
t715 = qJD(3) + qJD(4);
t671 = mrSges(5,1) * t715 - mrSges(5,3) * t680;
t712 = qJDD(3) + qJDD(4);
t661 = pkin(4) * t679 - qJ(5) * t680;
t711 = t715 ^ 2;
t612 = -pkin(4) * t711 + qJ(5) * t712 - t661 * t679 + t620;
t713 = t717 ^ 2;
t703 = g(1) * t723 - t726 * g(2);
t739 = qJDD(2) - t703;
t686 = (-pkin(1) - t754) * qJDD(1) + (-qJ(2) + (-t713 - t714) * pkin(7)) * t727 + t739;
t641 = -pkin(3) * t687 - pkin(8) * t697 + t699 * t693 + t686;
t647 = -t679 * qJD(4) + t721 * t687 + t688 * t755;
t615 = (t679 * t715 - t647) * qJ(5) + (t680 * t715 + t646) * pkin(4) + t641;
t716 = sin(pkin(11));
t718 = cos(pkin(11));
t667 = t680 * t718 + t715 * t716;
t607 = -0.2e1 * qJD(5) * t667 - t612 * t716 + t718 * t615;
t639 = t647 * t718 + t712 * t716;
t666 = -t680 * t716 + t715 * t718;
t605 = (t666 * t679 - t639) * pkin(9) + (t666 * t667 + t646) * pkin(5) + t607;
t608 = 0.2e1 * qJD(5) * t666 + t718 * t612 + t716 * t615;
t638 = -t647 * t716 + t712 * t718;
t655 = pkin(5) * t679 - pkin(9) * t667;
t665 = t666 ^ 2;
t606 = -pkin(5) * t665 + pkin(9) * t638 - t655 * t679 + t608;
t720 = sin(qJ(6));
t724 = cos(qJ(6));
t603 = t605 * t724 - t606 * t720;
t648 = t666 * t724 - t667 * t720;
t618 = qJD(6) * t648 + t638 * t720 + t639 * t724;
t649 = t666 * t720 + t667 * t724;
t625 = -mrSges(7,1) * t648 + mrSges(7,2) * t649;
t675 = qJD(6) + t679;
t631 = -mrSges(7,2) * t675 + mrSges(7,3) * t648;
t645 = qJDD(6) + t646;
t601 = m(7) * t603 + mrSges(7,1) * t645 - mrSges(7,3) * t618 - t625 * t649 + t631 * t675;
t604 = t605 * t720 + t606 * t724;
t617 = -qJD(6) * t649 + t638 * t724 - t639 * t720;
t632 = mrSges(7,1) * t675 - mrSges(7,3) * t649;
t602 = m(7) * t604 - mrSges(7,2) * t645 + mrSges(7,3) * t617 + t625 * t648 - t632 * t675;
t593 = t724 * t601 + t720 * t602;
t650 = -mrSges(6,1) * t666 + mrSges(6,2) * t667;
t653 = -mrSges(6,2) * t679 + mrSges(6,3) * t666;
t591 = m(6) * t607 + mrSges(6,1) * t646 - mrSges(6,3) * t639 - t650 * t667 + t653 * t679 + t593;
t654 = mrSges(6,1) * t679 - mrSges(6,3) * t667;
t740 = -t601 * t720 + t724 * t602;
t592 = m(6) * t608 - mrSges(6,2) * t646 + mrSges(6,3) * t638 + t650 * t666 - t654 * t679 + t740;
t741 = -t591 * t716 + t718 * t592;
t586 = m(5) * t620 - mrSges(5,2) * t712 - mrSges(5,3) * t646 - t662 * t679 - t671 * t715 + t741;
t619 = t627 * t755 - t721 * t630;
t670 = -mrSges(5,2) * t715 - mrSges(5,3) * t679;
t611 = -t712 * pkin(4) - t711 * qJ(5) + t680 * t661 + qJDD(5) - t619;
t609 = -t638 * pkin(5) - t665 * pkin(9) + t667 * t655 + t611;
t731 = m(7) * t609 - t617 * mrSges(7,1) + mrSges(7,2) * t618 - t648 * t631 + t632 * t649;
t729 = -m(6) * t611 + t638 * mrSges(6,1) - mrSges(6,2) * t639 + t666 * t653 - t654 * t667 - t731;
t597 = m(5) * t619 + mrSges(5,1) * t712 - mrSges(5,3) * t647 - t662 * t680 + t670 * t715 + t729;
t579 = t721 * t586 + t755 * t597;
t684 = -mrSges(4,1) * t698 + mrSges(4,2) * t699;
t691 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t698;
t577 = m(4) * t651 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t688 + qJD(3) * t691 - t684 * t699 + t579;
t692 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t699;
t742 = t755 * t586 - t597 * t721;
t578 = m(4) * t652 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t687 - qJD(3) * t692 + t684 * t698 + t742;
t572 = t725 * t577 + t722 * t578;
t689 = -t700 * t717 + t746;
t733 = mrSges(3,3) * qJDD(1) + t727 * (-mrSges(3,1) * t719 + t753);
t570 = m(3) * t689 - t717 * t733 + t572;
t743 = -t577 * t722 + t725 * t578;
t571 = m(3) * t690 + t719 * t733 + t743;
t744 = -t570 * t717 + t719 * t571;
t564 = m(2) * t704 - mrSges(2,1) * t727 - qJDD(1) * mrSges(2,2) + t744;
t696 = -qJDD(1) * pkin(1) - qJ(2) * t727 + t739;
t587 = t718 * t591 + t716 * t592;
t732 = m(5) * t641 + t646 * mrSges(5,1) + t647 * mrSges(5,2) + t679 * t670 + t680 * t671 + t587;
t730 = m(4) * t686 - t687 * mrSges(4,1) + mrSges(4,2) * t688 - t698 * t691 + t692 * t699 + t732;
t728 = -m(3) * t696 + mrSges(3,1) * t747 - t730 + (t713 * t727 + t752) * mrSges(3,3);
t583 = t728 + m(2) * t703 - mrSges(2,2) * t727 + (mrSges(2,1) - t753) * qJDD(1);
t751 = t723 * t564 + t726 * t583;
t565 = t719 * t570 + t717 * t571;
t736 = Ifges(3,5) * t717 + Ifges(3,6) * t719;
t750 = t727 * t736;
t745 = t726 * t564 - t583 * t723;
t738 = Ifges(3,1) * t717 + Ifges(3,4) * t719;
t737 = Ifges(3,4) * t717 + Ifges(3,2) * t719;
t678 = Ifges(4,1) * t699 + Ifges(4,4) * t698 + Ifges(4,5) * qJD(3);
t677 = Ifges(4,4) * t699 + Ifges(4,2) * t698 + Ifges(4,6) * qJD(3);
t676 = Ifges(4,5) * t699 + Ifges(4,6) * t698 + Ifges(4,3) * qJD(3);
t658 = Ifges(5,1) * t680 - Ifges(5,4) * t679 + Ifges(5,5) * t715;
t657 = Ifges(5,4) * t680 - Ifges(5,2) * t679 + Ifges(5,6) * t715;
t656 = Ifges(5,5) * t680 - Ifges(5,6) * t679 + Ifges(5,3) * t715;
t635 = Ifges(6,1) * t667 + Ifges(6,4) * t666 + Ifges(6,5) * t679;
t634 = Ifges(6,4) * t667 + Ifges(6,2) * t666 + Ifges(6,6) * t679;
t633 = Ifges(6,5) * t667 + Ifges(6,6) * t666 + Ifges(6,3) * t679;
t623 = Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t675;
t622 = Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t675;
t621 = Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t675;
t595 = mrSges(7,2) * t609 - mrSges(7,3) * t603 + Ifges(7,1) * t618 + Ifges(7,4) * t617 + Ifges(7,5) * t645 + t621 * t648 - t622 * t675;
t594 = -mrSges(7,1) * t609 + mrSges(7,3) * t604 + Ifges(7,4) * t618 + Ifges(7,2) * t617 + Ifges(7,6) * t645 - t621 * t649 + t623 * t675;
t581 = mrSges(6,2) * t611 - mrSges(6,3) * t607 + Ifges(6,1) * t639 + Ifges(6,4) * t638 + Ifges(6,5) * t646 - pkin(9) * t593 - t594 * t720 + t595 * t724 + t633 * t666 - t634 * t679;
t580 = -mrSges(6,1) * t611 + mrSges(6,3) * t608 + Ifges(6,4) * t639 + Ifges(6,2) * t638 + Ifges(6,6) * t646 - pkin(5) * t731 + pkin(9) * t740 + t724 * t594 + t720 * t595 - t667 * t633 + t679 * t635;
t573 = Ifges(5,4) * t647 + Ifges(5,6) * t712 - t680 * t656 + t715 * t658 - mrSges(5,1) * t641 + mrSges(5,3) * t620 - Ifges(6,5) * t639 - Ifges(6,6) * t638 - t667 * t634 + t666 * t635 - mrSges(6,1) * t607 + mrSges(6,2) * t608 - Ifges(7,5) * t618 - Ifges(7,6) * t617 - Ifges(7,3) * t645 - t649 * t622 + t648 * t623 - mrSges(7,1) * t603 + mrSges(7,2) * t604 - pkin(5) * t593 - pkin(4) * t587 + (-Ifges(5,2) - Ifges(6,3)) * t646;
t566 = mrSges(5,2) * t641 - mrSges(5,3) * t619 + Ifges(5,1) * t647 - Ifges(5,4) * t646 + Ifges(5,5) * t712 - qJ(5) * t587 - t580 * t716 + t581 * t718 - t656 * t679 - t657 * t715;
t561 = mrSges(4,2) * t686 - mrSges(4,3) * t651 + Ifges(4,1) * t688 + Ifges(4,4) * t687 + Ifges(4,5) * qJDD(3) - pkin(8) * t579 - qJD(3) * t677 + t566 * t755 - t721 * t573 + t698 * t676;
t560 = -mrSges(4,1) * t686 + mrSges(4,3) * t652 + Ifges(4,4) * t688 + Ifges(4,2) * t687 + Ifges(4,6) * qJDD(3) - pkin(3) * t732 + pkin(8) * t742 + qJD(3) * t678 + t721 * t566 + t573 * t755 - t699 * t676;
t559 = mrSges(2,1) * g(3) - pkin(2) * t572 - pkin(4) * t729 + (Ifges(2,6) - t736) * qJDD(1) - t716 * t581 - t718 * t580 + mrSges(2,3) * t704 - Ifges(5,3) * t712 + mrSges(3,2) * t690 + t698 * t678 - t699 * t677 - t680 * t657 - Ifges(4,6) * t687 - Ifges(4,5) * t688 - mrSges(3,1) * t689 - t679 * t658 - mrSges(4,1) * t651 + mrSges(4,2) * t652 + Ifges(5,6) * t646 - Ifges(5,5) * t647 + mrSges(5,2) * t620 - mrSges(5,1) * t619 - pkin(3) * t579 - qJ(5) * t741 + (-t717 * t737 + t719 * t738 + Ifges(2,5)) * t727 - pkin(1) * t565 - Ifges(4,3) * qJDD(3);
t558 = mrSges(3,2) * t696 - mrSges(3,3) * t689 - pkin(7) * t572 + qJDD(1) * t738 - t560 * t722 + t561 * t725 + t719 * t750;
t557 = -mrSges(3,1) * t696 + mrSges(3,3) * t690 - pkin(2) * t730 + pkin(7) * t743 + qJDD(1) * t737 + t725 * t560 + t722 * t561 - t717 * t750;
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t703 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t727 - qJ(2) * t565 - t557 * t717 + t558 * t719;
t1 = [-m(1) * g(1) + t745; -m(1) * g(2) + t751; (-m(1) - m(2)) * g(3) + t565; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t751 + t726 * t556 - t723 * t559; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t745 + t723 * t556 + t726 * t559; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t703 - mrSges(2,2) * t704 + t717 * t558 + t719 * t557 + pkin(1) * (-qJDD(1) * t753 + t728) + qJ(2) * t744;];
tauB  = t1;
