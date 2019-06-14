% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 13:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:02:18
% EndTime: 2019-05-06 13:02:32
% DurationCPUTime: 10.86s
% Computational Cost: add. (145967->367), mult. (342832->447), div. (0->0), fcn. (249246->10), ass. (0->142)
t778 = Ifges(5,1) + Ifges(6,2);
t777 = -Ifges(6,1) - Ifges(5,3);
t772 = Ifges(5,4) + Ifges(6,6);
t771 = Ifges(5,5) - Ifges(6,4);
t776 = Ifges(5,2) + Ifges(6,3);
t770 = Ifges(5,6) - Ifges(6,5);
t775 = -2 * qJD(5);
t774 = cos(qJ(4));
t746 = qJD(1) ^ 2;
t773 = pkin(2) * t746;
t737 = sin(pkin(10));
t738 = cos(pkin(10));
t741 = sin(qJ(2));
t744 = cos(qJ(2));
t713 = (-t737 * t741 + t738 * t744) * qJD(1);
t714 = (t737 * t744 + t738 * t741) * qJD(1);
t740 = sin(qJ(4));
t694 = -t774 * t713 + t714 * t740;
t735 = qJD(2) + qJD(4);
t769 = t694 * t735;
t742 = sin(qJ(1));
t745 = cos(qJ(1));
t730 = -g(1) * t745 - g(2) * t742;
t719 = -pkin(1) * t746 + qJDD(1) * pkin(7) + t730;
t768 = t719 * t741;
t760 = qJD(1) * qJD(2);
t724 = qJDD(1) * t741 + t744 * t760;
t676 = qJDD(2) * pkin(2) - qJ(3) * t724 - t768 + (qJ(3) * t760 + t741 * t773 - g(3)) * t744;
t704 = -g(3) * t741 + t744 * t719;
t725 = qJDD(1) * t744 - t741 * t760;
t762 = qJD(1) * t741;
t726 = qJD(2) * pkin(2) - qJ(3) * t762;
t736 = t744 ^ 2;
t677 = qJ(3) * t725 - qJD(2) * t726 - t736 * t773 + t704;
t645 = -0.2e1 * qJD(3) * t714 + t738 * t676 - t677 * t737;
t702 = t724 * t738 + t725 * t737;
t632 = (qJD(2) * t713 - t702) * pkin(8) + (t713 * t714 + qJDD(2)) * pkin(3) + t645;
t646 = 0.2e1 * qJD(3) * t713 + t737 * t676 + t738 * t677;
t701 = -t724 * t737 + t725 * t738;
t707 = qJD(2) * pkin(3) - pkin(8) * t714;
t712 = t713 ^ 2;
t635 = -pkin(3) * t712 + pkin(8) * t701 - qJD(2) * t707 + t646;
t629 = t774 * t632 - t740 * t635;
t656 = -t694 * qJD(4) + t740 * t701 + t774 * t702;
t695 = t740 * t713 + t774 * t714;
t671 = mrSges(5,1) * t694 + mrSges(5,2) * t695;
t682 = -mrSges(5,2) * t735 - mrSges(5,3) * t694;
t684 = mrSges(6,1) * t694 - mrSges(6,3) * t735;
t734 = qJDD(2) + qJDD(4);
t670 = pkin(4) * t694 - qJ(5) * t695;
t733 = t735 ^ 2;
t626 = -t734 * pkin(4) - t733 * qJ(5) + t695 * t670 + qJDD(5) - t629;
t621 = (t694 * t695 - t734) * pkin(9) + (t656 + t769) * pkin(5) + t626;
t655 = qJD(4) * t695 - t774 * t701 + t702 * t740;
t686 = pkin(5) * t695 - pkin(9) * t735;
t690 = t694 ^ 2;
t729 = t742 * g(1) - t745 * g(2);
t755 = -qJDD(1) * pkin(1) - t729;
t679 = -t725 * pkin(2) + qJDD(3) + t726 * t762 + (-qJ(3) * t736 - pkin(7)) * t746 + t755;
t644 = -t701 * pkin(3) - t712 * pkin(8) + t714 * t707 + t679;
t748 = (-t656 + t769) * qJ(5) + t644 + (pkin(4) * t735 + t775) * t695;
t624 = (pkin(4) + pkin(9)) * t655 - t695 * t686 - t690 * pkin(5) + t748;
t739 = sin(qJ(6));
t743 = cos(qJ(6));
t619 = t621 * t743 - t624 * t739;
t680 = t694 * t743 - t735 * t739;
t638 = qJD(6) * t680 + t655 * t739 + t734 * t743;
t654 = qJDD(6) + t656;
t681 = t694 * t739 + t735 * t743;
t657 = -mrSges(7,1) * t680 + mrSges(7,2) * t681;
t689 = qJD(6) + t695;
t658 = -mrSges(7,2) * t689 + mrSges(7,3) * t680;
t617 = m(7) * t619 + mrSges(7,1) * t654 - mrSges(7,3) * t638 - t657 * t681 + t658 * t689;
t620 = t621 * t739 + t624 * t743;
t637 = -qJD(6) * t681 + t655 * t743 - t734 * t739;
t659 = mrSges(7,1) * t689 - mrSges(7,3) * t681;
t618 = m(7) * t620 - mrSges(7,2) * t654 + mrSges(7,3) * t637 + t657 * t680 - t659 * t689;
t609 = t617 * t743 + t618 * t739;
t672 = -mrSges(6,2) * t694 - mrSges(6,3) * t695;
t753 = -m(6) * t626 - t656 * mrSges(6,1) - t695 * t672 - t609;
t607 = m(5) * t629 - mrSges(5,3) * t656 - t671 * t695 + (t682 - t684) * t735 + (mrSges(5,1) - mrSges(6,2)) * t734 + t753;
t630 = t740 * t632 + t774 * t635;
t683 = mrSges(5,1) * t735 - mrSges(5,3) * t695;
t752 = -pkin(4) * t733 + qJ(5) * t734 - t670 * t694 + t630;
t625 = t735 * t775 - t752;
t685 = mrSges(6,1) * t695 + mrSges(6,2) * t735;
t623 = -pkin(5) * t655 - pkin(9) * t690 + ((2 * qJD(5)) + t686) * t735 + t752;
t754 = -m(7) * t623 + mrSges(7,1) * t637 - t638 * mrSges(7,2) + t658 * t680 - t681 * t659;
t750 = -m(6) * t625 + t734 * mrSges(6,3) + t735 * t685 - t754;
t614 = m(5) * t630 - mrSges(5,2) * t734 - t683 * t735 + (-t671 - t672) * t694 + (-mrSges(5,3) - mrSges(6,1)) * t655 + t750;
t603 = t774 * t607 + t740 * t614;
t698 = -mrSges(4,1) * t713 + mrSges(4,2) * t714;
t705 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t713;
t601 = m(4) * t645 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t702 + qJD(2) * t705 - t698 * t714 + t603;
t706 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t714;
t756 = -t607 * t740 + t774 * t614;
t602 = m(4) * t646 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t701 - qJD(2) * t706 + t698 * t713 + t756;
t596 = t738 * t601 + t737 * t602;
t703 = -g(3) * t744 - t768;
t723 = (-mrSges(3,1) * t744 + mrSges(3,2) * t741) * qJD(1);
t761 = qJD(1) * t744;
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t761;
t594 = m(3) * t703 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t724 + qJD(2) * t728 - t723 * t762 + t596;
t727 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t762;
t757 = -t601 * t737 + t738 * t602;
t595 = m(3) * t704 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t725 - qJD(2) * t727 + t723 * t761 + t757;
t758 = -t594 * t741 + t744 * t595;
t588 = m(2) * t730 - mrSges(2,1) * t746 - qJDD(1) * mrSges(2,2) + t758;
t718 = -pkin(7) * t746 + t755;
t628 = t655 * pkin(4) + t748;
t766 = -t739 * t617 + t743 * t618;
t608 = m(6) * t628 - t655 * mrSges(6,2) - t656 * mrSges(6,3) - t694 * t684 - t695 * t685 + t766;
t751 = m(5) * t644 + t655 * mrSges(5,1) + t656 * mrSges(5,2) + t694 * t682 + t695 * t683 + t608;
t749 = m(4) * t679 - t701 * mrSges(4,1) + t702 * mrSges(4,2) - t713 * t705 + t714 * t706 + t751;
t747 = -m(3) * t718 + t725 * mrSges(3,1) - t724 * mrSges(3,2) - t727 * t762 + t728 * t761 - t749;
t605 = m(2) * t729 + qJDD(1) * mrSges(2,1) - t746 * mrSges(2,2) + t747;
t767 = t742 * t588 + t745 * t605;
t589 = t744 * t594 + t741 * t595;
t765 = t770 * t694 - t771 * t695 + t777 * t735;
t764 = t776 * t694 - t772 * t695 - t770 * t735;
t763 = -t772 * t694 + t778 * t695 + t771 * t735;
t759 = t745 * t588 - t605 * t742;
t717 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t741 + Ifges(3,4) * t744) * qJD(1);
t716 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t741 + Ifges(3,2) * t744) * qJD(1);
t715 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t741 + Ifges(3,6) * t744) * qJD(1);
t693 = Ifges(4,1) * t714 + Ifges(4,4) * t713 + Ifges(4,5) * qJD(2);
t692 = Ifges(4,4) * t714 + Ifges(4,2) * t713 + Ifges(4,6) * qJD(2);
t691 = Ifges(4,5) * t714 + Ifges(4,6) * t713 + Ifges(4,3) * qJD(2);
t641 = Ifges(7,1) * t681 + Ifges(7,4) * t680 + Ifges(7,5) * t689;
t640 = Ifges(7,4) * t681 + Ifges(7,2) * t680 + Ifges(7,6) * t689;
t639 = Ifges(7,5) * t681 + Ifges(7,6) * t680 + Ifges(7,3) * t689;
t611 = mrSges(7,2) * t623 - mrSges(7,3) * t619 + Ifges(7,1) * t638 + Ifges(7,4) * t637 + Ifges(7,5) * t654 + t639 * t680 - t640 * t689;
t610 = -mrSges(7,1) * t623 + mrSges(7,3) * t620 + Ifges(7,4) * t638 + Ifges(7,2) * t637 + Ifges(7,6) * t654 - t639 * t681 + t641 * t689;
t597 = mrSges(6,1) * t626 + mrSges(7,1) * t619 + mrSges(5,2) * t644 - mrSges(7,2) * t620 - mrSges(5,3) * t629 - mrSges(6,3) * t628 + Ifges(7,5) * t638 + Ifges(7,6) * t637 + Ifges(7,3) * t654 + pkin(5) * t609 - qJ(5) * t608 + t681 * t640 - t680 * t641 + t764 * t735 + t771 * t734 + t765 * t694 + t778 * t656 - t772 * t655;
t590 = -mrSges(5,1) * t644 - mrSges(6,1) * t625 + mrSges(6,2) * t628 + mrSges(5,3) * t630 - pkin(4) * t608 - pkin(5) * t754 - pkin(9) * t766 - t743 * t610 - t739 * t611 - t776 * t655 + t772 * t656 + t765 * t695 + t770 * t734 + t763 * t735;
t585 = mrSges(4,2) * t679 - mrSges(4,3) * t645 + Ifges(4,1) * t702 + Ifges(4,4) * t701 + Ifges(4,5) * qJDD(2) - pkin(8) * t603 - qJD(2) * t692 - t740 * t590 + t774 * t597 + t713 * t691;
t584 = -mrSges(4,1) * t679 + mrSges(4,3) * t646 + Ifges(4,4) * t702 + Ifges(4,2) * t701 + Ifges(4,6) * qJDD(2) - pkin(3) * t751 + pkin(8) * t756 + qJD(2) * t693 + t774 * t590 + t740 * t597 - t714 * t691;
t583 = t746 * Ifges(2,5) + t739 * t610 - t743 * t611 + mrSges(2,3) * t730 - t714 * t692 - Ifges(3,5) * t724 - Ifges(3,6) * t725 + mrSges(3,2) * t704 + t713 * t693 - Ifges(4,6) * t701 - Ifges(4,5) * t702 - mrSges(3,1) * t703 - mrSges(4,1) * t645 + mrSges(4,2) * t646 - mrSges(5,1) * t629 + mrSges(5,2) * t630 + mrSges(6,3) * t625 - mrSges(6,2) * t626 + pkin(9) * t609 + mrSges(2,1) * g(3) - pkin(3) * t603 - pkin(2) * t596 + (pkin(4) * mrSges(6,2) + t777) * t734 - qJ(5) * t750 + (qJ(5) * t672 - t763) * t694 + t764 * t695 - pkin(1) * t589 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (qJ(5) * mrSges(6,1) + t770) * t655 - t771 * t656 - pkin(4) * (-t684 * t735 + t753) + (-t741 * t716 + t744 * t717) * qJD(1);
t582 = mrSges(3,2) * t718 - mrSges(3,3) * t703 + Ifges(3,1) * t724 + Ifges(3,4) * t725 + Ifges(3,5) * qJDD(2) - qJ(3) * t596 - qJD(2) * t716 - t584 * t737 + t585 * t738 + t715 * t761;
t581 = -mrSges(3,1) * t718 + mrSges(3,3) * t704 + Ifges(3,4) * t724 + Ifges(3,2) * t725 + Ifges(3,6) * qJDD(2) - pkin(2) * t749 + qJ(3) * t757 + qJD(2) * t717 + t738 * t584 + t737 * t585 - t715 * t762;
t580 = -mrSges(2,2) * g(3) - mrSges(2,3) * t729 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t746 - pkin(7) * t589 - t581 * t741 + t582 * t744;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t767; (-m(1) - m(2)) * g(3) + t589; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t767 + t745 * t580 - t742 * t583; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t759 + t742 * t580 + t745 * t583; -mrSges(1,1) * g(2) + mrSges(2,1) * t729 + mrSges(1,2) * g(1) - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) + pkin(1) * t747 + pkin(7) * t758 + t744 * t581 + t741 * t582;];
tauB  = t1;
