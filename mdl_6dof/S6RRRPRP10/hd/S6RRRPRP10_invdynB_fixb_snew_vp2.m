% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 09:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:52:42
% EndTime: 2019-05-07 08:53:24
% DurationCPUTime: 26.83s
% Computational Cost: add. (423756->374), mult. (924086->475), div. (0->0), fcn. (732342->12), ass. (0->153)
t784 = Ifges(6,1) + Ifges(7,1);
t775 = Ifges(6,4) - Ifges(7,5);
t783 = -Ifges(6,5) - Ifges(7,4);
t782 = Ifges(6,2) + Ifges(7,3);
t773 = Ifges(6,6) - Ifges(7,6);
t781 = -Ifges(6,3) - Ifges(7,2);
t737 = sin(pkin(6));
t742 = sin(qJ(2));
t744 = cos(qJ(2));
t760 = qJD(1) * qJD(2);
t724 = (-qJDD(1) * t744 + t742 * t760) * t737;
t780 = cos(qJ(3));
t779 = cos(qJ(5));
t778 = pkin(8) * t737;
t739 = cos(pkin(6));
t777 = t739 * g(3);
t776 = -mrSges(6,3) - mrSges(7,2);
t772 = t737 * t742;
t771 = t737 * t744;
t770 = t739 * t742;
t769 = t739 * t744;
t743 = sin(qJ(1));
t745 = cos(qJ(1));
t729 = t743 * g(1) - g(2) * t745;
t746 = qJD(1) ^ 2;
t719 = qJDD(1) * pkin(1) + t746 * t778 + t729;
t730 = -g(1) * t745 - g(2) * t743;
t720 = -pkin(1) * t746 + qJDD(1) * t778 + t730;
t763 = t719 * t770 + t744 * t720;
t693 = -g(3) * t772 + t763;
t733 = qJD(1) * t739 + qJD(2);
t762 = qJD(1) * t737;
t758 = t742 * t762;
t717 = mrSges(3,1) * t733 - mrSges(3,3) * t758;
t721 = (-mrSges(3,1) * t744 + mrSges(3,2) * t742) * t762;
t732 = qJDD(1) * t739 + qJDD(2);
t722 = (-pkin(2) * t744 - pkin(9) * t742) * t762;
t731 = t733 ^ 2;
t761 = qJD(1) * t744;
t669 = -t731 * pkin(2) + t732 * pkin(9) + (-g(3) * t742 + t722 * t761) * t737 + t763;
t723 = (qJDD(1) * t742 + t744 * t760) * t737;
t670 = t724 * pkin(2) - t723 * pkin(9) - t777 + (-t719 + (pkin(2) * t742 - pkin(9) * t744) * t733 * qJD(1)) * t737;
t741 = sin(qJ(3));
t641 = t669 * t780 + t741 * t670;
t713 = t741 * t733 + t758 * t780;
t690 = qJD(3) * t713 + t723 * t741 - t732 * t780;
t712 = -t733 * t780 + t741 * t758;
t695 = mrSges(4,1) * t712 + mrSges(4,2) * t713;
t757 = t737 * t761;
t728 = qJD(3) - t757;
t703 = mrSges(4,1) * t728 - mrSges(4,3) * t713;
t716 = qJDD(3) + t724;
t694 = pkin(3) * t712 - qJ(4) * t713;
t727 = t728 ^ 2;
t631 = -pkin(3) * t727 + qJ(4) * t716 - t694 * t712 + t641;
t692 = -g(3) * t771 + t719 * t769 - t742 * t720;
t668 = -t732 * pkin(2) - t731 * pkin(9) + t722 * t758 - t692;
t691 = -t712 * qJD(3) + t723 * t780 + t741 * t732;
t634 = (t712 * t728 - t691) * qJ(4) + (t713 * t728 + t690) * pkin(3) + t668;
t736 = sin(pkin(11));
t738 = cos(pkin(11));
t701 = t713 * t738 + t728 * t736;
t626 = -0.2e1 * qJD(4) * t701 - t736 * t631 + t738 * t634;
t676 = t691 * t738 + t716 * t736;
t700 = -t713 * t736 + t728 * t738;
t623 = (t700 * t712 - t676) * pkin(10) + (t700 * t701 + t690) * pkin(4) + t626;
t627 = 0.2e1 * qJD(4) * t700 + t738 * t631 + t736 * t634;
t675 = -t691 * t736 + t716 * t738;
t681 = pkin(4) * t712 - pkin(10) * t701;
t699 = t700 ^ 2;
t625 = -pkin(4) * t699 + pkin(10) * t675 - t681 * t712 + t627;
t740 = sin(qJ(5));
t619 = t740 * t623 + t779 * t625;
t673 = t740 * t700 + t701 * t779;
t638 = t673 * qJD(5) - t675 * t779 + t740 * t676;
t711 = qJD(5) + t712;
t658 = mrSges(6,1) * t711 - mrSges(6,3) * t673;
t672 = -t700 * t779 + t740 * t701;
t688 = qJDD(5) + t690;
t651 = pkin(5) * t672 - qJ(6) * t673;
t710 = t711 ^ 2;
t616 = -pkin(5) * t710 + qJ(6) * t688 + 0.2e1 * qJD(6) * t711 - t651 * t672 + t619;
t659 = -mrSges(7,1) * t711 + mrSges(7,2) * t673;
t759 = m(7) * t616 + t688 * mrSges(7,3) + t711 * t659;
t652 = mrSges(7,1) * t672 - mrSges(7,3) * t673;
t764 = -mrSges(6,1) * t672 - mrSges(6,2) * t673 - t652;
t609 = m(6) * t619 - t688 * mrSges(6,2) + t638 * t776 - t711 * t658 + t672 * t764 + t759;
t618 = t623 * t779 - t740 * t625;
t639 = -t672 * qJD(5) + t740 * t675 + t676 * t779;
t657 = -mrSges(6,2) * t711 - mrSges(6,3) * t672;
t617 = -t688 * pkin(5) - t710 * qJ(6) + t673 * t651 + qJDD(6) - t618;
t656 = -mrSges(7,2) * t672 + mrSges(7,3) * t711;
t752 = -m(7) * t617 + t688 * mrSges(7,1) + t711 * t656;
t611 = m(6) * t618 + t688 * mrSges(6,1) + t639 * t776 + t711 * t657 + t673 * t764 + t752;
t604 = t740 * t609 + t779 * t611;
t677 = -mrSges(5,1) * t700 + mrSges(5,2) * t701;
t679 = -mrSges(5,2) * t712 + mrSges(5,3) * t700;
t602 = m(5) * t626 + mrSges(5,1) * t690 - mrSges(5,3) * t676 - t677 * t701 + t679 * t712 + t604;
t680 = mrSges(5,1) * t712 - mrSges(5,3) * t701;
t753 = t779 * t609 - t611 * t740;
t603 = m(5) * t627 - mrSges(5,2) * t690 + mrSges(5,3) * t675 + t677 * t700 - t680 * t712 + t753;
t754 = -t602 * t736 + t738 * t603;
t599 = m(4) * t641 - mrSges(4,2) * t716 - mrSges(4,3) * t690 - t695 * t712 - t703 * t728 + t754;
t640 = -t741 * t669 + t670 * t780;
t702 = -mrSges(4,2) * t728 - mrSges(4,3) * t712;
t630 = -t716 * pkin(3) - t727 * qJ(4) + t713 * t694 + qJDD(4) - t640;
t628 = -t675 * pkin(4) - t699 * pkin(10) + t701 * t681 + t630;
t621 = -0.2e1 * qJD(6) * t673 + (t672 * t711 - t639) * qJ(6) + (t673 * t711 + t638) * pkin(5) + t628;
t614 = m(7) * t621 + t638 * mrSges(7,1) - t639 * mrSges(7,3) + t672 * t656 - t673 * t659;
t749 = m(6) * t628 + t638 * mrSges(6,1) + mrSges(6,2) * t639 + t672 * t657 + t658 * t673 + t614;
t747 = -m(5) * t630 + t675 * mrSges(5,1) - mrSges(5,2) * t676 + t700 * t679 - t680 * t701 - t749;
t613 = m(4) * t640 + mrSges(4,1) * t716 - mrSges(4,3) * t691 - t695 * t713 + t702 * t728 + t747;
t755 = t599 * t780 - t613 * t741;
t590 = m(3) * t693 - mrSges(3,2) * t732 - mrSges(3,3) * t724 - t717 * t733 + t721 * t757 + t755;
t593 = t741 * t599 + t613 * t780;
t707 = -t737 * t719 - t777;
t718 = -mrSges(3,2) * t733 + mrSges(3,3) * t757;
t592 = m(3) * t707 + t724 * mrSges(3,1) + t723 * mrSges(3,2) + (t717 * t742 - t718 * t744) * t762 + t593;
t600 = t602 * t738 + t603 * t736;
t748 = -m(4) * t668 - t690 * mrSges(4,1) - mrSges(4,2) * t691 - t712 * t702 - t703 * t713 - t600;
t596 = m(3) * t692 + mrSges(3,1) * t732 - mrSges(3,3) * t723 + t718 * t733 - t721 * t758 + t748;
t579 = t590 * t770 - t592 * t737 + t596 * t769;
t577 = m(2) * t729 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t746 + t579;
t583 = t744 * t590 - t596 * t742;
t582 = m(2) * t730 - mrSges(2,1) * t746 - qJDD(1) * mrSges(2,2) + t583;
t768 = t745 * t577 + t743 * t582;
t767 = t782 * t672 - t775 * t673 - t773 * t711;
t766 = t773 * t672 + t783 * t673 + t781 * t711;
t765 = -t775 * t672 + t784 * t673 - t783 * t711;
t578 = t590 * t772 + t739 * t592 + t596 * t771;
t756 = -t577 * t743 + t745 * t582;
t605 = -mrSges(6,1) * t628 - mrSges(7,1) * t621 + mrSges(7,2) * t616 + mrSges(6,3) * t619 - pkin(5) * t614 - t782 * t638 + t775 * t639 + t766 * t673 + t773 * t688 + t765 * t711;
t606 = mrSges(6,2) * t628 + mrSges(7,2) * t617 - mrSges(6,3) * t618 - mrSges(7,3) * t621 - qJ(6) * t614 - t775 * t638 + t784 * t639 + t766 * t672 - t688 * t783 + t767 * t711;
t660 = Ifges(5,5) * t701 + Ifges(5,6) * t700 + Ifges(5,3) * t712;
t662 = Ifges(5,1) * t701 + Ifges(5,4) * t700 + Ifges(5,5) * t712;
t585 = -mrSges(5,1) * t630 + mrSges(5,3) * t627 + Ifges(5,4) * t676 + Ifges(5,2) * t675 + Ifges(5,6) * t690 - pkin(4) * t749 + pkin(10) * t753 + t605 * t779 + t740 * t606 - t701 * t660 + t712 * t662;
t661 = Ifges(5,4) * t701 + Ifges(5,2) * t700 + Ifges(5,6) * t712;
t586 = mrSges(5,2) * t630 - mrSges(5,3) * t626 + Ifges(5,1) * t676 + Ifges(5,4) * t675 + Ifges(5,5) * t690 - pkin(10) * t604 - t740 * t605 + t606 * t779 + t700 * t660 - t712 * t661;
t684 = Ifges(4,5) * t713 - Ifges(4,6) * t712 + Ifges(4,3) * t728;
t685 = Ifges(4,4) * t713 - Ifges(4,2) * t712 + Ifges(4,6) * t728;
t575 = mrSges(4,2) * t668 - mrSges(4,3) * t640 + Ifges(4,1) * t691 - Ifges(4,4) * t690 + Ifges(4,5) * t716 - qJ(4) * t600 - t585 * t736 + t586 * t738 - t684 * t712 - t685 * t728;
t686 = Ifges(4,1) * t713 - Ifges(4,4) * t712 + Ifges(4,5) * t728;
t584 = -qJ(6) * t759 + (mrSges(7,2) * qJ(6) + t773) * t638 + (-Ifges(5,3) - Ifges(4,2)) * t690 + t728 * t686 + Ifges(4,6) * t716 - t713 * t684 - pkin(5) * t752 - t701 * t661 + t700 * t662 + Ifges(4,4) * t691 - Ifges(5,6) * t675 - Ifges(5,5) * t676 - mrSges(4,1) * t668 + mrSges(4,3) * t641 - mrSges(5,1) * t626 + mrSges(5,2) * t627 + (qJ(6) * t652 - t765) * t672 + (pkin(5) * t652 + t767) * t673 - mrSges(6,1) * t618 + mrSges(6,2) * t619 + mrSges(7,1) * t617 - mrSges(7,3) * t616 - pkin(4) * t604 - pkin(3) * t600 + (mrSges(7,2) * pkin(5) + t783) * t639 + t781 * t688;
t704 = Ifges(3,3) * t733 + (Ifges(3,5) * t742 + Ifges(3,6) * t744) * t762;
t705 = Ifges(3,6) * t733 + (Ifges(3,4) * t742 + Ifges(3,2) * t744) * t762;
t573 = mrSges(3,2) * t707 - mrSges(3,3) * t692 + Ifges(3,1) * t723 - Ifges(3,4) * t724 + Ifges(3,5) * t732 - pkin(9) * t593 + t575 * t780 - t741 * t584 + t704 * t757 - t733 * t705;
t706 = Ifges(3,5) * t733 + (Ifges(3,1) * t742 + Ifges(3,4) * t744) * t762;
t574 = Ifges(3,4) * t723 - Ifges(3,2) * t724 + Ifges(3,6) * t732 - t704 * t758 + t733 * t706 - mrSges(3,1) * t707 + mrSges(3,3) * t693 - Ifges(4,5) * t691 + Ifges(4,6) * t690 - Ifges(4,3) * t716 - t713 * t685 - t712 * t686 - mrSges(4,1) * t640 + mrSges(4,2) * t641 - t736 * t586 - t738 * t585 - pkin(3) * t747 - qJ(4) * t754 - pkin(2) * t593;
t750 = pkin(8) * t583 + t573 * t742 + t574 * t744;
t572 = Ifges(3,5) * t723 - Ifges(3,6) * t724 + Ifges(3,3) * t732 + mrSges(3,1) * t692 - mrSges(3,2) * t693 + t741 * t575 + t780 * t584 + pkin(2) * t748 + pkin(9) * t755 + (t705 * t742 - t706 * t744) * t762;
t571 = -mrSges(2,2) * g(3) - mrSges(2,3) * t729 + Ifges(2,5) * qJDD(1) - t746 * Ifges(2,6) + t744 * t573 - t742 * t574 + (-t578 * t737 - t579 * t739) * pkin(8);
t570 = mrSges(2,1) * g(3) + mrSges(2,3) * t730 + t746 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t578 - t737 * t572 + t739 * t750;
t1 = [-m(1) * g(1) + t756; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t578; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t768 - t743 * t570 + t745 * t571; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t756 + t745 * t570 + t743 * t571; -mrSges(1,1) * g(2) + mrSges(2,1) * t729 + mrSges(1,2) * g(1) - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) + pkin(1) * t579 + t739 * t572 + t737 * t750;];
tauB  = t1;
