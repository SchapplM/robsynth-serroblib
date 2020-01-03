% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:57
% EndTime: 2019-12-31 22:09:07
% DurationCPUTime: 7.81s
% Computational Cost: add. (119766->305), mult. (256571->386), div. (0->0), fcn. (195375->10), ass. (0->131)
t815 = Ifges(5,1) + Ifges(6,1);
t808 = Ifges(5,4) + Ifges(6,4);
t807 = Ifges(5,5) + Ifges(6,5);
t814 = Ifges(5,2) + Ifges(6,2);
t806 = Ifges(5,6) + Ifges(6,6);
t813 = Ifges(5,3) + Ifges(6,3);
t769 = sin(pkin(5));
t773 = sin(qJ(2));
t777 = cos(qJ(2));
t792 = qJD(1) * qJD(2);
t757 = (-qJDD(1) * t777 + t773 * t792) * t769;
t770 = cos(pkin(5));
t765 = t770 * qJD(1) + qJD(2);
t772 = sin(qJ(3));
t776 = cos(qJ(3));
t794 = qJD(1) * t769;
t789 = t773 * t794;
t745 = t776 * t765 - t772 * t789;
t756 = (qJDD(1) * t773 + t777 * t792) * t769;
t764 = t770 * qJDD(1) + qJDD(2);
t728 = t745 * qJD(3) + t776 * t756 + t772 * t764;
t746 = t772 * t765 + t776 * t789;
t793 = qJD(1) * t777;
t788 = t769 * t793;
t760 = qJD(3) - t788;
t771 = sin(qJ(4));
t775 = cos(qJ(4));
t734 = -t771 * t746 + t775 * t760;
t749 = qJDD(3) + t757;
t695 = t734 * qJD(4) + t775 * t728 + t771 * t749;
t735 = t775 * t746 + t771 * t760;
t712 = -t734 * mrSges(6,1) + t735 * mrSges(6,2);
t755 = (-pkin(2) * t777 - pkin(8) * t773) * t794;
t763 = t765 ^ 2;
t774 = sin(qJ(1));
t778 = cos(qJ(1));
t761 = t774 * g(1) - t778 * g(2);
t779 = qJD(1) ^ 2;
t811 = pkin(7) * t769;
t752 = qJDD(1) * pkin(1) + t779 * t811 + t761;
t762 = -t778 * g(1) - t774 * g(2);
t753 = -t779 * pkin(1) + qJDD(1) * t811 + t762;
t802 = t770 * t773;
t795 = t752 * t802 + t777 * t753;
t709 = -t763 * pkin(2) + t764 * pkin(8) + (-g(3) * t773 + t755 * t793) * t769 + t795;
t810 = t770 * g(3);
t710 = t757 * pkin(2) - t756 * pkin(8) - t810 + (-t752 + (pkin(2) * t773 - pkin(8) * t777) * t765 * qJD(1)) * t769;
t690 = t776 * t709 + t772 * t710;
t732 = -t745 * pkin(3) - t746 * pkin(9);
t759 = t760 ^ 2;
t685 = -t759 * pkin(3) + t749 * pkin(9) + t745 * t732 + t690;
t801 = t770 * t777;
t803 = t769 * t777;
t729 = -g(3) * t803 + t752 * t801 - t773 * t753;
t708 = -t764 * pkin(2) - t763 * pkin(8) + t755 * t789 - t729;
t727 = -t746 * qJD(3) - t772 * t756 + t776 * t764;
t688 = (-t745 * t760 - t728) * pkin(9) + (t746 * t760 - t727) * pkin(3) + t708;
t680 = -t771 * t685 + t775 * t688;
t725 = qJDD(4) - t727;
t744 = qJD(4) - t745;
t677 = -0.2e1 * qJD(5) * t735 + (t734 * t744 - t695) * qJ(5) + (t734 * t735 + t725) * pkin(4) + t680;
t715 = -t744 * mrSges(6,2) + t734 * mrSges(6,3);
t791 = m(6) * t677 + t725 * mrSges(6,1) + t744 * t715;
t674 = -t695 * mrSges(6,3) - t735 * t712 + t791;
t681 = t775 * t685 + t771 * t688;
t694 = -t735 * qJD(4) - t771 * t728 + t775 * t749;
t717 = t744 * pkin(4) - t735 * qJ(5);
t733 = t734 ^ 2;
t679 = -t733 * pkin(4) + t694 * qJ(5) + 0.2e1 * qJD(5) * t734 - t744 * t717 + t681;
t797 = t808 * t734 + t815 * t735 + t807 * t744;
t798 = -t814 * t734 - t808 * t735 - t806 * t744;
t812 = mrSges(5,1) * t680 + mrSges(6,1) * t677 - mrSges(5,2) * t681 - mrSges(6,2) * t679 + pkin(4) * t674 + t806 * t694 + t807 * t695 + t813 * t725 - t797 * t734 - t798 * t735;
t809 = -mrSges(5,2) - mrSges(6,2);
t804 = t769 * t773;
t730 = -g(3) * t804 + t795;
t750 = t765 * mrSges(3,1) - mrSges(3,3) * t789;
t754 = (-mrSges(3,1) * t777 + mrSges(3,2) * t773) * t794;
t713 = -t734 * mrSges(5,1) + t735 * mrSges(5,2);
t716 = -t744 * mrSges(5,2) + t734 * mrSges(5,3);
t668 = m(5) * t680 + t725 * mrSges(5,1) + t744 * t716 + (-t712 - t713) * t735 + (-mrSges(5,3) - mrSges(6,3)) * t695 + t791;
t790 = m(6) * t679 + t694 * mrSges(6,3) + t734 * t712;
t718 = t744 * mrSges(6,1) - t735 * mrSges(6,3);
t796 = -t744 * mrSges(5,1) + t735 * mrSges(5,3) - t718;
t670 = m(5) * t681 + t694 * mrSges(5,3) + t734 * t713 + t809 * t725 + t796 * t744 + t790;
t667 = -t771 * t668 + t775 * t670;
t731 = -t745 * mrSges(4,1) + t746 * mrSges(4,2);
t737 = t760 * mrSges(4,1) - t746 * mrSges(4,3);
t664 = m(4) * t690 - t749 * mrSges(4,2) + t727 * mrSges(4,3) + t745 * t731 - t760 * t737 + t667;
t689 = -t772 * t709 + t776 * t710;
t684 = -t749 * pkin(3) - t759 * pkin(9) + t746 * t732 - t689;
t682 = -t694 * pkin(4) - t733 * qJ(5) + t735 * t717 + qJDD(5) + t684;
t785 = -m(6) * t682 + t694 * mrSges(6,1) + t734 * t715;
t673 = -m(5) * t684 + t694 * mrSges(5,1) + t809 * t695 + t734 * t716 + t796 * t735 + t785;
t736 = -t760 * mrSges(4,2) + t745 * mrSges(4,3);
t672 = m(4) * t689 + t749 * mrSges(4,1) - t728 * mrSges(4,3) - t746 * t731 + t760 * t736 + t673;
t786 = t776 * t664 - t772 * t672;
t654 = m(3) * t730 - t764 * mrSges(3,2) - t757 * mrSges(3,3) - t765 * t750 + t754 * t788 + t786;
t657 = t772 * t664 + t776 * t672;
t741 = -t769 * t752 - t810;
t751 = -t765 * mrSges(3,2) + mrSges(3,3) * t788;
t656 = m(3) * t741 + t757 * mrSges(3,1) + t756 * mrSges(3,2) + (t750 * t773 - t751 * t777) * t794 + t657;
t666 = t775 * t668 + t771 * t670;
t781 = -m(4) * t708 + t727 * mrSges(4,1) - t728 * mrSges(4,2) + t745 * t736 - t746 * t737 - t666;
t661 = m(3) * t729 + t764 * mrSges(3,1) - t756 * mrSges(3,3) + t765 * t751 - t754 * t789 + t781;
t643 = t654 * t802 - t769 * t656 + t661 * t801;
t640 = m(2) * t761 + qJDD(1) * mrSges(2,1) - t779 * mrSges(2,2) + t643;
t649 = t777 * t654 - t773 * t661;
t647 = m(2) * t762 - t779 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t649;
t800 = t778 * t640 + t774 * t647;
t799 = -t806 * t734 - t807 * t735 - t813 * t744;
t642 = t654 * t804 + t770 * t656 + t661 * t803;
t787 = -t774 * t640 + t778 * t647;
t675 = t695 * mrSges(6,2) + t735 * t718 - t785;
t658 = -mrSges(5,1) * t684 + mrSges(5,3) * t681 - mrSges(6,1) * t682 + mrSges(6,3) * t679 - pkin(4) * t675 + qJ(5) * t790 + (-qJ(5) * t718 + t797) * t744 + t799 * t735 + (-qJ(5) * mrSges(6,2) + t806) * t725 + t808 * t695 + t814 * t694;
t665 = mrSges(5,2) * t684 + mrSges(6,2) * t682 - mrSges(5,3) * t680 - mrSges(6,3) * t677 - qJ(5) * t674 + t808 * t694 + t815 * t695 + t807 * t725 - t799 * t734 + t798 * t744;
t721 = Ifges(4,5) * t746 + Ifges(4,6) * t745 + Ifges(4,3) * t760;
t722 = Ifges(4,4) * t746 + Ifges(4,2) * t745 + Ifges(4,6) * t760;
t644 = mrSges(4,2) * t708 - mrSges(4,3) * t689 + Ifges(4,1) * t728 + Ifges(4,4) * t727 + Ifges(4,5) * t749 - pkin(9) * t666 - t771 * t658 + t775 * t665 + t745 * t721 - t760 * t722;
t723 = Ifges(4,1) * t746 + Ifges(4,4) * t745 + Ifges(4,5) * t760;
t650 = -mrSges(4,1) * t708 + mrSges(4,3) * t690 + Ifges(4,4) * t728 + Ifges(4,2) * t727 + Ifges(4,6) * t749 - pkin(3) * t666 - t746 * t721 + t760 * t723 - t812;
t739 = Ifges(3,6) * t765 + (Ifges(3,4) * t773 + Ifges(3,2) * t777) * t794;
t740 = Ifges(3,5) * t765 + (Ifges(3,1) * t773 + Ifges(3,4) * t777) * t794;
t634 = Ifges(3,5) * t756 - Ifges(3,6) * t757 + Ifges(3,3) * t764 + mrSges(3,1) * t729 - mrSges(3,2) * t730 + t772 * t644 + t776 * t650 + pkin(2) * t781 + pkin(8) * t786 + (t739 * t773 - t740 * t777) * t794;
t738 = Ifges(3,3) * t765 + (Ifges(3,5) * t773 + Ifges(3,6) * t777) * t794;
t636 = mrSges(3,2) * t741 - mrSges(3,3) * t729 + Ifges(3,1) * t756 - Ifges(3,4) * t757 + Ifges(3,5) * t764 - pkin(8) * t657 + t776 * t644 - t772 * t650 + t738 * t788 - t765 * t739;
t780 = mrSges(4,1) * t689 - mrSges(4,2) * t690 + Ifges(4,5) * t728 + Ifges(4,6) * t727 + Ifges(4,3) * t749 + pkin(3) * t673 + pkin(9) * t667 + t775 * t658 + t771 * t665 + t746 * t722 - t745 * t723;
t638 = -mrSges(3,1) * t741 + mrSges(3,3) * t730 + Ifges(3,4) * t756 - Ifges(3,2) * t757 + Ifges(3,6) * t764 - pkin(2) * t657 - t738 * t789 + t765 * t740 - t780;
t783 = mrSges(2,1) * t761 - mrSges(2,2) * t762 + Ifges(2,3) * qJDD(1) + pkin(1) * t643 + t770 * t634 + t636 * t804 + t638 * t803 + t649 * t811;
t632 = -mrSges(2,2) * g(3) - mrSges(2,3) * t761 + Ifges(2,5) * qJDD(1) - t779 * Ifges(2,6) + t777 * t636 - t773 * t638 + (-t642 * t769 - t643 * t770) * pkin(7);
t631 = mrSges(2,1) * g(3) + mrSges(2,3) * t762 + t779 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t642 - t769 * t634 + (pkin(7) * t649 + t636 * t773 + t638 * t777) * t770;
t1 = [-m(1) * g(1) + t787; -m(1) * g(2) + t800; (-m(1) - m(2)) * g(3) + t642; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t800 - t774 * t631 + t778 * t632; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t787 + t778 * t631 + t774 * t632; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t783; t783; t634; t780; t812; t675;];
tauJB = t1;
