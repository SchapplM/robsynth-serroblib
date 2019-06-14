% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 00:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:38:28
% EndTime: 2019-05-07 00:38:48
% DurationCPUTime: 17.78s
% Computational Cost: add. (286339->383), mult. (639533->478), div. (0->0), fcn. (478957->12), ass. (0->161)
t838 = -2 * qJD(3);
t837 = Ifges(3,1) + Ifges(4,2);
t830 = Ifges(3,4) + Ifges(4,6);
t829 = Ifges(3,5) - Ifges(4,4);
t836 = Ifges(3,2) + Ifges(4,3);
t828 = Ifges(3,6) - Ifges(4,5);
t835 = Ifges(3,3) + Ifges(4,1);
t790 = sin(qJ(2));
t785 = sin(pkin(6));
t817 = qJD(1) * t785;
t778 = t790 * t817;
t786 = cos(pkin(6));
t780 = qJD(1) * t786 + qJD(2);
t834 = (pkin(2) * t780 + t838) * t778;
t791 = sin(qJ(1));
t796 = cos(qJ(1));
t774 = t791 * g(1) - g(2) * t796;
t797 = qJD(1) ^ 2;
t758 = pkin(8) * t785 * t797 + qJDD(1) * pkin(1) + t774;
t775 = -g(1) * t796 - g(2) * t791;
t814 = qJDD(1) * t785;
t759 = -pkin(1) * t797 + pkin(8) * t814 + t775;
t795 = cos(qJ(2));
t824 = t786 * t790;
t826 = t785 * t790;
t723 = -g(3) * t826 + t758 * t824 + t795 * t759;
t760 = (-pkin(2) * t795 - qJ(3) * t790) * t817;
t777 = t780 ^ 2;
t779 = qJDD(1) * t786 + qJDD(2);
t816 = qJD(1) * t795;
t811 = t785 * t816;
t705 = t777 * pkin(2) - t779 * qJ(3) - t760 * t811 + t780 * t838 - t723;
t833 = -pkin(2) - pkin(9);
t832 = t786 * g(3);
t831 = mrSges(3,1) - mrSges(4,2);
t827 = t785 ^ 2 * t797;
t825 = t785 * t795;
t823 = t786 * t795;
t737 = -t785 * t758 - t832;
t754 = mrSges(3,1) * t780 - mrSges(3,3) * t778;
t755 = -mrSges(3,2) * t780 + mrSges(3,3) * t811;
t757 = mrSges(4,1) * t778 + mrSges(4,2) * t780;
t764 = (qJD(2) * t816 + qJDD(1) * t790) * t785;
t765 = -qJD(2) * t778 + t795 * t814;
t706 = -t765 * pkin(2) + (-t780 * t811 - t764) * qJ(3) + t737 + t834;
t756 = -mrSges(4,1) * t811 - mrSges(4,3) * t780;
t763 = pkin(3) * t778 - pkin(9) * t780;
t813 = t795 ^ 2 * t827;
t693 = -pkin(3) * t813 - t832 - t764 * qJ(3) + t833 * t765 + (-t758 + (-qJ(3) * t780 * t795 - t763 * t790) * qJD(1)) * t785 + t834;
t818 = g(3) * t825 + t790 * t759;
t805 = -t777 * qJ(3) + t760 * t778 + qJDD(3) + t818;
t696 = t764 * pkin(3) + t833 * t779 + (-pkin(3) * t780 * t817 - pkin(9) * t790 * t827 - t758 * t786) * t795 + t805;
t789 = sin(qJ(4));
t794 = cos(qJ(4));
t675 = -t789 * t693 + t794 * t696;
t746 = -t780 * t789 - t794 * t811;
t721 = qJD(4) * t746 - t765 * t789 + t779 * t794;
t747 = t780 * t794 - t789 * t811;
t753 = qJDD(4) + t764;
t770 = t778 + qJD(4);
t672 = (t746 * t770 - t721) * pkin(10) + (t746 * t747 + t753) * pkin(4) + t675;
t676 = t794 * t693 + t789 * t696;
t720 = -qJD(4) * t747 - t765 * t794 - t779 * t789;
t730 = pkin(4) * t770 - pkin(10) * t747;
t745 = t746 ^ 2;
t674 = -pkin(4) * t745 + pkin(10) * t720 - t730 * t770 + t676;
t788 = sin(qJ(5));
t793 = cos(qJ(5));
t669 = t788 * t672 + t793 * t674;
t726 = t746 * t788 + t747 * t793;
t688 = -qJD(5) * t726 + t720 * t793 - t721 * t788;
t725 = t746 * t793 - t747 * t788;
t707 = -mrSges(6,1) * t725 + mrSges(6,2) * t726;
t768 = qJD(5) + t770;
t713 = mrSges(6,1) * t768 - mrSges(6,3) * t726;
t750 = qJDD(5) + t753;
t708 = -pkin(5) * t725 - pkin(11) * t726;
t767 = t768 ^ 2;
t667 = -pkin(5) * t767 + pkin(11) * t750 + t708 * t725 + t669;
t692 = t765 * pkin(3) - pkin(9) * t813 + t780 * t763 - t705;
t678 = -t720 * pkin(4) - t745 * pkin(10) + t747 * t730 + t692;
t689 = qJD(5) * t725 + t720 * t788 + t721 * t793;
t670 = t678 + (-t725 * t768 - t689) * pkin(11) + (t726 * t768 - t688) * pkin(5);
t787 = sin(qJ(6));
t792 = cos(qJ(6));
t664 = -t667 * t787 + t670 * t792;
t710 = -t726 * t787 + t768 * t792;
t681 = qJD(6) * t710 + t689 * t792 + t750 * t787;
t687 = qJDD(6) - t688;
t711 = t726 * t792 + t768 * t787;
t697 = -mrSges(7,1) * t710 + mrSges(7,2) * t711;
t724 = qJD(6) - t725;
t698 = -mrSges(7,2) * t724 + mrSges(7,3) * t710;
t662 = m(7) * t664 + mrSges(7,1) * t687 - mrSges(7,3) * t681 - t697 * t711 + t698 * t724;
t665 = t667 * t792 + t670 * t787;
t680 = -qJD(6) * t711 - t689 * t787 + t750 * t792;
t699 = mrSges(7,1) * t724 - mrSges(7,3) * t711;
t663 = m(7) * t665 - mrSges(7,2) * t687 + mrSges(7,3) * t680 + t697 * t710 - t699 * t724;
t807 = -t662 * t787 + t792 * t663;
t653 = m(6) * t669 - mrSges(6,2) * t750 + mrSges(6,3) * t688 + t707 * t725 - t713 * t768 + t807;
t668 = t672 * t793 - t674 * t788;
t712 = -mrSges(6,2) * t768 + mrSges(6,3) * t725;
t666 = -pkin(5) * t750 - pkin(11) * t767 + t708 * t726 - t668;
t802 = -m(7) * t666 + t680 * mrSges(7,1) - mrSges(7,2) * t681 + t710 * t698 - t699 * t711;
t658 = m(6) * t668 + mrSges(6,1) * t750 - mrSges(6,3) * t689 - t707 * t726 + t712 * t768 + t802;
t646 = t788 * t653 + t793 * t658;
t727 = -mrSges(5,1) * t746 + mrSges(5,2) * t747;
t728 = -mrSges(5,2) * t770 + mrSges(5,3) * t746;
t644 = m(5) * t675 + mrSges(5,1) * t753 - mrSges(5,3) * t721 - t727 * t747 + t728 * t770 + t646;
t729 = mrSges(5,1) * t770 - mrSges(5,3) * t747;
t808 = t793 * t653 - t658 * t788;
t645 = m(5) * t676 - mrSges(5,2) * t753 + mrSges(5,3) * t720 + t727 * t746 - t729 * t770 + t808;
t809 = -t789 * t644 + t794 * t645;
t806 = m(4) * t706 - t764 * mrSges(4,3) + t756 * t811 + t809;
t637 = m(3) * t737 + t764 * mrSges(3,2) - t831 * t765 + (-t755 * t795 + (t754 - t757) * t790) * t817 + t806;
t812 = t758 * t823;
t722 = t812 - t818;
t761 = (mrSges(4,2) * t795 - mrSges(4,3) * t790) * t817;
t762 = (-mrSges(3,1) * t795 + mrSges(3,2) * t790) * t817;
t640 = t794 * t644 + t789 * t645;
t709 = -t779 * pkin(2) + t805 - t812;
t803 = -m(4) * t709 - t764 * mrSges(4,1) - t640;
t638 = m(3) * t722 - t764 * mrSges(3,3) + (t755 - t756) * t780 + t831 * t779 + (-t761 - t762) * t778 + t803;
t654 = t792 * t662 + t787 * t663;
t800 = m(6) * t678 - t688 * mrSges(6,1) + t689 * mrSges(6,2) - t725 * t712 + t726 * t713 + t654;
t799 = -m(5) * t692 + t720 * mrSges(5,1) - t721 * mrSges(5,2) + t746 * t728 - t747 * t729 - t800;
t798 = -m(4) * t705 + t779 * mrSges(4,3) + t780 * t757 + t761 * t811 - t799;
t650 = t798 - t779 * mrSges(3,2) - t780 * t754 + t762 * t811 + m(3) * t723 + (mrSges(3,3) + mrSges(4,1)) * t765;
t627 = -t637 * t785 + t638 * t823 + t650 * t824;
t625 = m(2) * t774 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t797 + t627;
t633 = -t638 * t790 + t795 * t650;
t632 = m(2) * t775 - mrSges(2,1) * t797 - qJDD(1) * mrSges(2,2) + t633;
t822 = t796 * t625 + t791 * t632;
t821 = (t829 * t790 + t828 * t795) * t817 + t835 * t780;
t820 = (-t830 * t790 - t836 * t795) * t817 - t828 * t780;
t819 = (t837 * t790 + t830 * t795) * t817 + t829 * t780;
t626 = t786 * t637 + t638 * t825 + t650 * t826;
t810 = -t625 * t791 + t796 * t632;
t682 = Ifges(7,5) * t711 + Ifges(7,6) * t710 + Ifges(7,3) * t724;
t684 = Ifges(7,1) * t711 + Ifges(7,4) * t710 + Ifges(7,5) * t724;
t655 = -mrSges(7,1) * t666 + mrSges(7,3) * t665 + Ifges(7,4) * t681 + Ifges(7,2) * t680 + Ifges(7,6) * t687 - t682 * t711 + t684 * t724;
t683 = Ifges(7,4) * t711 + Ifges(7,2) * t710 + Ifges(7,6) * t724;
t656 = mrSges(7,2) * t666 - mrSges(7,3) * t664 + Ifges(7,1) * t681 + Ifges(7,4) * t680 + Ifges(7,5) * t687 + t682 * t710 - t683 * t724;
t700 = Ifges(6,5) * t726 + Ifges(6,6) * t725 + Ifges(6,3) * t768;
t701 = Ifges(6,4) * t726 + Ifges(6,2) * t725 + Ifges(6,6) * t768;
t641 = mrSges(6,2) * t678 - mrSges(6,3) * t668 + Ifges(6,1) * t689 + Ifges(6,4) * t688 + Ifges(6,5) * t750 - pkin(11) * t654 - t655 * t787 + t656 * t792 + t700 * t725 - t701 * t768;
t702 = Ifges(6,1) * t726 + Ifges(6,4) * t725 + Ifges(6,5) * t768;
t642 = -mrSges(6,1) * t678 - mrSges(7,1) * t664 + mrSges(7,2) * t665 + mrSges(6,3) * t669 + Ifges(6,4) * t689 - Ifges(7,5) * t681 + Ifges(6,2) * t688 + Ifges(6,6) * t750 - Ifges(7,6) * t680 - Ifges(7,3) * t687 - pkin(5) * t654 - t683 * t711 + t684 * t710 - t700 * t726 + t702 * t768;
t714 = Ifges(5,5) * t747 + Ifges(5,6) * t746 + Ifges(5,3) * t770;
t716 = Ifges(5,1) * t747 + Ifges(5,4) * t746 + Ifges(5,5) * t770;
t628 = -mrSges(5,1) * t692 + mrSges(5,3) * t676 + Ifges(5,4) * t721 + Ifges(5,2) * t720 + Ifges(5,6) * t753 - pkin(4) * t800 + pkin(10) * t808 + t788 * t641 + t793 * t642 - t747 * t714 + t770 * t716;
t715 = Ifges(5,4) * t747 + Ifges(5,2) * t746 + Ifges(5,6) * t770;
t629 = mrSges(5,2) * t692 - mrSges(5,3) * t675 + Ifges(5,1) * t721 + Ifges(5,4) * t720 + Ifges(5,5) * t753 - pkin(10) * t646 + t641 * t793 - t642 * t788 + t714 * t746 - t715 * t770;
t639 = t765 * mrSges(4,2) - t757 * t778 + t806;
t622 = -mrSges(3,1) * t737 - mrSges(4,1) * t705 + mrSges(4,2) * t706 + mrSges(3,3) * t723 - pkin(2) * t639 - pkin(3) * t799 - pkin(9) * t809 - t794 * t628 - t789 * t629 + t830 * t764 + t836 * t765 - t821 * t778 + t828 * t779 + t819 * t780;
t623 = t837 * t764 + t820 * t780 + t821 * t811 + t829 * t779 + t830 * t765 + pkin(5) * t802 + pkin(4) * t646 + pkin(3) * t640 - qJ(3) * t639 + pkin(11) * t807 + t787 * t656 + t792 * t655 + Ifges(6,3) * t750 + Ifges(5,3) * t753 - t746 * t716 + t747 * t715 + t726 * t701 + mrSges(3,2) * t737 + Ifges(5,6) * t720 + Ifges(5,5) * t721 - mrSges(3,3) * t722 - t725 * t702 - mrSges(4,3) * t706 + mrSges(4,1) * t709 + Ifges(6,6) * t688 + Ifges(6,5) * t689 + mrSges(5,1) * t675 - mrSges(5,2) * t676 + mrSges(6,1) * t668 - mrSges(6,2) * t669;
t804 = pkin(8) * t633 + t622 * t795 + t623 * t790;
t621 = mrSges(3,1) * t722 - mrSges(3,2) * t723 + mrSges(4,2) * t709 - mrSges(4,3) * t705 + t794 * t629 - t789 * t628 - pkin(9) * t640 + pkin(2) * (-t780 * t756 + t803) + qJ(3) * t798 + (-mrSges(4,2) * pkin(2) + t835) * t779 + (mrSges(4,1) * qJ(3) + t828) * t765 + t829 * t764 + (-t819 * t795 + (-pkin(2) * t761 - t820) * t790) * t817;
t620 = -mrSges(2,2) * g(3) - mrSges(2,3) * t774 + Ifges(2,5) * qJDD(1) - t797 * Ifges(2,6) - t790 * t622 + t795 * t623 + (-t626 * t785 - t627 * t786) * pkin(8);
t619 = mrSges(2,1) * g(3) + mrSges(2,3) * t775 + t797 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t626 - t785 * t621 + t786 * t804;
t1 = [-m(1) * g(1) + t810; -m(1) * g(2) + t822; (-m(1) - m(2)) * g(3) + t626; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t822 - t791 * t619 + t796 * t620; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t810 + t796 * t619 + t791 * t620; -mrSges(1,1) * g(2) + mrSges(2,1) * t774 + mrSges(1,2) * g(1) - mrSges(2,2) * t775 + Ifges(2,3) * qJDD(1) + pkin(1) * t627 + t786 * t621 + t785 * t804;];
tauB  = t1;
