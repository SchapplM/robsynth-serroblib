% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR13
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
% Datum: 2019-05-07 01:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 01:09:53
% EndTime: 2019-05-07 01:10:15
% DurationCPUTime: 18.50s
% Computational Cost: add. (292683->383), mult. (648924->478), div. (0->0), fcn. (481968->12), ass. (0->161)
t852 = -2 * qJD(3);
t851 = Ifges(3,1) + Ifges(4,2);
t844 = Ifges(3,4) + Ifges(4,6);
t843 = Ifges(3,5) - Ifges(4,4);
t850 = Ifges(3,2) + Ifges(4,3);
t842 = Ifges(3,6) - Ifges(4,5);
t849 = Ifges(3,3) + Ifges(4,1);
t799 = cos(pkin(6));
t793 = qJD(1) * t799 + qJD(2);
t803 = sin(qJ(2));
t798 = sin(pkin(6));
t831 = qJD(1) * t798;
t825 = t803 * t831;
t848 = (pkin(2) * t793 + t852) * t825;
t804 = sin(qJ(1));
t809 = cos(qJ(1));
t788 = g(1) * t804 - g(2) * t809;
t810 = qJD(1) ^ 2;
t770 = pkin(8) * t798 * t810 + qJDD(1) * pkin(1) + t788;
t789 = -g(1) * t809 - g(2) * t804;
t828 = qJDD(1) * t798;
t771 = -pkin(1) * t810 + pkin(8) * t828 + t789;
t808 = cos(qJ(2));
t838 = t799 * t803;
t840 = t798 * t803;
t734 = -g(3) * t840 + t770 * t838 + t771 * t808;
t772 = (-pkin(2) * t808 - qJ(3) * t803) * t831;
t791 = t793 ^ 2;
t792 = qJDD(1) * t799 + qJDD(2);
t830 = qJD(1) * t808;
t824 = t798 * t830;
t710 = t791 * pkin(2) - qJ(3) * t792 - t772 * t824 + t793 * t852 - t734;
t847 = -pkin(2) - pkin(9);
t846 = t799 * g(3);
t845 = mrSges(3,1) - mrSges(4,2);
t841 = t798 ^ 2 * t810;
t839 = t798 * t808;
t837 = t799 * t808;
t748 = -t798 * t770 - t846;
t766 = mrSges(3,1) * t793 - mrSges(3,3) * t825;
t767 = -mrSges(3,2) * t793 + mrSges(3,3) * t824;
t769 = mrSges(4,1) * t825 + mrSges(4,2) * t793;
t776 = (qJD(2) * t830 + qJDD(1) * t803) * t798;
t777 = -qJD(2) * t825 + t808 * t828;
t711 = -t777 * pkin(2) + (-t793 * t824 - t776) * qJ(3) + t748 + t848;
t768 = -mrSges(4,1) * t824 - mrSges(4,3) * t793;
t775 = pkin(3) * t825 - pkin(9) * t793;
t827 = t808 ^ 2 * t841;
t701 = -pkin(3) * t827 - t846 - t776 * qJ(3) + t847 * t777 + (-t770 + (-qJ(3) * t793 * t808 - t775 * t803) * qJD(1)) * t798 + t848;
t832 = g(3) * t839 + t771 * t803;
t818 = -t791 * qJ(3) + t772 * t825 + qJDD(3) + t832;
t703 = t776 * pkin(3) + t847 * t792 + (-pkin(3) * t793 * t831 - pkin(9) * t803 * t841 - t770 * t799) * t808 + t818;
t802 = sin(qJ(4));
t807 = cos(qJ(4));
t691 = t701 * t807 + t703 * t802;
t760 = t793 * t807 - t802 * t824;
t731 = -qJD(4) * t760 - t777 * t807 - t792 * t802;
t759 = -t793 * t802 - t807 * t824;
t735 = -mrSges(5,1) * t759 + mrSges(5,2) * t760;
t784 = qJD(4) + t825;
t741 = mrSges(5,1) * t784 - mrSges(5,3) * t760;
t765 = qJDD(4) + t776;
t736 = -pkin(4) * t759 - pkin(10) * t760;
t781 = t784 ^ 2;
t683 = -pkin(4) * t781 + pkin(10) * t765 + t736 * t759 + t691;
t700 = t777 * pkin(3) - pkin(9) * t827 + t775 * t793 - t710;
t732 = qJD(4) * t759 - t777 * t802 + t792 * t807;
t687 = (-t759 * t784 - t732) * pkin(10) + (t760 * t784 - t731) * pkin(4) + t700;
t801 = sin(qJ(5));
t806 = cos(qJ(5));
t678 = -t801 * t683 + t687 * t806;
t738 = -t760 * t801 + t784 * t806;
t706 = qJD(5) * t738 + t732 * t806 + t765 * t801;
t729 = qJDD(5) - t731;
t739 = t760 * t806 + t784 * t801;
t758 = qJD(5) - t759;
t676 = (t738 * t758 - t706) * pkin(11) + (t738 * t739 + t729) * pkin(5) + t678;
t679 = t683 * t806 + t687 * t801;
t705 = -qJD(5) * t739 - t732 * t801 + t765 * t806;
t723 = pkin(5) * t758 - pkin(11) * t739;
t737 = t738 ^ 2;
t677 = -pkin(5) * t737 + pkin(11) * t705 - t723 * t758 + t679;
t800 = sin(qJ(6));
t805 = cos(qJ(6));
t674 = t676 * t805 - t677 * t800;
t717 = t738 * t805 - t739 * t800;
t689 = qJD(6) * t717 + t705 * t800 + t706 * t805;
t718 = t738 * t800 + t739 * t805;
t696 = -mrSges(7,1) * t717 + mrSges(7,2) * t718;
t756 = qJD(6) + t758;
t707 = -mrSges(7,2) * t756 + mrSges(7,3) * t717;
t724 = qJDD(6) + t729;
t672 = m(7) * t674 + mrSges(7,1) * t724 - mrSges(7,3) * t689 - t696 * t718 + t707 * t756;
t675 = t676 * t800 + t677 * t805;
t688 = -qJD(6) * t718 + t705 * t805 - t706 * t800;
t708 = mrSges(7,1) * t756 - mrSges(7,3) * t718;
t673 = m(7) * t675 - mrSges(7,2) * t724 + mrSges(7,3) * t688 + t696 * t717 - t708 * t756;
t665 = t672 * t805 + t673 * t800;
t719 = -mrSges(6,1) * t738 + mrSges(6,2) * t739;
t721 = -mrSges(6,2) * t758 + mrSges(6,3) * t738;
t663 = m(6) * t678 + mrSges(6,1) * t729 - mrSges(6,3) * t706 - t719 * t739 + t721 * t758 + t665;
t722 = mrSges(6,1) * t758 - mrSges(6,3) * t739;
t820 = -t672 * t800 + t673 * t805;
t664 = m(6) * t679 - mrSges(6,2) * t729 + mrSges(6,3) * t705 + t719 * t738 - t722 * t758 + t820;
t821 = -t663 * t801 + t664 * t806;
t658 = m(5) * t691 - mrSges(5,2) * t765 + mrSges(5,3) * t731 + t735 * t759 - t741 * t784 + t821;
t690 = -t701 * t802 + t703 * t807;
t740 = -mrSges(5,2) * t784 + mrSges(5,3) * t759;
t682 = -pkin(4) * t765 - pkin(10) * t781 + t736 * t760 - t690;
t680 = -pkin(5) * t705 - pkin(11) * t737 + t723 * t739 + t682;
t815 = m(7) * t680 - mrSges(7,1) * t688 + mrSges(7,2) * t689 - t707 * t717 + t708 * t718;
t811 = -m(6) * t682 + mrSges(6,1) * t705 - mrSges(6,2) * t706 + t721 * t738 - t722 * t739 - t815;
t668 = m(5) * t690 + mrSges(5,1) * t765 - mrSges(5,3) * t732 - t735 * t760 + t740 * t784 + t811;
t822 = t658 * t807 - t802 * t668;
t819 = m(4) * t711 - mrSges(4,3) * t776 + t768 * t824 + t822;
t647 = m(3) * t748 + t776 * mrSges(3,2) - t845 * t777 + (-t767 * t808 + (t766 - t769) * t803) * t831 + t819;
t826 = t770 * t837;
t733 = t826 - t832;
t773 = (mrSges(4,2) * t808 - mrSges(4,3) * t803) * t831;
t774 = (-mrSges(3,1) * t808 + mrSges(3,2) * t803) * t831;
t650 = t802 * t658 + t807 * t668;
t716 = -t792 * pkin(2) + t818 - t826;
t816 = -m(4) * t716 - mrSges(4,1) * t776 - t650;
t648 = m(3) * t733 - t776 * mrSges(3,3) + (t767 - t768) * t793 + t845 * t792 + (-t773 - t774) * t825 + t816;
t659 = t663 * t806 + t664 * t801;
t813 = -m(5) * t700 + t731 * mrSges(5,1) - mrSges(5,2) * t732 + t759 * t740 - t741 * t760 - t659;
t812 = -m(4) * t710 + mrSges(4,3) * t792 + t769 * t793 + t773 * t824 - t813;
t656 = (mrSges(3,3) + mrSges(4,1)) * t777 + t812 + t774 * t824 - t793 * t766 - t792 * mrSges(3,2) + m(3) * t734;
t637 = -t647 * t798 + t648 * t837 + t656 * t838;
t635 = m(2) * t788 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t810 + t637;
t642 = -t648 * t803 + t656 * t808;
t641 = m(2) * t789 - mrSges(2,1) * t810 - qJDD(1) * mrSges(2,2) + t642;
t836 = t635 * t809 + t641 * t804;
t835 = (t803 * t843 + t808 * t842) * t831 + t849 * t793;
t834 = (-t803 * t844 - t808 * t850) * t831 - t842 * t793;
t833 = (t803 * t851 + t808 * t844) * t831 + t843 * t793;
t636 = t647 * t799 + t648 * t839 + t656 * t840;
t823 = -t635 * t804 + t641 * t809;
t692 = Ifges(7,5) * t718 + Ifges(7,6) * t717 + Ifges(7,3) * t756;
t694 = Ifges(7,1) * t718 + Ifges(7,4) * t717 + Ifges(7,5) * t756;
t666 = -mrSges(7,1) * t680 + mrSges(7,3) * t675 + Ifges(7,4) * t689 + Ifges(7,2) * t688 + Ifges(7,6) * t724 - t692 * t718 + t694 * t756;
t693 = Ifges(7,4) * t718 + Ifges(7,2) * t717 + Ifges(7,6) * t756;
t667 = mrSges(7,2) * t680 - mrSges(7,3) * t674 + Ifges(7,1) * t689 + Ifges(7,4) * t688 + Ifges(7,5) * t724 + t692 * t717 - t693 * t756;
t712 = Ifges(6,5) * t739 + Ifges(6,6) * t738 + Ifges(6,3) * t758;
t714 = Ifges(6,1) * t739 + Ifges(6,4) * t738 + Ifges(6,5) * t758;
t651 = -mrSges(6,1) * t682 + mrSges(6,3) * t679 + Ifges(6,4) * t706 + Ifges(6,2) * t705 + Ifges(6,6) * t729 - pkin(5) * t815 + pkin(11) * t820 + t805 * t666 + t800 * t667 - t739 * t712 + t758 * t714;
t713 = Ifges(6,4) * t739 + Ifges(6,2) * t738 + Ifges(6,6) * t758;
t652 = mrSges(6,2) * t682 - mrSges(6,3) * t678 + Ifges(6,1) * t706 + Ifges(6,4) * t705 + Ifges(6,5) * t729 - pkin(11) * t665 - t666 * t800 + t667 * t805 + t712 * t738 - t713 * t758;
t725 = Ifges(5,5) * t760 + Ifges(5,6) * t759 + Ifges(5,3) * t784;
t726 = Ifges(5,4) * t760 + Ifges(5,2) * t759 + Ifges(5,6) * t784;
t638 = mrSges(5,2) * t700 - mrSges(5,3) * t690 + Ifges(5,1) * t732 + Ifges(5,4) * t731 + Ifges(5,5) * t765 - pkin(10) * t659 - t651 * t801 + t652 * t806 + t725 * t759 - t726 * t784;
t727 = Ifges(5,1) * t760 + Ifges(5,4) * t759 + Ifges(5,5) * t784;
t643 = Ifges(5,4) * t732 + Ifges(5,2) * t731 + Ifges(5,6) * t765 - t760 * t725 + t784 * t727 - mrSges(5,1) * t700 + mrSges(5,3) * t691 - Ifges(6,5) * t706 - Ifges(6,6) * t705 - Ifges(6,3) * t729 - t739 * t713 + t738 * t714 - mrSges(6,1) * t678 + mrSges(6,2) * t679 - Ifges(7,5) * t689 - Ifges(7,6) * t688 - Ifges(7,3) * t724 - t718 * t693 + t717 * t694 - mrSges(7,1) * t674 + mrSges(7,2) * t675 - pkin(5) * t665 - pkin(4) * t659;
t649 = t777 * mrSges(4,2) - t769 * t825 + t819;
t632 = -mrSges(3,1) * t748 - mrSges(4,1) * t710 + mrSges(4,2) * t711 + mrSges(3,3) * t734 - pkin(2) * t649 - pkin(3) * t813 - pkin(9) * t822 - t802 * t638 - t807 * t643 + t776 * t844 + t777 * t850 + t792 * t842 + t793 * t833 - t825 * t835;
t633 = t806 * t651 + pkin(10) * t821 + t801 * t652 + Ifges(5,3) * t765 - t759 * t727 + t760 * t726 + pkin(4) * t811 + mrSges(3,2) * t748 + Ifges(5,6) * t731 + Ifges(5,5) * t732 - mrSges(3,3) * t733 - mrSges(4,3) * t711 + mrSges(4,1) * t716 + mrSges(5,1) * t690 - mrSges(5,2) * t691 + pkin(3) * t650 - qJ(3) * t649 + t834 * t793 + t843 * t792 + t844 * t777 + t851 * t776 + t835 * t824;
t817 = pkin(8) * t642 + t632 * t808 + t633 * t803;
t631 = mrSges(3,1) * t733 - mrSges(3,2) * t734 + mrSges(4,2) * t716 - mrSges(4,3) * t710 + t807 * t638 - t802 * t643 - pkin(9) * t650 + pkin(2) * (-t793 * t768 + t816) + qJ(3) * t812 + (-mrSges(4,2) * pkin(2) + t849) * t792 + (mrSges(4,1) * qJ(3) + t842) * t777 + t843 * t776 + (-t833 * t808 + (-pkin(2) * t773 - t834) * t803) * t831;
t630 = -mrSges(2,2) * g(3) - mrSges(2,3) * t788 + Ifges(2,5) * qJDD(1) - t810 * Ifges(2,6) - t803 * t632 + t808 * t633 + (-t636 * t798 - t637 * t799) * pkin(8);
t629 = mrSges(2,1) * g(3) + mrSges(2,3) * t789 + t810 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t636 - t798 * t631 + t799 * t817;
t1 = [-m(1) * g(1) + t823; -m(1) * g(2) + t836; (-m(1) - m(2)) * g(3) + t636; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t836 - t629 * t804 + t630 * t809; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t823 + t809 * t629 + t804 * t630; -mrSges(1,1) * g(2) + mrSges(2,1) * t788 + mrSges(1,2) * g(1) - mrSges(2,2) * t789 + Ifges(2,3) * qJDD(1) + pkin(1) * t637 + t799 * t631 + t798 * t817;];
tauB  = t1;
