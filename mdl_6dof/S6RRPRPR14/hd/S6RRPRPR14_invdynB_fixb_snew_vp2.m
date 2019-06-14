% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-05-06 17:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:55:21
% EndTime: 2019-05-06 16:55:33
% DurationCPUTime: 7.53s
% Computational Cost: add. (99011->367), mult. (222613->437), div. (0->0), fcn. (157211->10), ass. (0->159)
t868 = -2 * qJD(3);
t867 = Ifges(3,1) + Ifges(4,2);
t866 = Ifges(5,1) + Ifges(6,2);
t865 = -Ifges(6,1) - Ifges(5,3);
t853 = Ifges(3,4) + Ifges(4,6);
t852 = Ifges(5,4) + Ifges(6,6);
t851 = Ifges(3,5) - Ifges(4,4);
t850 = Ifges(5,5) - Ifges(6,4);
t864 = Ifges(3,2) + Ifges(4,3);
t863 = -Ifges(5,2) - Ifges(6,3);
t849 = Ifges(3,6) - Ifges(4,5);
t848 = Ifges(5,6) - Ifges(6,5);
t862 = Ifges(3,3) + Ifges(4,1);
t797 = cos(pkin(6));
t791 = qJD(1) * t797 + qJD(2);
t800 = sin(qJ(2));
t796 = sin(pkin(6));
t830 = qJD(1) * t796;
t824 = t800 * t830;
t861 = (pkin(2) * t791 + t868) * t824;
t801 = sin(qJ(1));
t805 = cos(qJ(1));
t786 = t801 * g(1) - g(2) * t805;
t806 = qJD(1) ^ 2;
t770 = pkin(8) * t796 * t806 + qJDD(1) * pkin(1) + t786;
t787 = -g(1) * t805 - g(2) * t801;
t827 = qJDD(1) * t796;
t771 = -pkin(1) * t806 + pkin(8) * t827 + t787;
t804 = cos(qJ(2));
t841 = t797 * t800;
t843 = t796 * t800;
t727 = -g(3) * t843 + t770 * t841 + t804 * t771;
t772 = (-pkin(2) * t804 - qJ(3) * t800) * t830;
t789 = t791 ^ 2;
t790 = qJDD(1) * t797 + qJDD(2);
t829 = qJD(1) * t804;
t823 = t796 * t829;
t700 = pkin(2) * t789 - t790 * qJ(3) - t772 * t823 + t791 * t868 - t727;
t799 = sin(qJ(4));
t803 = cos(qJ(4));
t758 = t791 * t803 - t799 * t823;
t777 = -qJD(2) * t824 + t804 * t827;
t724 = qJD(4) * t758 + t803 * t777 + t790 * t799;
t757 = t791 * t799 + t803 * t823;
t782 = qJD(4) + t824;
t734 = mrSges(6,1) * t757 - mrSges(6,3) * t782;
t860 = -t724 * mrSges(6,2) - t757 * t734;
t736 = -mrSges(5,2) * t782 - mrSges(5,3) * t757;
t835 = t734 - t736;
t854 = mrSges(5,1) - mrSges(6,2);
t859 = t854 * t724 - t835 * t757;
t858 = -2 * qJD(5);
t857 = -pkin(2) - pkin(9);
t856 = t797 * g(3);
t855 = mrSges(3,1) - mrSges(4,2);
t845 = t757 * t782;
t844 = t796 ^ 2 * t806;
t842 = t796 * t804;
t840 = t797 * t804;
t746 = -t796 * t770 - t856;
t766 = mrSges(3,1) * t791 - mrSges(3,3) * t824;
t767 = -mrSges(3,2) * t791 + mrSges(3,3) * t823;
t769 = mrSges(4,1) * t824 + mrSges(4,2) * t791;
t776 = (qJD(2) * t829 + qJDD(1) * t800) * t796;
t701 = -t777 * pkin(2) + (-t791 * t823 - t776) * qJ(3) + t746 + t861;
t768 = -mrSges(4,1) * t823 - mrSges(4,3) * t791;
t775 = pkin(3) * t824 - pkin(9) * t791;
t825 = t804 ^ 2 * t844;
t693 = -pkin(3) * t825 - t856 - t776 * qJ(3) + t857 * t777 + (-t770 + (-qJ(3) * t791 * t804 - t775 * t800) * qJD(1)) * t796 + t861;
t831 = g(3) * t842 + t800 * t771;
t817 = -qJ(3) * t789 + t772 * t824 + qJDD(3) + t831;
t695 = pkin(3) * t776 + t857 * t790 + (-pkin(3) * t791 * t830 - pkin(9) * t800 * t844 - t770 * t797) * t804 + t817;
t687 = -t799 * t693 + t695 * t803;
t725 = -qJD(4) * t757 - t777 * t799 + t790 * t803;
t729 = mrSges(5,1) * t757 + mrSges(5,2) * t758;
t765 = qJDD(4) + t776;
t728 = pkin(4) * t757 - qJ(5) * t758;
t779 = t782 ^ 2;
t684 = -pkin(4) * t765 - qJ(5) * t779 + t758 * t728 + qJDD(5) - t687;
t679 = (t757 * t758 - t765) * pkin(10) + (t725 + t845) * pkin(5) + t684;
t738 = pkin(5) * t758 - pkin(10) * t782;
t756 = t757 ^ 2;
t692 = pkin(3) * t777 - pkin(9) * t825 + t791 * t775 - t700;
t807 = (-t725 + t845) * qJ(5) + t692 + (t782 * pkin(4) + t858) * t758;
t682 = t807 + (pkin(4) + pkin(10)) * t724 - pkin(5) * t756 - t738 * t758;
t798 = sin(qJ(6));
t802 = cos(qJ(6));
t677 = t679 * t802 - t682 * t798;
t732 = t757 * t802 - t782 * t798;
t698 = qJD(6) * t732 + t724 * t798 + t765 * t802;
t733 = t757 * t798 + t782 * t802;
t707 = -mrSges(7,1) * t732 + mrSges(7,2) * t733;
t755 = qJD(6) + t758;
t710 = -mrSges(7,2) * t755 + mrSges(7,3) * t732;
t721 = qJDD(6) + t725;
t675 = m(7) * t677 + mrSges(7,1) * t721 - mrSges(7,3) * t698 - t707 * t733 + t710 * t755;
t678 = t679 * t798 + t682 * t802;
t697 = -qJD(6) * t733 + t724 * t802 - t765 * t798;
t711 = mrSges(7,1) * t755 - mrSges(7,3) * t733;
t676 = m(7) * t678 - mrSges(7,2) * t721 + mrSges(7,3) * t697 + t707 * t732 - t711 * t755;
t669 = t675 * t802 + t676 * t798;
t730 = -mrSges(6,2) * t757 - mrSges(6,3) * t758;
t813 = -m(6) * t684 - t725 * mrSges(6,1) - t758 * t730 - t669;
t667 = m(5) * t687 - mrSges(5,3) * t725 - t729 * t758 + t854 * t765 - t835 * t782 + t813;
t688 = t803 * t693 + t799 * t695;
t737 = mrSges(5,1) * t782 - mrSges(5,3) * t758;
t812 = -t779 * pkin(4) + t765 * qJ(5) - t757 * t728 + t688;
t683 = t782 * t858 - t812;
t735 = mrSges(6,1) * t758 + mrSges(6,2) * t782;
t681 = -t724 * pkin(5) - t756 * pkin(10) + ((2 * qJD(5)) + t738) * t782 + t812;
t815 = -m(7) * t681 + mrSges(7,1) * t697 - t698 * mrSges(7,2) + t710 * t732 - t733 * t711;
t809 = -m(6) * t683 + t765 * mrSges(6,3) + t782 * t735 - t815;
t673 = m(5) * t688 - mrSges(5,2) * t765 - t737 * t782 + (-t729 - t730) * t757 + (-mrSges(5,3) - mrSges(6,1)) * t724 + t809;
t821 = -t667 * t799 + t803 * t673;
t819 = m(4) * t701 - t776 * mrSges(4,3) + t768 * t823 + t821;
t659 = m(3) * t746 + mrSges(3,2) * t776 - t855 * t777 + (-t767 * t804 + (t766 - t769) * t800) * t830 + t819;
t826 = t770 * t840;
t726 = t826 - t831;
t773 = (mrSges(4,2) * t804 - mrSges(4,3) * t800) * t830;
t774 = (-mrSges(3,1) * t804 + mrSges(3,2) * t800) * t830;
t662 = t667 * t803 + t673 * t799;
t706 = -pkin(2) * t790 + t817 - t826;
t814 = -m(4) * t706 - t776 * mrSges(4,1) - t662;
t660 = m(3) * t726 - mrSges(3,3) * t776 + (t767 - t768) * t791 + t855 * t790 + (-t773 - t774) * t824 + t814;
t686 = pkin(4) * t724 + t807;
t820 = -t798 * t675 + t802 * t676;
t818 = -m(6) * t686 + t725 * mrSges(6,3) + t758 * t735 - t820;
t810 = -m(5) * t692 - t725 * mrSges(5,2) - t758 * t737 + t818;
t808 = -m(4) * t700 + t790 * mrSges(4,3) + t791 * t769 + t773 * t823 - t810;
t666 = t808 + t774 * t823 + m(3) * t727 - t790 * mrSges(3,2) - t791 * t766 + (mrSges(3,3) + mrSges(4,1)) * t777 + t859;
t649 = -t659 * t796 + t660 * t840 + t666 * t841;
t647 = m(2) * t786 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t806 + t649;
t655 = -t660 * t800 + t804 * t666;
t654 = m(2) * t787 - mrSges(2,1) * t806 - qJDD(1) * mrSges(2,2) + t655;
t839 = t805 * t647 + t801 * t654;
t838 = t757 * t848 - t758 * t850 + t782 * t865;
t837 = t757 * t863 + t758 * t852 + t782 * t848;
t836 = -t757 * t852 + t758 * t866 + t782 * t850;
t834 = (t800 * t851 + t804 * t849) * t830 + t862 * t791;
t833 = (-t800 * t853 - t804 * t864) * t830 - t849 * t791;
t832 = (t800 * t867 + t804 * t853) * t830 + t851 * t791;
t648 = t797 * t659 + t660 * t842 + t666 * t843;
t822 = -t647 * t801 + t805 * t654;
t668 = -t818 + t860;
t702 = Ifges(7,5) * t733 + Ifges(7,6) * t732 + Ifges(7,3) * t755;
t704 = Ifges(7,1) * t733 + Ifges(7,4) * t732 + Ifges(7,5) * t755;
t670 = -mrSges(7,1) * t681 + mrSges(7,3) * t678 + Ifges(7,4) * t698 + Ifges(7,2) * t697 + Ifges(7,6) * t721 - t702 * t733 + t704 * t755;
t703 = Ifges(7,4) * t733 + Ifges(7,2) * t732 + Ifges(7,6) * t755;
t671 = mrSges(7,2) * t681 - mrSges(7,3) * t677 + Ifges(7,1) * t698 + Ifges(7,4) * t697 + Ifges(7,5) * t721 + t702 * t732 - t703 * t755;
t650 = -mrSges(5,1) * t692 - mrSges(6,1) * t683 + mrSges(6,2) * t686 + mrSges(5,3) * t688 - pkin(4) * t668 - pkin(5) * t815 - pkin(10) * t820 - t802 * t670 - t798 * t671 + t724 * t863 + t852 * t725 + t838 * t758 + t848 * t765 + t836 * t782;
t651 = mrSges(6,1) * t684 + mrSges(7,1) * t677 + mrSges(5,2) * t692 - mrSges(7,2) * t678 - mrSges(5,3) * t687 - mrSges(6,3) * t686 + Ifges(7,5) * t698 + Ifges(7,6) * t697 + Ifges(7,3) * t721 + pkin(5) * t669 - qJ(5) * t668 + t733 * t703 - t732 * t704 - t837 * t782 + t850 * t765 + t838 * t757 + t866 * t725 - t852 * t724;
t661 = mrSges(4,2) * t777 - t769 * t824 + t819;
t644 = -mrSges(3,1) * t746 + mrSges(3,3) * t727 - mrSges(4,1) * t700 + mrSges(4,2) * t701 - t799 * t651 - t803 * t650 - pkin(3) * (t810 - t859) - pkin(9) * t821 - pkin(2) * t661 + t832 * t791 + t849 * t790 + t864 * t777 + t853 * t776 - t834 * t824;
t645 = pkin(3) * t662 - qJ(3) * t661 - t798 * t670 + t802 * t671 + (-qJ(5) * t730 + t836) * t757 + t837 * t758 + t833 * t791 + t834 * t823 + qJ(5) * t809 + pkin(4) * (-t734 * t782 + t813) - mrSges(6,3) * t683 + mrSges(6,2) * t684 - mrSges(3,3) * t726 + mrSges(5,1) * t687 - mrSges(5,2) * t688 - pkin(10) * t669 + mrSges(3,2) * t746 + (-mrSges(6,2) * pkin(4) - t865) * t765 + t867 * t776 - mrSges(4,3) * t701 + mrSges(4,1) * t706 + (-mrSges(6,1) * qJ(5) - t848) * t724 + t850 * t725 + t851 * t790 + t853 * t777;
t816 = pkin(8) * t655 + t644 * t804 + t645 * t800;
t643 = mrSges(3,1) * t726 - mrSges(3,2) * t727 + mrSges(4,2) * t706 - mrSges(4,3) * t700 + t803 * t651 - t799 * t650 - pkin(9) * t662 + pkin(2) * (-t768 * t791 + t814) + qJ(3) * (t724 * mrSges(5,1) + t757 * t736 + t808 + t860) + (-mrSges(4,2) * pkin(2) + t862) * t790 + (mrSges(4,1) * qJ(3) + t849) * t777 + t851 * t776 + (-t832 * t804 + (-pkin(2) * t773 - t833) * t800) * t830;
t642 = -mrSges(2,2) * g(3) - mrSges(2,3) * t786 + Ifges(2,5) * qJDD(1) - t806 * Ifges(2,6) - t800 * t644 + t804 * t645 + (-t648 * t796 - t649 * t797) * pkin(8);
t641 = mrSges(2,1) * g(3) + mrSges(2,3) * t787 + t806 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t648 - t796 * t643 + t816 * t797;
t1 = [-m(1) * g(1) + t822; -m(1) * g(2) + t839; (-m(1) - m(2)) * g(3) + t648; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t839 - t801 * t641 + t805 * t642; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t822 + t805 * t641 + t801 * t642; -mrSges(1,1) * g(2) + mrSges(2,1) * t786 + mrSges(1,2) * g(1) - mrSges(2,2) * t787 + Ifges(2,3) * qJDD(1) + pkin(1) * t649 + t643 * t797 + t816 * t796;];
tauB  = t1;
