% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:38:02
% EndTime: 2019-05-05 17:38:10
% DurationCPUTime: 7.24s
% Computational Cost: add. (112603->319), mult. (231458->387), div. (0->0), fcn. (149127->10), ass. (0->128)
t855 = Ifges(6,1) + Ifges(7,1);
t848 = Ifges(6,4) - Ifges(7,5);
t847 = -Ifges(6,5) - Ifges(7,4);
t854 = Ifges(6,2) + Ifges(7,3);
t846 = Ifges(6,6) - Ifges(7,6);
t853 = -Ifges(6,3) - Ifges(7,2);
t818 = sin(qJ(1));
t820 = cos(qJ(1));
t801 = t818 * g(1) - t820 * g(2);
t792 = qJDD(1) * pkin(1) + t801;
t802 = -t820 * g(1) - t818 * g(2);
t822 = qJD(1) ^ 2;
t795 = -t822 * pkin(1) + t802;
t813 = sin(pkin(9));
t815 = cos(pkin(9));
t766 = t813 * t792 + t815 * t795;
t756 = -t822 * pkin(2) + qJDD(1) * pkin(7) + t766;
t811 = -g(3) + qJDD(2);
t817 = sin(qJ(3));
t819 = cos(qJ(3));
t747 = -t817 * t756 + t819 * t811;
t793 = (-pkin(3) * t819 - qJ(4) * t817) * qJD(1);
t821 = qJD(3) ^ 2;
t839 = qJD(1) * t817;
t740 = -qJDD(3) * pkin(3) - t821 * qJ(4) + t793 * t839 + qJDD(4) - t747;
t837 = qJD(1) * qJD(3);
t835 = t819 * t837;
t796 = t817 * qJDD(1) + t835;
t812 = sin(pkin(10));
t814 = cos(pkin(10));
t771 = t814 * qJDD(3) - t812 * t796;
t789 = t812 * qJD(3) + t814 * t839;
t838 = t819 * qJD(1);
t773 = -pkin(4) * t838 - t789 * pkin(8);
t788 = t814 * qJD(3) - t812 * t839;
t787 = t788 ^ 2;
t723 = -t771 * pkin(4) - t787 * pkin(8) + t789 * t773 + t740;
t816 = sin(qJ(5));
t850 = cos(qJ(5));
t763 = t816 * t788 + t850 * t789;
t772 = t812 * qJDD(3) + t814 * t796;
t726 = t763 * qJD(5) - t850 * t771 + t816 * t772;
t762 = -t850 * t788 + t816 * t789;
t727 = -t762 * qJD(5) + t816 * t771 + t850 * t772;
t804 = qJD(5) - t838;
t716 = -0.2e1 * qJD(6) * t763 + (t762 * t804 - t727) * qJ(6) + (t763 * t804 + t726) * pkin(5) + t723;
t751 = -t804 * mrSges(7,1) + t763 * mrSges(7,2);
t752 = -t762 * mrSges(7,2) + t804 * mrSges(7,3);
t710 = m(7) * t716 + t726 * mrSges(7,1) - t727 * mrSges(7,3) - t763 * t751 + t762 * t752;
t765 = t815 * t792 - t813 * t795;
t755 = -qJDD(1) * pkin(2) - t822 * pkin(7) - t765;
t806 = t817 * t837;
t797 = t819 * qJDD(1) - t806;
t738 = (-t796 - t835) * qJ(4) + (-t797 + t806) * pkin(3) + t755;
t748 = t819 * t756 + t817 * t811;
t744 = -t821 * pkin(3) + qJDD(3) * qJ(4) + t793 * t838 + t748;
t721 = -0.2e1 * qJD(4) * t789 + t814 * t738 - t812 * t744;
t718 = (-t788 * t838 - t772) * pkin(8) + (t788 * t789 - t797) * pkin(4) + t721;
t722 = 0.2e1 * qJD(4) * t788 + t812 * t738 + t814 * t744;
t720 = -t787 * pkin(4) + t771 * pkin(8) + t773 * t838 + t722;
t715 = t816 * t718 + t850 * t720;
t741 = t762 * pkin(5) - t763 * qJ(6);
t791 = qJDD(5) - t797;
t803 = t804 ^ 2;
t712 = -t803 * pkin(5) + t791 * qJ(6) + 0.2e1 * qJD(6) * t804 - t762 * t741 + t715;
t841 = -t848 * t762 + t855 * t763 - t847 * t804;
t842 = t846 * t762 + t847 * t763 + t853 * t804;
t698 = -mrSges(6,1) * t723 - mrSges(7,1) * t716 + mrSges(7,2) * t712 + mrSges(6,3) * t715 - pkin(5) * t710 - t854 * t726 + t848 * t727 + t842 * t763 + t846 * t791 + t841 * t804;
t714 = t850 * t718 - t816 * t720;
t713 = -t791 * pkin(5) - t803 * qJ(6) + t763 * t741 + qJDD(6) - t714;
t843 = t854 * t762 - t848 * t763 - t846 * t804;
t699 = mrSges(6,2) * t723 + mrSges(7,2) * t713 - mrSges(6,3) * t714 - mrSges(7,3) * t716 - qJ(6) * t710 - t848 * t726 + t855 * t727 + t842 * t762 - t847 * t791 + t843 * t804;
t757 = Ifges(5,5) * t789 + Ifges(5,6) * t788 - Ifges(5,3) * t838;
t759 = Ifges(5,1) * t789 + Ifges(5,4) * t788 - Ifges(5,5) * t838;
t749 = -t804 * mrSges(6,2) - t762 * mrSges(6,3);
t750 = t804 * mrSges(6,1) - t763 * mrSges(6,3);
t824 = m(6) * t723 + t726 * mrSges(6,1) + t727 * mrSges(6,2) + t762 * t749 + t763 * t750 + t710;
t836 = m(7) * t712 + t791 * mrSges(7,3) + t804 * t751;
t742 = t762 * mrSges(7,1) - t763 * mrSges(7,3);
t840 = -t762 * mrSges(6,1) - t763 * mrSges(6,2) - t742;
t849 = -mrSges(6,3) - mrSges(7,2);
t702 = m(6) * t715 - t791 * mrSges(6,2) + t849 * t726 - t804 * t750 + t840 * t762 + t836;
t830 = -m(7) * t713 + t791 * mrSges(7,1) + t804 * t752;
t704 = m(6) * t714 + t791 * mrSges(6,1) + t849 * t727 + t804 * t749 + t840 * t763 + t830;
t831 = t850 * t702 - t816 * t704;
t678 = -mrSges(5,1) * t740 + mrSges(5,3) * t722 + Ifges(5,4) * t772 + Ifges(5,2) * t771 - Ifges(5,6) * t797 - pkin(4) * t824 + pkin(8) * t831 + t850 * t698 + t816 * t699 - t789 * t757 - t759 * t838;
t697 = t816 * t702 + t850 * t704;
t758 = Ifges(5,4) * t789 + Ifges(5,2) * t788 - Ifges(5,6) * t838;
t679 = mrSges(5,2) * t740 - mrSges(5,3) * t721 + Ifges(5,1) * t772 + Ifges(5,4) * t771 - Ifges(5,5) * t797 - pkin(8) * t697 - t816 * t698 + t850 * t699 + t788 * t757 + t758 * t838;
t767 = -t788 * mrSges(5,1) + t789 * mrSges(5,2);
t828 = mrSges(5,2) * t838 + t788 * mrSges(5,3);
t695 = m(5) * t721 - t797 * mrSges(5,1) - t772 * mrSges(5,3) - t789 * t767 - t828 * t838 + t697;
t770 = -mrSges(5,1) * t838 - t789 * mrSges(5,3);
t696 = m(5) * t722 + t797 * mrSges(5,2) + t771 * mrSges(5,3) + t788 * t767 + t770 * t838 + t831;
t693 = -t812 * t695 + t814 * t696;
t707 = m(5) * t740 - t771 * mrSges(5,1) + t772 * mrSges(5,2) + t789 * t770 - t788 * t828 + t824;
t782 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t817 + Ifges(4,2) * t819) * qJD(1);
t783 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t817 + Ifges(4,4) * t819) * qJD(1);
t852 = mrSges(4,1) * t747 - mrSges(4,2) * t748 + Ifges(4,5) * t796 + Ifges(4,6) * t797 + Ifges(4,3) * qJDD(3) - pkin(3) * t707 + qJ(4) * t693 + t814 * t678 + t812 * t679 + (t817 * t782 - t819 * t783) * qJD(1);
t709 = t727 * mrSges(7,2) + t763 * t742 - t830;
t851 = t846 * t726 + t847 * t727 - t841 * t762 + t843 * t763 + t853 * t791 - mrSges(6,1) * t714 + mrSges(7,1) * t713 + mrSges(6,2) * t715 - mrSges(7,3) * t712 + pkin(5) * t709 - qJ(6) * (-t726 * mrSges(7,2) - t762 * t742 + t836);
t794 = (-mrSges(4,1) * t819 + mrSges(4,2) * t817) * qJD(1);
t799 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t839;
t691 = m(4) * t748 - qJDD(3) * mrSges(4,2) + t797 * mrSges(4,3) - qJD(3) * t799 + t794 * t838 + t693;
t800 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t838;
t706 = m(4) * t747 + qJDD(3) * mrSges(4,1) - t796 * mrSges(4,3) + qJD(3) * t800 - t794 * t839 - t707;
t832 = t819 * t691 - t817 * t706;
t682 = m(3) * t766 - t822 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t832;
t692 = t814 * t695 + t812 * t696;
t825 = -m(4) * t755 + t797 * mrSges(4,1) - t796 * mrSges(4,2) - t799 * t839 + t800 * t838 - t692;
t687 = m(3) * t765 + qJDD(1) * mrSges(3,1) - t822 * mrSges(3,2) + t825;
t675 = t813 * t682 + t815 * t687;
t672 = m(2) * t801 + qJDD(1) * mrSges(2,1) - t822 * mrSges(2,2) + t675;
t833 = t815 * t682 - t813 * t687;
t673 = m(2) * t802 - t822 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t833;
t844 = t820 * t672 + t818 * t673;
t685 = t817 * t691 + t819 * t706;
t683 = m(3) * t811 + t685;
t834 = -t818 * t672 + t820 * t673;
t781 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t817 + Ifges(4,6) * t819) * qJD(1);
t668 = mrSges(4,2) * t755 - mrSges(4,3) * t747 + Ifges(4,1) * t796 + Ifges(4,4) * t797 + Ifges(4,5) * qJDD(3) - qJ(4) * t692 - qJD(3) * t782 - t812 * t678 + t814 * t679 + t781 * t838;
t677 = -t781 * t839 + t851 + Ifges(4,6) * qJDD(3) + (Ifges(4,2) + Ifges(5,3)) * t797 - t789 * t758 + Ifges(4,4) * t796 + qJD(3) * t783 + t788 * t759 - Ifges(5,6) * t771 - Ifges(5,5) * t772 + mrSges(4,3) * t748 - mrSges(4,1) * t755 - mrSges(5,1) * t721 + mrSges(5,2) * t722 - pkin(4) * t697 - pkin(3) * t692;
t827 = mrSges(2,1) * t801 + mrSges(3,1) * t765 - mrSges(2,2) * t802 - mrSges(3,2) * t766 + pkin(1) * t675 + pkin(2) * t825 + pkin(7) * t832 + t817 * t668 + t819 * t677 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t666 = -mrSges(3,1) * t811 + mrSges(3,3) * t766 + t822 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t685 - t852;
t665 = mrSges(3,2) * t811 - mrSges(3,3) * t765 + Ifges(3,5) * qJDD(1) - t822 * Ifges(3,6) - pkin(7) * t685 + t819 * t668 - t817 * t677;
t664 = -mrSges(2,2) * g(3) - mrSges(2,3) * t801 + Ifges(2,5) * qJDD(1) - t822 * Ifges(2,6) - qJ(2) * t675 + t815 * t665 - t813 * t666;
t663 = mrSges(2,1) * g(3) + mrSges(2,3) * t802 + t822 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t683 + qJ(2) * t833 + t813 * t665 + t815 * t666;
t1 = [-m(1) * g(1) + t834; -m(1) * g(2) + t844; (-m(1) - m(2)) * g(3) + t683; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t844 - t818 * t663 + t820 * t664; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t834 + t820 * t663 + t818 * t664; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t827; t827; t683; t852; t707; -t851; t709;];
tauJB  = t1;
