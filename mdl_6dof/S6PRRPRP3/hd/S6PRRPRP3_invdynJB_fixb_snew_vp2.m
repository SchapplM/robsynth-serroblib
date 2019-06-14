% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:52:35
% EndTime: 2019-05-05 03:52:46
% DurationCPUTime: 10.73s
% Computational Cost: add. (174160->317), mult. (357664->396), div. (0->0), fcn. (250965->12), ass. (0->135)
t886 = Ifges(6,1) + Ifges(7,1);
t879 = Ifges(6,4) - Ifges(7,5);
t878 = -Ifges(6,5) - Ifges(7,4);
t885 = Ifges(6,2) + Ifges(7,3);
t877 = Ifges(6,6) - Ifges(7,6);
t884 = -Ifges(6,3) - Ifges(7,2);
t838 = sin(pkin(10));
t841 = cos(pkin(10));
t827 = g(1) * t838 - g(2) * t841;
t828 = -g(1) * t841 - g(2) * t838;
t836 = -g(3) + qJDD(1);
t839 = sin(pkin(6));
t842 = cos(pkin(6));
t845 = sin(qJ(2));
t847 = cos(qJ(2));
t787 = -t845 * t828 + (t827 * t842 + t836 * t839) * t847;
t873 = t842 * t845;
t874 = t839 * t845;
t788 = t827 * t873 + t847 * t828 + t836 * t874;
t849 = qJD(2) ^ 2;
t782 = -pkin(2) * t849 + qJDD(2) * pkin(8) + t788;
t804 = -t827 * t839 + t836 * t842;
t844 = sin(qJ(3));
t846 = cos(qJ(3));
t772 = -t844 * t782 + t846 * t804;
t823 = (-pkin(3) * t846 - qJ(4) * t844) * qJD(2);
t848 = qJD(3) ^ 2;
t867 = qJD(2) * t844;
t755 = -qJDD(3) * pkin(3) - t848 * qJ(4) + t823 * t867 + qJDD(4) - t772;
t865 = qJD(2) * qJD(3);
t863 = t846 * t865;
t825 = qJDD(2) * t844 + t863;
t837 = sin(pkin(11));
t840 = cos(pkin(11));
t801 = qJDD(3) * t840 - t825 * t837;
t818 = qJD(3) * t837 + t840 * t867;
t866 = qJD(2) * t846;
t803 = -pkin(4) * t866 - pkin(9) * t818;
t817 = qJD(3) * t840 - t837 * t867;
t816 = t817 ^ 2;
t753 = -t801 * pkin(4) - t816 * pkin(9) + t818 * t803 + t755;
t843 = sin(qJ(5));
t881 = cos(qJ(5));
t795 = t843 * t817 + t881 * t818;
t802 = qJDD(3) * t837 + t825 * t840;
t762 = t795 * qJD(5) - t881 * t801 + t843 * t802;
t794 = -t881 * t817 + t843 * t818;
t763 = -t794 * qJD(5) + t843 * t801 + t881 * t802;
t833 = qJD(5) - t866;
t746 = -0.2e1 * qJD(6) * t795 + (t794 * t833 - t763) * qJ(6) + (t795 * t833 + t762) * pkin(5) + t753;
t785 = -mrSges(7,1) * t833 + mrSges(7,2) * t795;
t786 = -mrSges(7,2) * t794 + mrSges(7,3) * t833;
t740 = m(7) * t746 + t762 * mrSges(7,1) - t763 * mrSges(7,3) - t795 * t785 + t794 * t786;
t773 = t846 * t782 + t844 * t804;
t756 = -pkin(3) * t848 + qJDD(3) * qJ(4) + t823 * t866 + t773;
t781 = -qJDD(2) * pkin(2) - t849 * pkin(8) - t787;
t834 = t844 * t865;
t826 = qJDD(2) * t846 - t834;
t761 = (-t825 - t863) * qJ(4) + (-t826 + t834) * pkin(3) + t781;
t751 = -0.2e1 * qJD(4) * t818 - t837 * t756 + t840 * t761;
t748 = (-t817 * t866 - t802) * pkin(9) + (t817 * t818 - t826) * pkin(4) + t751;
t752 = 0.2e1 * qJD(4) * t817 + t840 * t756 + t837 * t761;
t750 = -pkin(4) * t816 + pkin(9) * t801 + t803 * t866 + t752;
t745 = t843 * t748 + t881 * t750;
t774 = pkin(5) * t794 - qJ(6) * t795;
t820 = qJDD(5) - t826;
t832 = t833 ^ 2;
t742 = -pkin(5) * t832 + qJ(6) * t820 + 0.2e1 * qJD(6) * t833 - t774 * t794 + t745;
t869 = t879 * t794 - t795 * t886 + t878 * t833;
t870 = t794 * t877 + t795 * t878 + t833 * t884;
t728 = -mrSges(6,1) * t753 - mrSges(7,1) * t746 + mrSges(7,2) * t742 + mrSges(6,3) * t745 - pkin(5) * t740 - t762 * t885 + t879 * t763 + t870 * t795 + t877 * t820 - t869 * t833;
t744 = t881 * t748 - t843 * t750;
t743 = -t820 * pkin(5) - t832 * qJ(6) + t795 * t774 + qJDD(6) - t744;
t871 = t794 * t885 - t795 * t879 - t833 * t877;
t729 = mrSges(6,2) * t753 + mrSges(7,2) * t743 - mrSges(6,3) * t744 - mrSges(7,3) * t746 - qJ(6) * t740 - t879 * t762 + t763 * t886 + t870 * t794 - t878 * t820 + t871 * t833;
t789 = Ifges(5,5) * t818 + Ifges(5,6) * t817 - Ifges(5,3) * t866;
t791 = Ifges(5,1) * t818 + Ifges(5,4) * t817 - Ifges(5,5) * t866;
t783 = -mrSges(6,2) * t833 - mrSges(6,3) * t794;
t784 = mrSges(6,1) * t833 - mrSges(6,3) * t795;
t851 = m(6) * t753 + t762 * mrSges(6,1) + t763 * mrSges(6,2) + t794 * t783 + t795 * t784 + t740;
t864 = m(7) * t742 + t820 * mrSges(7,3) + t833 * t785;
t775 = mrSges(7,1) * t794 - mrSges(7,3) * t795;
t868 = -mrSges(6,1) * t794 - mrSges(6,2) * t795 - t775;
t880 = -mrSges(6,3) - mrSges(7,2);
t732 = m(6) * t745 - t820 * mrSges(6,2) + t880 * t762 - t833 * t784 + t868 * t794 + t864;
t858 = -m(7) * t743 + t820 * mrSges(7,1) + t833 * t786;
t734 = m(6) * t744 + t820 * mrSges(6,1) + t880 * t763 + t833 * t783 + t868 * t795 + t858;
t860 = t881 * t732 - t734 * t843;
t707 = -mrSges(5,1) * t755 + mrSges(5,3) * t752 + Ifges(5,4) * t802 + Ifges(5,2) * t801 - Ifges(5,6) * t826 - pkin(4) * t851 + pkin(9) * t860 + t881 * t728 + t843 * t729 - t818 * t789 - t791 * t866;
t727 = t843 * t732 + t881 * t734;
t790 = Ifges(5,4) * t818 + Ifges(5,2) * t817 - Ifges(5,6) * t866;
t708 = mrSges(5,2) * t755 - mrSges(5,3) * t751 + Ifges(5,1) * t802 + Ifges(5,4) * t801 - Ifges(5,5) * t826 - pkin(9) * t727 - t843 * t728 + t881 * t729 + t817 * t789 + t790 * t866;
t796 = -mrSges(5,1) * t817 + mrSges(5,2) * t818;
t855 = mrSges(5,2) * t866 + mrSges(5,3) * t817;
t725 = m(5) * t751 - t826 * mrSges(5,1) - t802 * mrSges(5,3) - t818 * t796 - t855 * t866 + t727;
t800 = -mrSges(5,1) * t866 - mrSges(5,3) * t818;
t726 = m(5) * t752 + mrSges(5,2) * t826 + mrSges(5,3) * t801 + t796 * t817 + t800 * t866 + t860;
t723 = -t725 * t837 + t840 * t726;
t737 = m(5) * t755 - t801 * mrSges(5,1) + t802 * mrSges(5,2) + t818 * t800 - t817 * t855 + t851;
t811 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t844 + Ifges(4,2) * t846) * qJD(2);
t812 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t844 + Ifges(4,4) * t846) * qJD(2);
t883 = mrSges(4,1) * t772 - mrSges(4,2) * t773 + Ifges(4,5) * t825 + Ifges(4,6) * t826 + Ifges(4,3) * qJDD(3) - pkin(3) * t737 + qJ(4) * t723 + t840 * t707 + t837 * t708 + (t811 * t844 - t812 * t846) * qJD(2);
t739 = t763 * mrSges(7,2) + t795 * t775 - t858;
t882 = t877 * t762 + t878 * t763 + t869 * t794 + t871 * t795 + t884 * t820 - mrSges(6,1) * t744 + mrSges(7,1) * t743 + mrSges(6,2) * t745 - mrSges(7,3) * t742 + pkin(5) * t739 - qJ(6) * (-t762 * mrSges(7,2) - t794 * t775 + t864);
t722 = t725 * t840 + t726 * t837;
t829 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t867;
t830 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t866;
t852 = -m(4) * t781 + t826 * mrSges(4,1) - mrSges(4,2) * t825 - t829 * t867 + t830 * t866 - t722;
t718 = m(3) * t787 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t849 + t852;
t875 = t718 * t847;
t824 = (-mrSges(4,1) * t846 + mrSges(4,2) * t844) * qJD(2);
t721 = m(4) * t773 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t826 - qJD(3) * t829 + t824 * t866 + t723;
t736 = m(4) * t772 + qJDD(3) * mrSges(4,1) - t825 * mrSges(4,3) + qJD(3) * t830 - t824 * t867 - t737;
t861 = t846 * t721 - t736 * t844;
t712 = m(3) * t788 - mrSges(3,1) * t849 - qJDD(2) * mrSges(3,2) + t861;
t715 = t844 * t721 + t846 * t736;
t714 = m(3) * t804 + t715;
t701 = t712 * t873 - t714 * t839 + t842 * t875;
t699 = m(2) * t827 + t701;
t705 = t847 * t712 - t718 * t845;
t704 = m(2) * t828 + t705;
t872 = t841 * t699 + t838 * t704;
t700 = t712 * t874 + t842 * t714 + t839 * t875;
t862 = -t699 * t838 + t841 * t704;
t859 = m(2) * t836 + t700;
t810 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t844 + Ifges(4,6) * t846) * qJD(2);
t697 = mrSges(4,2) * t781 - mrSges(4,3) * t772 + Ifges(4,1) * t825 + Ifges(4,4) * t826 + Ifges(4,5) * qJDD(3) - qJ(4) * t722 - qJD(3) * t811 - t707 * t837 + t708 * t840 + t810 * t866;
t706 = -t810 * t867 + t882 + Ifges(4,6) * qJDD(3) + (Ifges(5,3) + Ifges(4,2)) * t826 + t817 * t791 - t818 * t790 + Ifges(4,4) * t825 - Ifges(5,6) * t801 - Ifges(5,5) * t802 + qJD(3) * t812 + mrSges(4,3) * t773 - mrSges(4,1) * t781 - mrSges(5,1) * t751 + mrSges(5,2) * t752 - pkin(4) * t727 - pkin(3) * t722;
t695 = mrSges(3,2) * t804 - mrSges(3,3) * t787 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t849 - pkin(8) * t715 + t697 * t846 - t706 * t844;
t696 = -mrSges(3,1) * t804 + mrSges(3,3) * t788 + t849 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t715 - t883;
t854 = pkin(7) * t705 + t695 * t845 + t696 * t847;
t694 = mrSges(3,1) * t787 - mrSges(3,2) * t788 + Ifges(3,3) * qJDD(2) + pkin(2) * t852 + pkin(8) * t861 + t844 * t697 + t846 * t706;
t693 = mrSges(2,2) * t836 - mrSges(2,3) * t827 + t847 * t695 - t845 * t696 + (-t700 * t839 - t701 * t842) * pkin(7);
t692 = -mrSges(2,1) * t836 + mrSges(2,3) * t828 - pkin(1) * t700 - t839 * t694 + t854 * t842;
t1 = [-m(1) * g(1) + t862; -m(1) * g(2) + t872; -m(1) * g(3) + t859; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t872 - t838 * t692 + t841 * t693; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t862 + t841 * t692 + t838 * t693; -mrSges(1,1) * g(2) + mrSges(2,1) * t827 + mrSges(1,2) * g(1) - mrSges(2,2) * t828 + pkin(1) * t701 + t842 * t694 + t854 * t839; t859; t694; t883; t737; -t882; t739;];
tauJB  = t1;
