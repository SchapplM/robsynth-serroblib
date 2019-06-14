% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:46:48
% EndTime: 2019-05-06 01:46:56
% DurationCPUTime: 5.79s
% Computational Cost: add. (69920->315), mult. (136595->374), div. (0->0), fcn. (89146->8), ass. (0->128)
t886 = Ifges(6,1) + Ifges(7,1);
t878 = Ifges(6,4) - Ifges(7,5);
t876 = -Ifges(6,5) - Ifges(7,4);
t885 = Ifges(6,2) + Ifges(7,3);
t874 = Ifges(6,6) - Ifges(7,6);
t884 = -Ifges(6,3) - Ifges(7,2);
t839 = sin(qJ(1));
t842 = cos(qJ(1));
t818 = -t842 * g(1) - t839 * g(2);
t854 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t818;
t837 = sin(qJ(4));
t838 = sin(qJ(3));
t840 = cos(qJ(4));
t841 = cos(qJ(3));
t802 = (t837 * t841 + t838 * t840) * qJD(1);
t865 = qJD(1) * qJD(3);
t811 = -qJDD(1) * t838 - t841 * t865;
t862 = t838 * t865;
t812 = qJDD(1) * t841 - t862;
t769 = -qJD(4) * t802 + t811 * t837 + t812 * t840;
t803 = (-t837 * t838 + t840 * t841) * qJD(1);
t827 = qJD(3) + qJD(4);
t836 = sin(qJ(5));
t881 = cos(qJ(5));
t783 = t836 * t803 - t827 * t881;
t826 = qJDD(3) + qJDD(4);
t739 = -t783 * qJD(5) + t769 * t881 + t836 * t826;
t784 = t803 * t881 + t836 * t827;
t756 = t783 * mrSges(7,1) - mrSges(7,3) * t784;
t866 = qJD(1) * t841;
t816 = (qJD(3) * pkin(3)) - pkin(8) * t866;
t833 = t838 ^ 2;
t843 = qJD(1) ^ 2;
t882 = -pkin(1) - pkin(7);
t760 = -t811 * pkin(3) + t816 * t866 + (-pkin(8) * t833 + t882) * t843 + t854;
t768 = -qJD(4) * t803 + t811 * t840 - t812 * t837;
t729 = (t802 * t827 - t769) * pkin(9) + (t803 * t827 - t768) * pkin(4) + t760;
t817 = t839 * g(1) - t842 * g(2);
t853 = -t843 * qJ(2) + qJDD(2) - t817;
t791 = qJDD(1) * t882 + t853;
t781 = t838 * g(3) + t841 * t791;
t753 = (-t812 - t862) * pkin(8) + (-t838 * t841 * t843 + qJDD(3)) * pkin(3) + t781;
t782 = -g(3) * t841 + t838 * t791;
t754 = -pkin(3) * t833 * t843 + pkin(8) * t811 - qJD(3) * t816 + t782;
t735 = t837 * t753 + t840 * t754;
t779 = pkin(4) * t802 - pkin(9) * t803;
t825 = t827 ^ 2;
t732 = -pkin(4) * t825 + pkin(9) * t826 - t779 * t802 + t735;
t726 = t729 * t881 - t836 * t732;
t755 = t783 * pkin(5) - qJ(6) * t784;
t767 = qJDD(5) - t768;
t798 = qJD(5) + t802;
t796 = t798 ^ 2;
t724 = -t767 * pkin(5) - t796 * qJ(6) + t784 * t755 + qJDD(6) - t726;
t770 = -t783 * mrSges(7,2) + mrSges(7,3) * t798;
t857 = -m(7) * t724 + t767 * mrSges(7,1) + t798 * t770;
t720 = t739 * mrSges(7,2) + t784 * t756 - t857;
t727 = t836 * t729 + t881 * t732;
t723 = -pkin(5) * t796 + t767 * qJ(6) + 0.2e1 * qJD(6) * t798 - t783 * t755 + t727;
t738 = t784 * qJD(5) + t836 * t769 - t826 * t881;
t773 = -mrSges(7,1) * t798 + mrSges(7,2) * t784;
t863 = m(7) * t723 + t767 * mrSges(7,3) + t798 * t773;
t869 = t878 * t783 - t886 * t784 + t876 * t798;
t870 = t885 * t783 - t878 * t784 - t874 * t798;
t883 = -t738 * t874 - t739 * t876 - t884 * t767 - t783 * t869 - t784 * t870 + mrSges(6,1) * t726 - mrSges(7,1) * t724 - mrSges(6,2) * t727 + mrSges(7,3) * t723 - pkin(5) * t720 + qJ(6) * (-t738 * mrSges(7,2) - t783 * t756 + t863);
t880 = mrSges(2,1) - mrSges(3,2);
t879 = -mrSges(6,3) - mrSges(7,2);
t877 = Ifges(2,5) - Ifges(3,4);
t875 = (-Ifges(2,6) + Ifges(3,5));
t778 = mrSges(5,1) * t802 + mrSges(5,2) * t803;
t789 = mrSges(5,1) * t827 - mrSges(5,3) * t803;
t772 = mrSges(6,1) * t798 - mrSges(6,3) * t784;
t868 = -t783 * mrSges(6,1) - mrSges(6,2) * t784 - t756;
t715 = m(6) * t727 - t767 * mrSges(6,2) + t738 * t879 - t798 * t772 + t783 * t868 + t863;
t771 = -mrSges(6,2) * t798 - t783 * mrSges(6,3);
t717 = m(6) * t726 + t767 * mrSges(6,1) + t739 * t879 + t798 * t771 + t784 * t868 + t857;
t858 = t881 * t715 - t717 * t836;
t704 = m(5) * t735 - mrSges(5,2) * t826 + t768 * mrSges(5,3) - t778 * t802 - t789 * t827 + t858;
t734 = t840 * t753 - t837 * t754;
t788 = -mrSges(5,2) * t827 - mrSges(5,3) * t802;
t731 = -t826 * pkin(4) - t825 * pkin(9) + t803 * t779 - t734;
t725 = -0.2e1 * qJD(6) * t784 + (t783 * t798 - t739) * qJ(6) + (t784 * t798 + t738) * pkin(5) + t731;
t721 = m(7) * t725 + t738 * mrSges(7,1) - t739 * mrSges(7,3) + t783 * t770 - t784 * t773;
t846 = -m(6) * t731 - t738 * mrSges(6,1) - t739 * mrSges(6,2) - t783 * t771 - t772 * t784 - t721;
t712 = m(5) * t734 + mrSges(5,1) * t826 - t769 * mrSges(5,3) - t778 * t803 + t788 * t827 + t846;
t697 = t837 * t704 + t840 * t712;
t810 = (mrSges(4,1) * t838 + mrSges(4,2) * t841) * qJD(1);
t867 = qJD(1) * t838;
t814 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t867;
t694 = m(4) * t781 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t812 + qJD(3) * t814 - t810 * t866 + t697;
t815 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t866;
t859 = t840 * t704 - t712 * t837;
t695 = m(4) * t782 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t811 - qJD(3) * t815 - t810 * t867 + t859;
t691 = t841 * t694 + t838 * t695;
t797 = -qJDD(1) * pkin(1) + t853;
t852 = -m(3) * t797 + (t843 * mrSges(3,3)) - t691;
t686 = m(2) * t817 - (t843 * mrSges(2,2)) + qJDD(1) * t880 + t852;
t794 = t843 * pkin(1) - t854;
t790 = t843 * t882 + t854;
t710 = t836 * t715 + t881 * t717;
t851 = m(5) * t760 - t768 * mrSges(5,1) + t769 * mrSges(5,2) + t788 * t802 + t803 * t789 + t710;
t848 = -m(4) * t790 + mrSges(4,1) * t811 - t812 * mrSges(4,2) - t814 * t867 - t815 * t866 - t851;
t845 = -m(3) * t794 + (t843 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t848;
t700 = m(2) * t818 - (mrSges(2,1) * t843) - qJDD(1) * mrSges(2,2) + t845;
t872 = t842 * t686 + t839 * t700;
t871 = t874 * t783 + t876 * t784 + t884 * t798;
t861 = -t686 * t839 + t842 * t700;
t860 = -t838 * t694 + t841 * t695;
t706 = -mrSges(6,1) * t731 - mrSges(7,1) * t725 + mrSges(7,2) * t723 + mrSges(6,3) * t727 - pkin(5) * t721 - t885 * t738 + t878 * t739 + t874 * t767 + t871 * t784 - t869 * t798;
t708 = mrSges(6,2) * t731 + mrSges(7,2) * t724 - mrSges(6,3) * t726 - mrSges(7,3) * t725 - qJ(6) * t721 - t878 * t738 + t886 * t739 - t876 * t767 + t871 * t783 + t870 * t798;
t775 = Ifges(5,4) * t803 - Ifges(5,2) * t802 + Ifges(5,6) * t827;
t776 = Ifges(5,1) * t803 - Ifges(5,4) * t802 + Ifges(5,5) * t827;
t850 = mrSges(5,1) * t734 - mrSges(5,2) * t735 + Ifges(5,5) * t769 + Ifges(5,6) * t768 + Ifges(5,3) * t826 + pkin(4) * t846 + pkin(9) * t858 + t881 * t706 + t836 * t708 + t803 * t775 + t802 * t776;
t774 = Ifges(5,5) * t803 - Ifges(5,6) * t802 + Ifges(5,3) * t827;
t687 = mrSges(5,2) * t760 - mrSges(5,3) * t734 + Ifges(5,1) * t769 + Ifges(5,4) * t768 + Ifges(5,5) * t826 - pkin(9) * t710 - t836 * t706 + t708 * t881 - t802 * t774 - t827 * t775;
t692 = -mrSges(5,1) * t760 + mrSges(5,3) * t735 + Ifges(5,4) * t769 + Ifges(5,2) * t768 + Ifges(5,6) * t826 - pkin(4) * t710 - t803 * t774 + t827 * t776 - t883;
t799 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t841 - Ifges(4,6) * t838) * qJD(1);
t801 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t841 - Ifges(4,4) * t838) * qJD(1);
t682 = -mrSges(4,1) * t790 + mrSges(4,3) * t782 + Ifges(4,4) * t812 + Ifges(4,2) * t811 + Ifges(4,6) * qJDD(3) - pkin(3) * t851 + pkin(8) * t859 + qJD(3) * t801 + t837 * t687 + t840 * t692 - t799 * t866;
t800 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t841 - Ifges(4,2) * t838) * qJD(1);
t684 = mrSges(4,2) * t790 - mrSges(4,3) * t781 + Ifges(4,1) * t812 + Ifges(4,4) * t811 + Ifges(4,5) * qJDD(3) - pkin(8) * t697 - qJD(3) * t800 + t687 * t840 - t692 * t837 - t799 * t867;
t689 = qJDD(1) * mrSges(3,2) - t852;
t847 = mrSges(2,1) * t817 - mrSges(2,2) * t818 + mrSges(3,2) * t797 - mrSges(3,3) * t794 - pkin(1) * t689 - pkin(7) * t691 + qJ(2) * t845 - t682 * t838 + t841 * t684 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t844 = mrSges(4,1) * t781 - mrSges(4,2) * t782 + Ifges(4,5) * t812 + Ifges(4,6) * t811 + Ifges(4,3) * qJDD(3) + pkin(3) * t697 + t800 * t866 + t801 * t867 + t850;
t690 = -m(3) * g(3) + t860;
t681 = pkin(2) * t691 - qJ(2) * t690 + t844 + t877 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t875 * t843) - mrSges(2,3) * t817 + mrSges(3,1) * t797;
t680 = -mrSges(3,1) * t794 + mrSges(2,3) * t818 - pkin(1) * t690 - pkin(2) * t848 - pkin(7) * t860 + g(3) * t880 - qJDD(1) * t875 - t841 * t682 - t838 * t684 + t843 * t877;
t1 = [-m(1) * g(1) + t861; -m(1) * g(2) + t872; (-m(1) - m(2) - m(3)) * g(3) + t860; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t872 - t839 * t680 + t842 * t681; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t861 + t842 * t680 + t839 * t681; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t847; t847; t689; t844; t850; t883; t720;];
tauJB  = t1;
