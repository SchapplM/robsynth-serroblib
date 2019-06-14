% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:50:15
% EndTime: 2019-05-05 06:50:24
% DurationCPUTime: 5.75s
% Computational Cost: add. (69460->296), mult. (130936->351), div. (0->0), fcn. (84882->10), ass. (0->129)
t889 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t869 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t868 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t888 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t867 = -Ifges(5,6) + Ifges(6,5) + Ifges(7,4);
t887 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t834 = sin(pkin(10));
t836 = cos(pkin(10));
t822 = g(1) * t834 - g(2) * t836;
t823 = -g(1) * t836 - g(2) * t834;
t833 = -g(3) + qJDD(1);
t842 = cos(qJ(2));
t837 = cos(pkin(6));
t840 = sin(qJ(2));
t875 = t837 * t840;
t835 = sin(pkin(6));
t876 = t835 * t840;
t763 = t822 * t875 + t842 * t823 + t833 * t876;
t844 = qJD(2) ^ 2;
t759 = -pkin(2) * t844 + qJDD(2) * pkin(8) + t763;
t800 = -t822 * t835 + t833 * t837;
t839 = sin(qJ(3));
t841 = cos(qJ(3));
t755 = t841 * t759 + t839 * t800;
t819 = (-pkin(3) * t841 - pkin(9) * t839) * qJD(2);
t843 = qJD(3) ^ 2;
t871 = qJD(2) * t841;
t751 = -pkin(3) * t843 + qJDD(3) * pkin(9) + t819 * t871 + t755;
t762 = -t840 * t823 + (t822 * t837 + t833 * t835) * t842;
t758 = -qJDD(2) * pkin(2) - t844 * pkin(8) - t762;
t870 = qJD(2) * qJD(3);
t860 = t841 * t870;
t820 = qJDD(2) * t839 + t860;
t861 = t839 * t870;
t821 = qJDD(2) * t841 - t861;
t753 = (-t820 - t860) * pkin(9) + (-t821 + t861) * pkin(3) + t758;
t838 = sin(qJ(4));
t881 = cos(qJ(4));
t746 = -t838 * t751 + t753 * t881;
t872 = qJD(2) * t839;
t816 = -qJD(3) * t881 + t838 * t872;
t817 = t838 * qJD(3) + t872 * t881;
t788 = pkin(4) * t816 - qJ(5) * t817;
t813 = qJDD(4) - t821;
t828 = -qJD(4) + t871;
t827 = t828 ^ 2;
t744 = -t813 * pkin(4) - t827 * qJ(5) + t817 * t788 + qJDD(5) - t746;
t783 = -t816 * qJD(4) + t838 * qJDD(3) + t820 * t881;
t790 = -mrSges(6,2) * t816 - mrSges(6,3) * t817;
t886 = -m(6) * t744 - t783 * mrSges(6,1) - t817 * t790;
t798 = mrSges(6,1) * t817 - mrSges(6,2) * t828;
t782 = qJD(4) * t817 - qJDD(3) * t881 + t820 * t838;
t754 = -t839 * t759 + t841 * t800;
t750 = -qJDD(3) * pkin(3) - t843 * pkin(9) + t819 * t872 - t754;
t877 = t816 * t828;
t883 = -2 * qJD(5);
t847 = (-t783 - t877) * qJ(5) + t750 + (-pkin(4) * t828 + t883) * t817;
t745 = t782 * pkin(4) + t847;
t796 = mrSges(6,1) * t816 + mrSges(6,3) * t828;
t794 = pkin(5) * t817 + qJ(6) * t828;
t812 = t816 ^ 2;
t882 = 2 * qJD(6);
t742 = -t812 * pkin(5) + t816 * t882 - t817 * t794 + (pkin(4) + qJ(6)) * t782 + t847;
t795 = mrSges(7,1) * t817 + mrSges(7,3) * t828;
t797 = -mrSges(7,1) * t816 - mrSges(7,2) * t828;
t855 = m(7) * t742 - t783 * mrSges(7,2) + t782 * mrSges(7,3) - t817 * t795 + t816 * t797;
t850 = -m(6) * t745 + t782 * mrSges(6,2) + t816 * t796 - t855;
t735 = -t783 * mrSges(6,3) - t817 * t798 - t850;
t787 = -mrSges(7,2) * t817 + mrSges(7,3) * t816;
t747 = t881 * t751 + t838 * t753;
t849 = -t827 * pkin(4) + t813 * qJ(5) - t816 * t788 + t747;
t740 = -t782 * pkin(5) - t812 * qJ(6) + qJDD(6) + (t883 - t794) * t828 + t849;
t865 = m(7) * t740 + t813 * mrSges(7,2) - t828 * t795;
t737 = -t782 * mrSges(7,1) - t816 * t787 + t865;
t743 = 0.2e1 * qJD(5) * t828 - t849;
t862 = -t869 * t816 + t889 * t817 - t868 * t828;
t864 = -t867 * t816 - t868 * t817 + t887 * t828;
t710 = -mrSges(5,1) * t750 - mrSges(6,1) * t743 + mrSges(7,1) * t740 + mrSges(6,2) * t745 + mrSges(5,3) * t747 - mrSges(7,3) * t742 - pkin(4) * t735 + pkin(5) * t737 - qJ(6) * t855 + t888 * t782 + t869 * t783 - t867 * t813 + t864 * t817 - t862 * t828;
t738 = t828 * t882 + (t816 * t817 - t813) * qJ(6) + (t783 - t877) * pkin(5) + t744;
t856 = -m(7) * t738 + t813 * mrSges(7,3) - t828 * t797;
t736 = t783 * mrSges(7,1) + t817 * t787 - t856;
t863 = t888 * t816 + t869 * t817 + t867 * t828;
t718 = mrSges(6,1) * t744 + mrSges(7,1) * t738 + mrSges(5,2) * t750 - mrSges(7,2) * t742 - mrSges(5,3) * t746 - mrSges(6,3) * t745 + pkin(5) * t736 - qJ(5) * t735 - t869 * t782 + t889 * t783 + t868 * t813 + t864 * t816 + t863 * t828;
t789 = mrSges(5,1) * t816 + mrSges(5,2) * t817;
t792 = mrSges(5,2) * t828 - mrSges(5,3) * t816;
t879 = -mrSges(7,1) - mrSges(5,3);
t729 = m(5) * t746 + (-t792 + t796) * t828 + (-t787 - t789) * t817 + (mrSges(5,1) - mrSges(6,2)) * t813 + t879 * t783 + t856 + t886;
t793 = -mrSges(5,1) * t828 - mrSges(5,3) * t817;
t852 = -m(6) * t743 + t813 * mrSges(6,3) - t828 * t798 + t865;
t873 = -t787 - t790;
t731 = m(5) * t747 - t813 * mrSges(5,2) + t828 * t793 + (-t789 + t873) * t816 + (-mrSges(6,1) + t879) * t782 + t852;
t726 = -t729 * t838 + t881 * t731;
t732 = -m(5) * t750 - t782 * mrSges(5,1) - t816 * t792 + (-t793 + t798) * t817 + (-mrSges(5,2) + mrSges(6,3)) * t783 + t850;
t805 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t839 + Ifges(4,2) * t841) * qJD(2);
t806 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t839 + Ifges(4,4) * t841) * qJD(2);
t885 = mrSges(4,1) * t754 - mrSges(4,2) * t755 + Ifges(4,5) * t820 + Ifges(4,6) * t821 + Ifges(4,3) * qJDD(3) + pkin(3) * t732 + pkin(9) * t726 + (t805 * t839 - t806 * t841) * qJD(2) + t710 * t881 + t838 * t718;
t733 = t813 * mrSges(6,2) - t828 * t796 + t736 - t886;
t884 = t782 * t867 + t783 * t868 + t887 * t813 + t816 * t862 + t817 * t863 + mrSges(5,1) * t746 - mrSges(5,2) * t747 + mrSges(6,2) * t744 + mrSges(7,2) * t740 - mrSges(6,3) * t743 - mrSges(7,3) * t738 - pkin(4) * t733 + qJ(5) * (t873 * t816 + (-mrSges(6,1) - mrSges(7,1)) * t782 + t852) - qJ(6) * t736;
t725 = t729 * t881 + t838 * t731;
t824 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t872;
t825 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t871;
t848 = -m(4) * t758 + t821 * mrSges(4,1) - t820 * mrSges(4,2) - t824 * t872 + t825 * t871 - t725;
t721 = m(3) * t762 + qJDD(2) * mrSges(3,1) - t844 * mrSges(3,2) + t848;
t878 = t721 * t842;
t818 = (-mrSges(4,1) * t841 + mrSges(4,2) * t839) * qJD(2);
t724 = m(4) * t755 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t821 - qJD(3) * t824 + t818 * t871 + t726;
t728 = m(4) * t754 + qJDD(3) * mrSges(4,1) - t820 * mrSges(4,3) + qJD(3) * t825 - t818 * t872 + t732;
t858 = t841 * t724 - t728 * t839;
t714 = m(3) * t763 - mrSges(3,1) * t844 - qJDD(2) * mrSges(3,2) + t858;
t717 = t839 * t724 + t841 * t728;
t716 = m(3) * t800 + t717;
t703 = t714 * t875 - t716 * t835 + t837 * t878;
t701 = m(2) * t822 + t703;
t709 = t842 * t714 - t721 * t840;
t708 = m(2) * t823 + t709;
t874 = t836 * t701 + t834 * t708;
t702 = t714 * t876 + t837 * t716 + t835 * t878;
t859 = -t701 * t834 + t836 * t708;
t857 = m(2) * t833 + t702;
t804 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t839 + Ifges(4,6) * t841) * qJD(2);
t704 = mrSges(4,2) * t758 - mrSges(4,3) * t754 + Ifges(4,1) * t820 + Ifges(4,4) * t821 + Ifges(4,5) * qJDD(3) - pkin(9) * t725 - qJD(3) * t805 - t838 * t710 + t718 * t881 + t804 * t871;
t705 = -mrSges(4,1) * t758 + mrSges(4,3) * t755 + Ifges(4,4) * t820 + Ifges(4,2) * t821 + Ifges(4,6) * qJDD(3) - pkin(3) * t725 + qJD(3) * t806 - t804 * t872 - t884;
t698 = mrSges(3,2) * t800 - mrSges(3,3) * t762 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t844 - pkin(8) * t717 + t704 * t841 - t705 * t839;
t699 = -mrSges(3,1) * t800 + mrSges(3,3) * t763 + t844 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t717 - t885;
t851 = pkin(7) * t709 + t698 * t840 + t699 * t842;
t697 = mrSges(3,1) * t762 - mrSges(3,2) * t763 + Ifges(3,3) * qJDD(2) + pkin(2) * t848 + pkin(8) * t858 + t839 * t704 + t841 * t705;
t696 = mrSges(2,2) * t833 - mrSges(2,3) * t822 + t842 * t698 - t840 * t699 + (-t702 * t835 - t703 * t837) * pkin(7);
t695 = -mrSges(2,1) * t833 + mrSges(2,3) * t823 - pkin(1) * t702 - t835 * t697 + t837 * t851;
t1 = [-m(1) * g(1) + t859; -m(1) * g(2) + t874; -m(1) * g(3) + t857; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t874 - t834 * t695 + t836 * t696; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t859 + t836 * t695 + t834 * t696; -mrSges(1,1) * g(2) + mrSges(2,1) * t822 + mrSges(1,2) * g(1) - mrSges(2,2) * t823 + pkin(1) * t703 + t837 * t697 + t835 * t851; t857; t697; t885; t884; t733; t737;];
tauJB  = t1;
