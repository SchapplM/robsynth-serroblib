% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP2
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
% Datum: 2019-05-05 17:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:33:08
% EndTime: 2019-05-05 17:33:15
% DurationCPUTime: 6.80s
% Computational Cost: add. (99199->319), mult. (212288->391), div. (0->0), fcn. (137415->10), ass. (0->129)
t883 = -2 * qJD(4);
t882 = Ifges(6,1) + Ifges(7,1);
t875 = Ifges(6,4) - Ifges(7,5);
t874 = -Ifges(6,5) - Ifges(7,4);
t881 = Ifges(6,2) + Ifges(7,3);
t873 = Ifges(6,6) - Ifges(7,6);
t880 = -Ifges(6,3) - Ifges(7,2);
t842 = sin(qJ(1));
t844 = cos(qJ(1));
t824 = t842 * g(1) - g(2) * t844;
t815 = qJDD(1) * pkin(1) + t824;
t825 = -g(1) * t844 - g(2) * t842;
t846 = qJD(1) ^ 2;
t817 = -pkin(1) * t846 + t825;
t837 = sin(pkin(9));
t839 = cos(pkin(9));
t793 = t837 * t815 + t839 * t817;
t783 = -pkin(2) * t846 + qJDD(1) * pkin(7) + t793;
t835 = -g(3) + qJDD(2);
t841 = sin(qJ(3));
t843 = cos(qJ(3));
t772 = -t841 * t783 + t843 * t835;
t862 = qJD(1) * qJD(3);
t860 = t843 * t862;
t818 = qJDD(1) * t841 + t860;
t751 = (-t818 + t860) * qJ(4) + (t841 * t843 * t846 + qJDD(3)) * pkin(3) + t772;
t773 = t843 * t783 + t841 * t835;
t819 = qJDD(1) * t843 - t841 * t862;
t865 = qJD(1) * t841;
t821 = qJD(3) * pkin(3) - qJ(4) * t865;
t834 = t843 ^ 2;
t752 = -pkin(3) * t834 * t846 + qJ(4) * t819 - qJD(3) * t821 + t773;
t836 = sin(pkin(10));
t838 = cos(pkin(10));
t805 = (t836 * t843 + t838 * t841) * qJD(1);
t744 = t838 * t751 - t836 * t752 + t805 * t883;
t804 = (t836 * t841 - t838 * t843) * qJD(1);
t745 = t836 * t751 + t838 * t752 + t804 * t883;
t786 = pkin(4) * t804 - pkin(8) * t805;
t845 = qJD(3) ^ 2;
t743 = -pkin(4) * t845 + qJDD(3) * pkin(8) - t786 * t804 + t745;
t792 = t839 * t815 - t837 * t817;
t851 = -qJDD(1) * pkin(2) - t792;
t753 = -t819 * pkin(3) + qJDD(4) + t821 * t865 + (-qJ(4) * t834 - pkin(7)) * t846 + t851;
t794 = -t818 * t836 + t819 * t838;
t795 = t818 * t838 + t819 * t836;
t747 = (qJD(3) * t804 - t795) * pkin(8) + (qJD(3) * t805 - t794) * pkin(4) + t753;
t840 = sin(qJ(5));
t877 = cos(qJ(5));
t740 = t877 * t743 + t840 * t747;
t797 = t840 * qJD(3) + t877 * t805;
t764 = t797 * qJD(5) - t877 * qJDD(3) + t840 * t795;
t803 = qJD(5) + t804;
t776 = mrSges(6,1) * t803 - mrSges(6,3) * t797;
t791 = qJDD(5) - t794;
t796 = -t877 * qJD(3) + t840 * t805;
t768 = pkin(5) * t796 - qJ(6) * t797;
t802 = t803 ^ 2;
t736 = -pkin(5) * t802 + qJ(6) * t791 + 0.2e1 * qJD(6) * t803 - t768 * t796 + t740;
t777 = -mrSges(7,1) * t803 + mrSges(7,2) * t797;
t861 = m(7) * t736 + t791 * mrSges(7,3) + t803 * t777;
t769 = mrSges(7,1) * t796 - mrSges(7,3) * t797;
t866 = -mrSges(6,1) * t796 - mrSges(6,2) * t797 - t769;
t876 = -mrSges(6,3) - mrSges(7,2);
t728 = m(6) * t740 - t791 * mrSges(6,2) + t876 * t764 - t803 * t776 + t866 * t796 + t861;
t739 = -t840 * t743 + t877 * t747;
t765 = -t796 * qJD(5) + t840 * qJDD(3) + t877 * t795;
t775 = -mrSges(6,2) * t803 - mrSges(6,3) * t796;
t737 = -t791 * pkin(5) - t802 * qJ(6) + t797 * t768 + qJDD(6) - t739;
t774 = -mrSges(7,2) * t796 + mrSges(7,3) * t803;
t854 = -m(7) * t737 + t791 * mrSges(7,1) + t803 * t774;
t730 = m(6) * t739 + t791 * mrSges(6,1) + t876 * t765 + t803 * t775 + t866 * t797 + t854;
t723 = t877 * t728 - t730 * t840;
t785 = mrSges(5,1) * t804 + mrSges(5,2) * t805;
t799 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t805;
t718 = m(5) * t745 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t794 - qJD(3) * t799 - t785 * t804 + t723;
t742 = -qJDD(3) * pkin(4) - t845 * pkin(8) + t805 * t786 - t744;
t738 = -0.2e1 * qJD(6) * t797 + (t796 * t803 - t765) * qJ(6) + (t797 * t803 + t764) * pkin(5) + t742;
t734 = m(7) * t738 + mrSges(7,1) * t764 - t765 * mrSges(7,3) + t774 * t796 - t797 * t777;
t731 = -m(6) * t742 - t764 * mrSges(6,1) - mrSges(6,2) * t765 - t796 * t775 - t776 * t797 - t734;
t798 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t804;
t725 = m(5) * t744 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t795 + qJD(3) * t798 - t785 * t805 + t731;
t712 = t836 * t718 + t838 * t725;
t867 = t875 * t796 - t882 * t797 + t874 * t803;
t869 = t873 * t796 + t874 * t797 + t880 * t803;
t719 = -mrSges(6,1) * t742 - mrSges(7,1) * t738 + mrSges(7,2) * t736 + mrSges(6,3) * t740 - pkin(5) * t734 - t881 * t764 + t875 * t765 + t873 * t791 + t869 * t797 - t867 * t803;
t868 = t881 * t796 - t875 * t797 - t873 * t803;
t720 = mrSges(6,2) * t742 + mrSges(7,2) * t737 - mrSges(6,3) * t739 - mrSges(7,3) * t738 - qJ(6) * t734 - t875 * t764 + t882 * t765 - t874 * t791 + t869 * t796 + t868 * t803;
t780 = Ifges(5,4) * t805 - Ifges(5,2) * t804 + Ifges(5,6) * qJD(3);
t781 = Ifges(5,1) * t805 - Ifges(5,4) * t804 + Ifges(5,5) * qJD(3);
t810 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t841 + Ifges(4,2) * t843) * qJD(1);
t811 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t841 + Ifges(4,4) * t843) * qJD(1);
t879 = (t810 * t841 - t811 * t843) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t772 + mrSges(5,1) * t744 - mrSges(4,2) * t773 - mrSges(5,2) * t745 + Ifges(4,5) * t818 + Ifges(5,5) * t795 + Ifges(4,6) * t819 + Ifges(5,6) * t794 + pkin(3) * t712 + pkin(4) * t731 + pkin(8) * t723 + t877 * t719 + t840 * t720 + t805 * t780 + t804 * t781;
t733 = t765 * mrSges(7,2) + t797 * t769 - t854;
t878 = -t873 * t764 - t874 * t765 - t880 * t791 - t867 * t796 - t868 * t797 + mrSges(6,1) * t739 - mrSges(7,1) * t737 - mrSges(6,2) * t740 + mrSges(7,3) * t736 - pkin(5) * t733 + qJ(6) * (-t764 * mrSges(7,2) - t796 * t769 + t861);
t816 = (-mrSges(4,1) * t843 + mrSges(4,2) * t841) * qJD(1);
t864 = qJD(1) * t843;
t823 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t864;
t710 = m(4) * t772 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t818 + qJD(3) * t823 - t816 * t865 + t712;
t822 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t865;
t856 = t838 * t718 - t725 * t836;
t711 = m(4) * t773 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t819 - qJD(3) * t822 + t816 * t864 + t856;
t857 = -t710 * t841 + t843 * t711;
t701 = m(3) * t793 - mrSges(3,1) * t846 - qJDD(1) * mrSges(3,2) + t857;
t722 = t840 * t728 + t877 * t730;
t721 = m(5) * t753 - t794 * mrSges(5,1) + mrSges(5,2) * t795 + t804 * t798 + t799 * t805 + t722;
t782 = -t846 * pkin(7) + t851;
t848 = -m(4) * t782 + t819 * mrSges(4,1) - mrSges(4,2) * t818 - t822 * t865 + t823 * t864 - t721;
t714 = m(3) * t792 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t846 + t848;
t698 = t837 * t701 + t839 * t714;
t695 = m(2) * t824 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t846 + t698;
t858 = t839 * t701 - t714 * t837;
t696 = m(2) * t825 - mrSges(2,1) * t846 - qJDD(1) * mrSges(2,2) + t858;
t870 = t844 * t695 + t842 * t696;
t704 = t843 * t710 + t841 * t711;
t702 = m(3) * t835 + t704;
t859 = -t695 * t842 + t844 * t696;
t779 = Ifges(5,5) * t805 - Ifges(5,6) * t804 + Ifges(5,3) * qJD(3);
t705 = mrSges(5,2) * t753 - mrSges(5,3) * t744 + Ifges(5,1) * t795 + Ifges(5,4) * t794 + Ifges(5,5) * qJDD(3) - pkin(8) * t722 - qJD(3) * t780 - t840 * t719 + t877 * t720 - t804 * t779;
t706 = -mrSges(5,1) * t753 + mrSges(5,3) * t745 + Ifges(5,4) * t795 + Ifges(5,2) * t794 + Ifges(5,6) * qJDD(3) - pkin(4) * t722 + qJD(3) * t781 - t805 * t779 - t878;
t809 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t841 + Ifges(4,6) * t843) * qJD(1);
t689 = -mrSges(4,1) * t782 + mrSges(4,3) * t773 + Ifges(4,4) * t818 + Ifges(4,2) * t819 + Ifges(4,6) * qJDD(3) - pkin(3) * t721 + qJ(4) * t856 + qJD(3) * t811 + t836 * t705 + t838 * t706 - t809 * t865;
t691 = mrSges(4,2) * t782 - mrSges(4,3) * t772 + Ifges(4,1) * t818 + Ifges(4,4) * t819 + Ifges(4,5) * qJDD(3) - qJ(4) * t712 - qJD(3) * t810 + t705 * t838 - t706 * t836 + t809 * t864;
t850 = mrSges(2,1) * t824 + mrSges(3,1) * t792 - mrSges(2,2) * t825 - mrSges(3,2) * t793 + pkin(1) * t698 + pkin(2) * t848 + pkin(7) * t857 + t843 * t689 + t841 * t691 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t687 = -mrSges(3,1) * t835 + mrSges(3,3) * t793 + t846 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t704 - t879;
t686 = mrSges(3,2) * t835 - mrSges(3,3) * t792 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t846 - pkin(7) * t704 - t689 * t841 + t691 * t843;
t685 = -mrSges(2,2) * g(3) - mrSges(2,3) * t824 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t846 - qJ(2) * t698 + t686 * t839 - t687 * t837;
t684 = mrSges(2,1) * g(3) + mrSges(2,3) * t825 + t846 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t702 + qJ(2) * t858 + t837 * t686 + t839 * t687;
t1 = [-m(1) * g(1) + t859; -m(1) * g(2) + t870; (-m(1) - m(2)) * g(3) + t702; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t870 - t842 * t684 + t844 * t685; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t859 + t844 * t684 + t842 * t685; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t850; t850; t702; t879; t721; t878; t733;];
tauJB  = t1;
