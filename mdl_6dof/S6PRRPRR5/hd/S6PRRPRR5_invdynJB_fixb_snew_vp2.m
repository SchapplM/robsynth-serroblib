% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 05:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:23:29
% EndTime: 2019-05-05 05:23:56
% DurationCPUTime: 25.97s
% Computational Cost: add. (441438->341), mult. (921591->438), div. (0->0), fcn. (674058->14), ass. (0->145)
t849 = sin(pkin(11));
t852 = cos(pkin(11));
t838 = g(1) * t849 - g(2) * t852;
t839 = -g(1) * t852 - g(2) * t849;
t847 = -g(3) + qJDD(1);
t850 = sin(pkin(6));
t853 = cos(pkin(6));
t857 = sin(qJ(2));
t861 = cos(qJ(2));
t799 = -t857 * t839 + (t838 * t853 + t847 * t850) * t861;
t884 = t853 * t857;
t885 = t850 * t857;
t800 = t838 * t884 + t861 * t839 + t847 * t885;
t863 = qJD(2) ^ 2;
t795 = -pkin(2) * t863 + qJDD(2) * pkin(8) + t800;
t816 = -t838 * t850 + t847 * t853;
t856 = sin(qJ(3));
t860 = cos(qJ(3));
t790 = t860 * t795 + t856 * t816;
t834 = (-pkin(3) * t860 - qJ(4) * t856) * qJD(2);
t862 = qJD(3) ^ 2;
t881 = qJD(2) * t860;
t774 = -pkin(3) * t862 + qJDD(3) * qJ(4) + t834 * t881 + t790;
t794 = -qJDD(2) * pkin(2) - t863 * pkin(8) - t799;
t880 = qJD(2) * qJD(3);
t879 = t860 * t880;
t836 = qJDD(2) * t856 + t879;
t845 = t856 * t880;
t837 = qJDD(2) * t860 - t845;
t780 = (-t836 - t879) * qJ(4) + (-t837 + t845) * pkin(3) + t794;
t848 = sin(pkin(12));
t851 = cos(pkin(12));
t882 = qJD(2) * t856;
t829 = qJD(3) * t848 + t851 * t882;
t764 = -0.2e1 * qJD(4) * t829 - t848 * t774 + t851 * t780;
t814 = qJDD(3) * t848 + t836 * t851;
t828 = qJD(3) * t851 - t848 * t882;
t755 = (-t828 * t881 - t814) * pkin(9) + (t828 * t829 - t837) * pkin(4) + t764;
t765 = 0.2e1 * qJD(4) * t828 + t851 * t774 + t848 * t780;
t813 = qJDD(3) * t851 - t836 * t848;
t815 = -pkin(4) * t881 - pkin(9) * t829;
t827 = t828 ^ 2;
t757 = -pkin(4) * t827 + pkin(9) * t813 + t815 * t881 + t765;
t855 = sin(qJ(5));
t859 = cos(qJ(5));
t750 = t859 * t755 - t855 * t757;
t807 = t828 * t859 - t829 * t855;
t782 = qJD(5) * t807 + t813 * t855 + t814 * t859;
t808 = t828 * t855 + t829 * t859;
t831 = qJDD(5) - t837;
t844 = qJD(5) - t881;
t748 = (t807 * t844 - t782) * pkin(10) + (t807 * t808 + t831) * pkin(5) + t750;
t751 = t855 * t755 + t859 * t757;
t781 = -qJD(5) * t808 + t813 * t859 - t814 * t855;
t798 = pkin(5) * t844 - pkin(10) * t808;
t806 = t807 ^ 2;
t749 = -pkin(5) * t806 + pkin(10) * t781 - t798 * t844 + t751;
t854 = sin(qJ(6));
t858 = cos(qJ(6));
t747 = t748 * t854 + t749 * t858;
t789 = -t856 * t795 + t816 * t860;
t773 = -qJDD(3) * pkin(3) - qJ(4) * t862 + t834 * t882 + qJDD(4) - t789;
t766 = -pkin(4) * t813 - pkin(9) * t827 + t829 * t815 + t773;
t752 = -pkin(5) * t781 - pkin(10) * t806 + t798 * t808 + t766;
t788 = t807 * t854 + t808 * t858;
t761 = -qJD(6) * t788 + t781 * t858 - t782 * t854;
t787 = t807 * t858 - t808 * t854;
t762 = qJD(6) * t787 + t781 * t854 + t782 * t858;
t843 = qJD(6) + t844;
t767 = Ifges(7,5) * t788 + Ifges(7,6) * t787 + Ifges(7,3) * t843;
t769 = Ifges(7,1) * t788 + Ifges(7,4) * t787 + Ifges(7,5) * t843;
t825 = qJDD(6) + t831;
t735 = -mrSges(7,1) * t752 + mrSges(7,3) * t747 + Ifges(7,4) * t762 + Ifges(7,2) * t761 + Ifges(7,6) * t825 - t767 * t788 + t769 * t843;
t746 = t748 * t858 - t749 * t854;
t768 = Ifges(7,4) * t788 + Ifges(7,2) * t787 + Ifges(7,6) * t843;
t736 = mrSges(7,2) * t752 - mrSges(7,3) * t746 + Ifges(7,1) * t762 + Ifges(7,4) * t761 + Ifges(7,5) * t825 + t767 * t787 - t768 * t843;
t783 = Ifges(6,5) * t808 + Ifges(6,6) * t807 + Ifges(6,3) * t844;
t785 = Ifges(6,1) * t808 + Ifges(6,4) * t807 + Ifges(6,5) * t844;
t777 = -mrSges(7,2) * t843 + mrSges(7,3) * t787;
t778 = mrSges(7,1) * t843 - mrSges(7,3) * t788;
t869 = m(7) * t752 - t761 * mrSges(7,1) + t762 * mrSges(7,2) - t787 * t777 + t788 * t778;
t771 = -mrSges(7,1) * t787 + mrSges(7,2) * t788;
t740 = m(7) * t746 + mrSges(7,1) * t825 - mrSges(7,3) * t762 - t771 * t788 + t777 * t843;
t741 = m(7) * t747 - mrSges(7,2) * t825 + mrSges(7,3) * t761 + t771 * t787 - t778 * t843;
t875 = -t740 * t854 + t858 * t741;
t723 = -mrSges(6,1) * t766 + mrSges(6,3) * t751 + Ifges(6,4) * t782 + Ifges(6,2) * t781 + Ifges(6,6) * t831 - pkin(5) * t869 + pkin(10) * t875 + t858 * t735 + t854 * t736 - t808 * t783 + t844 * t785;
t734 = t858 * t740 + t854 * t741;
t784 = Ifges(6,4) * t808 + Ifges(6,2) * t807 + Ifges(6,6) * t844;
t724 = mrSges(6,2) * t766 - mrSges(6,3) * t750 + Ifges(6,1) * t782 + Ifges(6,4) * t781 + Ifges(6,5) * t831 - pkin(10) * t734 - t735 * t854 + t736 * t858 + t783 * t807 - t784 * t844;
t801 = Ifges(5,5) * t829 + Ifges(5,6) * t828 - Ifges(5,3) * t881;
t803 = Ifges(5,1) * t829 + Ifges(5,4) * t828 - Ifges(5,5) * t881;
t796 = -mrSges(6,2) * t844 + mrSges(6,3) * t807;
t797 = mrSges(6,1) * t844 - mrSges(6,3) * t808;
t866 = m(6) * t766 - t781 * mrSges(6,1) + t782 * mrSges(6,2) - t807 * t796 + t808 * t797 + t869;
t791 = -mrSges(6,1) * t807 + mrSges(6,2) * t808;
t732 = m(6) * t750 + mrSges(6,1) * t831 - mrSges(6,3) * t782 - t791 * t808 + t796 * t844 + t734;
t733 = m(6) * t751 - mrSges(6,2) * t831 + mrSges(6,3) * t781 + t791 * t807 - t797 * t844 + t875;
t876 = -t732 * t855 + t859 * t733;
t706 = -mrSges(5,1) * t773 + mrSges(5,3) * t765 + Ifges(5,4) * t814 + Ifges(5,2) * t813 - Ifges(5,6) * t837 - pkin(4) * t866 + pkin(9) * t876 + t859 * t723 + t855 * t724 - t829 * t801 - t803 * t881;
t728 = t859 * t732 + t855 * t733;
t802 = Ifges(5,4) * t829 + Ifges(5,2) * t828 - Ifges(5,6) * t881;
t707 = mrSges(5,2) * t773 - mrSges(5,3) * t764 + Ifges(5,1) * t814 + Ifges(5,4) * t813 - Ifges(5,5) * t837 - pkin(9) * t728 - t723 * t855 + t724 * t859 + t801 * t828 + t802 * t881;
t809 = -mrSges(5,1) * t828 + mrSges(5,2) * t829;
t871 = mrSges(5,2) * t881 + mrSges(5,3) * t828;
t726 = m(5) * t764 - t837 * mrSges(5,1) - t814 * mrSges(5,3) - t829 * t809 - t871 * t881 + t728;
t812 = -mrSges(5,1) * t881 - mrSges(5,3) * t829;
t727 = m(5) * t765 + mrSges(5,2) * t837 + mrSges(5,3) * t813 + t809 * t828 + t812 * t881 + t876;
t722 = -t726 * t848 + t851 * t727;
t744 = m(5) * t773 - t813 * mrSges(5,1) + t814 * mrSges(5,2) + t829 * t812 - t828 * t871 + t866;
t823 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t856 + Ifges(4,2) * t860) * qJD(2);
t824 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t856 + Ifges(4,4) * t860) * qJD(2);
t887 = mrSges(4,1) * t789 - mrSges(4,2) * t790 + Ifges(4,5) * t836 + Ifges(4,6) * t837 + Ifges(4,3) * qJDD(3) - pkin(3) * t744 + qJ(4) * t722 + t851 * t706 + t848 * t707 + (t823 * t856 - t824 * t860) * qJD(2);
t721 = t726 * t851 + t727 * t848;
t840 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t882;
t841 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t881;
t867 = -m(4) * t794 + t837 * mrSges(4,1) - mrSges(4,2) * t836 - t840 * t882 + t841 * t881 - t721;
t717 = m(3) * t799 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t863 + t867;
t886 = t717 * t861;
t835 = (-mrSges(4,1) * t860 + mrSges(4,2) * t856) * qJD(2);
t720 = m(4) * t790 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t837 - qJD(3) * t840 + t835 * t881 + t722;
t743 = m(4) * t789 + qJDD(3) * mrSges(4,1) - t836 * mrSges(4,3) + qJD(3) * t841 - t835 * t882 - t744;
t877 = t860 * t720 - t743 * t856;
t711 = m(3) * t800 - mrSges(3,1) * t863 - qJDD(2) * mrSges(3,2) + t877;
t714 = t856 * t720 + t860 * t743;
t713 = m(3) * t816 + t714;
t700 = t711 * t884 - t713 * t850 + t853 * t886;
t698 = m(2) * t838 + t700;
t704 = t861 * t711 - t717 * t857;
t703 = m(2) * t839 + t704;
t883 = t852 * t698 + t849 * t703;
t699 = t711 * t885 + t853 * t713 + t850 * t886;
t878 = -t698 * t849 + t852 * t703;
t874 = m(2) * t847 + t699;
t822 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t856 + Ifges(4,6) * t860) * qJD(2);
t696 = mrSges(4,2) * t794 - mrSges(4,3) * t789 + Ifges(4,1) * t836 + Ifges(4,4) * t837 + Ifges(4,5) * qJDD(3) - qJ(4) * t721 - qJD(3) * t823 - t706 * t848 + t707 * t851 + t822 * t881;
t868 = -mrSges(7,1) * t746 + mrSges(7,2) * t747 - Ifges(7,5) * t762 - Ifges(7,6) * t761 - Ifges(7,3) * t825 - t788 * t768 + t787 * t769;
t864 = mrSges(6,1) * t750 - mrSges(6,2) * t751 + Ifges(6,5) * t782 + Ifges(6,6) * t781 + Ifges(6,3) * t831 + pkin(5) * t734 + t808 * t784 - t807 * t785 - t868;
t705 = Ifges(4,6) * qJDD(3) + (Ifges(4,2) + Ifges(5,3)) * t837 - t864 + Ifges(4,4) * t836 + t828 * t803 - t829 * t802 - Ifges(5,6) * t813 - Ifges(5,5) * t814 + qJD(3) * t824 + mrSges(4,3) * t790 - mrSges(4,1) * t794 - mrSges(5,1) * t764 + mrSges(5,2) * t765 - t822 * t882 - pkin(4) * t728 - pkin(3) * t721;
t694 = mrSges(3,2) * t816 - mrSges(3,3) * t799 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t863 - pkin(8) * t714 + t696 * t860 - t705 * t856;
t695 = -mrSges(3,1) * t816 + mrSges(3,3) * t800 + t863 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t714 - t887;
t870 = pkin(7) * t704 + t694 * t857 + t695 * t861;
t693 = mrSges(3,1) * t799 - mrSges(3,2) * t800 + Ifges(3,3) * qJDD(2) + pkin(2) * t867 + pkin(8) * t877 + t856 * t696 + t860 * t705;
t692 = mrSges(2,2) * t847 - mrSges(2,3) * t838 + t861 * t694 - t857 * t695 + (-t699 * t850 - t700 * t853) * pkin(7);
t691 = -mrSges(2,1) * t847 + mrSges(2,3) * t839 - pkin(1) * t699 - t850 * t693 + t853 * t870;
t1 = [-m(1) * g(1) + t878; -m(1) * g(2) + t883; -m(1) * g(3) + t874; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t883 - t849 * t691 + t852 * t692; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t878 + t852 * t691 + t849 * t692; -mrSges(1,1) * g(2) + mrSges(2,1) * t838 + mrSges(1,2) * g(1) - mrSges(2,2) * t839 + pkin(1) * t700 + t853 * t693 + t850 * t870; t874; t693; t887; t744; t864; -t868;];
tauJB  = t1;
