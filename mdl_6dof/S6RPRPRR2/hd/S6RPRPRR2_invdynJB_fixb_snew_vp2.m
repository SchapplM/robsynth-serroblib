% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:23:14
% EndTime: 2019-05-05 18:23:28
% DurationCPUTime: 13.81s
% Computational Cost: add. (224687->344), mult. (487688->433), div. (0->0), fcn. (331778->12), ass. (0->138)
t855 = sin(qJ(1));
t859 = cos(qJ(1));
t837 = t855 * g(1) - g(2) * t859;
t828 = qJDD(1) * pkin(1) + t837;
t838 = -g(1) * t859 - g(2) * t855;
t861 = qJD(1) ^ 2;
t830 = -pkin(1) * t861 + t838;
t849 = sin(pkin(10));
t851 = cos(pkin(10));
t804 = t849 * t828 + t851 * t830;
t794 = -pkin(2) * t861 + qJDD(1) * pkin(7) + t804;
t847 = -g(3) + qJDD(2);
t854 = sin(qJ(3));
t858 = cos(qJ(3));
t784 = -t854 * t794 + t858 * t847;
t877 = qJD(1) * qJD(3);
t876 = t858 * t877;
t831 = qJDD(1) * t854 + t876;
t768 = (-t831 + t876) * qJ(4) + (t854 * t858 * t861 + qJDD(3)) * pkin(3) + t784;
t785 = t858 * t794 + t854 * t847;
t832 = qJDD(1) * t858 - t854 * t877;
t880 = qJD(1) * t854;
t834 = qJD(3) * pkin(3) - qJ(4) * t880;
t846 = t858 ^ 2;
t771 = -pkin(3) * t846 * t861 + qJ(4) * t832 - qJD(3) * t834 + t785;
t848 = sin(pkin(11));
t850 = cos(pkin(11));
t817 = (t848 * t858 + t850 * t854) * qJD(1);
t755 = -0.2e1 * qJD(4) * t817 + t768 * t850 - t848 * t771;
t879 = qJD(1) * t858;
t816 = -t848 * t880 + t850 * t879;
t756 = 0.2e1 * qJD(4) * t816 + t848 * t768 + t850 * t771;
t798 = -pkin(4) * t816 - pkin(8) * t817;
t860 = qJD(3) ^ 2;
t748 = -pkin(4) * t860 + qJDD(3) * pkin(8) + t798 * t816 + t756;
t803 = t851 * t828 - t849 * t830;
t868 = -qJDD(1) * pkin(2) - t803;
t772 = -t832 * pkin(3) + qJDD(4) + t834 * t880 + (-qJ(4) * t846 - pkin(7)) * t861 + t868;
t805 = -t848 * t831 + t832 * t850;
t806 = t831 * t850 + t832 * t848;
t759 = (-qJD(3) * t816 - t806) * pkin(8) + (qJD(3) * t817 - t805) * pkin(4) + t772;
t853 = sin(qJ(5));
t857 = cos(qJ(5));
t743 = -t853 * t748 + t857 * t759;
t808 = qJD(3) * t857 - t817 * t853;
t779 = qJD(5) * t808 + qJDD(3) * t853 + t806 * t857;
t802 = qJDD(5) - t805;
t809 = qJD(3) * t853 + t817 * t857;
t815 = qJD(5) - t816;
t741 = (t808 * t815 - t779) * pkin(9) + (t808 * t809 + t802) * pkin(5) + t743;
t744 = t857 * t748 + t853 * t759;
t778 = -qJD(5) * t809 + qJDD(3) * t857 - t806 * t853;
t788 = pkin(5) * t815 - pkin(9) * t809;
t807 = t808 ^ 2;
t742 = -pkin(5) * t807 + pkin(9) * t778 - t788 * t815 + t744;
t852 = sin(qJ(6));
t856 = cos(qJ(6));
t739 = t741 * t856 - t742 * t852;
t780 = t808 * t856 - t809 * t852;
t753 = qJD(6) * t780 + t778 * t852 + t779 * t856;
t781 = t808 * t852 + t809 * t856;
t764 = -mrSges(7,1) * t780 + mrSges(7,2) * t781;
t812 = qJD(6) + t815;
t769 = -mrSges(7,2) * t812 + mrSges(7,3) * t780;
t799 = qJDD(6) + t802;
t735 = m(7) * t739 + mrSges(7,1) * t799 - mrSges(7,3) * t753 - t764 * t781 + t769 * t812;
t740 = t741 * t852 + t742 * t856;
t752 = -qJD(6) * t781 + t778 * t856 - t779 * t852;
t770 = mrSges(7,1) * t812 - mrSges(7,3) * t781;
t736 = m(7) * t740 - mrSges(7,2) * t799 + mrSges(7,3) * t752 + t764 * t780 - t770 * t812;
t727 = t856 * t735 + t852 * t736;
t782 = -mrSges(6,1) * t808 + mrSges(6,2) * t809;
t786 = -mrSges(6,2) * t815 + mrSges(6,3) * t808;
t725 = m(6) * t743 + mrSges(6,1) * t802 - mrSges(6,3) * t779 - t782 * t809 + t786 * t815 + t727;
t787 = mrSges(6,1) * t815 - mrSges(6,3) * t809;
t871 = -t735 * t852 + t856 * t736;
t726 = m(6) * t744 - mrSges(6,2) * t802 + mrSges(6,3) * t778 + t782 * t808 - t787 * t815 + t871;
t721 = -t725 * t853 + t857 * t726;
t796 = -mrSges(5,1) * t816 + mrSges(5,2) * t817;
t811 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t817;
t718 = m(5) * t756 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t805 - qJD(3) * t811 + t796 * t816 + t721;
t747 = -qJDD(3) * pkin(4) - pkin(8) * t860 + t817 * t798 - t755;
t745 = -pkin(5) * t778 - pkin(9) * t807 + t788 * t809 + t747;
t867 = m(7) * t745 - t752 * mrSges(7,1) + mrSges(7,2) * t753 - t780 * t769 + t770 * t781;
t737 = -m(6) * t747 + t778 * mrSges(6,1) - mrSges(6,2) * t779 + t808 * t786 - t787 * t809 - t867;
t810 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t816;
t731 = m(5) * t755 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t806 + qJD(3) * t810 - t796 * t817 + t737;
t710 = t848 * t718 + t850 * t731;
t760 = Ifges(7,5) * t781 + Ifges(7,6) * t780 + Ifges(7,3) * t812;
t762 = Ifges(7,1) * t781 + Ifges(7,4) * t780 + Ifges(7,5) * t812;
t728 = -mrSges(7,1) * t745 + mrSges(7,3) * t740 + Ifges(7,4) * t753 + Ifges(7,2) * t752 + Ifges(7,6) * t799 - t760 * t781 + t762 * t812;
t761 = Ifges(7,4) * t781 + Ifges(7,2) * t780 + Ifges(7,6) * t812;
t729 = mrSges(7,2) * t745 - mrSges(7,3) * t739 + Ifges(7,1) * t753 + Ifges(7,4) * t752 + Ifges(7,5) * t799 + t760 * t780 - t761 * t812;
t773 = Ifges(6,5) * t809 + Ifges(6,6) * t808 + Ifges(6,3) * t815;
t775 = Ifges(6,1) * t809 + Ifges(6,4) * t808 + Ifges(6,5) * t815;
t711 = -mrSges(6,1) * t747 + mrSges(6,3) * t744 + Ifges(6,4) * t779 + Ifges(6,2) * t778 + Ifges(6,6) * t802 - pkin(5) * t867 + pkin(9) * t871 + t856 * t728 + t852 * t729 - t809 * t773 + t815 * t775;
t774 = Ifges(6,4) * t809 + Ifges(6,2) * t808 + Ifges(6,6) * t815;
t712 = mrSges(6,2) * t747 - mrSges(6,3) * t743 + Ifges(6,1) * t779 + Ifges(6,4) * t778 + Ifges(6,5) * t802 - pkin(9) * t727 - t728 * t852 + t729 * t856 + t773 * t808 - t774 * t815;
t791 = Ifges(5,4) * t817 + Ifges(5,2) * t816 + Ifges(5,6) * qJD(3);
t792 = Ifges(5,1) * t817 + Ifges(5,4) * t816 + Ifges(5,5) * qJD(3);
t822 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t854 + Ifges(4,2) * t858) * qJD(1);
t823 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t854 + Ifges(4,4) * t858) * qJD(1);
t883 = (t822 * t854 - t823 * t858) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t784 + mrSges(5,1) * t755 - mrSges(4,2) * t785 - mrSges(5,2) * t756 + Ifges(4,5) * t831 + Ifges(5,5) * t806 + Ifges(4,6) * t832 + Ifges(5,6) * t805 + pkin(3) * t710 + pkin(4) * t737 + pkin(8) * t721 + t857 * t711 + t853 * t712 + t817 * t791 - t816 * t792;
t829 = (-mrSges(4,1) * t858 + mrSges(4,2) * t854) * qJD(1);
t836 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t879;
t708 = m(4) * t784 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t831 + qJD(3) * t836 - t829 * t880 + t710;
t835 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t880;
t872 = t850 * t718 - t731 * t848;
t709 = m(4) * t785 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t832 - qJD(3) * t835 + t829 * t879 + t872;
t873 = -t708 * t854 + t858 * t709;
t700 = m(3) * t804 - mrSges(3,1) * t861 - qJDD(1) * mrSges(3,2) + t873;
t720 = t857 * t725 + t853 * t726;
t719 = m(5) * t772 - t805 * mrSges(5,1) + mrSges(5,2) * t806 - t816 * t810 + t811 * t817 + t720;
t793 = -t861 * pkin(7) + t868;
t864 = -m(4) * t793 + t832 * mrSges(4,1) - mrSges(4,2) * t831 - t835 * t880 + t836 * t879 - t719;
t714 = m(3) * t803 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t861 + t864;
t696 = t849 * t700 + t851 * t714;
t693 = m(2) * t837 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t861 + t696;
t874 = t851 * t700 - t714 * t849;
t694 = m(2) * t838 - mrSges(2,1) * t861 - qJDD(1) * mrSges(2,2) + t874;
t881 = t859 * t693 + t855 * t694;
t703 = t858 * t708 + t854 * t709;
t701 = m(3) * t847 + t703;
t875 = -t693 * t855 + t859 * t694;
t866 = -mrSges(7,1) * t739 + mrSges(7,2) * t740 - Ifges(7,5) * t753 - Ifges(7,6) * t752 - Ifges(7,3) * t799 - t781 * t761 + t780 * t762;
t790 = Ifges(5,5) * t817 + Ifges(5,6) * t816 + Ifges(5,3) * qJD(3);
t697 = mrSges(5,2) * t772 - mrSges(5,3) * t755 + Ifges(5,1) * t806 + Ifges(5,4) * t805 + Ifges(5,5) * qJDD(3) - pkin(8) * t720 - qJD(3) * t791 - t711 * t853 + t712 * t857 + t790 * t816;
t863 = mrSges(6,1) * t743 - mrSges(6,2) * t744 + Ifges(6,5) * t779 + Ifges(6,6) * t778 + Ifges(6,3) * t802 + pkin(5) * t727 + t809 * t774 - t808 * t775 - t866;
t704 = -mrSges(5,1) * t772 + mrSges(5,3) * t756 + Ifges(5,4) * t806 + Ifges(5,2) * t805 + Ifges(5,6) * qJDD(3) - pkin(4) * t720 + qJD(3) * t792 - t817 * t790 - t863;
t821 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t854 + Ifges(4,6) * t858) * qJD(1);
t686 = -mrSges(4,1) * t793 + mrSges(4,3) * t785 + Ifges(4,4) * t831 + Ifges(4,2) * t832 + Ifges(4,6) * qJDD(3) - pkin(3) * t719 + qJ(4) * t872 + qJD(3) * t823 + t848 * t697 + t850 * t704 - t821 * t880;
t689 = mrSges(4,2) * t793 - mrSges(4,3) * t784 + Ifges(4,1) * t831 + Ifges(4,4) * t832 + Ifges(4,5) * qJDD(3) - qJ(4) * t710 - qJD(3) * t822 + t697 * t850 - t704 * t848 + t821 * t879;
t865 = mrSges(2,1) * t837 + mrSges(3,1) * t803 - mrSges(2,2) * t838 - mrSges(3,2) * t804 + pkin(1) * t696 + pkin(2) * t864 + pkin(7) * t873 + t858 * t686 + t854 * t689 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t687 = -mrSges(3,1) * t847 + mrSges(3,3) * t804 + t861 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t703 - t883;
t684 = mrSges(3,2) * t847 - mrSges(3,3) * t803 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t861 - pkin(7) * t703 - t686 * t854 + t689 * t858;
t683 = -mrSges(2,2) * g(3) - mrSges(2,3) * t837 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t861 - qJ(2) * t696 + t684 * t851 - t687 * t849;
t682 = mrSges(2,1) * g(3) + mrSges(2,3) * t838 + t861 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t701 + qJ(2) * t874 + t849 * t684 + t851 * t687;
t1 = [-m(1) * g(1) + t875; -m(1) * g(2) + t881; (-m(1) - m(2)) * g(3) + t701; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t881 - t855 * t682 + t859 * t683; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t875 + t859 * t682 + t855 * t683; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t865; t865; t701; t883; t719; t863; -t866;];
tauJB  = t1;
