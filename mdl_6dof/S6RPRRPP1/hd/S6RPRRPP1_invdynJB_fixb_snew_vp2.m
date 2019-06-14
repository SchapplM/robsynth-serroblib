% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:15:56
% EndTime: 2019-05-05 21:16:05
% DurationCPUTime: 7.73s
% Computational Cost: add. (120333->319), mult. (238135->387), div. (0->0), fcn. (153893->10), ass. (0->130)
t861 = Ifges(6,1) + Ifges(7,1);
t854 = Ifges(6,4) - Ifges(7,5);
t853 = Ifges(6,5) + Ifges(7,4);
t860 = -Ifges(6,2) - Ifges(7,3);
t859 = -Ifges(7,2) - Ifges(6,3);
t852 = Ifges(6,6) - Ifges(7,6);
t822 = sin(qJ(1));
t825 = cos(qJ(1));
t807 = t822 * g(1) - t825 * g(2);
t798 = qJDD(1) * pkin(1) + t807;
t808 = -t825 * g(1) - t822 * g(2);
t827 = qJD(1) ^ 2;
t800 = -t827 * pkin(1) + t808;
t818 = sin(pkin(9));
t819 = cos(pkin(9));
t777 = t818 * t798 + t819 * t800;
t763 = -t827 * pkin(2) + qJDD(1) * pkin(7) + t777;
t816 = -g(3) + qJDD(2);
t821 = sin(qJ(3));
t824 = cos(qJ(3));
t753 = -t821 * t763 + t824 * t816;
t801 = (-pkin(3) * t824 - pkin(8) * t821) * qJD(1);
t826 = qJD(3) ^ 2;
t845 = qJD(1) * t821;
t746 = -qJDD(3) * pkin(3) - t826 * pkin(8) + t801 * t845 - t753;
t820 = sin(qJ(4));
t823 = cos(qJ(4));
t797 = t820 * qJD(3) + t823 * t845;
t843 = qJD(1) * qJD(3);
t839 = t824 * t843;
t802 = t821 * qJDD(1) + t839;
t771 = -t797 * qJD(4) + t823 * qJDD(3) - t820 * t802;
t844 = t824 * qJD(1);
t810 = qJD(4) - t844;
t780 = t810 * pkin(4) - t797 * qJ(5);
t796 = t823 * qJD(3) - t820 * t845;
t794 = t796 ^ 2;
t725 = -t771 * pkin(4) - t794 * qJ(5) + t797 * t780 + qJDD(5) + t746;
t772 = t796 * qJD(4) + t820 * qJDD(3) + t823 * t802;
t817 = sin(pkin(10));
t851 = cos(pkin(10));
t740 = -t851 * t771 + t817 * t772;
t741 = t817 * t771 + t851 * t772;
t774 = -t851 * t796 + t817 * t797;
t775 = t817 * t796 + t851 * t797;
t720 = -0.2e1 * qJD(6) * t775 + (t774 * t810 - t741) * qJ(6) + (t775 * t810 + t740) * pkin(5) + t725;
t757 = -t810 * mrSges(7,1) + t775 * mrSges(7,2);
t758 = -t774 * mrSges(7,2) + t810 * mrSges(7,3);
t713 = m(7) * t720 + t740 * mrSges(7,1) - t741 * mrSges(7,3) - t775 * t757 + t774 * t758;
t776 = t819 * t798 - t818 * t800;
t762 = -qJDD(1) * pkin(2) - t827 * pkin(7) - t776;
t840 = t821 * t843;
t803 = t824 * qJDD(1) - t840;
t739 = (-t802 - t839) * pkin(8) + (-t803 + t840) * pkin(3) + t762;
t754 = t824 * t763 + t821 * t816;
t747 = -t826 * pkin(3) + qJDD(3) * pkin(8) + t801 * t844 + t754;
t726 = t823 * t739 - t820 * t747;
t795 = qJDD(4) - t803;
t722 = (t796 * t810 - t772) * qJ(5) + (t796 * t797 + t795) * pkin(4) + t726;
t727 = t820 * t739 + t823 * t747;
t724 = -t794 * pkin(4) + t771 * qJ(5) - t810 * t780 + t727;
t856 = -2 * qJD(5);
t718 = t817 * t722 + t851 * t724 + t774 * t856;
t748 = t774 * pkin(5) - t775 * qJ(6);
t809 = t810 ^ 2;
t715 = -t809 * pkin(5) + t795 * qJ(6) + 0.2e1 * qJD(6) * t810 - t774 * t748 + t718;
t847 = -t854 * t774 + t861 * t775 + t853 * t810;
t849 = t852 * t774 - t853 * t775 + t859 * t810;
t700 = -mrSges(6,1) * t725 - mrSges(7,1) * t720 + mrSges(7,2) * t715 + mrSges(6,3) * t718 - pkin(5) * t713 + t860 * t740 + t854 * t741 + t849 * t775 + t852 * t795 + t847 * t810;
t832 = t851 * t722 - t817 * t724;
t716 = -t795 * pkin(5) - t809 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t748) * t775 - t832;
t717 = t775 * t856 + t832;
t848 = t860 * t774 + t854 * t775 + t852 * t810;
t701 = mrSges(6,2) * t725 + mrSges(7,2) * t716 - mrSges(6,3) * t717 - mrSges(7,3) * t720 - qJ(6) * t713 - t854 * t740 + t861 * t741 + t849 * t774 + t853 * t795 - t848 * t810;
t755 = -t810 * mrSges(6,2) - t774 * mrSges(6,3);
t756 = t810 * mrSges(6,1) - t775 * mrSges(6,3);
t710 = m(6) * t725 + t740 * mrSges(6,1) + t741 * mrSges(6,2) + t774 * t755 + t775 * t756 + t713;
t764 = Ifges(5,5) * t797 + Ifges(5,6) * t796 + Ifges(5,3) * t810;
t766 = Ifges(5,1) * t797 + Ifges(5,4) * t796 + Ifges(5,5) * t810;
t841 = m(7) * t715 + t795 * mrSges(7,3) + t810 * t757;
t749 = t774 * mrSges(7,1) - t775 * mrSges(7,3);
t846 = -t774 * mrSges(6,1) - t775 * mrSges(6,2) - t749;
t855 = -mrSges(6,3) - mrSges(7,2);
t704 = m(6) * t718 - t795 * mrSges(6,2) + t855 * t740 - t810 * t756 + t846 * t774 + t841;
t834 = -m(7) * t716 + t795 * mrSges(7,1) + t810 * t758;
t706 = m(6) * t717 + t795 * mrSges(6,1) + t855 * t741 + t810 * t755 + t846 * t775 + t834;
t835 = t851 * t704 - t817 * t706;
t680 = -mrSges(5,1) * t746 + mrSges(5,3) * t727 + Ifges(5,4) * t772 + Ifges(5,2) * t771 + Ifges(5,6) * t795 - pkin(4) * t710 + qJ(5) * t835 + t851 * t700 + t817 * t701 - t797 * t764 + t810 * t766;
t699 = t817 * t704 + t851 * t706;
t765 = Ifges(5,4) * t797 + Ifges(5,2) * t796 + Ifges(5,6) * t810;
t681 = mrSges(5,2) * t746 - mrSges(5,3) * t726 + Ifges(5,1) * t772 + Ifges(5,4) * t771 + Ifges(5,5) * t795 - qJ(5) * t699 - t817 * t700 + t851 * t701 + t796 * t764 - t810 * t765;
t778 = -t796 * mrSges(5,1) + t797 * mrSges(5,2);
t779 = -t810 * mrSges(5,2) + t796 * mrSges(5,3);
t697 = m(5) * t726 + t795 * mrSges(5,1) - t772 * mrSges(5,3) - t797 * t778 + t810 * t779 + t699;
t781 = t810 * mrSges(5,1) - t797 * mrSges(5,3);
t698 = m(5) * t727 - t795 * mrSges(5,2) + t771 * mrSges(5,3) + t796 * t778 - t810 * t781 + t835;
t695 = -t820 * t697 + t823 * t698;
t709 = -m(5) * t746 + t771 * mrSges(5,1) - t772 * mrSges(5,2) + t796 * t779 - t797 * t781 - t710;
t788 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t821 + Ifges(4,2) * t824) * qJD(1);
t789 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t821 + Ifges(4,4) * t824) * qJD(1);
t858 = mrSges(4,1) * t753 - mrSges(4,2) * t754 + Ifges(4,5) * t802 + Ifges(4,6) * t803 + Ifges(4,3) * qJDD(3) + pkin(3) * t709 + pkin(8) * t695 + t823 * t680 + t820 * t681 + (t821 * t788 - t824 * t789) * qJD(1);
t712 = t741 * mrSges(7,2) + t775 * t749 - t834;
t857 = -t852 * t740 + t853 * t741 + t847 * t774 + t848 * t775 + (Ifges(5,3) - t859) * t795 + mrSges(5,1) * t726 + mrSges(6,1) * t717 - mrSges(7,1) * t716 - mrSges(5,2) * t727 - mrSges(6,2) * t718 + mrSges(7,3) * t715 + Ifges(5,5) * t772 + Ifges(5,6) * t771 + pkin(4) * t699 - pkin(5) * t712 + qJ(6) * (-t740 * mrSges(7,2) - t774 * t749 + t841) + t797 * t765 - t796 * t766;
t799 = (-mrSges(4,1) * t824 + mrSges(4,2) * t821) * qJD(1);
t805 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t845;
t693 = m(4) * t754 - qJDD(3) * mrSges(4,2) + t803 * mrSges(4,3) - qJD(3) * t805 + t799 * t844 + t695;
t806 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t844;
t708 = m(4) * t753 + qJDD(3) * mrSges(4,1) - t802 * mrSges(4,3) + qJD(3) * t806 - t799 * t845 + t709;
t836 = t824 * t693 - t821 * t708;
t684 = m(3) * t777 - t827 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t836;
t694 = t823 * t697 + t820 * t698;
t830 = -m(4) * t762 + t803 * mrSges(4,1) - t802 * mrSges(4,2) - t805 * t845 + t806 * t844 - t694;
t689 = m(3) * t776 + qJDD(1) * mrSges(3,1) - t827 * mrSges(3,2) + t830;
t677 = t818 * t684 + t819 * t689;
t674 = m(2) * t807 + qJDD(1) * mrSges(2,1) - t827 * mrSges(2,2) + t677;
t837 = t819 * t684 - t818 * t689;
t675 = m(2) * t808 - t827 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t837;
t850 = t825 * t674 + t822 * t675;
t687 = t821 * t693 + t824 * t708;
t685 = m(3) * t816 + t687;
t838 = -t822 * t674 + t825 * t675;
t787 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t821 + Ifges(4,6) * t824) * qJD(1);
t670 = mrSges(4,2) * t762 - mrSges(4,3) * t753 + Ifges(4,1) * t802 + Ifges(4,4) * t803 + Ifges(4,5) * qJDD(3) - pkin(8) * t694 - qJD(3) * t788 - t820 * t680 + t823 * t681 + t787 * t844;
t679 = -mrSges(4,1) * t762 + mrSges(4,3) * t754 + Ifges(4,4) * t802 + Ifges(4,2) * t803 + Ifges(4,6) * qJDD(3) - pkin(3) * t694 + qJD(3) * t789 - t787 * t845 - t857;
t831 = mrSges(2,1) * t807 + mrSges(3,1) * t776 - mrSges(2,2) * t808 - mrSges(3,2) * t777 + pkin(1) * t677 + pkin(2) * t830 + pkin(7) * t836 + t821 * t670 + t824 * t679 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t668 = -mrSges(3,1) * t816 + mrSges(3,3) * t777 + t827 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t687 - t858;
t667 = mrSges(3,2) * t816 - mrSges(3,3) * t776 + Ifges(3,5) * qJDD(1) - t827 * Ifges(3,6) - pkin(7) * t687 + t824 * t670 - t821 * t679;
t666 = -mrSges(2,2) * g(3) - mrSges(2,3) * t807 + Ifges(2,5) * qJDD(1) - t827 * Ifges(2,6) - qJ(2) * t677 + t819 * t667 - t818 * t668;
t665 = mrSges(2,1) * g(3) + mrSges(2,3) * t808 + t827 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t685 + qJ(2) * t837 + t818 * t667 + t819 * t668;
t1 = [-m(1) * g(1) + t838; -m(1) * g(2) + t850; (-m(1) - m(2)) * g(3) + t685; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t850 - t822 * t665 + t825 * t666; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t838 + t825 * t665 + t822 * t666; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t831; t831; t685; t858; t857; t710; t712;];
tauJB  = t1;
