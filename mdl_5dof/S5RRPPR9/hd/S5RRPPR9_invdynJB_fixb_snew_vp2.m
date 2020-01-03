% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:53
% EndTime: 2019-12-31 19:40:55
% DurationCPUTime: 1.95s
% Computational Cost: add. (12613->275), mult. (26133->322), div. (0->0), fcn. (12107->6), ass. (0->114)
t873 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t849 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t848 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t872 = Ifges(3,2) + Ifges(5,1) + Ifges(4,3);
t847 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t871 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t820 = sin(qJ(2));
t823 = cos(qJ(2));
t782 = (-mrSges(4,1) * t823 - mrSges(4,3) * t820) * qJD(1);
t851 = qJD(1) * qJD(2);
t841 = t823 * t851;
t786 = t820 * qJDD(1) + t841;
t853 = qJD(1) * t823;
t796 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t853;
t842 = t820 * t851;
t787 = t823 * qJDD(1) - t842;
t852 = t820 * qJD(1);
t792 = -qJD(2) * pkin(3) - qJ(4) * t852;
t821 = sin(qJ(1));
t824 = cos(qJ(1));
t799 = t821 * g(1) - t824 * g(2);
t826 = qJD(1) ^ 2;
t762 = -qJDD(1) * pkin(1) - t826 * pkin(6) - t799;
t836 = -t787 * pkin(2) + t762 + (-t786 - t841) * qJ(3);
t856 = t823 ^ 2 * t826;
t863 = 2 * qJD(3);
t830 = -qJ(4) * t856 + qJDD(4) - t836 + (t792 + t863) * t852;
t861 = pkin(3) + pkin(7);
t862 = -pkin(2) - pkin(7);
t720 = t861 * t787 + t786 * pkin(4) + (pkin(4) * t823 + t862 * t820) * t851 + t830;
t785 = (pkin(4) * t820 + pkin(7) * t823) * qJD(1);
t825 = qJD(2) ^ 2;
t800 = -t824 * g(1) - t821 * g(2);
t763 = -t826 * pkin(1) + qJDD(1) * pkin(6) + t800;
t745 = -t823 * g(3) - t820 * t763;
t781 = (-pkin(2) * t823 - qJ(3) * t820) * qJD(1);
t838 = t781 * t852 + qJDD(3) - t745;
t855 = t823 * t826;
t850 = qJD(1) * qJD(4);
t868 = -0.2e1 * t820 * t850 + (-t786 + t841) * qJ(4);
t723 = (-pkin(4) - qJ(3)) * t825 + (-pkin(3) * t855 - qJD(1) * t785) * t820 + (-pkin(2) - t861) * qJDD(2) + t838 + t868;
t819 = sin(qJ(5));
t822 = cos(qJ(5));
t718 = t822 * t720 - t819 * t723;
t779 = -t822 * qJD(2) + t819 * t853;
t741 = t779 * qJD(5) - t819 * qJDD(2) - t822 * t787;
t780 = -t819 * qJD(2) - t822 * t853;
t742 = -t779 * mrSges(6,1) + t780 * mrSges(6,2);
t803 = qJD(5) + t852;
t743 = -t803 * mrSges(6,2) + t779 * mrSges(6,3);
t777 = qJDD(5) + t786;
t714 = m(6) * t718 + t777 * mrSges(6,1) - t741 * mrSges(6,3) - t780 * t742 + t803 * t743;
t719 = t819 * t720 + t822 * t723;
t740 = -t780 * qJD(5) - t822 * qJDD(2) + t819 * t787;
t744 = t803 * mrSges(6,1) - t780 * mrSges(6,3);
t715 = m(6) * t719 - t777 * mrSges(6,2) + t740 * mrSges(6,3) + t779 * t742 - t803 * t744;
t704 = -t819 * t714 + t822 * t715;
t732 = -qJDD(2) * pkin(2) - t825 * qJ(3) + t838;
t727 = (-t820 * t855 - qJDD(2)) * pkin(3) + t732 + t868;
t784 = (mrSges(5,1) * t820 - mrSges(5,2) * t823) * qJD(1);
t837 = -m(5) * t727 + t784 * t852 - t704;
t832 = qJDD(2) * mrSges(5,2) + qJD(2) * t796 - t837;
t798 = mrSges(4,2) * t853 + qJD(2) * mrSges(4,3);
t867 = m(4) * t732 - qJDD(2) * mrSges(4,1) - qJD(2) * t798;
t700 = t782 * t852 + (mrSges(4,2) - mrSges(5,3)) * t786 + t832 + t867;
t702 = -t786 * mrSges(5,3) + t832;
t807 = qJD(2) * t863;
t746 = -t820 * g(3) + t823 * t763;
t869 = qJDD(2) * qJ(3) + t781 * t853 + t746;
t833 = pkin(3) * t856 + t787 * qJ(4) - t869;
t722 = qJDD(2) * pkin(4) + qJD(2) * t792 + t807 + t862 * t825 + (-0.2e1 * qJD(4) - t785) * t853 - t833;
t733 = Ifges(6,5) * t780 + Ifges(6,6) * t779 + Ifges(6,3) * t803;
t735 = Ifges(6,1) * t780 + Ifges(6,4) * t779 + Ifges(6,5) * t803;
t708 = -mrSges(6,1) * t722 + mrSges(6,3) * t719 + Ifges(6,4) * t741 + Ifges(6,2) * t740 + Ifges(6,6) * t777 - t780 * t733 + t803 * t735;
t734 = Ifges(6,4) * t780 + Ifges(6,2) * t779 + Ifges(6,6) * t803;
t709 = mrSges(6,2) * t722 - mrSges(6,3) * t718 + Ifges(6,1) * t741 + Ifges(6,4) * t740 + Ifges(6,5) * t777 + t779 * t733 - t803 * t734;
t716 = -m(6) * t722 + t740 * mrSges(6,1) - t741 * mrSges(6,2) + t779 * t743 - t780 * t744;
t859 = t825 * pkin(2);
t864 = -2 * qJD(3);
t726 = 0.2e1 * t823 * t850 + t859 + (t864 - t792) * qJD(2) + t833;
t731 = t807 - t859 + t869;
t795 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t852;
t793 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t852;
t831 = -m(5) * t726 + qJDD(2) * mrSges(5,1) - t787 * mrSges(5,3) + qJD(2) * t793 - t716;
t829 = m(4) * t731 + qJDD(2) * mrSges(4,3) + qJD(2) * t795 + t782 * t853 + t831;
t843 = -t848 * qJD(2) + (-t873 * t820 - t849 * t823) * qJD(1);
t844 = -t847 * qJD(2) + (-t849 * t820 - t872 * t823) * qJD(1);
t870 = -(t844 * t820 - t843 * t823) * qJD(1) + t871 * qJDD(2) + t848 * t786 + t847 * t787 + mrSges(3,1) * t745 - mrSges(4,1) * t732 - mrSges(5,1) * t726 - mrSges(3,2) * t746 + mrSges(5,2) * t727 + mrSges(4,3) * t731 - pkin(2) * t700 - pkin(3) * t702 - pkin(4) * t716 - pkin(7) * t704 + qJ(3) * (t787 * mrSges(4,2) - t784 * t853 + t829) - t822 * t708 - t819 * t709;
t858 = mrSges(3,3) + mrSges(4,2);
t783 = (-mrSges(3,1) * t823 + mrSges(3,2) * t820) * qJD(1);
t797 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t853;
t698 = m(3) * t745 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t796 + t797) * qJD(2) + (-t782 - t783) * t852 + (mrSges(5,3) - t858) * t786 + t837 - t867;
t794 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t852;
t707 = -qJD(2) * t794 + m(3) * t746 + t829 + t858 * t787 + (t783 - t784) * t853 - qJDD(2) * mrSges(3,2);
t839 = -t820 * t698 + t823 * t707;
t691 = m(2) * t800 - t826 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t839;
t703 = t822 * t714 + t819 * t715;
t725 = -pkin(2) * t842 + t787 * pkin(3) + t830;
t701 = m(5) * t725 + t786 * mrSges(5,1) - t787 * mrSges(5,2) + t793 * t852 - t796 * t853 + t703;
t728 = (pkin(2) * qJD(2) + t864) * t852 + t836;
t699 = m(4) * t728 - t787 * mrSges(4,1) - t786 * mrSges(4,3) - t795 * t852 - t798 * t853 - t701;
t828 = -m(3) * t762 + t787 * mrSges(3,1) - t786 * mrSges(3,2) - t794 * t852 + t797 * t853 - t699;
t695 = m(2) * t799 + qJDD(1) * mrSges(2,1) - t826 * mrSges(2,2) + t828;
t854 = t821 * t691 + t824 * t695;
t693 = t823 * t698 + t820 * t707;
t845 = -t871 * qJD(2) + (-t848 * t820 - t847 * t823) * qJD(1);
t840 = t824 * t691 - t821 * t695;
t686 = t819 * t708 - t822 * t709 - qJ(4) * t831 - mrSges(3,1) * t762 + mrSges(3,3) * t746 - mrSges(5,2) * t725 + mrSges(5,3) * t726 - mrSges(4,1) * t728 + mrSges(4,2) * t731 + pkin(7) * t703 + pkin(3) * t701 - pkin(2) * t699 + t872 * t787 + t849 * t786 + t847 * qJDD(2) - t843 * qJD(2) + (qJ(4) * t784 * t823 + t845 * t820) * qJD(1);
t834 = mrSges(6,1) * t718 - mrSges(6,2) * t719 + Ifges(6,5) * t741 + Ifges(6,6) * t740 + Ifges(6,3) * t777 + t780 * t734 - t779 * t735;
t688 = mrSges(5,1) * t725 + mrSges(3,2) * t762 + mrSges(4,2) * t732 - mrSges(3,3) * t745 - mrSges(4,3) * t728 - mrSges(5,3) * t727 + pkin(4) * t703 - qJ(3) * t699 - qJ(4) * t702 + t844 * qJD(2) + t848 * qJDD(2) + t873 * t786 + t849 * t787 - t845 * t853 + t834;
t835 = mrSges(2,1) * t799 - mrSges(2,2) * t800 + Ifges(2,3) * qJDD(1) + pkin(1) * t828 + pkin(6) * t839 + t823 * t686 + t820 * t688;
t684 = mrSges(2,1) * g(3) + mrSges(2,3) * t800 + t826 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t693 - t870;
t683 = -mrSges(2,2) * g(3) - mrSges(2,3) * t799 + Ifges(2,5) * qJDD(1) - t826 * Ifges(2,6) - pkin(6) * t693 - t820 * t686 + t823 * t688;
t1 = [-m(1) * g(1) + t840; -m(1) * g(2) + t854; (-m(1) - m(2)) * g(3) + t693; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t854 + t824 * t683 - t821 * t684; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t840 + t821 * t683 + t824 * t684; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t835; t835; t870; t700; t701; t834;];
tauJB = t1;
