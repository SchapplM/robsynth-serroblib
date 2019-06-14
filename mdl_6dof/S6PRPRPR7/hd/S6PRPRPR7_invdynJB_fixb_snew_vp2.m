% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-05-04 23:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:19:14
% EndTime: 2019-05-04 23:19:18
% DurationCPUTime: 3.49s
% Computational Cost: add. (36713->282), mult. (68910->341), div. (0->0), fcn. (38853->10), ass. (0->126)
t843 = Ifges(5,5) - Ifges(6,4);
t842 = Ifges(5,6) - Ifges(6,5);
t841 = (Ifges(5,3) + Ifges(6,1));
t787 = sin(pkin(10));
t789 = cos(pkin(10));
t763 = t787 * g(1) - t789 * g(2);
t764 = -t789 * g(1) - t787 * g(2);
t784 = -g(3) + qJDD(1);
t796 = cos(qJ(2));
t790 = cos(pkin(6));
t793 = sin(qJ(2));
t829 = t790 * t793;
t788 = sin(pkin(6));
t830 = t788 * t793;
t715 = t763 * t829 + t796 * t764 + t784 * t830;
t814 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t715;
t798 = qJD(2) ^ 2;
t837 = (-pkin(2) - pkin(8));
t823 = t837 * t798;
t709 = t823 + t814;
t792 = sin(qJ(4));
t795 = cos(qJ(4));
t825 = qJD(2) * qJD(4);
t822 = t795 * t825;
t760 = t792 * qJDD(2) + t822;
t775 = t792 * t825;
t761 = t795 * qJDD(2) - t775;
t827 = qJD(2) * t792;
t765 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t827;
t826 = t795 * qJD(2);
t766 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t826;
t838 = -2 * qJD(5);
t805 = pkin(4) * t822 + t826 * t838 + t814 + (-t761 + t775) * qJ(5);
t703 = t760 * pkin(4) + t805 + t823;
t767 = mrSges(6,1) * t827 - (qJD(4) * mrSges(6,3));
t768 = mrSges(6,1) * t826 + (qJD(4) * mrSges(6,2));
t714 = -t793 * t764 + (t763 * t790 + t784 * t788) * t796;
t803 = -t798 * qJ(3) + qJDD(3) - t714;
t711 = t837 * qJDD(2) + t803;
t731 = -t788 * t763 + t790 * t784;
t705 = t795 * t711 - t792 * t731;
t757 = (pkin(4) * t792 - qJ(5) * t795) * qJD(2);
t797 = qJD(4) ^ 2;
t701 = -qJDD(4) * pkin(4) - t797 * qJ(5) + t757 * t826 + qJDD(5) - t705;
t697 = (t792 * t795 * t798 - qJDD(4)) * pkin(9) + (t761 + t775) * pkin(5) + t701;
t770 = pkin(5) * t826 - qJD(4) * pkin(9);
t783 = t792 ^ 2;
t698 = -t770 * t826 + (pkin(4) + pkin(9)) * t760 + (-pkin(5) * t783 + t837) * t798 + t805;
t791 = sin(qJ(6));
t794 = cos(qJ(6));
t693 = t794 * t697 - t791 * t698;
t755 = -t791 * qJD(4) + t794 * t827;
t724 = t755 * qJD(6) + t794 * qJDD(4) + t791 * t760;
t756 = t794 * qJD(4) + t791 * t827;
t725 = -t755 * mrSges(7,1) + t756 * mrSges(7,2);
t773 = qJD(6) + t826;
t728 = -t773 * mrSges(7,2) + t755 * mrSges(7,3);
t752 = qJDD(6) + t761;
t690 = m(7) * t693 + t752 * mrSges(7,1) - t724 * mrSges(7,3) - t756 * t725 + t773 * t728;
t694 = t791 * t697 + t794 * t698;
t723 = -t756 * qJD(6) - t791 * qJDD(4) + t794 * t760;
t729 = t773 * mrSges(7,1) - t756 * mrSges(7,3);
t691 = m(7) * t694 - t752 * mrSges(7,2) + t723 * mrSges(7,3) + t755 * t725 - t773 * t729;
t819 = -t791 * t690 + t794 * t691;
t801 = m(6) * t703 - t761 * mrSges(6,3) - (t767 * t792 + t768 * t795) * qJD(2) + t819;
t835 = mrSges(5,1) - mrSges(6,2);
t840 = -m(5) * t709 - t761 * mrSges(5,2) - t835 * t760 - t765 * t827 - t766 * t826 - t801;
t836 = mrSges(3,1) - mrSges(4,2);
t834 = Ifges(5,4) + Ifges(6,6);
t833 = (Ifges(3,5) - Ifges(4,4));
t832 = Ifges(3,6) - Ifges(4,5);
t682 = t794 * t690 + t791 * t691;
t807 = -m(6) * t701 - t761 * mrSges(6,1) - t682;
t758 = (-mrSges(6,2) * t792 - mrSges(6,3) * t795) * qJD(2);
t817 = qJD(2) * (-t758 - (mrSges(5,1) * t792 + mrSges(5,2) * t795) * qJD(2));
t678 = m(5) * t705 - t761 * mrSges(5,3) + t835 * qJDD(4) + (t765 - t767) * qJD(4) + t795 * t817 + t807;
t706 = t792 * t711 + t795 * t731;
t804 = -t797 * pkin(4) + qJDD(4) * qJ(5) - t757 * t827 + t706;
t699 = (qJD(4) * t838) - t804;
t696 = -t783 * t798 * pkin(9) - t760 * pkin(5) + ((2 * qJD(5)) + t770) * qJD(4) + t804;
t809 = -m(7) * t696 + t723 * mrSges(7,1) - t724 * mrSges(7,2) + t755 * t728 - t756 * t729;
t802 = -m(6) * t699 + qJDD(4) * mrSges(6,3) + qJD(4) * t768 - t809;
t687 = m(5) * t706 - qJDD(4) * mrSges(5,2) - qJD(4) * t766 + (-mrSges(5,3) - mrSges(6,1)) * t760 + t792 * t817 + t802;
t673 = t795 * t678 + t792 * t687;
t713 = -qJDD(2) * pkin(2) + t803;
t808 = -m(4) * t713 + (t798 * mrSges(4,3)) - t673;
t668 = m(3) * t714 - (t798 * mrSges(3,2)) + t836 * qJDD(2) + t808;
t831 = t668 * t796;
t820 = -t792 * t678 + t795 * t687;
t672 = m(4) * t731 + t820;
t671 = m(3) * t731 + t672;
t712 = t798 * pkin(2) - t814;
t800 = -m(4) * t712 + (t798 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t840;
t677 = m(3) * t715 - (t798 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t800;
t659 = -t788 * t671 + t677 * t829 + t790 * t831;
t657 = m(2) * t763 + t659;
t665 = -t793 * t668 + t796 * t677;
t664 = m(2) * t764 + t665;
t828 = t789 * t657 + t787 * t664;
t658 = t790 * t671 + t677 * t830 + t788 * t831;
t821 = -t787 * t657 + t789 * t664;
t818 = qJD(2) * (-(t841 * qJD(4)) + (t842 * t792 - t795 * t843) * qJD(2));
t815 = m(2) * t784 + t658;
t679 = -t760 * mrSges(6,2) + t801;
t716 = Ifges(7,5) * t756 + Ifges(7,6) * t755 + Ifges(7,3) * t773;
t718 = Ifges(7,1) * t756 + Ifges(7,4) * t755 + Ifges(7,5) * t773;
t684 = -mrSges(7,1) * t696 + mrSges(7,3) * t694 + Ifges(7,4) * t724 + Ifges(7,2) * t723 + Ifges(7,6) * t752 - t756 * t716 + t773 * t718;
t717 = Ifges(7,4) * t756 + Ifges(7,2) * t755 + Ifges(7,6) * t773;
t685 = mrSges(7,2) * t696 - mrSges(7,3) * t693 + Ifges(7,1) * t724 + Ifges(7,4) * t723 + Ifges(7,5) * t752 + t755 * t716 - t773 * t717;
t738 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t795 - Ifges(5,4) * t792) * qJD(2);
t740 = (Ifges(6,4) * qJD(4)) + (-Ifges(6,2) * t795 + Ifges(6,6) * t792) * qJD(2);
t660 = -mrSges(5,1) * t709 + mrSges(5,3) * t706 - mrSges(6,1) * t699 + mrSges(6,2) * t703 - t791 * t685 - t794 * t684 - pkin(5) * t809 - pkin(9) * t819 - pkin(4) * t679 + t834 * t761 + (-Ifges(5,2) - Ifges(6,3)) * t760 + t842 * qJDD(4) + (t738 - t740) * qJD(4) + t795 * t818;
t737 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t795 - Ifges(5,2) * t792) * qJD(2);
t739 = (Ifges(6,5) * qJD(4)) + (-Ifges(6,6) * t795 + Ifges(6,3) * t792) * qJD(2);
t806 = mrSges(7,1) * t693 - mrSges(7,2) * t694 + Ifges(7,5) * t724 + Ifges(7,6) * t723 + Ifges(7,3) * t752 + t756 * t717 - t755 * t718;
t661 = t792 * t818 + (Ifges(5,1) + Ifges(6,2)) * t761 - t834 * t760 + (-t737 + t739) * qJD(4) + t806 + mrSges(5,2) * t709 + mrSges(6,1) * t701 - mrSges(6,3) * t703 - mrSges(5,3) * t705 + pkin(5) * t682 - qJ(5) * t679 + t843 * qJDD(4);
t654 = -mrSges(4,1) * t712 + mrSges(3,3) * t715 - pkin(2) * t672 - pkin(3) * t840 - pkin(8) * t820 + t832 * qJDD(2) - t795 * t660 - t792 * t661 - t836 * t731 + (t833 * t798);
t681 = qJDD(4) * mrSges(6,2) + qJD(4) * t767 + t758 * t826 - t807;
t799 = -mrSges(5,2) * t706 - mrSges(6,3) * t699 - pkin(9) * t682 - t791 * t684 - pkin(4) * t681 + t794 * t685 + qJ(5) * (-t758 * t827 + t802) + mrSges(6,2) * t701 + mrSges(5,1) * t705 + t738 * t827 + t737 * t826 + (-t795 * t739 - t792 * t740) * qJD(2) + t843 * t761 + (-qJ(5) * mrSges(6,1) - t842) * t760 + t841 * qJDD(4);
t655 = -t832 * t798 + (mrSges(3,2) - mrSges(4,3)) * t731 + mrSges(4,1) * t713 - mrSges(3,3) * t714 + pkin(3) * t673 - qJ(3) * t672 + t799 + t833 * qJDD(2);
t810 = pkin(7) * t665 + t654 * t796 + t655 * t793;
t669 = qJDD(2) * mrSges(4,2) - t808;
t653 = mrSges(3,1) * t714 - mrSges(3,2) * t715 + mrSges(4,2) * t713 - mrSges(4,3) * t712 + t795 * t661 - t792 * t660 - pkin(8) * t673 - pkin(2) * t669 + qJ(3) * t800 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t652 = mrSges(2,2) * t784 - mrSges(2,3) * t763 - t793 * t654 + t796 * t655 + (-t658 * t788 - t659 * t790) * pkin(7);
t651 = -mrSges(2,1) * t784 + mrSges(2,3) * t764 - pkin(1) * t658 - t788 * t653 + t810 * t790;
t1 = [-m(1) * g(1) + t821; -m(1) * g(2) + t828; -m(1) * g(3) + t815; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t828 - t787 * t651 + t789 * t652; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t821 + t789 * t651 + t787 * t652; -mrSges(1,1) * g(2) + mrSges(2,1) * t763 + mrSges(1,2) * g(1) - mrSges(2,2) * t764 + pkin(1) * t659 + t790 * t653 + t810 * t788; t815; t653; t669; t799; t681; t806;];
tauJB  = t1;
