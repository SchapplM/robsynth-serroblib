% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-05-04 21:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:40:26
% EndTime: 2019-05-04 21:40:36
% DurationCPUTime: 10.38s
% Computational Cost: add. (176048->277), mult. (357223->354), div. (0->0), fcn. (263631->14), ass. (0->135)
t794 = qJD(2) ^ 2;
t779 = sin(pkin(12));
t783 = cos(pkin(12));
t788 = sin(qJ(5));
t791 = cos(qJ(5));
t801 = t779 * t788 - t783 * t791;
t758 = t801 * qJD(2);
t802 = t779 * t791 + t783 * t788;
t759 = t802 * qJD(2);
t814 = t759 * qJD(5);
t746 = -qJDD(2) * t801 - t814;
t825 = pkin(4) * t783;
t824 = mrSges(5,2) * t779;
t777 = t783 ^ 2;
t823 = t777 * t794;
t782 = sin(pkin(6));
t789 = sin(qJ(2));
t822 = t782 * t789;
t792 = cos(qJ(2));
t821 = t782 * t792;
t786 = cos(pkin(6));
t820 = t786 * t789;
t819 = t786 * t792;
t781 = sin(pkin(10));
t785 = cos(pkin(10));
t765 = g(1) * t781 - g(2) * t785;
t766 = -g(1) * t785 - g(2) * t781;
t778 = -g(3) + qJDD(1);
t739 = t765 * t819 - t766 * t789 + t778 * t821;
t734 = qJDD(2) * pkin(2) + t739;
t740 = t765 * t820 + t792 * t766 + t778 * t822;
t735 = -pkin(2) * t794 + t740;
t780 = sin(pkin(11));
t784 = cos(pkin(11));
t720 = t780 * t734 + t784 * t735;
t718 = -pkin(3) * t794 + qJDD(2) * qJ(4) + t720;
t756 = -t765 * t782 + t786 * t778;
t753 = qJDD(3) + t756;
t813 = qJD(2) * qJD(4);
t817 = t783 * t753 - 0.2e1 * t779 * t813;
t711 = (-pkin(8) * qJDD(2) + t794 * t825 - t718) * t779 + t817;
t714 = t779 * t753 + (t718 + 0.2e1 * t813) * t783;
t812 = qJDD(2) * t783;
t712 = -pkin(4) * t823 + pkin(8) * t812 + t714;
t707 = t788 * t711 + t791 * t712;
t745 = pkin(5) * t758 - pkin(9) * t759;
t793 = qJD(5) ^ 2;
t705 = -pkin(5) * t793 + qJDD(5) * pkin(9) - t745 * t758 + t707;
t776 = t779 ^ 2;
t719 = t784 * t734 - t780 * t735;
t803 = qJDD(4) - t719;
t715 = (-pkin(3) - t825) * qJDD(2) + (-qJ(4) + (-t776 - t777) * pkin(8)) * t794 + t803;
t815 = t758 * qJD(5);
t747 = qJDD(2) * t802 - t815;
t708 = (-t747 + t815) * pkin(9) + (-t746 + t814) * pkin(5) + t715;
t787 = sin(qJ(6));
t790 = cos(qJ(6));
t702 = -t705 * t787 + t708 * t790;
t750 = qJD(5) * t790 - t759 * t787;
t727 = qJD(6) * t750 + qJDD(5) * t787 + t747 * t790;
t751 = qJD(5) * t787 + t759 * t790;
t728 = -mrSges(7,1) * t750 + mrSges(7,2) * t751;
t757 = qJD(6) + t758;
t729 = -mrSges(7,2) * t757 + mrSges(7,3) * t750;
t744 = qJDD(6) - t746;
t699 = m(7) * t702 + mrSges(7,1) * t744 - mrSges(7,3) * t727 - t728 * t751 + t729 * t757;
t703 = t705 * t790 + t708 * t787;
t726 = -qJD(6) * t751 + qJDD(5) * t790 - t747 * t787;
t730 = mrSges(7,1) * t757 - mrSges(7,3) * t751;
t700 = m(7) * t703 - mrSges(7,2) * t744 + mrSges(7,3) * t726 + t728 * t750 - t730 * t757;
t691 = -t699 * t787 + t790 * t700;
t742 = mrSges(6,1) * t758 + mrSges(6,2) * t759;
t755 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t759;
t688 = m(6) * t707 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t746 - qJD(5) * t755 - t742 * t758 + t691;
t706 = t711 * t791 - t712 * t788;
t704 = -qJDD(5) * pkin(5) - pkin(9) * t793 + t745 * t759 - t706;
t701 = -m(7) * t704 + t726 * mrSges(7,1) - mrSges(7,2) * t727 + t750 * t729 - t730 * t751;
t754 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t758;
t695 = m(6) * t706 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t747 + qJD(5) * t754 - t742 * t759 + t701;
t683 = t788 * t688 + t791 * t695;
t713 = -t718 * t779 + t817;
t800 = mrSges(5,3) * qJDD(2) + t794 * (-mrSges(5,1) * t783 + t824);
t681 = m(5) * t713 - t779 * t800 + t683;
t808 = t791 * t688 - t788 * t695;
t682 = m(5) * t714 + t783 * t800 + t808;
t809 = -t681 * t779 + t783 * t682;
t671 = m(4) * t720 - mrSges(4,1) * t794 - qJDD(2) * mrSges(4,2) + t809;
t717 = -qJDD(2) * pkin(3) - t794 * qJ(4) + t803;
t690 = t790 * t699 + t787 * t700;
t798 = m(6) * t715 - t746 * mrSges(6,1) + t747 * mrSges(6,2) + t758 * t754 + t759 * t755 + t690;
t797 = -m(5) * t717 + mrSges(5,1) * t812 - t798 + (t776 * t794 + t823) * mrSges(5,3);
t685 = t797 + (mrSges(4,1) - t824) * qJDD(2) - t794 * mrSges(4,2) + m(4) * t719;
t668 = t780 * t671 + t784 * t685;
t666 = m(3) * t739 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t794 + t668;
t810 = t784 * t671 - t685 * t780;
t667 = m(3) * t740 - mrSges(3,1) * t794 - qJDD(2) * mrSges(3,2) + t810;
t675 = t783 * t681 + t779 * t682;
t674 = m(4) * t753 + t675;
t673 = m(3) * t756 + t674;
t653 = t666 * t819 + t667 * t820 - t673 * t782;
t651 = m(2) * t765 + t653;
t657 = -t666 * t789 + t792 * t667;
t656 = m(2) * t766 + t657;
t818 = t785 * t651 + t781 * t656;
t804 = Ifges(5,5) * t779 + Ifges(5,6) * t783;
t816 = t794 * t804;
t652 = t666 * t821 + t667 * t822 + t786 * t673;
t811 = -t651 * t781 + t785 * t656;
t807 = m(2) * t778 + t652;
t806 = Ifges(5,1) * t779 + Ifges(5,4) * t783;
t805 = Ifges(5,4) * t779 + Ifges(5,2) * t783;
t721 = Ifges(7,5) * t751 + Ifges(7,6) * t750 + Ifges(7,3) * t757;
t723 = Ifges(7,1) * t751 + Ifges(7,4) * t750 + Ifges(7,5) * t757;
t692 = -mrSges(7,1) * t704 + mrSges(7,3) * t703 + Ifges(7,4) * t727 + Ifges(7,2) * t726 + Ifges(7,6) * t744 - t721 * t751 + t723 * t757;
t722 = Ifges(7,4) * t751 + Ifges(7,2) * t750 + Ifges(7,6) * t757;
t693 = mrSges(7,2) * t704 - mrSges(7,3) * t702 + Ifges(7,1) * t727 + Ifges(7,4) * t726 + Ifges(7,5) * t744 + t721 * t750 - t722 * t757;
t736 = Ifges(6,5) * t759 - Ifges(6,6) * t758 + Ifges(6,3) * qJD(5);
t737 = Ifges(6,4) * t759 - Ifges(6,2) * t758 + Ifges(6,6) * qJD(5);
t676 = mrSges(6,2) * t715 - mrSges(6,3) * t706 + Ifges(6,1) * t747 + Ifges(6,4) * t746 + Ifges(6,5) * qJDD(5) - pkin(9) * t690 - qJD(5) * t737 - t692 * t787 + t693 * t790 - t736 * t758;
t738 = Ifges(6,1) * t759 - Ifges(6,4) * t758 + Ifges(6,5) * qJD(5);
t796 = mrSges(7,1) * t702 - mrSges(7,2) * t703 + Ifges(7,5) * t727 + Ifges(7,6) * t726 + Ifges(7,3) * t744 + t722 * t751 - t723 * t750;
t677 = -mrSges(6,1) * t715 + mrSges(6,3) * t707 + Ifges(6,4) * t747 + Ifges(6,2) * t746 + Ifges(6,6) * qJDD(5) - pkin(5) * t690 + qJD(5) * t738 - t736 * t759 - t796;
t659 = -mrSges(5,1) * t717 + mrSges(5,3) * t714 - pkin(4) * t798 + pkin(8) * t808 + qJDD(2) * t805 + t788 * t676 + t791 * t677 - t779 * t816;
t660 = mrSges(5,2) * t717 - mrSges(5,3) * t713 - pkin(8) * t683 + qJDD(2) * t806 + t791 * t676 - t788 * t677 + t783 * t816;
t649 = mrSges(4,2) * t753 - mrSges(4,3) * t719 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t794 - qJ(4) * t675 - t659 * t779 + t660 * t783;
t795 = mrSges(6,1) * t706 - mrSges(6,2) * t707 + Ifges(6,5) * t747 + Ifges(6,6) * t746 + Ifges(6,3) * qJDD(5) + pkin(5) * t701 + pkin(9) * t691 + t790 * t692 + t787 * t693 + t759 * t737 + t758 * t738;
t658 = -t795 + (Ifges(4,6) - t804) * qJDD(2) - mrSges(4,1) * t753 - mrSges(5,1) * t713 + mrSges(5,2) * t714 + mrSges(4,3) * t720 - pkin(4) * t683 - pkin(3) * t675 + (-t779 * t805 + t783 * t806 + Ifges(4,5)) * t794;
t646 = -mrSges(3,1) * t756 + mrSges(3,3) * t740 + t794 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t674 + qJ(3) * t810 + t780 * t649 + t784 * t658;
t647 = mrSges(3,2) * t756 - mrSges(3,3) * t739 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t794 - qJ(3) * t668 + t649 * t784 - t658 * t780;
t799 = pkin(7) * t657 + t646 * t792 + t647 * t789;
t689 = qJDD(2) * t824 - t797;
t648 = mrSges(3,1) * t739 - mrSges(3,2) * t740 + mrSges(4,1) * t719 - mrSges(4,2) * t720 + t779 * t660 + t783 * t659 - pkin(3) * t689 + qJ(4) * t809 + pkin(2) * t668 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t645 = mrSges(2,2) * t778 - mrSges(2,3) * t765 - t789 * t646 + t792 * t647 + (-t652 * t782 - t653 * t786) * pkin(7);
t644 = -mrSges(2,1) * t778 + mrSges(2,3) * t766 - pkin(1) * t652 - t782 * t648 + t786 * t799;
t1 = [-m(1) * g(1) + t811; -m(1) * g(2) + t818; -m(1) * g(3) + t807; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t818 - t781 * t644 + t785 * t645; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t811 + t785 * t644 + t781 * t645; -mrSges(1,1) * g(2) + mrSges(2,1) * t765 + mrSges(1,2) * g(1) - mrSges(2,2) * t766 + pkin(1) * t653 + t786 * t648 + t782 * t799; t807; t648; t674; t689; t795; t796;];
tauJB  = t1;
