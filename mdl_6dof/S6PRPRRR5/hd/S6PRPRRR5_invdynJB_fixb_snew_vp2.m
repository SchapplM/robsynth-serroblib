% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:16:10
% EndTime: 2019-05-05 01:16:17
% DurationCPUTime: 6.47s
% Computational Cost: add. (109072->294), mult. (206353->371), div. (0->0), fcn. (140398->12), ass. (0->131)
t780 = sin(pkin(11));
t782 = cos(pkin(11));
t761 = t780 * g(1) - t782 * g(2);
t762 = -t782 * g(1) - t780 * g(2);
t777 = -g(3) + qJDD(1);
t791 = cos(qJ(2));
t783 = cos(pkin(6));
t787 = sin(qJ(2));
t817 = t783 * t787;
t781 = sin(pkin(6));
t818 = t781 * t787;
t730 = t761 * t817 + t791 * t762 + t777 * t818;
t805 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t730;
t785 = sin(qJ(5));
t786 = sin(qJ(4));
t789 = cos(qJ(5));
t790 = cos(qJ(4));
t748 = (t785 * t790 + t786 * t789) * qJD(2);
t729 = -t787 * t762 + (t761 * t783 + t777 * t781) * t791;
t824 = -pkin(2) - pkin(8);
t792 = qJD(2) ^ 2;
t823 = pkin(4) * t792;
t822 = mrSges(3,1) - mrSges(4,2);
t821 = (Ifges(3,5) - Ifges(4,4));
t820 = -Ifges(3,6) + Ifges(4,5);
t798 = -t792 * qJ(3) + qJDD(3) - t729;
t715 = qJDD(2) * t824 + t798;
t712 = t790 * t715;
t740 = -t781 * t761 + t783 * t777;
t813 = qJD(2) * qJD(4);
t760 = t790 * qJDD(2) - t786 * t813;
t695 = (qJDD(4) * pkin(4)) - t760 * pkin(9) + t712 + (-pkin(9) * t813 - t790 * t823 - t740) * t786;
t704 = t786 * t715 + t790 * t740;
t759 = -t786 * qJDD(2) - t790 * t813;
t814 = qJD(2) * t790;
t766 = (qJD(4) * pkin(4)) - pkin(9) * t814;
t776 = t786 ^ 2;
t696 = t759 * pkin(9) - qJD(4) * t766 - t776 * t823 + t704;
t691 = t785 * t695 + t789 * t696;
t749 = (-t785 * t786 + t789 * t790) * qJD(2);
t721 = -t749 * qJD(5) + t789 * t759 - t785 * t760;
t732 = t748 * mrSges(6,1) + t749 * mrSges(6,2);
t772 = qJD(4) + qJD(5);
t738 = t772 * mrSges(6,1) - t749 * mrSges(6,3);
t771 = qJDD(4) + qJDD(5);
t733 = t748 * pkin(5) - t749 * pkin(10);
t770 = t772 ^ 2;
t688 = -t770 * pkin(5) + t771 * pkin(10) - t748 * t733 + t691;
t701 = -t759 * pkin(4) + t766 * t814 + (-pkin(9) * t776 + t824) * t792 + t805;
t722 = -t748 * qJD(5) + t785 * t759 + t789 * t760;
t692 = (t748 * t772 - t722) * pkin(10) + (t749 * t772 - t721) * pkin(5) + t701;
t784 = sin(qJ(6));
t788 = cos(qJ(6));
t685 = -t784 * t688 + t788 * t692;
t734 = -t784 * t749 + t788 * t772;
t699 = t734 * qJD(6) + t788 * t722 + t784 * t771;
t735 = t788 * t749 + t784 * t772;
t709 = -t734 * mrSges(7,1) + t735 * mrSges(7,2);
t719 = qJDD(6) - t721;
t743 = qJD(6) + t748;
t724 = -t743 * mrSges(7,2) + t734 * mrSges(7,3);
t682 = m(7) * t685 + t719 * mrSges(7,1) - t699 * mrSges(7,3) - t735 * t709 + t743 * t724;
t686 = t788 * t688 + t784 * t692;
t698 = -t735 * qJD(6) - t784 * t722 + t788 * t771;
t725 = t743 * mrSges(7,1) - t735 * mrSges(7,3);
t683 = m(7) * t686 - t719 * mrSges(7,2) + t698 * mrSges(7,3) + t734 * t709 - t743 * t725;
t808 = -t784 * t682 + t788 * t683;
t670 = m(6) * t691 - t771 * mrSges(6,2) + t721 * mrSges(6,3) - t748 * t732 - t772 * t738 + t808;
t690 = t789 * t695 - t785 * t696;
t737 = -t772 * mrSges(6,2) - t748 * mrSges(6,3);
t687 = -t771 * pkin(5) - t770 * pkin(10) + t749 * t733 - t690;
t800 = -m(7) * t687 + t698 * mrSges(7,1) - t699 * mrSges(7,2) + t734 * t724 - t735 * t725;
t678 = m(6) * t690 + t771 * mrSges(6,1) - t722 * mrSges(6,3) - t749 * t732 + t772 * t737 + t800;
t663 = t785 * t670 + t789 * t678;
t703 = -t786 * t740 + t712;
t758 = (mrSges(5,1) * t786 + mrSges(5,2) * t790) * qJD(2);
t815 = qJD(2) * t786;
t763 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t815;
t660 = m(5) * t703 + qJDD(4) * mrSges(5,1) - t760 * mrSges(5,3) + qJD(4) * t763 - t758 * t814 + t663;
t764 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t814;
t809 = t789 * t670 - t785 * t678;
t661 = m(5) * t704 - qJDD(4) * mrSges(5,2) + t759 * mrSges(5,3) - qJD(4) * t764 - t758 * t815 + t809;
t656 = t790 * t660 + t786 * t661;
t723 = -qJDD(2) * pkin(2) + t798;
t801 = -m(4) * t723 + (t792 * mrSges(4,3)) - t656;
t651 = m(3) * t729 - (t792 * mrSges(3,2)) + qJDD(2) * t822 + t801;
t819 = t651 * t791;
t810 = -t786 * t660 + t790 * t661;
t655 = m(4) * t740 + t810;
t654 = m(3) * t740 + t655;
t720 = t792 * pkin(2) - t805;
t714 = t792 * t824 + t805;
t672 = t788 * t682 + t784 * t683;
t799 = m(6) * t701 - t721 * mrSges(6,1) + t722 * mrSges(6,2) + t748 * t737 + t749 * t738 + t672;
t796 = -m(5) * t714 + t759 * mrSges(5,1) - t760 * mrSges(5,2) - t763 * t815 - t764 * t814 - t799;
t794 = -m(4) * t720 + (t792 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t796;
t667 = m(3) * t730 - (t792 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t794;
t642 = -t781 * t654 + t667 * t817 + t783 * t819;
t640 = m(2) * t761 + t642;
t648 = -t787 * t651 + t791 * t667;
t647 = m(2) * t762 + t648;
t816 = t782 * t640 + t780 * t647;
t641 = t783 * t654 + t667 * t818 + t781 * t819;
t811 = -t780 * t640 + t782 * t647;
t806 = m(2) * t777 + t641;
t705 = Ifges(7,5) * t735 + Ifges(7,6) * t734 + Ifges(7,3) * t743;
t707 = Ifges(7,1) * t735 + Ifges(7,4) * t734 + Ifges(7,5) * t743;
t675 = -mrSges(7,1) * t687 + mrSges(7,3) * t686 + Ifges(7,4) * t699 + Ifges(7,2) * t698 + Ifges(7,6) * t719 - t735 * t705 + t743 * t707;
t706 = Ifges(7,4) * t735 + Ifges(7,2) * t734 + Ifges(7,6) * t743;
t676 = mrSges(7,2) * t687 - mrSges(7,3) * t685 + Ifges(7,1) * t699 + Ifges(7,4) * t698 + Ifges(7,5) * t719 + t734 * t705 - t743 * t706;
t726 = Ifges(6,5) * t749 - Ifges(6,6) * t748 + Ifges(6,3) * t772;
t727 = Ifges(6,4) * t749 - Ifges(6,2) * t748 + Ifges(6,6) * t772;
t657 = mrSges(6,2) * t701 - mrSges(6,3) * t690 + Ifges(6,1) * t722 + Ifges(6,4) * t721 + Ifges(6,5) * t771 - pkin(10) * t672 - t784 * t675 + t788 * t676 - t748 * t726 - t772 * t727;
t728 = Ifges(6,1) * t749 - Ifges(6,4) * t748 + Ifges(6,5) * t772;
t795 = mrSges(7,1) * t685 - mrSges(7,2) * t686 + Ifges(7,5) * t699 + Ifges(7,6) * t698 + Ifges(7,3) * t719 + t735 * t706 - t734 * t707;
t658 = -mrSges(6,1) * t701 + mrSges(6,3) * t691 + Ifges(6,4) * t722 + Ifges(6,2) * t721 + Ifges(6,6) * t771 - pkin(5) * t672 - t749 * t726 + t772 * t728 - t795;
t745 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t790 - Ifges(5,6) * t786) * qJD(2);
t747 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t790 - Ifges(5,4) * t786) * qJD(2);
t643 = -mrSges(5,1) * t714 + mrSges(5,3) * t704 + Ifges(5,4) * t760 + Ifges(5,2) * t759 + Ifges(5,6) * qJDD(4) - pkin(4) * t799 + pkin(9) * t809 + qJD(4) * t747 + t785 * t657 + t789 * t658 - t745 * t814;
t746 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t790 - Ifges(5,2) * t786) * qJD(2);
t644 = mrSges(5,2) * t714 - mrSges(5,3) * t703 + Ifges(5,1) * t760 + Ifges(5,4) * t759 + Ifges(5,5) * qJDD(4) - pkin(9) * t663 - qJD(4) * t746 + t789 * t657 - t785 * t658 - t745 * t815;
t637 = -mrSges(4,1) * t720 + mrSges(3,3) * t730 - pkin(2) * t655 - pkin(3) * t796 - pkin(8) * t810 - qJDD(2) * t820 - t790 * t643 - t786 * t644 - t740 * t822 + (t792 * t821);
t797 = mrSges(6,1) * t690 - mrSges(6,2) * t691 + Ifges(6,5) * t722 + Ifges(6,6) * t721 + Ifges(6,3) * t771 + pkin(5) * t800 + pkin(10) * t808 + t788 * t675 + t784 * t676 + t749 * t727 + t748 * t728;
t793 = mrSges(5,1) * t703 - mrSges(5,2) * t704 + Ifges(5,5) * t760 + Ifges(5,6) * t759 + (Ifges(5,3) * qJDD(4)) + pkin(4) * t663 + t746 * t814 + t747 * t815 + t797;
t638 = t820 * t792 + (mrSges(3,2) - mrSges(4,3)) * t740 - mrSges(3,3) * t729 + mrSges(4,1) * t723 + t821 * qJDD(2) + t793 + pkin(3) * t656 - qJ(3) * t655;
t802 = pkin(7) * t648 + t637 * t791 + t638 * t787;
t652 = qJDD(2) * mrSges(4,2) - t801;
t636 = mrSges(3,1) * t729 - mrSges(3,2) * t730 + mrSges(4,2) * t723 - mrSges(4,3) * t720 + t790 * t644 - t786 * t643 - pkin(8) * t656 - pkin(2) * t652 + qJ(3) * t794 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t635 = mrSges(2,2) * t777 - mrSges(2,3) * t761 - t787 * t637 + t791 * t638 + (-t641 * t781 - t642 * t783) * pkin(7);
t634 = -mrSges(2,1) * t777 + mrSges(2,3) * t762 - pkin(1) * t641 - t781 * t636 + t783 * t802;
t1 = [-m(1) * g(1) + t811; -m(1) * g(2) + t816; -m(1) * g(3) + t806; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t816 - t780 * t634 + t782 * t635; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t811 + t782 * t634 + t780 * t635; -mrSges(1,1) * g(2) + mrSges(2,1) * t761 + mrSges(1,2) * g(1) - mrSges(2,2) * t762 + pkin(1) * t642 + t783 * t636 + t781 * t802; t806; t636; t652; t793; t797; t795;];
tauJB  = t1;
