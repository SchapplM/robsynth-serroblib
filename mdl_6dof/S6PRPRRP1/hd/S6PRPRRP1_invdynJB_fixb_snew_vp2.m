% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:29:38
% EndTime: 2019-05-04 23:29:47
% DurationCPUTime: 7.71s
% Computational Cost: add. (101178->276), mult. (185230->340), div. (0->0), fcn. (125869->12), ass. (0->125)
t810 = Ifges(6,1) + Ifges(7,1);
t804 = Ifges(6,4) + Ifges(7,4);
t803 = Ifges(6,5) + Ifges(7,5);
t809 = Ifges(6,2) + Ifges(7,2);
t802 = Ifges(6,6) + Ifges(7,6);
t808 = Ifges(6,3) + Ifges(7,3);
t762 = sin(pkin(10));
t765 = cos(pkin(10));
t751 = g(1) * t762 - g(2) * t765;
t752 = -g(1) * t765 - g(2) * t762;
t760 = -g(3) + qJDD(1);
t769 = sin(qJ(2));
t766 = cos(pkin(6));
t772 = cos(qJ(2));
t797 = t766 * t772;
t763 = sin(pkin(6));
t799 = t763 * t772;
t706 = t751 * t797 - t752 * t769 + t760 * t799;
t704 = qJDD(2) * pkin(2) + t706;
t798 = t766 * t769;
t800 = t763 * t769;
t707 = t751 * t798 + t772 * t752 + t760 * t800;
t774 = qJD(2) ^ 2;
t705 = -pkin(2) * t774 + t707;
t761 = sin(pkin(11));
t764 = cos(pkin(11));
t699 = t761 * t704 + t764 * t705;
t697 = -pkin(3) * t774 + qJDD(2) * pkin(8) + t699;
t732 = -t751 * t763 + t766 * t760;
t731 = qJDD(3) + t732;
t768 = sin(qJ(4));
t771 = cos(qJ(4));
t693 = t771 * t697 + t768 * t731;
t748 = (-pkin(4) * t771 - pkin(9) * t768) * qJD(2);
t773 = qJD(4) ^ 2;
t790 = qJD(2) * t771;
t688 = -pkin(4) * t773 + qJDD(4) * pkin(9) + t748 * t790 + t693;
t698 = t764 * t704 - t761 * t705;
t696 = -qJDD(2) * pkin(3) - t774 * pkin(8) - t698;
t789 = qJD(2) * qJD(4);
t785 = t771 * t789;
t749 = qJDD(2) * t768 + t785;
t786 = t768 * t789;
t750 = qJDD(2) * t771 - t786;
t691 = (-t749 - t785) * pkin(9) + (-t750 + t786) * pkin(4) + t696;
t767 = sin(qJ(5));
t770 = cos(qJ(5));
t684 = t770 * t688 + t767 * t691;
t791 = qJD(2) * t768;
t746 = qJD(4) * t767 + t770 * t791;
t719 = -qJD(5) * t746 + qJDD(4) * t770 - t749 * t767;
t757 = qJD(5) - t790;
t728 = pkin(5) * t757 - qJ(6) * t746;
t745 = qJD(4) * t770 - t767 * t791;
t742 = t745 ^ 2;
t681 = -pkin(5) * t742 + qJ(6) * t719 + 0.2e1 * qJD(6) * t745 - t728 * t757 + t684;
t720 = qJD(5) * t745 + qJDD(4) * t767 + t749 * t770;
t729 = mrSges(7,1) * t757 - mrSges(7,3) * t746;
t692 = -t768 * t697 + t731 * t771;
t687 = -qJDD(4) * pkin(4) - pkin(9) * t773 + t748 * t791 - t692;
t685 = -pkin(5) * t719 - qJ(6) * t742 + t728 * t746 + qJDD(6) + t687;
t726 = -mrSges(7,2) * t757 + mrSges(7,3) * t745;
t780 = -m(7) * t685 + t719 * mrSges(7,1) + t745 * t726;
t682 = t720 * mrSges(7,2) + t746 * t729 - t780;
t743 = qJDD(5) - t750;
t722 = -mrSges(7,1) * t745 + mrSges(7,2) * t746;
t787 = m(7) * t681 + t719 * mrSges(7,3) + t745 * t722;
t793 = t804 * t745 + t810 * t746 + t803 * t757;
t795 = -t802 * t745 - t803 * t746 - t808 * t757;
t662 = -mrSges(6,1) * t687 + mrSges(6,3) * t684 - mrSges(7,1) * t685 + mrSges(7,3) * t681 - pkin(5) * t682 + qJ(6) * t787 + (-qJ(6) * t729 + t793) * t757 + t795 * t746 + (-mrSges(7,2) * qJ(6) + t802) * t743 + t804 * t720 + t809 * t719;
t683 = -t767 * t688 + t770 * t691;
t679 = -0.2e1 * qJD(6) * t746 + (t745 * t757 - t720) * qJ(6) + (t745 * t746 + t743) * pkin(5) + t683;
t788 = m(7) * t679 + t743 * mrSges(7,1) + t757 * t726;
t677 = -t720 * mrSges(7,3) - t746 * t722 + t788;
t794 = -t809 * t745 - t804 * t746 - t802 * t757;
t668 = mrSges(6,2) * t687 + mrSges(7,2) * t685 - mrSges(6,3) * t683 - mrSges(7,3) * t679 - qJ(6) * t677 + t804 * t719 + t810 * t720 + t803 * t743 - t795 * t745 + t794 * t757;
t723 = -mrSges(6,1) * t745 + mrSges(6,2) * t746;
t727 = -mrSges(6,2) * t757 + mrSges(6,3) * t745;
t671 = m(6) * t683 + t743 * mrSges(6,1) + t757 * t727 + (-t722 - t723) * t746 + (-mrSges(6,3) - mrSges(7,3)) * t720 + t788;
t792 = -mrSges(6,1) * t757 + mrSges(6,3) * t746 - t729;
t805 = -mrSges(6,2) - mrSges(7,2);
t673 = m(6) * t684 + t719 * mrSges(6,3) + t745 * t723 + t743 * t805 + t792 * t757 + t787;
t670 = -t671 * t767 + t770 * t673;
t676 = -m(6) * t687 + t719 * mrSges(6,1) + t720 * t805 + t745 * t727 + t792 * t746 + t780;
t737 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t768 + Ifges(5,2) * t771) * qJD(2);
t738 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t768 + Ifges(5,4) * t771) * qJD(2);
t807 = mrSges(5,1) * t692 - mrSges(5,2) * t693 + Ifges(5,5) * t749 + Ifges(5,6) * t750 + Ifges(5,3) * qJDD(4) + pkin(4) * t676 + pkin(9) * t670 + t770 * t662 + t767 * t668 + (t737 * t768 - t738 * t771) * qJD(2);
t806 = mrSges(6,1) * t683 + mrSges(7,1) * t679 - mrSges(6,2) * t684 - mrSges(7,2) * t681 + pkin(5) * t677 + t719 * t802 + t720 * t803 + t808 * t743 - t793 * t745 - t794 * t746;
t747 = (-mrSges(5,1) * t771 + mrSges(5,2) * t768) * qJD(2);
t753 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t791;
t667 = m(5) * t693 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t750 - qJD(4) * t753 + t747 * t790 + t670;
t754 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t790;
t675 = m(5) * t692 + qJDD(4) * mrSges(5,1) - t749 * mrSges(5,3) + qJD(4) * t754 - t747 * t791 + t676;
t782 = t771 * t667 - t675 * t768;
t657 = m(4) * t699 - mrSges(4,1) * t774 - qJDD(2) * mrSges(4,2) + t782;
t669 = t671 * t770 + t673 * t767;
t776 = -m(5) * t696 + t750 * mrSges(5,1) - mrSges(5,2) * t749 - t753 * t791 + t754 * t790 - t669;
t664 = m(4) * t698 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t774 + t776;
t653 = t761 * t657 + t764 * t664;
t651 = m(3) * t706 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t774 + t653;
t783 = t764 * t657 - t664 * t761;
t652 = m(3) * t707 - mrSges(3,1) * t774 - qJDD(2) * mrSges(3,2) + t783;
t661 = t768 * t667 + t771 * t675;
t660 = m(4) * t731 + t661;
t659 = m(3) * t732 + t660;
t639 = t651 * t797 + t652 * t798 - t659 * t763;
t637 = m(2) * t751 + t639;
t643 = -t651 * t769 + t772 * t652;
t642 = m(2) * t752 + t643;
t796 = t765 * t637 + t762 * t642;
t638 = t651 * t799 + t652 * t800 + t766 * t659;
t784 = -t637 * t762 + t765 * t642;
t781 = m(2) * t760 + t638;
t736 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t768 + Ifges(5,6) * t771) * qJD(2);
t645 = mrSges(5,2) * t696 - mrSges(5,3) * t692 + Ifges(5,1) * t749 + Ifges(5,4) * t750 + Ifges(5,5) * qJDD(4) - pkin(9) * t669 - qJD(4) * t737 - t662 * t767 + t668 * t770 + t736 * t790;
t654 = -mrSges(5,1) * t696 + mrSges(5,3) * t693 + Ifges(5,4) * t749 + Ifges(5,2) * t750 + Ifges(5,6) * qJDD(4) - pkin(4) * t669 + qJD(4) * t738 - t736 * t791 - t806;
t635 = mrSges(4,2) * t731 - mrSges(4,3) * t698 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t774 - pkin(8) * t661 + t645 * t771 - t654 * t768;
t644 = -mrSges(4,1) * t731 + mrSges(4,3) * t699 + t774 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t661 - t807;
t632 = -mrSges(3,1) * t732 + mrSges(3,3) * t707 + t774 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t660 + qJ(3) * t783 + t761 * t635 + t764 * t644;
t633 = mrSges(3,2) * t732 - mrSges(3,3) * t706 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t774 - qJ(3) * t653 + t635 * t764 - t644 * t761;
t778 = pkin(7) * t643 + t632 * t772 + t633 * t769;
t634 = mrSges(3,1) * t706 - mrSges(3,2) * t707 + mrSges(4,1) * t698 - mrSges(4,2) * t699 + t768 * t645 + t771 * t654 + pkin(3) * t776 + pkin(8) * t782 + pkin(2) * t653 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t631 = mrSges(2,2) * t760 - mrSges(2,3) * t751 - t769 * t632 + t772 * t633 + (-t638 * t763 - t639 * t766) * pkin(7);
t630 = -mrSges(2,1) * t760 + mrSges(2,3) * t752 - pkin(1) * t638 - t763 * t634 + t766 * t778;
t1 = [-m(1) * g(1) + t784; -m(1) * g(2) + t796; -m(1) * g(3) + t781; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t796 - t762 * t630 + t765 * t631; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t784 + t765 * t630 + t762 * t631; -mrSges(1,1) * g(2) + mrSges(2,1) * t751 + mrSges(1,2) * g(1) - mrSges(2,2) * t752 + pkin(1) * t639 + t766 * t634 + t763 * t778; t781; t634; t660; t807; t806; t682;];
tauJB  = t1;
