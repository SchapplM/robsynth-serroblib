% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:08:25
% EndTime: 2019-05-04 22:08:36
% DurationCPUTime: 10.97s
% Computational Cost: add. (189957->299), mult. (378757->388), div. (0->0), fcn. (267734->14), ass. (0->133)
t775 = sin(pkin(12));
t779 = cos(pkin(12));
t784 = sin(qJ(4));
t787 = cos(qJ(4));
t750 = (t775 * t784 - t779 * t787) * qJD(2);
t777 = sin(pkin(10));
t781 = cos(pkin(10));
t764 = g(1) * t777 - g(2) * t781;
t765 = -g(1) * t781 - g(2) * t777;
t774 = -g(3) + qJDD(1);
t785 = sin(qJ(2));
t782 = cos(pkin(6));
t788 = cos(qJ(2));
t809 = t782 * t788;
t778 = sin(pkin(6));
t811 = t778 * t788;
t730 = t764 * t809 - t765 * t785 + t774 * t811;
t725 = qJDD(2) * pkin(2) + t730;
t810 = t782 * t785;
t812 = t778 * t785;
t731 = t764 * t810 + t788 * t765 + t774 * t812;
t790 = qJD(2) ^ 2;
t726 = -pkin(2) * t790 + t731;
t776 = sin(pkin(11));
t780 = cos(pkin(11));
t711 = t776 * t725 + t780 * t726;
t709 = -pkin(3) * t790 + qJDD(2) * pkin(8) + t711;
t747 = -t764 * t778 + t782 * t774;
t744 = qJDD(3) + t747;
t705 = -t784 * t709 + t787 * t744;
t805 = qJD(2) * qJD(4);
t804 = t787 * t805;
t762 = qJDD(2) * t784 + t804;
t702 = (-t762 + t804) * qJ(5) + (t784 * t787 * t790 + qJDD(4)) * pkin(4) + t705;
t706 = t787 * t709 + t784 * t744;
t763 = qJDD(2) * t787 - t784 * t805;
t807 = qJD(2) * t784;
t766 = qJD(4) * pkin(4) - qJ(5) * t807;
t773 = t787 ^ 2;
t703 = -pkin(4) * t773 * t790 + qJ(5) * t763 - qJD(4) * t766 + t706;
t814 = 2 * qJD(5);
t698 = t775 * t702 + t779 * t703 - t750 * t814;
t751 = (t775 * t787 + t779 * t784) * qJD(2);
t734 = pkin(5) * t750 - pkin(9) * t751;
t789 = qJD(4) ^ 2;
t696 = -pkin(5) * t789 + qJDD(4) * pkin(9) - t734 * t750 + t698;
t710 = t780 * t725 - t776 * t726;
t795 = -qJDD(2) * pkin(3) - t710;
t704 = -t763 * pkin(4) + qJDD(5) + t766 * t807 + (-qJ(5) * t773 - pkin(8)) * t790 + t795;
t737 = -t762 * t775 + t763 * t779;
t738 = t762 * t779 + t763 * t775;
t699 = (qJD(4) * t750 - t738) * pkin(9) + (qJD(4) * t751 - t737) * pkin(5) + t704;
t783 = sin(qJ(6));
t786 = cos(qJ(6));
t693 = -t696 * t783 + t699 * t786;
t741 = qJD(4) * t786 - t751 * t783;
t718 = qJD(6) * t741 + qJDD(4) * t783 + t738 * t786;
t742 = qJD(4) * t783 + t751 * t786;
t719 = -mrSges(7,1) * t741 + mrSges(7,2) * t742;
t749 = qJD(6) + t750;
t720 = -mrSges(7,2) * t749 + mrSges(7,3) * t741;
t736 = qJDD(6) - t737;
t690 = m(7) * t693 + mrSges(7,1) * t736 - mrSges(7,3) * t718 - t719 * t742 + t720 * t749;
t694 = t696 * t786 + t699 * t783;
t717 = -qJD(6) * t742 + qJDD(4) * t786 - t738 * t783;
t721 = mrSges(7,1) * t749 - mrSges(7,3) * t742;
t691 = m(7) * t694 - mrSges(7,2) * t736 + mrSges(7,3) * t717 + t719 * t741 - t721 * t749;
t682 = -t690 * t783 + t786 * t691;
t733 = mrSges(6,1) * t750 + mrSges(6,2) * t751;
t746 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t751;
t679 = m(6) * t698 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t737 - qJD(4) * t746 - t733 * t750 + t682;
t798 = -t779 * t702 + t775 * t703;
t695 = -qJDD(4) * pkin(5) - t789 * pkin(9) + (t814 + t734) * t751 + t798;
t692 = -m(7) * t695 + t717 * mrSges(7,1) - mrSges(7,2) * t718 + t741 * t720 - t721 * t742;
t697 = -0.2e1 * qJD(5) * t751 - t798;
t745 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t750;
t686 = m(6) * t697 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t738 + qJD(4) * t745 - t733 * t751 + t692;
t674 = t775 * t679 + t779 * t686;
t712 = Ifges(7,5) * t742 + Ifges(7,6) * t741 + Ifges(7,3) * t749;
t714 = Ifges(7,1) * t742 + Ifges(7,4) * t741 + Ifges(7,5) * t749;
t683 = -mrSges(7,1) * t695 + mrSges(7,3) * t694 + Ifges(7,4) * t718 + Ifges(7,2) * t717 + Ifges(7,6) * t736 - t712 * t742 + t714 * t749;
t713 = Ifges(7,4) * t742 + Ifges(7,2) * t741 + Ifges(7,6) * t749;
t684 = mrSges(7,2) * t695 - mrSges(7,3) * t693 + Ifges(7,1) * t718 + Ifges(7,4) * t717 + Ifges(7,5) * t736 + t712 * t741 - t713 * t749;
t728 = Ifges(6,4) * t751 - Ifges(6,2) * t750 + Ifges(6,6) * qJD(4);
t729 = Ifges(6,1) * t751 - Ifges(6,4) * t750 + Ifges(6,5) * qJD(4);
t755 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t784 + Ifges(5,2) * t787) * qJD(2);
t756 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t784 + Ifges(5,4) * t787) * qJD(2);
t815 = (t755 * t784 - t756 * t787) * qJD(2) + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + mrSges(5,1) * t705 + mrSges(6,1) * t697 - mrSges(5,2) * t706 - mrSges(6,2) * t698 + Ifges(5,5) * t762 + Ifges(6,5) * t738 + Ifges(5,6) * t763 + Ifges(6,6) * t737 + pkin(4) * t674 + pkin(5) * t692 + pkin(9) * t682 + t786 * t683 + t783 * t684 + t751 * t728 + t750 * t729;
t761 = (-mrSges(5,1) * t787 + mrSges(5,2) * t784) * qJD(2);
t806 = qJD(2) * t787;
t768 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t806;
t672 = m(5) * t705 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t762 + qJD(4) * t768 - t761 * t807 + t674;
t767 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t807;
t800 = t779 * t679 - t686 * t775;
t673 = m(5) * t706 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t763 - qJD(4) * t767 + t761 * t806 + t800;
t801 = -t672 * t784 + t787 * t673;
t662 = m(4) * t711 - mrSges(4,1) * t790 - qJDD(2) * mrSges(4,2) + t801;
t681 = t786 * t690 + t783 * t691;
t680 = m(6) * t704 - t737 * mrSges(6,1) + mrSges(6,2) * t738 + t750 * t745 + t746 * t751 + t681;
t708 = -t790 * pkin(8) + t795;
t792 = -m(5) * t708 + t763 * mrSges(5,1) - mrSges(5,2) * t762 - t767 * t807 + t768 * t806 - t680;
t676 = m(4) * t710 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t790 + t792;
t659 = t776 * t662 + t780 * t676;
t657 = m(3) * t730 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t790 + t659;
t802 = t780 * t662 - t676 * t776;
t658 = m(3) * t731 - mrSges(3,1) * t790 - qJDD(2) * mrSges(3,2) + t802;
t666 = t787 * t672 + t784 * t673;
t665 = m(4) * t744 + t666;
t664 = m(3) * t747 + t665;
t644 = t657 * t809 + t658 * t810 - t664 * t778;
t642 = m(2) * t764 + t644;
t648 = -t657 * t785 + t788 * t658;
t647 = m(2) * t765 + t648;
t808 = t781 * t642 + t777 * t647;
t643 = t657 * t811 + t658 * t812 + t782 * t664;
t803 = -t642 * t777 + t781 * t647;
t799 = m(2) * t774 + t643;
t727 = Ifges(6,5) * t751 - Ifges(6,6) * t750 + Ifges(6,3) * qJD(4);
t667 = mrSges(6,2) * t704 - mrSges(6,3) * t697 + Ifges(6,1) * t738 + Ifges(6,4) * t737 + Ifges(6,5) * qJDD(4) - pkin(9) * t681 - qJD(4) * t728 - t683 * t783 + t684 * t786 - t727 * t750;
t793 = mrSges(7,1) * t693 - mrSges(7,2) * t694 + Ifges(7,5) * t718 + Ifges(7,6) * t717 + Ifges(7,3) * t736 + t713 * t742 - t714 * t741;
t668 = -mrSges(6,1) * t704 + mrSges(6,3) * t698 + Ifges(6,4) * t738 + Ifges(6,2) * t737 + Ifges(6,6) * qJDD(4) - pkin(5) * t681 + qJD(4) * t729 - t727 * t751 - t793;
t754 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t784 + Ifges(5,6) * t787) * qJD(2);
t650 = -mrSges(5,1) * t708 + mrSges(5,3) * t706 + Ifges(5,4) * t762 + Ifges(5,2) * t763 + Ifges(5,6) * qJDD(4) - pkin(4) * t680 + qJ(5) * t800 + qJD(4) * t756 + t775 * t667 + t779 * t668 - t754 * t807;
t651 = mrSges(5,2) * t708 - mrSges(5,3) * t705 + Ifges(5,1) * t762 + Ifges(5,4) * t763 + Ifges(5,5) * qJDD(4) - qJ(5) * t674 - qJD(4) * t755 + t667 * t779 - t668 * t775 + t754 * t806;
t640 = mrSges(4,2) * t744 - mrSges(4,3) * t710 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t790 - pkin(8) * t666 - t650 * t784 + t651 * t787;
t649 = -mrSges(4,1) * t744 + mrSges(4,3) * t711 + t790 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t666 - t815;
t637 = -mrSges(3,1) * t747 + mrSges(3,3) * t731 + t790 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t665 + qJ(3) * t802 + t776 * t640 + t780 * t649;
t638 = mrSges(3,2) * t747 - mrSges(3,3) * t730 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t790 - qJ(3) * t659 + t640 * t780 - t649 * t776;
t794 = pkin(7) * t648 + t637 * t788 + t638 * t785;
t639 = mrSges(3,1) * t730 - mrSges(3,2) * t731 + mrSges(4,1) * t710 - mrSges(4,2) * t711 + t784 * t651 + t787 * t650 + pkin(3) * t792 + pkin(8) * t801 + pkin(2) * t659 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t636 = mrSges(2,2) * t774 - mrSges(2,3) * t764 - t785 * t637 + t788 * t638 + (-t643 * t778 - t644 * t782) * pkin(7);
t635 = -mrSges(2,1) * t774 + mrSges(2,3) * t765 - pkin(1) * t643 - t778 * t639 + t782 * t794;
t1 = [-m(1) * g(1) + t803; -m(1) * g(2) + t808; -m(1) * g(3) + t799; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t808 - t777 * t635 + t781 * t636; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t803 + t781 * t635 + t777 * t636; -mrSges(1,1) * g(2) + mrSges(2,1) * t764 + mrSges(1,2) * g(1) - mrSges(2,2) * t765 + pkin(1) * t644 + t782 * t639 + t778 * t794; t799; t639; t665; t815; t680; t793;];
tauJB  = t1;
