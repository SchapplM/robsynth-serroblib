% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRP2
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
% Datum: 2019-05-04 23:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:35:30
% EndTime: 2019-05-04 23:35:38
% DurationCPUTime: 7.21s
% Computational Cost: add. (99128->274), mult. (180206->340), div. (0->0), fcn. (122136->12), ass. (0->124)
t808 = Ifges(6,1) + Ifges(7,1);
t801 = Ifges(6,4) - Ifges(7,5);
t800 = -Ifges(6,5) - Ifges(7,4);
t807 = Ifges(6,2) + Ifges(7,3);
t799 = Ifges(6,6) - Ifges(7,6);
t806 = -Ifges(6,3) - Ifges(7,2);
t761 = sin(pkin(10));
t764 = cos(pkin(10));
t747 = g(1) * t761 - g(2) * t764;
t748 = -g(1) * t764 - g(2) * t761;
t759 = -g(3) + qJDD(1);
t768 = sin(qJ(2));
t765 = cos(pkin(6));
t770 = cos(qJ(2));
t794 = t765 * t770;
t762 = sin(pkin(6));
t796 = t762 * t770;
t702 = t747 * t794 - t748 * t768 + t759 * t796;
t700 = qJDD(2) * pkin(2) + t702;
t795 = t765 * t768;
t797 = t762 * t768;
t703 = t747 * t795 + t770 * t748 + t759 * t797;
t772 = qJD(2) ^ 2;
t701 = -pkin(2) * t772 + t703;
t760 = sin(pkin(11));
t763 = cos(pkin(11));
t696 = t760 * t700 + t763 * t701;
t694 = -pkin(3) * t772 + qJDD(2) * pkin(8) + t696;
t728 = -t747 * t762 + t765 * t759;
t727 = qJDD(3) + t728;
t767 = sin(qJ(4));
t769 = cos(qJ(4));
t689 = -t767 * t694 + t769 * t727;
t744 = (-pkin(4) * t769 - pkin(9) * t767) * qJD(2);
t771 = qJD(4) ^ 2;
t788 = qJD(2) * t767;
t685 = -qJDD(4) * pkin(4) - t771 * pkin(9) + t744 * t788 - t689;
t766 = sin(qJ(5));
t803 = cos(qJ(5));
t742 = t766 * qJD(4) + t803 * t788;
t786 = qJD(2) * qJD(4);
t783 = t769 * t786;
t745 = qJDD(2) * t767 + t783;
t714 = t742 * qJD(5) - t803 * qJDD(4) + t766 * t745;
t741 = -t803 * qJD(4) + t766 * t788;
t715 = -t741 * qJD(5) + t766 * qJDD(4) + t803 * t745;
t787 = qJD(2) * t769;
t754 = qJD(5) - t787;
t683 = -0.2e1 * qJD(6) * t742 + (t741 * t754 - t715) * qJ(6) + (t742 * t754 + t714) * pkin(5) + t685;
t725 = -mrSges(7,1) * t754 + mrSges(7,2) * t742;
t726 = -mrSges(7,2) * t741 + mrSges(7,3) * t754;
t677 = m(7) * t683 + mrSges(7,1) * t714 - t715 * mrSges(7,3) - t742 * t725 + t726 * t741;
t690 = t769 * t694 + t767 * t727;
t686 = -pkin(4) * t771 + qJDD(4) * pkin(9) + t744 * t787 + t690;
t695 = t763 * t700 - t760 * t701;
t693 = -qJDD(2) * pkin(3) - t772 * pkin(8) - t695;
t784 = t767 * t786;
t746 = qJDD(2) * t769 - t784;
t688 = (-t745 - t783) * pkin(9) + (-t746 + t784) * pkin(4) + t693;
t682 = t803 * t686 + t766 * t688;
t718 = pkin(5) * t741 - qJ(6) * t742;
t739 = qJDD(5) - t746;
t753 = t754 ^ 2;
t679 = -pkin(5) * t753 + qJ(6) * t739 + 0.2e1 * qJD(6) * t754 - t718 * t741 + t682;
t790 = t801 * t741 - t808 * t742 + t800 * t754;
t792 = t799 * t741 + t800 * t742 + t806 * t754;
t665 = -mrSges(6,1) * t685 - mrSges(7,1) * t683 + mrSges(7,2) * t679 + mrSges(6,3) * t682 - pkin(5) * t677 - t807 * t714 + t801 * t715 + t799 * t739 + t792 * t742 - t790 * t754;
t681 = -t766 * t686 + t803 * t688;
t680 = -t739 * pkin(5) - t753 * qJ(6) + t742 * t718 + qJDD(6) - t681;
t791 = t807 * t741 - t801 * t742 - t799 * t754;
t666 = mrSges(6,2) * t685 + mrSges(7,2) * t680 - mrSges(6,3) * t681 - mrSges(7,3) * t683 - qJ(6) * t677 - t801 * t714 + t808 * t715 - t800 * t739 + t792 * t741 + t791 * t754;
t724 = mrSges(6,1) * t754 - mrSges(6,3) * t742;
t785 = m(7) * t679 + t739 * mrSges(7,3) + t754 * t725;
t719 = mrSges(7,1) * t741 - mrSges(7,3) * t742;
t789 = -mrSges(6,1) * t741 - mrSges(6,2) * t742 - t719;
t802 = -mrSges(6,3) - mrSges(7,2);
t670 = m(6) * t682 - t739 * mrSges(6,2) + t802 * t714 - t754 * t724 + t789 * t741 + t785;
t723 = -mrSges(6,2) * t754 - mrSges(6,3) * t741;
t778 = -m(7) * t680 + t739 * mrSges(7,1) + t754 * t726;
t671 = m(6) * t681 + t739 * mrSges(6,1) + t802 * t715 + t754 * t723 + t789 * t742 + t778;
t668 = t803 * t670 - t671 * t766;
t674 = -m(6) * t685 - t714 * mrSges(6,1) - mrSges(6,2) * t715 - t741 * t723 - t724 * t742 - t677;
t733 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t767 + Ifges(5,2) * t769) * qJD(2);
t734 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t767 + Ifges(5,4) * t769) * qJD(2);
t805 = mrSges(5,1) * t689 - mrSges(5,2) * t690 + Ifges(5,5) * t745 + Ifges(5,6) * t746 + Ifges(5,3) * qJDD(4) + pkin(4) * t674 + pkin(9) * t668 + (t733 * t767 - t734 * t769) * qJD(2) + t803 * t665 + t766 * t666;
t676 = t715 * mrSges(7,2) + t742 * t719 - t778;
t804 = -t799 * t714 - t800 * t715 - t806 * t739 - t790 * t741 - t791 * t742 + mrSges(6,1) * t681 - mrSges(7,1) * t680 - mrSges(6,2) * t682 + mrSges(7,3) * t679 - pkin(5) * t676 + qJ(6) * (-t714 * mrSges(7,2) - t741 * t719 + t785);
t743 = (-mrSges(5,1) * t769 + mrSges(5,2) * t767) * qJD(2);
t749 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t788;
t664 = m(5) * t690 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t746 - qJD(4) * t749 + t743 * t787 + t668;
t750 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t787;
t673 = m(5) * t689 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t745 + qJD(4) * t750 - t743 * t788 + t674;
t780 = t769 * t664 - t673 * t767;
t655 = m(4) * t696 - mrSges(4,1) * t772 - qJDD(2) * mrSges(4,2) + t780;
t667 = t766 * t670 + t803 * t671;
t774 = -m(5) * t693 + t746 * mrSges(5,1) - t745 * mrSges(5,2) - t749 * t788 + t750 * t787 - t667;
t661 = m(4) * t695 + qJDD(2) * mrSges(4,1) - t772 * mrSges(4,2) + t774;
t651 = t760 * t655 + t763 * t661;
t649 = m(3) * t702 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t772 + t651;
t781 = t763 * t655 - t661 * t760;
t650 = m(3) * t703 - mrSges(3,1) * t772 - qJDD(2) * mrSges(3,2) + t781;
t659 = t767 * t664 + t769 * t673;
t658 = m(4) * t727 + t659;
t657 = m(3) * t728 + t658;
t637 = t649 * t794 + t650 * t795 - t657 * t762;
t635 = m(2) * t747 + t637;
t641 = -t649 * t768 + t770 * t650;
t640 = m(2) * t748 + t641;
t793 = t764 * t635 + t761 * t640;
t636 = t649 * t796 + t650 * t797 + t765 * t657;
t782 = -t635 * t761 + t764 * t640;
t779 = m(2) * t759 + t636;
t732 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t767 + Ifges(5,6) * t769) * qJD(2);
t643 = mrSges(5,2) * t693 - mrSges(5,3) * t689 + Ifges(5,1) * t745 + Ifges(5,4) * t746 + Ifges(5,5) * qJDD(4) - pkin(9) * t667 - qJD(4) * t733 - t766 * t665 + t803 * t666 + t732 * t787;
t652 = -mrSges(5,1) * t693 + mrSges(5,3) * t690 + Ifges(5,4) * t745 + Ifges(5,2) * t746 + Ifges(5,6) * qJDD(4) - pkin(4) * t667 + qJD(4) * t734 - t732 * t788 - t804;
t633 = mrSges(4,2) * t727 - mrSges(4,3) * t695 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t772 - pkin(8) * t659 + t643 * t769 - t652 * t767;
t642 = -mrSges(4,1) * t727 + mrSges(4,3) * t696 + t772 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t659 - t805;
t630 = -mrSges(3,1) * t728 + mrSges(3,3) * t703 + t772 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t658 + qJ(3) * t781 + t760 * t633 + t763 * t642;
t631 = mrSges(3,2) * t728 - mrSges(3,3) * t702 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t772 - qJ(3) * t651 + t633 * t763 - t642 * t760;
t776 = pkin(7) * t641 + t630 * t770 + t631 * t768;
t632 = mrSges(3,1) * t702 - mrSges(3,2) * t703 + mrSges(4,1) * t695 - mrSges(4,2) * t696 + t767 * t643 + t769 * t652 + pkin(3) * t774 + pkin(8) * t780 + pkin(2) * t651 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t629 = mrSges(2,2) * t759 - mrSges(2,3) * t747 - t768 * t630 + t770 * t631 + (-t636 * t762 - t637 * t765) * pkin(7);
t628 = -mrSges(2,1) * t759 + mrSges(2,3) * t748 - pkin(1) * t636 - t762 * t632 + t776 * t765;
t1 = [-m(1) * g(1) + t782; -m(1) * g(2) + t793; -m(1) * g(3) + t779; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t793 - t761 * t628 + t764 * t629; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t782 + t764 * t628 + t761 * t629; -mrSges(1,1) * g(2) + mrSges(2,1) * t747 + mrSges(1,2) * g(1) - mrSges(2,2) * t748 + pkin(1) * t637 + t765 * t632 + t776 * t762; t779; t632; t658; t805; t804; t676;];
tauJB  = t1;
