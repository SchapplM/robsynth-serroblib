% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:27:37
% EndTime: 2019-05-04 20:27:50
% DurationCPUTime: 10.53s
% Computational Cost: add. (189526->272), mult. (340497->344), div. (0->0), fcn. (259840->14), ass. (0->129)
t809 = Ifges(6,1) + Ifges(7,1);
t802 = Ifges(6,4) + Ifges(7,4);
t801 = Ifges(6,5) + Ifges(7,5);
t808 = Ifges(6,2) + Ifges(7,2);
t800 = Ifges(6,6) + Ifges(7,6);
t807 = Ifges(6,3) + Ifges(7,3);
t755 = sin(pkin(11));
t759 = cos(pkin(11));
t748 = -g(1) * t759 - g(2) * t755;
t754 = sin(pkin(12));
t758 = cos(pkin(12));
t747 = g(1) * t755 - g(2) * t759;
t753 = -g(3) + qJDD(1);
t757 = sin(pkin(6));
t761 = cos(pkin(6));
t775 = t747 * t761 + t753 * t757;
t707 = -t754 * t748 + t775 * t758;
t708 = t758 * t748 + t775 * t754;
t731 = -t747 * t757 + t753 * t761 + qJDD(2);
t767 = cos(qJ(3));
t760 = cos(pkin(7));
t764 = sin(qJ(3));
t794 = t760 * t764;
t756 = sin(pkin(7));
t795 = t756 * t764;
t700 = t707 * t794 + t767 * t708 + t731 * t795;
t769 = qJD(3) ^ 2;
t698 = -pkin(3) * t769 + qJDD(3) * pkin(9) + t700;
t702 = -t707 * t756 + t731 * t760;
t763 = sin(qJ(4));
t766 = cos(qJ(4));
t691 = t766 * t698 + t763 * t702;
t744 = (-pkin(4) * t766 - pkin(10) * t763) * qJD(3);
t768 = qJD(4) ^ 2;
t787 = qJD(3) * t766;
t689 = -pkin(4) * t768 + qJDD(4) * pkin(10) + t744 * t787 + t691;
t699 = -t764 * t708 + (t707 * t760 + t731 * t756) * t767;
t697 = -qJDD(3) * pkin(3) - t769 * pkin(9) - t699;
t786 = qJD(3) * qJD(4);
t782 = t766 * t786;
t745 = qJDD(3) * t763 + t782;
t783 = t763 * t786;
t746 = qJDD(3) * t766 - t783;
t694 = (-t745 - t782) * pkin(10) + (-t746 + t783) * pkin(4) + t697;
t762 = sin(qJ(5));
t765 = cos(qJ(5));
t684 = -t762 * t689 + t765 * t694;
t788 = qJD(3) * t763;
t741 = qJD(4) * t765 - t762 * t788;
t721 = qJD(5) * t741 + qJDD(4) * t762 + t745 * t765;
t742 = qJD(4) * t762 + t765 * t788;
t723 = -mrSges(7,1) * t741 + mrSges(7,2) * t742;
t724 = -mrSges(6,1) * t741 + mrSges(6,2) * t742;
t751 = qJD(5) - t787;
t727 = -mrSges(6,2) * t751 + mrSges(6,3) * t741;
t740 = qJDD(5) - t746;
t680 = -0.2e1 * qJD(6) * t742 + (t741 * t751 - t721) * qJ(6) + (t741 * t742 + t740) * pkin(5) + t684;
t726 = -mrSges(7,2) * t751 + mrSges(7,3) * t741;
t785 = m(7) * t680 + t740 * mrSges(7,1) + t751 * t726;
t672 = m(6) * t684 + t740 * mrSges(6,1) + t751 * t727 + (-t723 - t724) * t742 + (-mrSges(6,3) - mrSges(7,3)) * t721 + t785;
t685 = t765 * t689 + t762 * t694;
t720 = -qJD(5) * t742 + qJDD(4) * t765 - t745 * t762;
t728 = pkin(5) * t751 - qJ(6) * t742;
t739 = t741 ^ 2;
t682 = -pkin(5) * t739 + qJ(6) * t720 + 0.2e1 * qJD(6) * t741 - t728 * t751 + t685;
t784 = m(7) * t682 + t720 * mrSges(7,3) + t741 * t723;
t729 = mrSges(7,1) * t751 - mrSges(7,3) * t742;
t789 = -mrSges(6,1) * t751 + mrSges(6,3) * t742 - t729;
t803 = -mrSges(6,2) - mrSges(7,2);
t674 = m(6) * t685 + t720 * mrSges(6,3) + t741 * t724 + t803 * t740 + t789 * t751 + t784;
t671 = -t672 * t762 + t765 * t674;
t743 = (-mrSges(5,1) * t766 + mrSges(5,2) * t763) * qJD(3);
t749 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t788;
t668 = m(5) * t691 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t746 - qJD(4) * t749 + t743 * t787 + t671;
t690 = -t763 * t698 + t702 * t766;
t688 = -qJDD(4) * pkin(4) - pkin(10) * t768 + t744 * t788 - t690;
t686 = -pkin(5) * t720 - qJ(6) * t739 + t728 * t742 + qJDD(6) + t688;
t778 = -m(7) * t686 + t720 * mrSges(7,1) + t741 * t726;
t677 = -m(6) * t688 + t720 * mrSges(6,1) + t803 * t721 + t741 * t727 + t789 * t742 + t778;
t750 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t787;
t676 = m(5) * t690 + qJDD(4) * mrSges(5,1) - t745 * mrSges(5,3) + qJD(4) * t750 - t743 * t788 + t677;
t780 = t766 * t668 - t676 * t763;
t658 = m(4) * t700 - mrSges(4,1) * t769 - qJDD(3) * mrSges(4,2) + t780;
t661 = t763 * t668 + t766 * t676;
t660 = m(4) * t702 + t661;
t670 = t672 * t765 + t674 * t762;
t771 = -m(5) * t697 + t746 * mrSges(5,1) - mrSges(5,2) * t745 - t749 * t788 + t750 * t787 - t670;
t665 = m(4) * t699 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t769 + t771;
t796 = t665 * t767;
t647 = t658 * t794 - t660 * t756 + t760 * t796;
t643 = m(3) * t707 + t647;
t653 = t767 * t658 - t665 * t764;
t652 = m(3) * t708 + t653;
t806 = t643 * t758 + t652 * t754;
t683 = t721 * mrSges(7,2) + t742 * t729 - t778;
t790 = t802 * t741 + t809 * t742 + t801 * t751;
t792 = -t800 * t741 - t801 * t742 - t807 * t751;
t662 = -mrSges(6,1) * t688 + mrSges(6,3) * t685 - mrSges(7,1) * t686 + mrSges(7,3) * t682 - pkin(5) * t683 + qJ(6) * t784 + (-qJ(6) * t729 + t790) * t751 + t792 * t742 + (-mrSges(7,2) * qJ(6) + t800) * t740 + t802 * t721 + t808 * t720;
t678 = -t721 * mrSges(7,3) - t742 * t723 + t785;
t791 = -t808 * t741 - t802 * t742 - t800 * t751;
t669 = mrSges(6,2) * t688 + mrSges(7,2) * t686 - mrSges(6,3) * t684 - mrSges(7,3) * t680 - qJ(6) * t678 + t802 * t720 + t809 * t721 + t801 * t740 - t792 * t741 + t791 * t751;
t734 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t763 + Ifges(5,2) * t766) * qJD(3);
t735 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t763 + Ifges(5,4) * t766) * qJD(3);
t805 = mrSges(5,1) * t690 - mrSges(5,2) * t691 + Ifges(5,5) * t745 + Ifges(5,6) * t746 + Ifges(5,3) * qJDD(4) + pkin(4) * t677 + pkin(10) * t671 + t765 * t662 + t762 * t669 + (t734 * t763 - t735 * t766) * qJD(3);
t804 = mrSges(6,1) * t684 + mrSges(7,1) * t680 - mrSges(6,2) * t685 - mrSges(7,2) * t682 + pkin(5) * t678 + t800 * t720 + t801 * t721 + t807 * t740 - t790 * t741 - t791 * t742;
t646 = t658 * t795 + t760 * t660 + t756 * t796;
t645 = m(3) * t731 + t646;
t633 = -t645 * t757 + t806 * t761;
t631 = m(2) * t747 + t633;
t639 = -t643 * t754 + t758 * t652;
t638 = m(2) * t748 + t639;
t793 = t759 * t631 + t755 * t638;
t632 = t761 * t645 + t806 * t757;
t781 = -t631 * t755 + t759 * t638;
t779 = m(2) * t753 + t632;
t733 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t763 + Ifges(5,6) * t766) * qJD(3);
t648 = mrSges(5,2) * t697 - mrSges(5,3) * t690 + Ifges(5,1) * t745 + Ifges(5,4) * t746 + Ifges(5,5) * qJDD(4) - pkin(10) * t670 - qJD(4) * t734 - t662 * t762 + t669 * t765 + t733 * t787;
t654 = -mrSges(5,1) * t697 + mrSges(5,3) * t691 + Ifges(5,4) * t745 + Ifges(5,2) * t746 + Ifges(5,6) * qJDD(4) - pkin(4) * t670 + qJD(4) * t735 - t733 * t788 - t804;
t635 = mrSges(4,2) * t702 - mrSges(4,3) * t699 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t769 - pkin(9) * t661 + t648 * t766 - t654 * t763;
t640 = -mrSges(4,1) * t702 + mrSges(4,3) * t700 + t769 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t661 - t805;
t774 = pkin(8) * t653 + t635 * t764 + t640 * t767;
t634 = mrSges(4,1) * t699 - mrSges(4,2) * t700 + Ifges(4,3) * qJDD(3) + pkin(3) * t771 + pkin(9) * t780 + t763 * t648 + t766 * t654;
t628 = -mrSges(3,1) * t731 + mrSges(3,3) * t708 - pkin(2) * t646 - t756 * t634 + t774 * t760;
t629 = mrSges(3,2) * t731 - mrSges(3,3) * t707 + t767 * t635 - t764 * t640 + (-t646 * t756 - t647 * t760) * pkin(8);
t773 = qJ(2) * t639 + t628 * t758 + t629 * t754;
t627 = mrSges(3,1) * t707 - mrSges(3,2) * t708 + pkin(2) * t647 + t760 * t634 + t774 * t756;
t626 = mrSges(2,2) * t753 - mrSges(2,3) * t747 - t754 * t628 + t758 * t629 + (-t632 * t757 - t633 * t761) * qJ(2);
t625 = -mrSges(2,1) * t753 + mrSges(2,3) * t748 - pkin(1) * t632 - t757 * t627 + t773 * t761;
t1 = [-m(1) * g(1) + t781; -m(1) * g(2) + t793; -m(1) * g(3) + t779; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t793 - t755 * t625 + t759 * t626; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t781 + t759 * t625 + t755 * t626; -mrSges(1,1) * g(2) + mrSges(2,1) * t747 + mrSges(1,2) * g(1) - mrSges(2,2) * t748 + pkin(1) * t633 + t761 * t627 + t773 * t757; t779; t645; t634; t805; t804; t683;];
tauJB  = t1;
