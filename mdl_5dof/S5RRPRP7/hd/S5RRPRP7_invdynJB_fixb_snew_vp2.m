% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:00:19
% EndTime: 2019-12-31 20:00:25
% DurationCPUTime: 4.52s
% Computational Cost: add. (45328->289), mult. (102104->357), div. (0->0), fcn. (66993->8), ass. (0->116)
t775 = -2 * qJD(3);
t774 = Ifges(5,1) + Ifges(6,1);
t766 = Ifges(5,4) - Ifges(6,5);
t765 = -Ifges(5,5) - Ifges(6,4);
t773 = Ifges(5,2) + Ifges(6,3);
t764 = Ifges(5,6) - Ifges(6,6);
t772 = -Ifges(5,3) - Ifges(6,2);
t733 = sin(qJ(2));
t735 = cos(qJ(2));
t752 = qJD(1) * qJD(2);
t717 = t733 * qJDD(1) + t735 * t752;
t734 = sin(qJ(1));
t736 = cos(qJ(1));
t724 = -t736 * g(1) - t734 * g(2);
t738 = qJD(1) ^ 2;
t712 = -t738 * pkin(1) + qJDD(1) * pkin(6) + t724;
t761 = t733 * t712;
t768 = pkin(2) * t738;
t672 = qJDD(2) * pkin(2) - t717 * qJ(3) - t761 + (qJ(3) * t752 + t733 * t768 - g(3)) * t735;
t698 = -t733 * g(3) + t735 * t712;
t718 = t735 * qJDD(1) - t733 * t752;
t755 = qJD(1) * t733;
t720 = qJD(2) * pkin(2) - qJ(3) * t755;
t729 = t735 ^ 2;
t673 = t718 * qJ(3) - qJD(2) * t720 - t729 * t768 + t698;
t730 = sin(pkin(8));
t731 = cos(pkin(8));
t707 = (t730 * t735 + t731 * t733) * qJD(1);
t653 = t731 * t672 - t730 * t673 + t707 * t775;
t706 = (t730 * t733 - t731 * t735) * qJD(1);
t654 = t730 * t672 + t731 * t673 + t706 * t775;
t688 = t706 * pkin(3) - t707 * pkin(7);
t737 = qJD(2) ^ 2;
t650 = -t737 * pkin(3) + qJDD(2) * pkin(7) - t706 * t688 + t654;
t723 = t734 * g(1) - t736 * g(2);
t743 = -qJDD(1) * pkin(1) - t723;
t678 = -t718 * pkin(2) + qJDD(3) + t720 * t755 + (-qJ(3) * t729 - pkin(6)) * t738 + t743;
t693 = -t730 * t717 + t731 * t718;
t694 = t731 * t717 + t730 * t718;
t652 = (qJD(2) * t706 - t694) * pkin(7) + (qJD(2) * t707 - t693) * pkin(3) + t678;
t732 = sin(qJ(4));
t769 = cos(qJ(4));
t647 = t769 * t650 + t732 * t652;
t696 = t732 * qJD(2) + t769 * t707;
t665 = t696 * qJD(4) - t769 * qJDD(2) + t732 * t694;
t705 = qJD(4) + t706;
t681 = t705 * mrSges(5,1) - t696 * mrSges(5,3);
t692 = qJDD(4) - t693;
t695 = -t769 * qJD(2) + t732 * t707;
t674 = t695 * pkin(4) - t696 * qJ(5);
t704 = t705 ^ 2;
t643 = -t704 * pkin(4) + t692 * qJ(5) + 0.2e1 * qJD(5) * t705 - t695 * t674 + t647;
t682 = -t705 * mrSges(6,1) + t696 * mrSges(6,2);
t751 = m(6) * t643 + t692 * mrSges(6,3) + t705 * t682;
t675 = t695 * mrSges(6,1) - t696 * mrSges(6,3);
t756 = -t695 * mrSges(5,1) - t696 * mrSges(5,2) - t675;
t767 = -mrSges(5,3) - mrSges(6,2);
t635 = m(5) * t647 - t692 * mrSges(5,2) + t767 * t665 - t705 * t681 + t756 * t695 + t751;
t646 = -t732 * t650 + t769 * t652;
t666 = -t695 * qJD(4) + t732 * qJDD(2) + t769 * t694;
t680 = -t705 * mrSges(5,2) - t695 * mrSges(5,3);
t644 = -t692 * pkin(4) - t704 * qJ(5) + t696 * t674 + qJDD(5) - t646;
t679 = -t695 * mrSges(6,2) + t705 * mrSges(6,3);
t746 = -m(6) * t644 + t692 * mrSges(6,1) + t705 * t679;
t637 = m(5) * t646 + t692 * mrSges(5,1) + t767 * t666 + t705 * t680 + t756 * t696 + t746;
t630 = t769 * t635 - t732 * t637;
t687 = t706 * mrSges(4,1) + t707 * mrSges(4,2);
t700 = qJD(2) * mrSges(4,1) - t707 * mrSges(4,3);
t626 = m(4) * t654 - qJDD(2) * mrSges(4,2) + t693 * mrSges(4,3) - qJD(2) * t700 - t706 * t687 + t630;
t649 = -qJDD(2) * pkin(3) - t737 * pkin(7) + t707 * t688 - t653;
t645 = -0.2e1 * qJD(5) * t696 + (t695 * t705 - t666) * qJ(5) + (t696 * t705 + t665) * pkin(4) + t649;
t641 = m(6) * t645 + t665 * mrSges(6,1) - t666 * mrSges(6,3) + t695 * t679 - t696 * t682;
t638 = -m(5) * t649 - t665 * mrSges(5,1) - t666 * mrSges(5,2) - t695 * t680 - t696 * t681 - t641;
t699 = -qJD(2) * mrSges(4,2) - t706 * mrSges(4,3);
t632 = m(4) * t653 + qJDD(2) * mrSges(4,1) - t694 * mrSges(4,3) + qJD(2) * t699 - t707 * t687 + t638;
t619 = t730 * t626 + t731 * t632;
t757 = t766 * t695 - t774 * t696 + t765 * t705;
t759 = t764 * t695 + t765 * t696 + t772 * t705;
t623 = -mrSges(5,1) * t649 - mrSges(6,1) * t645 + mrSges(6,2) * t643 + mrSges(5,3) * t647 - pkin(4) * t641 - t773 * t665 + t766 * t666 + t764 * t692 + t759 * t696 - t757 * t705;
t758 = t773 * t695 - t766 * t696 - t764 * t705;
t627 = mrSges(5,2) * t649 + mrSges(6,2) * t644 - mrSges(5,3) * t646 - mrSges(6,3) * t645 - qJ(5) * t641 - t766 * t665 + t774 * t666 - t765 * t692 + t759 * t695 + t758 * t705;
t684 = Ifges(4,4) * t707 - Ifges(4,2) * t706 + Ifges(4,6) * qJD(2);
t685 = Ifges(4,1) * t707 - Ifges(4,4) * t706 + Ifges(4,5) * qJD(2);
t697 = -t735 * g(3) - t761;
t709 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t733 + Ifges(3,2) * t735) * qJD(1);
t710 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t733 + Ifges(3,4) * t735) * qJD(1);
t771 = (t733 * t709 - t735 * t710) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t697 + mrSges(4,1) * t653 - mrSges(3,2) * t698 - mrSges(4,2) * t654 + Ifges(3,5) * t717 + Ifges(4,5) * t694 + Ifges(3,6) * t718 + Ifges(4,6) * t693 + pkin(2) * t619 + pkin(3) * t638 + pkin(7) * t630 + t769 * t623 + t732 * t627 + t707 * t684 + t706 * t685;
t640 = t666 * mrSges(6,2) + t696 * t675 - t746;
t770 = -t764 * t665 - t765 * t666 - t772 * t692 - t757 * t695 - t758 * t696 + mrSges(5,1) * t646 - mrSges(6,1) * t644 - mrSges(5,2) * t647 + mrSges(6,3) * t643 - pkin(4) * t640 + qJ(5) * (-t665 * mrSges(6,2) - t695 * t675 + t751);
t716 = (-mrSges(3,1) * t735 + mrSges(3,2) * t733) * qJD(1);
t754 = qJD(1) * t735;
t722 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t754;
t617 = m(3) * t697 + qJDD(2) * mrSges(3,1) - t717 * mrSges(3,3) + qJD(2) * t722 - t716 * t755 + t619;
t721 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t755;
t748 = t731 * t626 - t730 * t632;
t618 = m(3) * t698 - qJDD(2) * mrSges(3,2) + t718 * mrSges(3,3) - qJD(2) * t721 + t716 * t754 + t748;
t749 = -t733 * t617 + t735 * t618;
t609 = m(2) * t724 - t738 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t749;
t629 = t732 * t635 + t769 * t637;
t628 = m(4) * t678 - t693 * mrSges(4,1) + t694 * mrSges(4,2) + t706 * t699 + t707 * t700 + t629;
t711 = -t738 * pkin(6) + t743;
t740 = -m(3) * t711 + t718 * mrSges(3,1) - t717 * mrSges(3,2) - t721 * t755 + t722 * t754 - t628;
t621 = m(2) * t723 + qJDD(1) * mrSges(2,1) - t738 * mrSges(2,2) + t740;
t760 = t734 * t609 + t736 * t621;
t611 = t735 * t617 + t733 * t618;
t750 = t736 * t609 - t734 * t621;
t683 = Ifges(4,5) * t707 - Ifges(4,6) * t706 + Ifges(4,3) * qJD(2);
t612 = mrSges(4,2) * t678 - mrSges(4,3) * t653 + Ifges(4,1) * t694 + Ifges(4,4) * t693 + Ifges(4,5) * qJDD(2) - pkin(7) * t629 - qJD(2) * t684 - t732 * t623 + t769 * t627 - t706 * t683;
t613 = -mrSges(4,1) * t678 + mrSges(4,3) * t654 + Ifges(4,4) * t694 + Ifges(4,2) * t693 + Ifges(4,6) * qJDD(2) - pkin(3) * t629 + qJD(2) * t685 - t707 * t683 - t770;
t708 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t733 + Ifges(3,6) * t735) * qJD(1);
t604 = -mrSges(3,1) * t711 + mrSges(3,3) * t698 + Ifges(3,4) * t717 + Ifges(3,2) * t718 + Ifges(3,6) * qJDD(2) - pkin(2) * t628 + qJ(3) * t748 + qJD(2) * t710 + t730 * t612 + t731 * t613 - t708 * t755;
t606 = mrSges(3,2) * t711 - mrSges(3,3) * t697 + Ifges(3,1) * t717 + Ifges(3,4) * t718 + Ifges(3,5) * qJDD(2) - qJ(3) * t619 - qJD(2) * t709 + t731 * t612 - t730 * t613 + t708 * t754;
t742 = mrSges(2,1) * t723 - mrSges(2,2) * t724 + Ifges(2,3) * qJDD(1) + pkin(1) * t740 + pkin(6) * t749 + t735 * t604 + t733 * t606;
t602 = mrSges(2,1) * g(3) + mrSges(2,3) * t724 + t738 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t611 - t771;
t601 = -mrSges(2,2) * g(3) - mrSges(2,3) * t723 + Ifges(2,5) * qJDD(1) - t738 * Ifges(2,6) - pkin(6) * t611 - t733 * t604 + t735 * t606;
t1 = [-m(1) * g(1) + t750; -m(1) * g(2) + t760; (-m(1) - m(2)) * g(3) + t611; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t760 + t736 * t601 - t734 * t602; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t750 + t734 * t601 + t736 * t602; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t742; t742; t771; t628; t770; t640;];
tauJB = t1;
