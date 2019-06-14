% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-05-05 14:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:18
% EndTime: 2019-05-05 14:21:22
% DurationCPUTime: 3.27s
% Computational Cost: add. (38108->288), mult. (76460->341), div. (0->0), fcn. (43015->8), ass. (0->112)
t733 = sin(qJ(1));
t736 = cos(qJ(1));
t708 = t733 * g(1) - t736 * g(2);
t738 = qJD(1) ^ 2;
t684 = -qJDD(1) * pkin(1) - t738 * qJ(2) + qJDD(2) - t708;
t672 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t684;
t709 = -t736 * g(1) - t733 * g(2);
t766 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t709;
t765 = -m(3) - m(4);
t764 = mrSges(2,1) - mrSges(3,2);
t763 = t738 * mrSges(4,3);
t682 = t738 * pkin(1) - t766;
t673 = qJDD(3) + (-pkin(1) - qJ(3)) * t738 + t766;
t669 = -qJDD(1) * pkin(7) + t673;
t732 = sin(qJ(4));
t735 = cos(qJ(4));
t662 = -t735 * g(3) + t732 * t669;
t702 = (mrSges(5,1) * t732 + mrSges(5,2) * t735) * qJD(1);
t759 = qJD(1) * qJD(4);
t753 = t735 * t759;
t703 = t732 * qJDD(1) + t753;
t761 = qJD(1) * t735;
t707 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t761;
t668 = -t738 * pkin(7) - t672;
t754 = t732 * t759;
t704 = t735 * qJDD(1) - t754;
t650 = (-t704 + t754) * qJ(5) + (t703 + t753) * pkin(4) + t668;
t701 = (pkin(4) * t732 - qJ(5) * t735) * qJD(1);
t737 = qJD(4) ^ 2;
t760 = t732 * qJD(1);
t653 = -t737 * pkin(4) + qJDD(4) * qJ(5) - t701 * t760 + t662;
t729 = sin(pkin(9));
t730 = cos(pkin(9));
t695 = t729 * qJD(4) + t730 * t761;
t637 = -0.2e1 * qJD(5) * t695 + t730 * t650 - t729 * t653;
t678 = t729 * qJDD(4) + t730 * t704;
t694 = t730 * qJD(4) - t729 * t761;
t635 = (t694 * t760 - t678) * pkin(8) + (t694 * t695 + t703) * pkin(5) + t637;
t638 = 0.2e1 * qJD(5) * t694 + t729 * t650 + t730 * t653;
t677 = t730 * qJDD(4) - t729 * t704;
t679 = pkin(5) * t760 - t695 * pkin(8);
t693 = t694 ^ 2;
t636 = -t693 * pkin(5) + t677 * pkin(8) - t679 * t760 + t638;
t731 = sin(qJ(6));
t734 = cos(qJ(6));
t633 = t734 * t635 - t731 * t636;
t663 = t734 * t694 - t731 * t695;
t642 = t663 * qJD(6) + t731 * t677 + t734 * t678;
t664 = t731 * t694 + t734 * t695;
t647 = -t663 * mrSges(7,1) + t664 * mrSges(7,2);
t710 = qJD(6) + t760;
t654 = -t710 * mrSges(7,2) + t663 * mrSges(7,3);
t700 = qJDD(6) + t703;
t629 = m(7) * t633 + t700 * mrSges(7,1) - t642 * mrSges(7,3) - t664 * t647 + t710 * t654;
t634 = t731 * t635 + t734 * t636;
t641 = -t664 * qJD(6) + t734 * t677 - t731 * t678;
t655 = t710 * mrSges(7,1) - t664 * mrSges(7,3);
t630 = m(7) * t634 - t700 * mrSges(7,2) + t641 * mrSges(7,3) + t663 * t647 - t710 * t655;
t621 = t734 * t629 + t731 * t630;
t665 = -t694 * mrSges(6,1) + t695 * mrSges(6,2);
t675 = -mrSges(6,2) * t760 + t694 * mrSges(6,3);
t619 = m(6) * t637 + t703 * mrSges(6,1) - t678 * mrSges(6,3) - t695 * t665 + t675 * t760 + t621;
t676 = mrSges(6,1) * t760 - t695 * mrSges(6,3);
t749 = -t731 * t629 + t734 * t630;
t620 = m(6) * t638 - t703 * mrSges(6,2) + t677 * mrSges(6,3) + t694 * t665 - t676 * t760 + t749;
t750 = -t729 * t619 + t730 * t620;
t613 = m(5) * t662 - qJDD(4) * mrSges(5,2) - t703 * mrSges(5,3) - qJD(4) * t707 - t702 * t760 + t750;
t661 = t732 * g(3) + t735 * t669;
t652 = -qJDD(4) * pkin(4) - t737 * qJ(5) + t701 * t761 + qJDD(5) - t661;
t639 = -t677 * pkin(5) - t693 * pkin(8) + t695 * t679 + t652;
t743 = m(7) * t639 - t641 * mrSges(7,1) + t642 * mrSges(7,2) - t663 * t654 + t664 * t655;
t632 = m(6) * t652 - t677 * mrSges(6,1) + t678 * mrSges(6,2) - t694 * t675 + t695 * t676 + t743;
t706 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t760;
t625 = m(5) * t661 + qJDD(4) * mrSges(5,1) - t704 * mrSges(5,3) + qJD(4) * t706 - t702 * t761 - t632;
t603 = t732 * t613 + t735 * t625;
t748 = m(4) * t673 + qJDD(1) * mrSges(4,2) + t603;
t744 = -m(3) * t682 + t738 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t748;
t597 = m(2) * t709 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t738 + t744;
t615 = t730 * t619 + t729 * t620;
t745 = -m(5) * t668 - t703 * mrSges(5,1) - t704 * mrSges(5,2) - t706 * t760 - t707 * t761 - t615;
t610 = m(4) * t672 - t738 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t745;
t741 = -m(3) * t684 + t738 * mrSges(3,3) - t610;
t607 = m(2) * t708 - t738 * mrSges(2,2) + t764 * qJDD(1) + t741;
t762 = t733 * t597 + t736 * t607;
t756 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t755 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t752 = t736 * t597 - t733 * t607;
t751 = t735 * t613 - t732 * t625;
t643 = Ifges(7,5) * t664 + Ifges(7,6) * t663 + Ifges(7,3) * t710;
t645 = Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t710;
t622 = -mrSges(7,1) * t639 + mrSges(7,3) * t634 + Ifges(7,4) * t642 + Ifges(7,2) * t641 + Ifges(7,6) * t700 - t664 * t643 + t710 * t645;
t644 = Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t710;
t623 = mrSges(7,2) * t639 - mrSges(7,3) * t633 + Ifges(7,1) * t642 + Ifges(7,4) * t641 + Ifges(7,5) * t700 + t663 * t643 - t710 * t644;
t657 = Ifges(6,5) * t695 + Ifges(6,6) * t694 + Ifges(6,3) * t760;
t659 = Ifges(6,1) * t695 + Ifges(6,4) * t694 + Ifges(6,5) * t760;
t601 = -mrSges(6,1) * t652 + mrSges(6,3) * t638 + Ifges(6,4) * t678 + Ifges(6,2) * t677 + Ifges(6,6) * t703 - pkin(5) * t743 + pkin(8) * t749 + t734 * t622 + t731 * t623 - t695 * t657 + t659 * t760;
t658 = Ifges(6,4) * t695 + Ifges(6,2) * t694 + Ifges(6,6) * t760;
t605 = mrSges(6,2) * t652 - mrSges(6,3) * t637 + Ifges(6,1) * t678 + Ifges(6,4) * t677 + Ifges(6,5) * t703 - pkin(8) * t621 - t731 * t622 + t734 * t623 + t694 * t657 - t658 * t760;
t689 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t735 - Ifges(5,2) * t732) * qJD(1);
t690 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t735 - Ifges(5,4) * t732) * qJD(1);
t742 = mrSges(5,1) * t661 - mrSges(5,2) * t662 + Ifges(5,5) * t704 - Ifges(5,6) * t703 + Ifges(5,3) * qJDD(4) - pkin(4) * t632 + qJ(5) * t750 + t730 * t601 + t729 * t605 + t689 * t761 + t690 * t760;
t740 = mrSges(7,1) * t633 - mrSges(7,2) * t634 + Ifges(7,5) * t642 + Ifges(7,6) * t641 + Ifges(7,3) * t700 + t664 * t644 - t663 * t645;
t688 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t735 - Ifges(5,6) * t732) * qJD(1);
t593 = mrSges(5,2) * t668 - mrSges(5,3) * t661 + Ifges(5,1) * t704 - Ifges(5,4) * t703 + Ifges(5,5) * qJDD(4) - qJ(5) * t615 - qJD(4) * t689 - t729 * t601 + t730 * t605 - t688 * t760;
t594 = -t688 * t761 + (-Ifges(5,2) - Ifges(6,3)) * t703 + Ifges(5,4) * t704 + t694 * t659 - t695 * t658 - Ifges(6,6) * t677 - Ifges(6,5) * t678 + qJD(4) * t690 + mrSges(5,3) * t662 - mrSges(5,1) * t668 - mrSges(6,1) * t637 + mrSges(6,2) * t638 - pkin(5) * t621 - pkin(4) * t615 + Ifges(5,6) * qJDD(4) - t740;
t609 = qJDD(1) * mrSges(3,2) - t741;
t739 = -mrSges(2,2) * t709 - mrSges(3,3) * t682 - mrSges(4,3) * t672 - pkin(7) * t603 - qJ(3) * t610 + t735 * t593 - t732 * t594 + qJ(2) * (t744 - t763) - pkin(1) * t609 + mrSges(4,2) * t673 + mrSges(3,2) * t684 + mrSges(2,1) * t708 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t602 = t765 * g(3) + t751;
t599 = t748 - t763;
t591 = t742 + t756 * t738 - qJ(3) * t751 + t755 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t764) * g(3) + mrSges(2,3) * t709 - mrSges(3,1) * t682 + mrSges(4,1) * t673 + pkin(3) * t603 - pkin(1) * t602 + pkin(2) * t599;
t590 = -qJ(2) * t602 - mrSges(2,3) * t708 + pkin(2) * t610 + mrSges(3,1) * t684 + pkin(7) * t751 + t732 * t593 + t735 * t594 + pkin(3) * t745 + mrSges(4,1) * t672 - t755 * t738 + t756 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t752; -m(1) * g(2) + t762; (-m(1) - m(2) + t765) * g(3) + t751; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t762 + t736 * t590 - t733 * t591; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t752 + t733 * t590 + t736 * t591; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t739; t739; t609; t599; t742; t632; t740;];
tauJB  = t1;
