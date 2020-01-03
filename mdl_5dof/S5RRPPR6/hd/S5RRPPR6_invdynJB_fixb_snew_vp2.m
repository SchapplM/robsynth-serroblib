% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:32:11
% EndTime: 2019-12-31 19:32:19
% DurationCPUTime: 7.27s
% Computational Cost: add. (98270->312), mult. (231357->399), div. (0->0), fcn. (157898->10), ass. (0->123)
t768 = -2 * qJD(3);
t738 = sin(qJ(2));
t741 = cos(qJ(2));
t758 = qJD(1) * qJD(2);
t722 = t738 * qJDD(1) + t741 * t758;
t739 = sin(qJ(1));
t742 = cos(qJ(1));
t729 = -t742 * g(1) - t739 * g(2);
t744 = qJD(1) ^ 2;
t717 = -t744 * pkin(1) + qJDD(1) * pkin(6) + t729;
t763 = t738 * t717;
t766 = pkin(2) * t744;
t677 = qJDD(2) * pkin(2) - t722 * qJ(3) - t763 + (qJ(3) * t758 + t738 * t766 - g(3)) * t741;
t703 = -t738 * g(3) + t741 * t717;
t723 = t741 * qJDD(1) - t738 * t758;
t761 = qJD(1) * t738;
t725 = qJD(2) * pkin(2) - qJ(3) * t761;
t733 = t741 ^ 2;
t679 = t723 * qJ(3) - qJD(2) * t725 - t733 * t766 + t703;
t735 = sin(pkin(8));
t764 = cos(pkin(8));
t711 = (t735 * t741 + t764 * t738) * qJD(1);
t663 = t764 * t677 - t735 * t679 + t711 * t768;
t760 = qJD(1) * t741;
t710 = t735 * t761 - t764 * t760;
t664 = t735 * t677 + t764 * t679 + t710 * t768;
t691 = t710 * pkin(3) - t711 * qJ(4);
t743 = qJD(2) ^ 2;
t652 = -t743 * pkin(3) + qJDD(2) * qJ(4) - t710 * t691 + t664;
t728 = t739 * g(1) - t742 * g(2);
t750 = -qJDD(1) * pkin(1) - t728;
t681 = -t723 * pkin(2) + qJDD(3) + t725 * t761 + (-qJ(3) * t733 - pkin(6)) * t744 + t750;
t695 = t735 * t722 - t764 * t723;
t696 = t764 * t722 + t735 * t723;
t655 = (qJD(2) * t710 - t696) * qJ(4) + (qJD(2) * t711 + t695) * pkin(3) + t681;
t734 = sin(pkin(9));
t736 = cos(pkin(9));
t701 = t734 * qJD(2) + t736 * t711;
t647 = -0.2e1 * qJD(4) * t701 - t734 * t652 + t736 * t655;
t686 = t734 * qJDD(2) + t736 * t696;
t700 = t736 * qJD(2) - t734 * t711;
t645 = (t710 * t700 - t686) * pkin(7) + (t700 * t701 + t695) * pkin(4) + t647;
t648 = 0.2e1 * qJD(4) * t700 + t736 * t652 + t734 * t655;
t683 = t710 * pkin(4) - t701 * pkin(7);
t685 = t736 * qJDD(2) - t734 * t696;
t699 = t700 ^ 2;
t646 = -t699 * pkin(4) + t685 * pkin(7) - t710 * t683 + t648;
t737 = sin(qJ(5));
t740 = cos(qJ(5));
t643 = t740 * t645 - t737 * t646;
t673 = t740 * t700 - t737 * t701;
t658 = t673 * qJD(5) + t737 * t685 + t740 * t686;
t674 = t737 * t700 + t740 * t701;
t665 = -t673 * mrSges(6,1) + t674 * mrSges(6,2);
t709 = qJD(5) + t710;
t666 = -t709 * mrSges(6,2) + t673 * mrSges(6,3);
t694 = qJDD(5) + t695;
t640 = m(6) * t643 + t694 * mrSges(6,1) - t658 * mrSges(6,3) - t674 * t665 + t709 * t666;
t644 = t737 * t645 + t740 * t646;
t657 = -t674 * qJD(5) + t740 * t685 - t737 * t686;
t667 = t709 * mrSges(6,1) - t674 * mrSges(6,3);
t641 = m(6) * t644 - t694 * mrSges(6,2) + t657 * mrSges(6,3) + t673 * t665 - t709 * t667;
t632 = t740 * t640 + t737 * t641;
t678 = -t700 * mrSges(5,1) + t701 * mrSges(5,2);
t753 = -t710 * mrSges(5,2) + t700 * mrSges(5,3);
t630 = m(5) * t647 + t695 * mrSges(5,1) - t686 * mrSges(5,3) - t701 * t678 + t710 * t753 + t632;
t682 = t710 * mrSges(5,1) - t701 * mrSges(5,3);
t754 = -t737 * t640 + t740 * t641;
t631 = m(5) * t648 - t695 * mrSges(5,2) + t685 * mrSges(5,3) + t700 * t678 - t710 * t682 + t754;
t626 = -t734 * t630 + t736 * t631;
t692 = t710 * mrSges(4,1) + t711 * mrSges(4,2);
t705 = qJD(2) * mrSges(4,1) - t711 * mrSges(4,3);
t623 = m(4) * t664 - qJDD(2) * mrSges(4,2) - t695 * mrSges(4,3) - qJD(2) * t705 - t710 * t692 + t626;
t651 = -qJDD(2) * pkin(3) - t743 * qJ(4) + t711 * t691 + qJDD(4) - t663;
t649 = -t685 * pkin(4) - t699 * pkin(7) + t701 * t683 + t651;
t748 = m(6) * t649 - t657 * mrSges(6,1) + t658 * mrSges(6,2) - t673 * t666 + t674 * t667;
t642 = m(5) * t651 - t685 * mrSges(5,1) + t686 * mrSges(5,2) + t701 * t682 - t700 * t753 + t748;
t704 = -qJD(2) * mrSges(4,2) - t710 * mrSges(4,3);
t636 = m(4) * t663 + qJDD(2) * mrSges(4,1) - t696 * mrSges(4,3) + qJD(2) * t704 - t711 * t692 - t642;
t615 = t735 * t623 + t764 * t636;
t659 = Ifges(6,5) * t674 + Ifges(6,6) * t673 + Ifges(6,3) * t709;
t661 = Ifges(6,1) * t674 + Ifges(6,4) * t673 + Ifges(6,5) * t709;
t633 = -mrSges(6,1) * t649 + mrSges(6,3) * t644 + Ifges(6,4) * t658 + Ifges(6,2) * t657 + Ifges(6,6) * t694 - t674 * t659 + t709 * t661;
t660 = Ifges(6,4) * t674 + Ifges(6,2) * t673 + Ifges(6,6) * t709;
t634 = mrSges(6,2) * t649 - mrSges(6,3) * t643 + Ifges(6,1) * t658 + Ifges(6,4) * t657 + Ifges(6,5) * t694 + t673 * t659 - t709 * t660;
t668 = Ifges(5,5) * t701 + Ifges(5,6) * t700 + Ifges(5,3) * t710;
t670 = Ifges(5,1) * t701 + Ifges(5,4) * t700 + Ifges(5,5) * t710;
t616 = -mrSges(5,1) * t651 + mrSges(5,3) * t648 + Ifges(5,4) * t686 + Ifges(5,2) * t685 + Ifges(5,6) * t695 - pkin(4) * t748 + pkin(7) * t754 + t740 * t633 + t737 * t634 - t701 * t668 + t710 * t670;
t669 = Ifges(5,4) * t701 + Ifges(5,2) * t700 + Ifges(5,6) * t710;
t617 = mrSges(5,2) * t651 - mrSges(5,3) * t647 + Ifges(5,1) * t686 + Ifges(5,4) * t685 + Ifges(5,5) * t695 - pkin(7) * t632 - t737 * t633 + t740 * t634 + t700 * t668 - t710 * t669;
t688 = Ifges(4,4) * t711 - Ifges(4,2) * t710 + Ifges(4,6) * qJD(2);
t689 = Ifges(4,1) * t711 - Ifges(4,4) * t710 + Ifges(4,5) * qJD(2);
t702 = -t741 * g(3) - t763;
t713 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t738 + Ifges(3,2) * t741) * qJD(1);
t714 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t738 + Ifges(3,4) * t741) * qJD(1);
t767 = (t738 * t713 - t741 * t714) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t702 + mrSges(4,1) * t663 - mrSges(3,2) * t703 - mrSges(4,2) * t664 + Ifges(3,5) * t722 + Ifges(4,5) * t696 + Ifges(3,6) * t723 - Ifges(4,6) * t695 + pkin(2) * t615 - pkin(3) * t642 + qJ(4) * t626 + t736 * t616 + t734 * t617 + t711 * t688 + t710 * t689;
t721 = (-mrSges(3,1) * t741 + mrSges(3,2) * t738) * qJD(1);
t727 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t760;
t613 = m(3) * t702 + qJDD(2) * mrSges(3,1) - t722 * mrSges(3,3) + qJD(2) * t727 - t721 * t761 + t615;
t726 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t761;
t755 = t764 * t623 - t735 * t636;
t614 = m(3) * t703 - qJDD(2) * mrSges(3,2) + t723 * mrSges(3,3) - qJD(2) * t726 + t721 * t760 + t755;
t756 = -t738 * t613 + t741 * t614;
t606 = m(2) * t729 - t744 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t756;
t625 = t736 * t630 + t734 * t631;
t624 = m(4) * t681 + t695 * mrSges(4,1) + t696 * mrSges(4,2) + t710 * t704 + t711 * t705 + t625;
t716 = -t744 * pkin(6) + t750;
t746 = -m(3) * t716 + t723 * mrSges(3,1) - t722 * mrSges(3,2) - t726 * t761 + t727 * t760 - t624;
t619 = m(2) * t728 + qJDD(1) * mrSges(2,1) - t744 * mrSges(2,2) + t746;
t762 = t739 * t606 + t742 * t619;
t608 = t741 * t613 + t738 * t614;
t757 = t742 * t606 - t739 * t619;
t687 = Ifges(4,5) * t711 - Ifges(4,6) * t710 + Ifges(4,3) * qJD(2);
t603 = mrSges(4,2) * t681 - mrSges(4,3) * t663 + Ifges(4,1) * t696 - Ifges(4,4) * t695 + Ifges(4,5) * qJDD(2) - qJ(4) * t625 - qJD(2) * t688 - t734 * t616 + t736 * t617 - t710 * t687;
t747 = mrSges(6,1) * t643 - mrSges(6,2) * t644 + Ifges(6,5) * t658 + Ifges(6,6) * t657 + Ifges(6,3) * t694 + t674 * t660 - t673 * t661;
t609 = (-Ifges(4,2) - Ifges(5,3)) * t695 - t711 * t687 + Ifges(4,4) * t696 + t700 * t670 - t701 * t669 - Ifges(5,6) * t685 - Ifges(5,5) * t686 + qJD(2) * t689 - mrSges(4,1) * t681 + mrSges(4,3) * t664 + Ifges(4,6) * qJDD(2) + mrSges(5,2) * t648 - mrSges(5,1) * t647 - pkin(4) * t632 - t747 - pkin(3) * t625;
t712 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t738 + Ifges(3,6) * t741) * qJD(1);
t599 = -mrSges(3,1) * t716 + mrSges(3,3) * t703 + Ifges(3,4) * t722 + Ifges(3,2) * t723 + Ifges(3,6) * qJDD(2) - pkin(2) * t624 + qJ(3) * t755 + qJD(2) * t714 + t735 * t603 + t764 * t609 - t712 * t761;
t602 = mrSges(3,2) * t716 - mrSges(3,3) * t702 + Ifges(3,1) * t722 + Ifges(3,4) * t723 + Ifges(3,5) * qJDD(2) - qJ(3) * t615 - qJD(2) * t713 + t764 * t603 - t735 * t609 + t712 * t760;
t749 = mrSges(2,1) * t728 - mrSges(2,2) * t729 + Ifges(2,3) * qJDD(1) + pkin(1) * t746 + pkin(6) * t756 + t741 * t599 + t738 * t602;
t600 = mrSges(2,1) * g(3) + mrSges(2,3) * t729 + t744 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t608 - t767;
t597 = -mrSges(2,2) * g(3) - mrSges(2,3) * t728 + Ifges(2,5) * qJDD(1) - t744 * Ifges(2,6) - pkin(6) * t608 - t738 * t599 + t741 * t602;
t1 = [-m(1) * g(1) + t757; -m(1) * g(2) + t762; (-m(1) - m(2)) * g(3) + t608; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t762 + t742 * t597 - t739 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t757 + t739 * t597 + t742 * t600; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t749; t749; t767; t624; t642; t747;];
tauJB = t1;
