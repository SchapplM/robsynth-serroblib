% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:23:04
% EndTime: 2019-12-31 21:23:14
% DurationCPUTime: 9.23s
% Computational Cost: add. (143116->312), mult. (298124->394), div. (0->0), fcn. (206448->10), ass. (0->124)
t738 = sin(qJ(1));
t742 = cos(qJ(1));
t725 = t738 * g(1) - t742 * g(2);
t744 = qJD(1) ^ 2;
t709 = -qJDD(1) * pkin(1) - t744 * pkin(6) - t725;
t737 = sin(qJ(2));
t741 = cos(qJ(2));
t757 = qJD(1) * qJD(2);
t756 = t741 * t757;
t720 = t737 * qJDD(1) + t756;
t729 = t737 * t757;
t721 = t741 * qJDD(1) - t729;
t676 = (-t720 - t756) * pkin(7) + (-t721 + t729) * pkin(2) + t709;
t726 = -t742 * g(1) - t738 * g(2);
t710 = -t744 * pkin(1) + qJDD(1) * pkin(6) + t726;
t701 = -t737 * g(3) + t741 * t710;
t719 = (-pkin(2) * t741 - pkin(7) * t737) * qJD(1);
t743 = qJD(2) ^ 2;
t758 = t741 * qJD(1);
t679 = -t743 * pkin(2) + qJDD(2) * pkin(7) + t719 * t758 + t701;
t736 = sin(qJ(3));
t740 = cos(qJ(3));
t660 = t740 * t676 - t736 * t679;
t759 = qJD(1) * t737;
t716 = t740 * qJD(2) - t736 * t759;
t692 = t716 * qJD(3) + t736 * qJDD(2) + t740 * t720;
t715 = qJDD(3) - t721;
t717 = t736 * qJD(2) + t740 * t759;
t728 = qJD(3) - t758;
t650 = (t716 * t728 - t692) * qJ(4) + (t716 * t717 + t715) * pkin(3) + t660;
t661 = t736 * t676 + t740 * t679;
t691 = -t717 * qJD(3) + t740 * qJDD(2) - t736 * t720;
t698 = t728 * pkin(3) - t717 * qJ(4);
t714 = t716 ^ 2;
t652 = -t714 * pkin(3) + t691 * qJ(4) - t728 * t698 + t661;
t733 = sin(pkin(9));
t734 = cos(pkin(9));
t695 = t733 * t716 + t734 * t717;
t637 = -0.2e1 * qJD(4) * t695 + t734 * t650 - t733 * t652;
t669 = t733 * t691 + t734 * t692;
t694 = t734 * t716 - t733 * t717;
t635 = (t694 * t728 - t669) * pkin(8) + (t694 * t695 + t715) * pkin(4) + t637;
t638 = 0.2e1 * qJD(4) * t694 + t733 * t650 + t734 * t652;
t668 = t734 * t691 - t733 * t692;
t682 = t728 * pkin(4) - t695 * pkin(8);
t693 = t694 ^ 2;
t636 = -t693 * pkin(4) + t668 * pkin(8) - t728 * t682 + t638;
t735 = sin(qJ(5));
t739 = cos(qJ(5));
t633 = t739 * t635 - t735 * t636;
t671 = t739 * t694 - t735 * t695;
t646 = t671 * qJD(5) + t735 * t668 + t739 * t669;
t672 = t735 * t694 + t739 * t695;
t658 = -t671 * mrSges(6,1) + t672 * mrSges(6,2);
t727 = qJD(5) + t728;
t662 = -t727 * mrSges(6,2) + t671 * mrSges(6,3);
t711 = qJDD(5) + t715;
t628 = m(6) * t633 + t711 * mrSges(6,1) - t646 * mrSges(6,3) - t672 * t658 + t727 * t662;
t634 = t735 * t635 + t739 * t636;
t645 = -t672 * qJD(5) + t739 * t668 - t735 * t669;
t663 = t727 * mrSges(6,1) - t672 * mrSges(6,3);
t629 = m(6) * t634 - t711 * mrSges(6,2) + t645 * mrSges(6,3) + t671 * t658 - t727 * t663;
t620 = t739 * t628 + t735 * t629;
t673 = -t694 * mrSges(5,1) + t695 * mrSges(5,2);
t680 = -t728 * mrSges(5,2) + t694 * mrSges(5,3);
t618 = m(5) * t637 + t715 * mrSges(5,1) - t669 * mrSges(5,3) - t695 * t673 + t728 * t680 + t620;
t681 = t728 * mrSges(5,1) - t695 * mrSges(5,3);
t752 = -t735 * t628 + t739 * t629;
t619 = m(5) * t638 - t715 * mrSges(5,2) + t668 * mrSges(5,3) + t694 * t673 - t728 * t681 + t752;
t614 = t734 * t618 + t733 * t619;
t666 = Ifges(5,4) * t695 + Ifges(5,2) * t694 + Ifges(5,6) * t728;
t667 = Ifges(5,1) * t695 + Ifges(5,4) * t694 + Ifges(5,5) * t728;
t684 = Ifges(4,4) * t717 + Ifges(4,2) * t716 + Ifges(4,6) * t728;
t685 = Ifges(4,1) * t717 + Ifges(4,4) * t716 + Ifges(4,5) * t728;
t654 = Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t727;
t655 = Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t727;
t748 = -mrSges(6,1) * t633 + mrSges(6,2) * t634 - Ifges(6,5) * t646 - Ifges(6,6) * t645 - Ifges(6,3) * t711 - t672 * t654 + t671 * t655;
t763 = mrSges(4,1) * t660 + mrSges(5,1) * t637 - mrSges(4,2) * t661 - mrSges(5,2) * t638 + Ifges(4,5) * t692 + Ifges(5,5) * t669 + Ifges(4,6) * t691 + Ifges(5,6) * t668 + pkin(3) * t614 + pkin(4) * t620 + t695 * t666 - t694 * t667 + t717 * t684 - t716 * t685 + (Ifges(4,3) + Ifges(5,3)) * t715 - t748;
t700 = -t741 * g(3) - t737 * t710;
t678 = -qJDD(2) * pkin(2) - t743 * pkin(7) + t719 * t759 - t700;
t659 = -t691 * pkin(3) - t714 * qJ(4) + t717 * t698 + qJDD(4) + t678;
t640 = -t668 * pkin(4) - t693 * pkin(8) + t695 * t682 + t659;
t653 = Ifges(6,5) * t672 + Ifges(6,6) * t671 + Ifges(6,3) * t727;
t621 = -mrSges(6,1) * t640 + mrSges(6,3) * t634 + Ifges(6,4) * t646 + Ifges(6,2) * t645 + Ifges(6,6) * t711 - t672 * t653 + t727 * t655;
t622 = mrSges(6,2) * t640 - mrSges(6,3) * t633 + Ifges(6,1) * t646 + Ifges(6,4) * t645 + Ifges(6,5) * t711 + t671 * t653 - t727 * t654;
t665 = Ifges(5,5) * t695 + Ifges(5,6) * t694 + Ifges(5,3) * t728;
t751 = m(6) * t640 - t645 * mrSges(6,1) + t646 * mrSges(6,2) - t671 * t662 + t672 * t663;
t609 = -mrSges(5,1) * t659 + mrSges(5,3) * t638 + Ifges(5,4) * t669 + Ifges(5,2) * t668 + Ifges(5,6) * t715 - pkin(4) * t751 + pkin(8) * t752 + t739 * t621 + t735 * t622 - t695 * t665 + t728 * t667;
t610 = mrSges(5,2) * t659 - mrSges(5,3) * t637 + Ifges(5,1) * t669 + Ifges(5,4) * t668 + Ifges(5,5) * t715 - pkin(8) * t620 - t735 * t621 + t739 * t622 + t694 * t665 - t728 * t666;
t631 = m(5) * t659 - t668 * mrSges(5,1) + t669 * mrSges(5,2) - t694 * t680 + t695 * t681 + t751;
t683 = Ifges(4,5) * t717 + Ifges(4,6) * t716 + Ifges(4,3) * t728;
t753 = -t733 * t618 + t734 * t619;
t594 = -mrSges(4,1) * t678 + mrSges(4,3) * t661 + Ifges(4,4) * t692 + Ifges(4,2) * t691 + Ifges(4,6) * t715 - pkin(3) * t631 + qJ(4) * t753 + t734 * t609 + t733 * t610 - t717 * t683 + t728 * t685;
t595 = mrSges(4,2) * t678 - mrSges(4,3) * t660 + Ifges(4,1) * t692 + Ifges(4,4) * t691 + Ifges(4,5) * t715 - qJ(4) * t614 - t733 * t609 + t734 * t610 + t716 * t683 - t728 * t684;
t696 = -t716 * mrSges(4,1) + t717 * mrSges(4,2);
t697 = -t728 * mrSges(4,2) + t716 * mrSges(4,3);
t612 = m(4) * t660 + t715 * mrSges(4,1) - t692 * mrSges(4,3) - t717 * t696 + t728 * t697 + t614;
t699 = t728 * mrSges(4,1) - t717 * mrSges(4,3);
t613 = m(4) * t661 - t715 * mrSges(4,2) + t691 * mrSges(4,3) + t716 * t696 - t728 * t699 + t753;
t608 = -t736 * t612 + t740 * t613;
t630 = -m(4) * t678 + t691 * mrSges(4,1) - t692 * mrSges(4,2) + t716 * t697 - t717 * t699 - t631;
t707 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t737 + Ifges(3,2) * t741) * qJD(1);
t708 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t737 + Ifges(3,4) * t741) * qJD(1);
t762 = mrSges(3,1) * t700 - mrSges(3,2) * t701 + Ifges(3,5) * t720 + Ifges(3,6) * t721 + Ifges(3,3) * qJDD(2) + pkin(2) * t630 + pkin(7) * t608 + t740 * t594 + t736 * t595 + (t737 * t707 - t741 * t708) * qJD(1);
t718 = (-mrSges(3,1) * t741 + mrSges(3,2) * t737) * qJD(1);
t723 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t759;
t606 = m(3) * t701 - qJDD(2) * mrSges(3,2) + t721 * mrSges(3,3) - qJD(2) * t723 + t718 * t758 + t608;
t724 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t758;
t624 = m(3) * t700 + qJDD(2) * mrSges(3,1) - t720 * mrSges(3,3) + qJD(2) * t724 - t718 * t759 + t630;
t754 = t741 * t606 - t737 * t624;
t598 = m(2) * t726 - t744 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t754;
t607 = t740 * t612 + t736 * t613;
t747 = -m(3) * t709 + t721 * mrSges(3,1) - t720 * mrSges(3,2) - t723 * t759 + t724 * t758 - t607;
t602 = m(2) * t725 + qJDD(1) * mrSges(2,1) - t744 * mrSges(2,2) + t747;
t760 = t738 * t598 + t742 * t602;
t600 = t737 * t606 + t741 * t624;
t755 = t742 * t598 - t738 * t602;
t706 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t737 + Ifges(3,6) * t741) * qJD(1);
t591 = mrSges(3,2) * t709 - mrSges(3,3) * t700 + Ifges(3,1) * t720 + Ifges(3,4) * t721 + Ifges(3,5) * qJDD(2) - pkin(7) * t607 - qJD(2) * t707 - t736 * t594 + t740 * t595 + t706 * t758;
t593 = -mrSges(3,1) * t709 + mrSges(3,3) * t701 + Ifges(3,4) * t720 + Ifges(3,2) * t721 + Ifges(3,6) * qJDD(2) - pkin(2) * t607 + qJD(2) * t708 - t706 * t759 - t763;
t749 = mrSges(2,1) * t725 - mrSges(2,2) * t726 + Ifges(2,3) * qJDD(1) + pkin(1) * t747 + pkin(6) * t754 + t737 * t591 + t741 * t593;
t589 = mrSges(2,1) * g(3) + mrSges(2,3) * t726 + t744 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t600 - t762;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t725 + Ifges(2,5) * qJDD(1) - t744 * Ifges(2,6) - pkin(6) * t600 + t741 * t591 - t737 * t593;
t1 = [-m(1) * g(1) + t755; -m(1) * g(2) + t760; (-m(1) - m(2)) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t760 + t742 * t588 - t738 * t589; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t755 + t738 * t588 + t742 * t589; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t749; t749; t762; t763; t631; -t748;];
tauJB = t1;
