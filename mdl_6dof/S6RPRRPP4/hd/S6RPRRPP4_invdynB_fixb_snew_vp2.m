% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:30:51
% EndTime: 2019-05-05 21:31:02
% DurationCPUTime: 10.36s
% Computational Cost: add. (150666->341), mult. (360956->414), div. (0->0), fcn. (270558->10), ass. (0->141)
t749 = Ifges(6,1) + Ifges(7,1);
t743 = Ifges(6,4) - Ifges(7,5);
t742 = Ifges(6,5) + Ifges(7,4);
t748 = Ifges(6,2) + Ifges(7,3);
t747 = -Ifges(7,2) - Ifges(6,3);
t741 = Ifges(6,6) - Ifges(7,6);
t707 = qJD(1) ^ 2;
t698 = sin(pkin(9));
t699 = cos(pkin(9));
t701 = sin(qJ(3));
t704 = cos(qJ(3));
t714 = t698 * t701 - t699 * t704;
t682 = t714 * qJD(1);
t715 = t698 * t704 + t699 * t701;
t683 = t715 * qJD(1);
t730 = t683 * qJD(3);
t668 = -t714 * qJDD(1) - t730;
t746 = -2 * qJD(5);
t745 = pkin(2) * t699;
t744 = -mrSges(6,3) - mrSges(7,2);
t740 = mrSges(3,2) * t698;
t739 = cos(pkin(10));
t696 = t699 ^ 2;
t738 = t696 * t707;
t702 = sin(qJ(1));
t705 = cos(qJ(1));
t688 = -t705 * g(1) - t702 * g(2);
t684 = -t707 * pkin(1) + qJDD(1) * qJ(2) + t688;
t729 = qJD(1) * qJD(2);
t726 = -t699 * g(3) - 0.2e1 * t698 * t729;
t655 = (-pkin(7) * qJDD(1) + t707 * t745 - t684) * t698 + t726;
t671 = -t698 * g(3) + (t684 + 0.2e1 * t729) * t699;
t728 = qJDD(1) * t699;
t656 = -pkin(2) * t738 + pkin(7) * t728 + t671;
t630 = t701 * t655 + t704 * t656;
t663 = t682 * mrSges(4,1) + t683 * mrSges(4,2);
t677 = qJD(3) * mrSges(4,1) - t683 * mrSges(4,3);
t666 = t682 * pkin(3) - t683 * pkin(8);
t706 = qJD(3) ^ 2;
t606 = -t706 * pkin(3) + qJDD(3) * pkin(8) - t682 * t666 + t630;
t695 = t698 ^ 2;
t687 = t702 * g(1) - t705 * g(2);
t719 = qJDD(2) - t687;
t667 = (-pkin(1) - t745) * qJDD(1) + (-qJ(2) + (-t695 - t696) * pkin(7)) * t707 + t719;
t731 = t682 * qJD(3);
t669 = t715 * qJDD(1) - t731;
t614 = (-t669 + t731) * pkin(8) + (-t668 + t730) * pkin(3) + t667;
t700 = sin(qJ(4));
t703 = cos(qJ(4));
t602 = -t700 * t606 + t703 * t614;
t674 = t703 * qJD(3) - t700 * t683;
t643 = t674 * qJD(4) + t700 * qJDD(3) + t703 * t669;
t665 = qJDD(4) - t668;
t675 = t700 * qJD(3) + t703 * t683;
t680 = qJD(4) + t682;
t598 = (t674 * t680 - t643) * qJ(5) + (t674 * t675 + t665) * pkin(4) + t602;
t603 = t703 * t606 + t700 * t614;
t642 = -t675 * qJD(4) + t703 * qJDD(3) - t700 * t669;
t651 = t680 * pkin(4) - t675 * qJ(5);
t673 = t674 ^ 2;
t600 = -t673 * pkin(4) + t642 * qJ(5) - t680 * t651 + t603;
t697 = sin(pkin(10));
t645 = -t739 * t674 + t697 * t675;
t594 = t697 * t598 + t739 * t600 + t645 * t746;
t610 = -t739 * t642 + t697 * t643;
t646 = t697 * t674 + t739 * t675;
t633 = t680 * mrSges(6,1) - t646 * mrSges(6,3);
t624 = t645 * pkin(5) - t646 * qJ(6);
t679 = t680 ^ 2;
t591 = -t679 * pkin(5) + t665 * qJ(6) + 0.2e1 * qJD(6) * t680 - t645 * t624 + t594;
t634 = -t680 * mrSges(7,1) + t646 * mrSges(7,2);
t727 = m(7) * t591 + t665 * mrSges(7,3) + t680 * t634;
t625 = t645 * mrSges(7,1) - t646 * mrSges(7,3);
t733 = -t645 * mrSges(6,1) - t646 * mrSges(6,2) - t625;
t584 = m(6) * t594 - t665 * mrSges(6,2) + t744 * t610 - t680 * t633 + t733 * t645 + t727;
t712 = t739 * t598 - t697 * t600;
t593 = t646 * t746 + t712;
t611 = t697 * t642 + t739 * t643;
t632 = -t680 * mrSges(6,2) - t645 * mrSges(6,3);
t592 = -t665 * pkin(5) - t679 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t624) * t646 - t712;
t631 = -t645 * mrSges(7,2) + t680 * mrSges(7,3);
t720 = -m(7) * t592 + t665 * mrSges(7,1) + t680 * t631;
t586 = m(6) * t593 + t665 * mrSges(6,1) + t744 * t611 + t680 * t632 + t733 * t646 + t720;
t579 = t697 * t584 + t739 * t586;
t647 = -t674 * mrSges(5,1) + t675 * mrSges(5,2);
t650 = -t680 * mrSges(5,2) + t674 * mrSges(5,3);
t577 = m(5) * t602 + t665 * mrSges(5,1) - t643 * mrSges(5,3) - t675 * t647 + t680 * t650 + t579;
t652 = t680 * mrSges(5,1) - t675 * mrSges(5,3);
t721 = t739 * t584 - t697 * t586;
t578 = m(5) * t603 - t665 * mrSges(5,2) + t642 * mrSges(5,3) + t674 * t647 - t680 * t652 + t721;
t722 = -t700 * t577 + t703 * t578;
t572 = m(4) * t630 - qJDD(3) * mrSges(4,2) + t668 * mrSges(4,3) - qJD(3) * t677 - t682 * t663 + t722;
t629 = t704 * t655 - t701 * t656;
t676 = -qJD(3) * mrSges(4,2) - t682 * mrSges(4,3);
t605 = -qJDD(3) * pkin(3) - t706 * pkin(8) + t683 * t666 - t629;
t601 = -t642 * pkin(4) - t673 * qJ(5) + t675 * t651 + qJDD(5) + t605;
t596 = -0.2e1 * qJD(6) * t646 + (t645 * t680 - t611) * qJ(6) + (t646 * t680 + t610) * pkin(5) + t601;
t589 = m(7) * t596 + t610 * mrSges(7,1) - t611 * mrSges(7,3) + t645 * t631 - t646 * t634;
t710 = m(6) * t601 + t610 * mrSges(6,1) + t611 * mrSges(6,2) + t645 * t632 + t646 * t633 + t589;
t708 = -m(5) * t605 + t642 * mrSges(5,1) - t643 * mrSges(5,2) + t674 * t650 - t675 * t652 - t710;
t588 = m(4) * t629 + qJDD(3) * mrSges(4,1) - t669 * mrSges(4,3) + qJD(3) * t676 - t683 * t663 + t708;
t567 = t701 * t572 + t704 * t588;
t670 = -t698 * t684 + t726;
t713 = mrSges(3,3) * qJDD(1) + t707 * (-mrSges(3,1) * t699 + t740);
t565 = m(3) * t670 - t713 * t698 + t567;
t723 = t704 * t572 - t701 * t588;
t566 = m(3) * t671 + t713 * t699 + t723;
t724 = -t698 * t565 + t699 * t566;
t557 = m(2) * t688 - t707 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t724;
t681 = -qJDD(1) * pkin(1) - t707 * qJ(2) + t719;
t573 = t703 * t577 + t700 * t578;
t711 = m(4) * t667 - t668 * mrSges(4,1) + t669 * mrSges(4,2) + t682 * t676 + t683 * t677 + t573;
t709 = -m(3) * t681 + mrSges(3,1) * t728 - t711 + (t695 * t707 + t738) * mrSges(3,3);
t569 = t709 - t707 * mrSges(2,2) + m(2) * t687 + (mrSges(2,1) - t740) * qJDD(1);
t737 = t702 * t557 + t705 * t569;
t558 = t699 * t565 + t698 * t566;
t736 = t645 * t748 - t646 * t743 - t680 * t741;
t735 = t645 * t741 - t646 * t742 + t680 * t747;
t734 = -t743 * t645 + t646 * t749 + t742 * t680;
t716 = Ifges(3,5) * t698 + Ifges(3,6) * t699;
t732 = t707 * t716;
t725 = t705 * t557 - t702 * t569;
t718 = Ifges(3,1) * t698 + Ifges(3,4) * t699;
t717 = Ifges(3,4) * t698 + Ifges(3,2) * t699;
t659 = Ifges(4,1) * t683 - Ifges(4,4) * t682 + Ifges(4,5) * qJD(3);
t658 = Ifges(4,4) * t683 - Ifges(4,2) * t682 + Ifges(4,6) * qJD(3);
t657 = Ifges(4,5) * t683 - Ifges(4,6) * t682 + Ifges(4,3) * qJD(3);
t638 = Ifges(5,1) * t675 + Ifges(5,4) * t674 + Ifges(5,5) * t680;
t637 = Ifges(5,4) * t675 + Ifges(5,2) * t674 + Ifges(5,6) * t680;
t636 = Ifges(5,5) * t675 + Ifges(5,6) * t674 + Ifges(5,3) * t680;
t581 = mrSges(6,2) * t601 + mrSges(7,2) * t592 - mrSges(6,3) * t593 - mrSges(7,3) * t596 - qJ(6) * t589 - t743 * t610 + t611 * t749 + t735 * t645 + t742 * t665 + t736 * t680;
t580 = -mrSges(6,1) * t601 - mrSges(7,1) * t596 + mrSges(7,2) * t591 + mrSges(6,3) * t594 - pkin(5) * t589 - t610 * t748 + t743 * t611 + t735 * t646 + t741 * t665 + t734 * t680;
t561 = mrSges(5,2) * t605 - mrSges(5,3) * t602 + Ifges(5,1) * t643 + Ifges(5,4) * t642 + Ifges(5,5) * t665 - qJ(5) * t579 - t697 * t580 + t739 * t581 + t674 * t636 - t680 * t637;
t560 = -mrSges(5,1) * t605 + mrSges(5,3) * t603 + Ifges(5,4) * t643 + Ifges(5,2) * t642 + Ifges(5,6) * t665 - pkin(4) * t710 + qJ(5) * t721 + t739 * t580 + t697 * t581 - t675 * t636 + t680 * t638;
t559 = Ifges(4,6) * qJDD(3) - pkin(5) * t720 - t683 * t657 + t674 * t638 - t675 * t637 - mrSges(4,1) * t667 + Ifges(4,2) * t668 + Ifges(4,4) * t669 + qJD(3) * t659 - Ifges(5,6) * t642 - Ifges(5,5) * t643 + mrSges(4,3) * t630 + mrSges(5,2) * t603 - mrSges(5,1) * t602 + mrSges(7,1) * t592 - mrSges(6,1) * t593 + mrSges(6,2) * t594 - mrSges(7,3) * t591 + (-Ifges(5,3) + t747) * t665 - qJ(6) * t727 + (qJ(6) * mrSges(7,2) + t741) * t610 + (pkin(5) * mrSges(7,2) - t742) * t611 - pkin(4) * t579 + (qJ(6) * t625 - t734) * t645 + (pkin(5) * t625 + t736) * t646 - pkin(3) * t573;
t554 = mrSges(4,2) * t667 - mrSges(4,3) * t629 + Ifges(4,1) * t669 + Ifges(4,4) * t668 + Ifges(4,5) * qJDD(3) - pkin(8) * t573 - qJD(3) * t658 - t700 * t560 + t703 * t561 - t682 * t657;
t553 = mrSges(3,2) * t681 - mrSges(3,3) * t670 - pkin(7) * t567 + t718 * qJDD(1) + t704 * t554 - t701 * t559 + t699 * t732;
t552 = mrSges(2,1) * g(3) - pkin(1) * t558 + mrSges(2,3) * t688 - pkin(2) * t567 - mrSges(3,1) * t670 + mrSges(3,2) * t671 - t700 * t561 - t703 * t560 - pkin(3) * t708 - pkin(8) * t722 - mrSges(4,1) * t629 + mrSges(4,2) * t630 - Ifges(4,5) * t669 - Ifges(4,6) * t668 - Ifges(4,3) * qJDD(3) - t683 * t658 - t682 * t659 + (Ifges(2,6) - t716) * qJDD(1) + (-t698 * t717 + t699 * t718 + Ifges(2,5)) * t707;
t551 = -mrSges(3,1) * t681 + mrSges(3,3) * t671 - pkin(2) * t711 + pkin(7) * t723 + t717 * qJDD(1) + t701 * t554 + t704 * t559 - t698 * t732;
t550 = -mrSges(2,2) * g(3) - mrSges(2,3) * t687 + Ifges(2,5) * qJDD(1) - t707 * Ifges(2,6) - qJ(2) * t558 - t698 * t551 + t699 * t553;
t1 = [-m(1) * g(1) + t725; -m(1) * g(2) + t737; (-m(1) - m(2)) * g(3) + t558; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t737 + t705 * t550 - t702 * t552; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t725 + t702 * t550 + t705 * t552; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t687 - mrSges(2,2) * t688 + t698 * t553 + t699 * t551 + pkin(1) * (-qJDD(1) * t740 + t709) + qJ(2) * t724;];
tauB  = t1;
