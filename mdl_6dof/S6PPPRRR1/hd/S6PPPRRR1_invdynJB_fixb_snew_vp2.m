% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-05-04 19:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPPRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:33:25
% EndTime: 2019-05-04 19:33:56
% DurationCPUTime: 31.27s
% Computational Cost: add. (601543->248), mult. (1018761->336), div. (0->0), fcn. (904666->18), ass. (0->129)
t667 = sin(pkin(12));
t673 = cos(pkin(12));
t659 = -g(1) * t673 - g(2) * t667;
t666 = sin(pkin(13));
t672 = cos(pkin(13));
t658 = g(1) * t667 - g(2) * t673;
t664 = -g(3) + qJDD(1);
t670 = sin(pkin(6));
t676 = cos(pkin(6));
t691 = t658 * t676 + t664 * t670;
t634 = t672 * t659 + t691 * t666;
t665 = sin(pkin(14));
t671 = cos(pkin(14));
t633 = -t666 * t659 + t691 * t672;
t645 = -t658 * t670 + t664 * t676 + qJDD(2);
t669 = sin(pkin(7));
t675 = cos(pkin(7));
t693 = t633 * t675 + t645 * t669;
t629 = -t665 * t634 + t693 * t671;
t630 = t671 * t634 + t693 * t665;
t632 = -t633 * t669 + t645 * t675 + qJDD(3);
t682 = cos(qJ(4));
t674 = cos(pkin(8));
t679 = sin(qJ(4));
t705 = t674 * t679;
t668 = sin(pkin(8));
t706 = t668 * t679;
t623 = t629 * t705 + t682 * t630 + t632 * t706;
t684 = qJD(4) ^ 2;
t621 = -pkin(4) * t684 + qJDD(4) * pkin(10) + t623;
t625 = -t629 * t668 + t632 * t674;
t678 = sin(qJ(5));
t681 = cos(qJ(5));
t617 = t681 * t621 + t678 * t625;
t655 = (-pkin(5) * t681 - pkin(11) * t678) * qJD(4);
t683 = qJD(5) ^ 2;
t701 = qJD(4) * t681;
t615 = -pkin(5) * t683 + qJDD(5) * pkin(11) + t655 * t701 + t617;
t622 = -t679 * t630 + (t629 * t674 + t632 * t668) * t682;
t620 = -qJDD(4) * pkin(4) - t684 * pkin(10) - t622;
t700 = qJD(4) * qJD(5);
t698 = t681 * t700;
t656 = qJDD(4) * t678 + t698;
t699 = t678 * t700;
t657 = qJDD(4) * t681 - t699;
t618 = (-t656 - t698) * pkin(11) + (-t657 + t699) * pkin(5) + t620;
t677 = sin(qJ(6));
t680 = cos(qJ(6));
t611 = -t615 * t677 + t618 * t680;
t702 = qJD(4) * t678;
t652 = qJD(5) * t680 - t677 * t702;
t641 = qJD(6) * t652 + qJDD(5) * t677 + t656 * t680;
t653 = qJD(5) * t677 + t680 * t702;
t642 = -mrSges(7,1) * t652 + mrSges(7,2) * t653;
t662 = qJD(6) - t701;
t643 = -mrSges(7,2) * t662 + mrSges(7,3) * t652;
t651 = qJDD(6) - t657;
t609 = m(7) * t611 + mrSges(7,1) * t651 - mrSges(7,3) * t641 - t642 * t653 + t643 * t662;
t612 = t615 * t680 + t618 * t677;
t640 = -qJD(6) * t653 + qJDD(5) * t680 - t656 * t677;
t644 = mrSges(7,1) * t662 - mrSges(7,3) * t653;
t610 = m(7) * t612 - mrSges(7,2) * t651 + mrSges(7,3) * t640 + t642 * t652 - t644 * t662;
t603 = -t609 * t677 + t680 * t610;
t654 = (-mrSges(6,1) * t681 + mrSges(6,2) * t678) * qJD(4);
t660 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t702;
t601 = m(6) * t617 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t657 - qJD(5) * t660 + t654 * t701 + t603;
t704 = t681 * t625;
t614 = -qJDD(5) * pkin(5) - t683 * pkin(11) - t704 + (qJD(4) * t655 + t621) * t678;
t613 = -m(7) * t614 + t640 * mrSges(7,1) - mrSges(7,2) * t641 + t652 * t643 - t644 * t653;
t616 = -t621 * t678 + t704;
t661 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t701;
t607 = m(6) * t616 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t656 + qJD(5) * t661 - t654 * t702 + t613;
t696 = t681 * t601 - t607 * t678;
t592 = m(5) * t623 - mrSges(5,1) * t684 - qJDD(4) * mrSges(5,2) + t696;
t595 = t678 * t601 + t681 * t607;
t594 = m(5) * t625 + t595;
t602 = t609 * t680 + t610 * t677;
t687 = -m(6) * t620 + t657 * mrSges(6,1) - mrSges(6,2) * t656 - t660 * t702 + t661 * t701 - t602;
t598 = m(5) * t622 + qJDD(4) * mrSges(5,1) - mrSges(5,2) * t684 + t687;
t707 = t598 * t682;
t581 = t592 * t705 - t594 * t668 + t674 * t707;
t577 = m(4) * t629 + t581;
t586 = t682 * t592 - t598 * t679;
t585 = m(4) * t630 + t586;
t714 = t577 * t671 + t585 * t665;
t580 = t592 * t706 + t674 * t594 + t668 * t707;
t579 = m(4) * t632 + t580;
t566 = -t579 * t669 + t675 * t714;
t562 = m(3) * t633 + t566;
t571 = -t577 * t665 + t671 * t585;
t570 = m(3) * t634 + t571;
t713 = t562 * t672 + t570 * t666;
t635 = Ifges(7,5) * t653 + Ifges(7,6) * t652 + Ifges(7,3) * t662;
t637 = Ifges(7,1) * t653 + Ifges(7,4) * t652 + Ifges(7,5) * t662;
t604 = -mrSges(7,1) * t614 + mrSges(7,3) * t612 + Ifges(7,4) * t641 + Ifges(7,2) * t640 + Ifges(7,6) * t651 - t635 * t653 + t637 * t662;
t636 = Ifges(7,4) * t653 + Ifges(7,2) * t652 + Ifges(7,6) * t662;
t605 = mrSges(7,2) * t614 - mrSges(7,3) * t611 + Ifges(7,1) * t641 + Ifges(7,4) * t640 + Ifges(7,5) * t651 + t635 * t652 - t636 * t662;
t647 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t678 + Ifges(6,2) * t681) * qJD(4);
t648 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t678 + Ifges(6,4) * t681) * qJD(4);
t712 = mrSges(6,1) * t616 - mrSges(6,2) * t617 + Ifges(6,5) * t656 + Ifges(6,6) * t657 + Ifges(6,3) * qJDD(5) + pkin(5) * t613 + pkin(11) * t603 + t680 * t604 + t677 * t605 + (t647 * t678 - t648 * t681) * qJD(4);
t565 = t675 * t579 + t669 * t714;
t564 = m(3) * t645 + t565;
t552 = -t564 * t670 + t676 * t713;
t550 = m(2) * t658 + t552;
t559 = -t562 * t666 + t672 * t570;
t558 = m(2) * t659 + t559;
t703 = t673 * t550 + t667 * t558;
t551 = t676 * t564 + t670 * t713;
t697 = -t550 * t667 + t673 * t558;
t695 = m(2) * t664 + t551;
t646 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t678 + Ifges(6,6) * t681) * qJD(4);
t587 = mrSges(6,2) * t620 - mrSges(6,3) * t616 + Ifges(6,1) * t656 + Ifges(6,4) * t657 + Ifges(6,5) * qJDD(5) - pkin(11) * t602 - qJD(5) * t647 - t604 * t677 + t605 * t680 + t646 * t701;
t686 = mrSges(7,1) * t611 - mrSges(7,2) * t612 + Ifges(7,5) * t641 + Ifges(7,6) * t640 + Ifges(7,3) * t651 + t636 * t653 - t637 * t652;
t588 = -mrSges(6,1) * t620 + mrSges(6,3) * t617 + Ifges(6,4) * t656 + Ifges(6,2) * t657 + Ifges(6,6) * qJDD(5) - pkin(5) * t602 + qJD(5) * t648 - t646 * t702 - t686;
t573 = mrSges(5,2) * t625 - mrSges(5,3) * t622 + Ifges(5,5) * qJDD(4) - Ifges(5,6) * t684 - pkin(10) * t595 + t587 * t681 - t588 * t678;
t574 = -mrSges(5,1) * t625 + mrSges(5,3) * t623 + t684 * Ifges(5,5) + Ifges(5,6) * qJDD(4) - pkin(4) * t595 - t712;
t690 = pkin(9) * t586 + t573 * t679 + t574 * t682;
t572 = mrSges(5,1) * t622 - mrSges(5,2) * t623 + Ifges(5,3) * qJDD(4) + pkin(4) * t687 + pkin(10) * t696 + t678 * t587 + t681 * t588;
t553 = mrSges(4,1) * t629 - mrSges(4,2) * t630 + pkin(3) * t581 + t674 * t572 + t690 * t668;
t554 = -mrSges(4,1) * t632 + mrSges(4,3) * t630 - pkin(3) * t580 - t668 * t572 + t690 * t674;
t555 = mrSges(4,2) * t632 - mrSges(4,3) * t629 + t682 * t573 - t679 * t574 + (-t580 * t668 - t581 * t674) * pkin(9);
t688 = qJ(3) * t571 + t554 * t671 + t555 * t665;
t547 = -mrSges(3,1) * t645 + mrSges(3,3) * t634 - pkin(2) * t565 - t669 * t553 + t688 * t675;
t548 = mrSges(3,2) * t645 - mrSges(3,3) * t633 - t665 * t554 + t671 * t555 + (-t565 * t669 - t566 * t675) * qJ(3);
t689 = qJ(2) * t559 + t547 * t672 + t548 * t666;
t546 = mrSges(3,1) * t633 - mrSges(3,2) * t634 + pkin(2) * t566 + t675 * t553 + t688 * t669;
t545 = mrSges(2,2) * t664 - mrSges(2,3) * t658 - t666 * t547 + t672 * t548 + (-t551 * t670 - t552 * t676) * qJ(2);
t544 = -mrSges(2,1) * t664 + mrSges(2,3) * t659 - pkin(1) * t551 - t670 * t546 + t689 * t676;
t1 = [-m(1) * g(1) + t697; -m(1) * g(2) + t703; -m(1) * g(3) + t695; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t703 - t667 * t544 + t673 * t545; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t697 + t673 * t544 + t667 * t545; -mrSges(1,1) * g(2) + mrSges(2,1) * t658 + mrSges(1,2) * g(1) - mrSges(2,2) * t659 + pkin(1) * t552 + t676 * t546 + t689 * t670; t695; t564; t579; t572; t712; t686;];
tauJB  = t1;
