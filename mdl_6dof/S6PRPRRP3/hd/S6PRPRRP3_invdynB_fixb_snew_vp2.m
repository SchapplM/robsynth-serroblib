% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRP3
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:42:06
% EndTime: 2019-05-04 23:42:15
% DurationCPUTime: 8.95s
% Computational Cost: add. (133498->298), mult. (294289->367), div. (0->0), fcn. (218327->12), ass. (0->135)
t714 = Ifges(6,1) + Ifges(7,1);
t708 = Ifges(6,4) + Ifges(7,4);
t707 = Ifges(6,5) + Ifges(7,5);
t713 = Ifges(6,2) + Ifges(7,2);
t712 = Ifges(6,6) + Ifges(7,6);
t711 = Ifges(6,3) + Ifges(7,3);
t670 = qJD(2) ^ 2;
t657 = sin(pkin(11));
t660 = cos(pkin(11));
t664 = sin(qJ(4));
t667 = cos(qJ(4));
t677 = t657 * t664 - t660 * t667;
t639 = t677 * qJD(2);
t658 = sin(pkin(10));
t661 = cos(pkin(10));
t646 = g(1) * t658 - g(2) * t661;
t647 = -g(1) * t661 - g(2) * t658;
t656 = -g(3) + qJDD(1);
t659 = sin(pkin(6));
t662 = cos(pkin(6));
t665 = sin(qJ(2));
t668 = cos(qJ(2));
t620 = -t665 * t647 + (t646 * t662 + t656 * t659) * t668;
t678 = t657 * t667 + t660 * t664;
t640 = t678 * qJD(2);
t692 = qJD(4) * t640;
t628 = -qJDD(2) * t677 - t692;
t710 = pkin(3) * t660;
t709 = -mrSges(6,2) - mrSges(7,2);
t705 = mrSges(4,2) * t657;
t674 = qJDD(3) - t620;
t610 = -qJDD(2) * pkin(2) - qJ(3) * t670 + t674;
t654 = t657 ^ 2;
t701 = t662 * t665;
t702 = t659 * t665;
t621 = t646 * t701 + t647 * t668 + t656 * t702;
t616 = -pkin(2) * t670 + qJDD(2) * qJ(3) + t621;
t637 = -t646 * t659 + t656 * t662;
t691 = qJD(2) * qJD(3);
t695 = t637 * t660 - 0.2e1 * t657 * t691;
t586 = (-pkin(8) * qJDD(2) + t670 * t710 - t616) * t657 + t695;
t589 = t657 * t637 + (t616 + 0.2e1 * t691) * t660;
t690 = qJDD(2) * t660;
t655 = t660 ^ 2;
t703 = t655 * t670;
t587 = -pkin(3) * t703 + pkin(8) * t690 + t589;
t579 = t586 * t664 + t587 * t667;
t627 = pkin(4) * t639 - pkin(9) * t640;
t669 = qJD(4) ^ 2;
t577 = -pkin(4) * t669 + qJDD(4) * pkin(9) - t627 * t639 + t579;
t604 = (-pkin(2) - t710) * qJDD(2) + (-qJ(3) + (-t654 - t655) * pkin(8)) * t670 + t674;
t693 = qJD(4) * t639;
t629 = qJDD(2) * t678 - t693;
t582 = (-t629 + t693) * pkin(9) + (-t628 + t692) * pkin(4) + t604;
t663 = sin(qJ(5));
t666 = cos(qJ(5));
t572 = -t577 * t663 + t582 * t666;
t631 = qJD(4) * t666 - t640 * t663;
t603 = qJD(5) * t631 + qJDD(4) * t663 + t629 * t666;
t632 = qJD(4) * t663 + t640 * t666;
t606 = -mrSges(7,1) * t631 + mrSges(7,2) * t632;
t607 = -mrSges(6,1) * t631 + mrSges(6,2) * t632;
t638 = qJD(5) + t639;
t612 = -mrSges(6,2) * t638 + mrSges(6,3) * t631;
t626 = qJDD(5) - t628;
t569 = -0.2e1 * qJD(6) * t632 + (t631 * t638 - t603) * qJ(6) + (t631 * t632 + t626) * pkin(5) + t572;
t611 = -mrSges(7,2) * t638 + mrSges(7,3) * t631;
t689 = m(7) * t569 + mrSges(7,1) * t626 + t611 * t638;
t561 = m(6) * t572 + mrSges(6,1) * t626 + t612 * t638 + (-t606 - t607) * t632 + (-mrSges(6,3) - mrSges(7,3)) * t603 + t689;
t573 = t577 * t666 + t582 * t663;
t602 = -qJD(5) * t632 + qJDD(4) * t666 - t629 * t663;
t613 = pkin(5) * t638 - qJ(6) * t632;
t630 = t631 ^ 2;
t571 = -pkin(5) * t630 + qJ(6) * t602 + 0.2e1 * qJD(6) * t631 - t613 * t638 + t573;
t688 = m(7) * t571 + mrSges(7,3) * t602 + t606 * t631;
t614 = mrSges(7,1) * t638 - mrSges(7,3) * t632;
t696 = -mrSges(6,1) * t638 + mrSges(6,3) * t632 - t614;
t565 = m(6) * t573 + mrSges(6,3) * t602 + t607 * t631 + t626 * t709 + t638 * t696 + t688;
t559 = t561 * t666 + t565 * t663;
t635 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t639;
t636 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t640;
t673 = m(5) * t604 - mrSges(5,1) * t628 + t629 * mrSges(5,2) + t635 * t639 + t640 * t636 + t559;
t671 = -m(4) * t610 + mrSges(4,1) * t690 - t673 + (t654 * t670 + t703) * mrSges(4,3);
t554 = t671 - mrSges(3,2) * t670 + (mrSges(3,1) - t705) * qJDD(2) + m(3) * t620;
t704 = t554 * t668;
t624 = mrSges(5,1) * t639 + mrSges(5,2) * t640;
t684 = -t561 * t663 + t565 * t666;
t557 = m(5) * t579 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t628 - qJD(4) * t636 - t624 * t639 + t684;
t578 = t586 * t667 - t587 * t664;
t576 = -qJDD(4) * pkin(4) - pkin(9) * t669 + t627 * t640 - t578;
t574 = -pkin(5) * t602 - qJ(6) * t630 + t613 * t632 + qJDD(6) + t576;
t683 = m(7) * t574 - mrSges(7,1) * t602 - t611 * t631;
t672 = -m(6) * t576 + mrSges(6,1) * t602 + t603 * t709 + t612 * t631 + t632 * t696 - t683;
t566 = m(5) * t578 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t629 + qJD(4) * t635 - t624 * t640 + t672;
t550 = t557 * t664 + t566 * t667;
t588 = -t616 * t657 + t695;
t676 = mrSges(4,3) * qJDD(2) + t670 * (-mrSges(4,1) * t660 + t705);
t548 = m(4) * t588 - t657 * t676 + t550;
t685 = t557 * t667 - t664 * t566;
t549 = m(4) * t589 + t660 * t676 + t685;
t686 = -t548 * t657 + t549 * t660;
t540 = m(3) * t621 - mrSges(3,1) * t670 - qJDD(2) * mrSges(3,2) + t686;
t543 = t548 * t660 + t549 * t657;
t542 = m(3) * t637 + t543;
t530 = t540 * t701 - t542 * t659 + t662 * t704;
t528 = m(2) * t646 + t530;
t535 = t540 * t668 - t554 * t665;
t534 = m(2) * t647 + t535;
t700 = t528 * t661 + t534 * t658;
t699 = t631 * t712 + t632 * t707 + t638 * t711;
t698 = -t631 * t713 - t632 * t708 - t638 * t712;
t697 = t631 * t708 + t632 * t714 + t638 * t707;
t680 = Ifges(4,5) * t657 + Ifges(4,6) * t660;
t694 = t670 * t680;
t529 = t540 * t702 + t542 * t662 + t659 * t704;
t687 = -t528 * t658 + t534 * t661;
t682 = Ifges(4,1) * t657 + Ifges(4,4) * t660;
t681 = Ifges(4,4) * t657 + Ifges(4,2) * t660;
t551 = -mrSges(6,1) * t576 + mrSges(6,3) * t573 - mrSges(7,1) * t574 + mrSges(7,3) * t571 - pkin(5) * t683 + qJ(6) * t688 + (-qJ(6) * t614 + t697) * t638 + (-pkin(5) * t614 - t699) * t632 + (-mrSges(7,2) * qJ(6) + t712) * t626 + (-mrSges(7,2) * pkin(5) + t708) * t603 + t713 * t602;
t567 = -mrSges(7,3) * t603 - t606 * t632 + t689;
t558 = mrSges(6,2) * t576 + mrSges(7,2) * t574 - mrSges(6,3) * t572 - mrSges(7,3) * t569 - qJ(6) * t567 + t602 * t708 + t603 * t714 + t626 * t707 + t631 * t699 + t638 * t698;
t617 = Ifges(5,5) * t640 - Ifges(5,6) * t639 + Ifges(5,3) * qJD(4);
t618 = Ifges(5,4) * t640 - Ifges(5,2) * t639 + Ifges(5,6) * qJD(4);
t536 = mrSges(5,2) * t604 - mrSges(5,3) * t578 + Ifges(5,1) * t629 + Ifges(5,4) * t628 + Ifges(5,5) * qJDD(4) - pkin(9) * t559 - qJD(4) * t618 - t551 * t663 + t558 * t666 - t617 * t639;
t619 = Ifges(5,1) * t640 - Ifges(5,4) * t639 + Ifges(5,5) * qJD(4);
t544 = -mrSges(5,1) * t604 - mrSges(6,1) * t572 - mrSges(7,1) * t569 + mrSges(6,2) * t573 + mrSges(7,2) * t571 + mrSges(5,3) * t579 + Ifges(5,4) * t629 + Ifges(5,2) * t628 + Ifges(5,6) * qJDD(4) - pkin(4) * t559 - pkin(5) * t567 + qJD(4) * t619 - t640 * t617 + t698 * t632 + t697 * t631 - t711 * t626 - t707 * t603 - t712 * t602;
t526 = -mrSges(4,1) * t610 + mrSges(4,3) * t589 - pkin(3) * t673 + pkin(8) * t685 + qJDD(2) * t681 + t664 * t536 + t667 * t544 - t657 * t694;
t531 = mrSges(4,2) * t610 - mrSges(4,3) * t588 - pkin(8) * t550 + qJDD(2) * t682 + t536 * t667 - t544 * t664 + t660 * t694;
t524 = mrSges(3,2) * t637 - mrSges(3,3) * t620 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t670 - qJ(3) * t543 - t526 * t657 + t531 * t660;
t525 = -pkin(4) * t672 - pkin(9) * t684 - t663 * t558 - t666 * t551 - t639 * t619 - mrSges(5,1) * t578 + mrSges(5,2) * t579 - Ifges(5,5) * t629 - Ifges(5,6) * t628 - Ifges(5,3) * qJDD(4) - t640 * t618 - mrSges(4,1) * t588 + mrSges(4,2) * t589 - pkin(3) * t550 - pkin(2) * t543 + mrSges(3,3) * t621 - mrSges(3,1) * t637 + (Ifges(3,6) - t680) * qJDD(2) + (-t657 * t681 + t660 * t682 + Ifges(3,5)) * t670;
t675 = pkin(7) * t535 + t524 * t665 + t525 * t668;
t523 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t620 - mrSges(3,2) * t621 + t657 * t531 + t660 * t526 + pkin(2) * (-qJDD(2) * t705 + t671) + qJ(3) * t686;
t522 = mrSges(2,2) * t656 - mrSges(2,3) * t646 + t524 * t668 - t525 * t665 + (-t529 * t659 - t530 * t662) * pkin(7);
t521 = -mrSges(2,1) * t656 + mrSges(2,3) * t647 - pkin(1) * t529 - t523 * t659 + t662 * t675;
t1 = [-m(1) * g(1) + t687; -m(1) * g(2) + t700; -m(1) * g(3) + m(2) * t656 + t529; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t700 - t521 * t658 + t522 * t661; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t687 + t661 * t521 + t658 * t522; -mrSges(1,1) * g(2) + mrSges(2,1) * t646 + mrSges(1,2) * g(1) - mrSges(2,2) * t647 + pkin(1) * t530 + t523 * t662 + t659 * t675;];
tauB  = t1;
