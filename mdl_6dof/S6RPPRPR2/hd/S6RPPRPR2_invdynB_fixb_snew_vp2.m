% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-05-05 14:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:02:05
% EndTime: 2019-05-05 14:02:11
% DurationCPUTime: 5.48s
% Computational Cost: add. (59650->303), mult. (134011->361), div. (0->0), fcn. (88714->10), ass. (0->130)
t699 = Ifges(5,1) + Ifges(6,2);
t693 = Ifges(5,4) + Ifges(6,6);
t692 = Ifges(5,5) - Ifges(6,4);
t698 = Ifges(5,2) + Ifges(6,3);
t691 = Ifges(5,6) - Ifges(6,5);
t697 = -Ifges(5,3) - Ifges(6,1);
t658 = qJD(1) ^ 2;
t696 = -2 * qJD(5);
t695 = cos(qJ(4));
t650 = cos(pkin(10));
t694 = pkin(3) * t650;
t648 = sin(pkin(10));
t690 = mrSges(4,2) * t648;
t645 = t650 ^ 2;
t689 = t645 * t658;
t654 = sin(qJ(1));
t656 = cos(qJ(1));
t631 = t654 * g(1) - g(2) * t656;
t629 = qJDD(1) * pkin(1) + t631;
t632 = -g(1) * t656 - g(2) * t654;
t630 = -pkin(1) * t658 + t632;
t649 = sin(pkin(9));
t651 = cos(pkin(9));
t609 = t649 * t629 + t651 * t630;
t594 = -pkin(2) * t658 + qJDD(1) * qJ(3) + t609;
t647 = -g(3) + qJDD(2);
t679 = qJD(1) * qJD(3);
t683 = t650 * t647 - 0.2e1 * t648 * t679;
t574 = (-pkin(7) * qJDD(1) + t658 * t694 - t594) * t648 + t683;
t580 = t648 * t647 + (t594 + 0.2e1 * t679) * t650;
t678 = qJDD(1) * t650;
t577 = -pkin(3) * t689 + pkin(7) * t678 + t580;
t653 = sin(qJ(4));
t564 = t574 * t695 - t653 * t577;
t676 = t650 * t695;
t682 = qJD(1) * t648;
t622 = -qJD(1) * t676 + t653 * t682;
t666 = t648 * t695 + t650 * t653;
t623 = t666 * qJD(1);
t599 = mrSges(5,1) * t622 + mrSges(5,2) * t623;
t680 = t622 * qJD(4);
t607 = qJDD(1) * t666 - t680;
t613 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t622;
t615 = mrSges(6,1) * t622 - qJD(4) * mrSges(6,3);
t598 = pkin(4) * t622 - qJ(5) * t623;
t657 = qJD(4) ^ 2;
t561 = -qJDD(4) * pkin(4) - t657 * qJ(5) + t623 * t598 + qJDD(5) - t564;
t556 = (t622 * t623 - qJDD(4)) * pkin(8) + (t607 + t680) * pkin(5) + t561;
t681 = qJD(4) * t623;
t606 = t681 + (t648 * t653 - t676) * qJDD(1);
t617 = pkin(5) * t623 - qJD(4) * pkin(8);
t621 = t622 ^ 2;
t644 = t648 ^ 2;
t608 = t651 * t629 - t649 * t630;
t668 = qJDD(3) - t608;
t578 = (-pkin(2) - t694) * qJDD(1) + (-qJ(3) + (-t644 - t645) * pkin(7)) * t658 + t668;
t659 = pkin(4) * t681 + t623 * t696 + (-t607 + t680) * qJ(5) + t578;
t559 = -t621 * pkin(5) - t623 * t617 + (pkin(4) + pkin(8)) * t606 + t659;
t652 = sin(qJ(6));
t655 = cos(qJ(6));
t554 = t556 * t655 - t559 * t652;
t610 = -qJD(4) * t652 + t622 * t655;
t576 = qJD(6) * t610 + qJDD(4) * t655 + t606 * t652;
t611 = qJD(4) * t655 + t622 * t652;
t581 = -mrSges(7,1) * t610 + mrSges(7,2) * t611;
t620 = qJD(6) + t623;
t584 = -mrSges(7,2) * t620 + mrSges(7,3) * t610;
t605 = qJDD(6) + t607;
t552 = m(7) * t554 + mrSges(7,1) * t605 - mrSges(7,3) * t576 - t581 * t611 + t584 * t620;
t555 = t556 * t652 + t559 * t655;
t575 = -qJD(6) * t611 - qJDD(4) * t652 + t606 * t655;
t585 = mrSges(7,1) * t620 - mrSges(7,3) * t611;
t553 = m(7) * t555 - mrSges(7,2) * t605 + mrSges(7,3) * t575 + t581 * t610 - t585 * t620;
t544 = t655 * t552 + t652 * t553;
t600 = -mrSges(6,2) * t622 - mrSges(6,3) * t623;
t664 = -m(6) * t561 - t607 * mrSges(6,1) - t623 * t600 - t544;
t542 = m(5) * t564 - t607 * mrSges(5,3) - t623 * t599 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t613 - t615) * qJD(4) + t664;
t565 = t653 * t574 + t695 * t577;
t614 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t623;
t663 = -t657 * pkin(4) + qJDD(4) * qJ(5) - t622 * t598 + t565;
t560 = qJD(4) * t696 - t663;
t616 = mrSges(6,1) * t623 + qJD(4) * mrSges(6,2);
t558 = -t606 * pkin(5) - t621 * pkin(8) + ((2 * qJD(5)) + t617) * qJD(4) + t663;
t665 = -m(7) * t558 + t575 * mrSges(7,1) - t576 * mrSges(7,2) + t610 * t584 - t611 * t585;
t662 = -m(6) * t560 + qJDD(4) * mrSges(6,3) + qJD(4) * t616 - t665;
t549 = m(5) * t565 - qJDD(4) * mrSges(5,2) - qJD(4) * t614 + (-t599 - t600) * t622 + (-mrSges(5,3) - mrSges(6,1)) * t606 + t662;
t538 = t695 * t542 + t653 * t549;
t579 = -t594 * t648 + t683;
t667 = mrSges(4,3) * qJDD(1) + t658 * (-mrSges(4,1) * t650 + t690);
t536 = m(4) * t579 - t648 * t667 + t538;
t672 = -t653 * t542 + t695 * t549;
t537 = m(4) * t580 + t650 * t667 + t672;
t673 = -t536 * t648 + t650 * t537;
t530 = m(3) * t609 - mrSges(3,1) * t658 - qJDD(1) * mrSges(3,2) + t673;
t587 = -qJDD(1) * pkin(2) - t658 * qJ(3) + t668;
t563 = t606 * pkin(4) + t659;
t687 = -t652 * t552 + t655 * t553;
t543 = m(6) * t563 - t606 * mrSges(6,2) - t607 * mrSges(6,3) - t622 * t615 - t623 * t616 + t687;
t661 = m(5) * t578 + t606 * mrSges(5,1) + t607 * mrSges(5,2) + t622 * t613 + t623 * t614 + t543;
t660 = -m(4) * t587 + mrSges(4,1) * t678 - t661 + (t644 * t658 + t689) * mrSges(4,3);
t540 = (mrSges(3,1) - t690) * qJDD(1) + t660 - t658 * mrSges(3,2) + m(3) * t608;
t526 = t649 * t530 + t651 * t540;
t524 = m(2) * t631 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t658 + t526;
t674 = t651 * t530 - t540 * t649;
t525 = m(2) * t632 - mrSges(2,1) * t658 - qJDD(1) * mrSges(2,2) + t674;
t688 = t656 * t524 + t654 * t525;
t531 = t650 * t536 + t648 * t537;
t686 = t697 * qJD(4) + t691 * t622 - t692 * t623;
t685 = -t691 * qJD(4) + t698 * t622 - t693 * t623;
t684 = t692 * qJD(4) - t693 * t622 + t699 * t623;
t677 = m(3) * t647 + t531;
t675 = -t524 * t654 + t656 * t525;
t671 = Ifges(4,1) * t648 + Ifges(4,4) * t650;
t670 = Ifges(4,4) * t648 + Ifges(4,2) * t650;
t669 = Ifges(4,5) * t648 + Ifges(4,6) * t650;
t628 = t669 * qJD(1);
t569 = Ifges(7,1) * t611 + Ifges(7,4) * t610 + Ifges(7,5) * t620;
t568 = Ifges(7,4) * t611 + Ifges(7,2) * t610 + Ifges(7,6) * t620;
t567 = Ifges(7,5) * t611 + Ifges(7,6) * t610 + Ifges(7,3) * t620;
t546 = mrSges(7,2) * t558 - mrSges(7,3) * t554 + Ifges(7,1) * t576 + Ifges(7,4) * t575 + Ifges(7,5) * t605 + t567 * t610 - t568 * t620;
t545 = -mrSges(7,1) * t558 + mrSges(7,3) * t555 + Ifges(7,4) * t576 + Ifges(7,2) * t575 + Ifges(7,6) * t605 - t567 * t611 + t569 * t620;
t532 = mrSges(6,1) * t561 + mrSges(7,1) * t554 + mrSges(5,2) * t578 - mrSges(7,2) * t555 - mrSges(5,3) * t564 - mrSges(6,3) * t563 + Ifges(7,5) * t576 + Ifges(7,6) * t575 + Ifges(7,3) * t605 + pkin(5) * t544 - qJ(5) * t543 + t611 * t568 - t610 * t569 + t686 * t622 + t699 * t607 - t693 * t606 + t692 * qJDD(4) + t685 * qJD(4);
t527 = -mrSges(5,1) * t578 - mrSges(6,1) * t560 + mrSges(6,2) * t563 + mrSges(5,3) * t565 - pkin(4) * t543 - pkin(5) * t665 - pkin(8) * t687 + t684 * qJD(4) + t691 * qJDD(4) - t655 * t545 - t652 * t546 - t698 * t606 + t693 * t607 + t686 * t623;
t520 = t650 * qJD(1) * t628 + mrSges(4,2) * t587 - mrSges(4,3) * t579 - pkin(7) * t538 + qJDD(1) * t671 - t653 * t527 + t532 * t695;
t519 = t652 * t545 - pkin(4) * (-qJD(4) * t615 + t664) - t655 * t546 - qJ(5) * t662 - mrSges(3,1) * t647 + mrSges(3,3) * t609 - mrSges(4,1) * t579 + mrSges(4,2) * t580 + mrSges(5,2) * t565 - mrSges(6,2) * t561 - mrSges(5,1) * t564 + mrSges(6,3) * t560 + pkin(8) * t544 - pkin(3) * t538 - pkin(2) * t531 + t685 * t623 + (qJ(5) * t600 - t684) * t622 - t692 * t607 + (mrSges(6,1) * qJ(5) + t691) * t606 + (mrSges(6,2) * pkin(4) + t697) * qJDD(4) + (Ifges(3,6) - t669) * qJDD(1) + (-t648 * t670 + t650 * t671 + Ifges(3,5)) * t658;
t518 = -mrSges(4,1) * t587 + mrSges(4,3) * t580 - pkin(3) * t661 + pkin(7) * t672 + qJDD(1) * t670 + t527 * t695 + t653 * t532 - t628 * t682;
t517 = mrSges(3,2) * t647 - mrSges(3,3) * t608 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t658 - qJ(3) * t531 - t518 * t648 + t520 * t650;
t516 = -mrSges(2,2) * g(3) - mrSges(2,3) * t631 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t658 - qJ(2) * t526 + t517 * t651 - t519 * t649;
t515 = mrSges(2,1) * g(3) + mrSges(2,3) * t632 + t658 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t677 + qJ(2) * t674 + t649 * t517 + t651 * t519;
t1 = [-m(1) * g(1) + t675; -m(1) * g(2) + t688; (-m(1) - m(2)) * g(3) + t677; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t688 - t654 * t515 + t656 * t516; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t675 + t656 * t515 + t654 * t516; pkin(1) * t526 + mrSges(2,1) * t631 - mrSges(2,2) * t632 + qJ(3) * t673 + t648 * t520 + t650 * t518 + pkin(2) * t660 + mrSges(3,1) * t608 - mrSges(3,2) * t609 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t690 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
