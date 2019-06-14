% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-05-05 13:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:54:59
% EndTime: 2019-05-05 13:55:10
% DurationCPUTime: 10.83s
% Computational Cost: add. (164832->322), mult. (381014->404), div. (0->0), fcn. (270247->12), ass. (0->136)
t658 = qJD(1) ^ 2;
t687 = cos(qJ(4));
t650 = cos(pkin(10));
t686 = pkin(3) * t650;
t647 = sin(pkin(10));
t685 = mrSges(4,2) * t647;
t644 = t650 ^ 2;
t684 = t644 * t658;
t654 = sin(qJ(1));
t656 = cos(qJ(1));
t630 = t654 * g(1) - g(2) * t656;
t628 = qJDD(1) * pkin(1) + t630;
t631 = -g(1) * t656 - g(2) * t654;
t629 = -pkin(1) * t658 + t631;
t648 = sin(pkin(9));
t651 = cos(pkin(9));
t612 = t648 * t628 + t651 * t629;
t603 = -pkin(2) * t658 + qJDD(1) * qJ(3) + t612;
t645 = -g(3) + qJDD(2);
t678 = qJD(1) * qJD(3);
t682 = t650 * t645 - 0.2e1 * t647 * t678;
t583 = (-pkin(7) * qJDD(1) + t658 * t686 - t603) * t647 + t682;
t589 = t647 * t645 + (t603 + 0.2e1 * t678) * t650;
t677 = qJDD(1) * t650;
t584 = -pkin(3) * t684 + pkin(7) * t677 + t589;
t653 = sin(qJ(4));
t568 = t653 * t583 + t687 * t584;
t675 = t650 * t687;
t681 = qJD(1) * t647;
t621 = -qJD(1) * t675 + t653 * t681;
t663 = t647 * t687 + t650 * t653;
t622 = t663 * qJD(1);
t606 = mrSges(5,1) * t621 + mrSges(5,2) * t622;
t680 = qJD(4) * t622;
t609 = t680 + (t647 * t653 - t675) * qJDD(1);
t619 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t622;
t605 = pkin(4) * t621 - qJ(5) * t622;
t657 = qJD(4) ^ 2;
t560 = -pkin(4) * t657 + qJDD(4) * qJ(5) - t605 * t621 + t568;
t643 = t647 ^ 2;
t611 = t651 * t628 - t648 * t629;
t665 = qJDD(3) - t611;
t585 = (-pkin(2) - t686) * qJDD(1) + (-qJ(3) + (-t643 - t644) * pkin(7)) * t658 + t665;
t679 = t621 * qJD(4);
t610 = qJDD(1) * t663 - t679;
t563 = (-t610 + t679) * qJ(5) + (t609 + t680) * pkin(4) + t585;
t646 = sin(pkin(11));
t649 = cos(pkin(11));
t617 = qJD(4) * t646 + t622 * t649;
t555 = -0.2e1 * qJD(5) * t617 - t646 * t560 + t649 * t563;
t597 = qJDD(4) * t646 + t610 * t649;
t616 = qJD(4) * t649 - t622 * t646;
t553 = (t616 * t621 - t597) * pkin(8) + (t616 * t617 + t609) * pkin(5) + t555;
t556 = 0.2e1 * qJD(5) * t616 + t649 * t560 + t646 * t563;
t595 = pkin(5) * t621 - pkin(8) * t617;
t596 = qJDD(4) * t649 - t610 * t646;
t615 = t616 ^ 2;
t554 = -pkin(5) * t615 + pkin(8) * t596 - t595 * t621 + t556;
t652 = sin(qJ(6));
t655 = cos(qJ(6));
t551 = t553 * t655 - t554 * t652;
t586 = t616 * t655 - t617 * t652;
t566 = qJD(6) * t586 + t596 * t652 + t597 * t655;
t587 = t616 * t652 + t617 * t655;
t573 = -mrSges(7,1) * t586 + mrSges(7,2) * t587;
t620 = qJD(6) + t621;
t574 = -mrSges(7,2) * t620 + mrSges(7,3) * t586;
t608 = qJDD(6) + t609;
t549 = m(7) * t551 + mrSges(7,1) * t608 - mrSges(7,3) * t566 - t573 * t587 + t574 * t620;
t552 = t553 * t652 + t554 * t655;
t565 = -qJD(6) * t587 + t596 * t655 - t597 * t652;
t575 = mrSges(7,1) * t620 - mrSges(7,3) * t587;
t550 = m(7) * t552 - mrSges(7,2) * t608 + mrSges(7,3) * t565 + t573 * t586 - t575 * t620;
t541 = t655 * t549 + t652 * t550;
t590 = -mrSges(6,1) * t616 + mrSges(6,2) * t617;
t593 = -mrSges(6,2) * t621 + mrSges(6,3) * t616;
t539 = m(6) * t555 + mrSges(6,1) * t609 - mrSges(6,3) * t597 - t590 * t617 + t593 * t621 + t541;
t594 = mrSges(6,1) * t621 - mrSges(6,3) * t617;
t669 = -t549 * t652 + t655 * t550;
t540 = m(6) * t556 - mrSges(6,2) * t609 + mrSges(6,3) * t596 + t590 * t616 - t594 * t621 + t669;
t670 = -t539 * t646 + t649 * t540;
t534 = m(5) * t568 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t609 - qJD(4) * t619 - t606 * t621 + t670;
t567 = t583 * t687 - t653 * t584;
t618 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t621;
t559 = -qJDD(4) * pkin(4) - t657 * qJ(5) + t622 * t605 + qJDD(5) - t567;
t557 = -t596 * pkin(5) - t615 * pkin(8) + t617 * t595 + t559;
t662 = m(7) * t557 - t565 * mrSges(7,1) + mrSges(7,2) * t566 - t586 * t574 + t575 * t587;
t659 = -m(6) * t559 + t596 * mrSges(6,1) - mrSges(6,2) * t597 + t616 * t593 - t594 * t617 - t662;
t545 = m(5) * t567 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t610 + qJD(4) * t618 - t606 * t622 + t659;
t527 = t653 * t534 + t687 * t545;
t588 = -t603 * t647 + t682;
t664 = mrSges(4,3) * qJDD(1) + t658 * (-mrSges(4,1) * t650 + t685);
t525 = m(4) * t588 - t647 * t664 + t527;
t671 = t687 * t534 - t653 * t545;
t526 = m(4) * t589 + t650 * t664 + t671;
t672 = -t525 * t647 + t650 * t526;
t519 = m(3) * t612 - mrSges(3,1) * t658 - qJDD(1) * mrSges(3,2) + t672;
t599 = -qJDD(1) * pkin(2) - t658 * qJ(3) + t665;
t535 = t649 * t539 + t646 * t540;
t661 = m(5) * t585 + t609 * mrSges(5,1) + t610 * mrSges(5,2) + t621 * t618 + t622 * t619 + t535;
t660 = -m(4) * t599 + mrSges(4,1) * t677 - t661 + (t643 * t658 + t684) * mrSges(4,3);
t531 = (mrSges(3,1) - t685) * qJDD(1) + t660 - t658 * mrSges(3,2) + m(3) * t611;
t515 = t648 * t519 + t651 * t531;
t513 = m(2) * t630 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t658 + t515;
t673 = t651 * t519 - t531 * t648;
t514 = m(2) * t631 - mrSges(2,1) * t658 - qJDD(1) * mrSges(2,2) + t673;
t683 = t656 * t513 + t654 * t514;
t520 = t650 * t525 + t647 * t526;
t676 = m(3) * t645 + t520;
t674 = -t513 * t654 + t656 * t514;
t668 = Ifges(4,1) * t647 + Ifges(4,4) * t650;
t667 = Ifges(4,4) * t647 + Ifges(4,2) * t650;
t666 = Ifges(4,5) * t647 + Ifges(4,6) * t650;
t627 = t666 * qJD(1);
t602 = Ifges(5,1) * t622 - Ifges(5,4) * t621 + Ifges(5,5) * qJD(4);
t601 = Ifges(5,4) * t622 - Ifges(5,2) * t621 + Ifges(5,6) * qJD(4);
t600 = Ifges(5,5) * t622 - Ifges(5,6) * t621 + Ifges(5,3) * qJD(4);
t579 = Ifges(6,1) * t617 + Ifges(6,4) * t616 + Ifges(6,5) * t621;
t578 = Ifges(6,4) * t617 + Ifges(6,2) * t616 + Ifges(6,6) * t621;
t577 = Ifges(6,5) * t617 + Ifges(6,6) * t616 + Ifges(6,3) * t621;
t571 = Ifges(7,1) * t587 + Ifges(7,4) * t586 + Ifges(7,5) * t620;
t570 = Ifges(7,4) * t587 + Ifges(7,2) * t586 + Ifges(7,6) * t620;
t569 = Ifges(7,5) * t587 + Ifges(7,6) * t586 + Ifges(7,3) * t620;
t543 = mrSges(7,2) * t557 - mrSges(7,3) * t551 + Ifges(7,1) * t566 + Ifges(7,4) * t565 + Ifges(7,5) * t608 + t569 * t586 - t570 * t620;
t542 = -mrSges(7,1) * t557 + mrSges(7,3) * t552 + Ifges(7,4) * t566 + Ifges(7,2) * t565 + Ifges(7,6) * t608 - t569 * t587 + t571 * t620;
t529 = mrSges(6,2) * t559 - mrSges(6,3) * t555 + Ifges(6,1) * t597 + Ifges(6,4) * t596 + Ifges(6,5) * t609 - pkin(8) * t541 - t542 * t652 + t543 * t655 + t577 * t616 - t578 * t621;
t528 = -mrSges(6,1) * t559 + mrSges(6,3) * t556 + Ifges(6,4) * t597 + Ifges(6,2) * t596 + Ifges(6,6) * t609 - pkin(5) * t662 + pkin(8) * t669 + t655 * t542 + t652 * t543 - t617 * t577 + t621 * t579;
t521 = Ifges(5,4) * t610 + Ifges(5,6) * qJDD(4) - t622 * t600 + qJD(4) * t602 - mrSges(5,1) * t585 + mrSges(5,3) * t568 - Ifges(6,5) * t597 - Ifges(6,6) * t596 - t617 * t578 + t616 * t579 - mrSges(6,1) * t555 + mrSges(6,2) * t556 - Ifges(7,5) * t566 - Ifges(7,6) * t565 - Ifges(7,3) * t608 - t587 * t570 + t586 * t571 - mrSges(7,1) * t551 + mrSges(7,2) * t552 - pkin(5) * t541 - pkin(4) * t535 + (-Ifges(5,2) - Ifges(6,3)) * t609;
t516 = mrSges(5,2) * t585 - mrSges(5,3) * t567 + Ifges(5,1) * t610 - Ifges(5,4) * t609 + Ifges(5,5) * qJDD(4) - qJ(5) * t535 - qJD(4) * t601 - t528 * t646 + t529 * t649 - t600 * t621;
t509 = t650 * qJD(1) * t627 + mrSges(4,2) * t599 - mrSges(4,3) * t588 - pkin(7) * t527 + qJDD(1) * t668 + t516 * t687 - t653 * t521;
t508 = -pkin(2) * t520 - mrSges(3,1) * t645 + mrSges(3,3) * t612 - pkin(3) * t527 - mrSges(4,1) * t588 + mrSges(4,2) * t589 - t646 * t529 - t649 * t528 - pkin(4) * t659 - qJ(5) * t670 - Ifges(5,5) * t610 + Ifges(5,6) * t609 - Ifges(5,3) * qJDD(4) - t622 * t601 - t621 * t602 - mrSges(5,1) * t567 + mrSges(5,2) * t568 + (Ifges(3,6) - t666) * qJDD(1) + (-t647 * t667 + t650 * t668 + Ifges(3,5)) * t658;
t507 = -mrSges(4,1) * t599 + mrSges(4,3) * t589 - pkin(3) * t661 + pkin(7) * t671 + qJDD(1) * t667 + t653 * t516 + t521 * t687 - t627 * t681;
t506 = mrSges(3,2) * t645 - mrSges(3,3) * t611 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t658 - qJ(3) * t520 - t507 * t647 + t509 * t650;
t505 = -mrSges(2,2) * g(3) - mrSges(2,3) * t630 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t658 - qJ(2) * t515 + t506 * t651 - t508 * t648;
t504 = mrSges(2,1) * g(3) + mrSges(2,3) * t631 + t658 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t676 + qJ(2) * t673 + t648 * t506 + t651 * t508;
t1 = [-m(1) * g(1) + t674; -m(1) * g(2) + t683; (-m(1) - m(2)) * g(3) + t676; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t683 - t654 * t504 + t656 * t505; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t674 + t656 * t504 + t654 * t505; pkin(1) * t515 + mrSges(2,1) * t630 - mrSges(2,2) * t631 + t647 * t509 + t650 * t507 + pkin(2) * t660 + qJ(3) * t672 + mrSges(3,1) * t611 - mrSges(3,2) * t612 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t685 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
