% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-05-05 21:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:41:46
% EndTime: 2019-05-05 21:41:54
% DurationCPUTime: 4.78s
% Computational Cost: add. (52870->314), mult. (105285->371), div. (0->0), fcn. (66216->8), ass. (0->123)
t669 = Ifges(6,1) + Ifges(7,1);
t660 = Ifges(6,4) - Ifges(7,5);
t658 = Ifges(6,5) + Ifges(7,4);
t668 = Ifges(6,2) + Ifges(7,3);
t667 = -Ifges(7,2) - Ifges(6,3);
t657 = Ifges(6,6) - Ifges(7,6);
t625 = sin(qJ(1));
t628 = cos(qJ(1));
t613 = -t628 * g(1) - t625 * g(2);
t666 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t613;
t665 = -2 * qJD(5);
t664 = -pkin(1) - pkin(7);
t663 = mrSges(2,1) - mrSges(3,2);
t662 = -mrSges(6,3) - mrSges(7,2);
t661 = -Ifges(3,4) + Ifges(2,5);
t659 = (Ifges(3,5) - Ifges(2,6));
t656 = cos(pkin(9));
t612 = t625 * g(1) - g(2) * t628;
t630 = qJD(1) ^ 2;
t636 = -t630 * qJ(2) + qJDD(2) - t612;
t589 = qJDD(1) * t664 + t636;
t624 = sin(qJ(3));
t627 = cos(qJ(3));
t582 = -g(3) * t627 + t589 * t624;
t606 = (mrSges(4,1) * t624 + mrSges(4,2) * t627) * qJD(1);
t648 = qJD(1) * qJD(3);
t644 = t627 * t648;
t608 = -qJDD(1) * t624 - t644;
t650 = qJD(1) * t627;
t611 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t650;
t588 = t630 * t664 - t666;
t645 = t624 * t648;
t609 = qJDD(1) * t627 - t645;
t558 = (-t609 + t645) * pkin(8) + (-t608 + t644) * pkin(3) + t588;
t607 = (pkin(3) * t624 - pkin(8) * t627) * qJD(1);
t629 = qJD(3) ^ 2;
t649 = t624 * qJD(1);
t563 = -pkin(3) * t629 + qJDD(3) * pkin(8) - t607 * t649 + t582;
t623 = sin(qJ(4));
t626 = cos(qJ(4));
t537 = t558 * t626 - t623 * t563;
t604 = qJD(3) * t626 - t623 * t650;
t577 = qJD(4) * t604 + qJDD(3) * t623 + t609 * t626;
t603 = qJDD(4) - t608;
t605 = qJD(3) * t623 + t626 * t650;
t615 = qJD(4) + t649;
t533 = (t604 * t615 - t577) * qJ(5) + (t604 * t605 + t603) * pkin(4) + t537;
t538 = t558 * t623 + t563 * t626;
t576 = -qJD(4) * t605 + qJDD(3) * t626 - t609 * t623;
t584 = pkin(4) * t615 - qJ(5) * t605;
t602 = t604 ^ 2;
t535 = -pkin(4) * t602 + qJ(5) * t576 - t584 * t615 + t538;
t622 = sin(pkin(9));
t578 = -t604 * t656 + t605 * t622;
t529 = t533 * t622 + t535 * t656 + t578 * t665;
t548 = -t576 * t656 + t577 * t622;
t579 = t604 * t622 + t605 * t656;
t566 = mrSges(6,1) * t615 - mrSges(6,3) * t579;
t553 = pkin(5) * t578 - qJ(6) * t579;
t614 = t615 ^ 2;
t526 = -pkin(5) * t614 + qJ(6) * t603 + 0.2e1 * qJD(6) * t615 - t553 * t578 + t529;
t567 = -mrSges(7,1) * t615 + mrSges(7,2) * t579;
t646 = m(7) * t526 + mrSges(7,3) * t603 + t567 * t615;
t554 = mrSges(7,1) * t578 - mrSges(7,3) * t579;
t651 = -mrSges(6,1) * t578 - mrSges(6,2) * t579 - t554;
t520 = m(6) * t529 - mrSges(6,2) * t603 + t548 * t662 - t566 * t615 + t578 * t651 + t646;
t637 = t533 * t656 - t622 * t535;
t528 = t579 * t665 + t637;
t549 = t576 * t622 + t577 * t656;
t564 = -mrSges(6,2) * t615 - mrSges(6,3) * t578;
t527 = -t603 * pkin(5) - t614 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t553) * t579 - t637;
t565 = -mrSges(7,2) * t578 + mrSges(7,3) * t615;
t639 = -m(7) * t527 + mrSges(7,1) * t603 + t565 * t615;
t522 = m(6) * t528 + mrSges(6,1) * t603 + t549 * t662 + t564 * t615 + t579 * t651 + t639;
t515 = t520 * t622 + t522 * t656;
t580 = -mrSges(5,1) * t604 + mrSges(5,2) * t605;
t583 = -mrSges(5,2) * t615 + mrSges(5,3) * t604;
t513 = m(5) * t537 + mrSges(5,1) * t603 - mrSges(5,3) * t577 - t580 * t605 + t583 * t615 + t515;
t585 = mrSges(5,1) * t615 - mrSges(5,3) * t605;
t640 = t520 * t656 - t522 * t622;
t514 = m(5) * t538 - mrSges(5,2) * t603 + mrSges(5,3) * t576 + t580 * t604 - t585 * t615 + t640;
t641 = -t513 * t623 + t514 * t626;
t508 = m(4) * t582 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t608 - qJD(3) * t611 - t606 * t649 + t641;
t581 = t624 * g(3) + t627 * t589;
t610 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t649;
t562 = -qJDD(3) * pkin(3) - t629 * pkin(8) + t607 * t650 - t581;
t536 = -t576 * pkin(4) - t602 * qJ(5) + t584 * t605 + qJDD(5) + t562;
t531 = -0.2e1 * qJD(6) * t579 + (t578 * t615 - t549) * qJ(6) + (t579 * t615 + t548) * pkin(5) + t536;
t524 = m(7) * t531 + mrSges(7,1) * t548 - mrSges(7,3) * t549 + t565 * t578 - t567 * t579;
t632 = m(6) * t536 + mrSges(6,1) * t548 + mrSges(6,2) * t549 + t564 * t578 + t566 * t579 + t524;
t631 = -m(5) * t562 + mrSges(5,1) * t576 - mrSges(5,2) * t577 + t583 * t604 - t585 * t605 - t632;
t523 = m(4) * t581 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t609 + qJD(3) * t610 - t606 * t650 + t631;
t503 = t624 * t508 + t627 * t523;
t591 = -qJDD(1) * pkin(1) + t636;
t635 = -m(3) * t591 + (mrSges(3,3) * t630) - t503;
t501 = m(2) * t612 - (t630 * mrSges(2,2)) + qJDD(1) * t663 + t635;
t590 = t630 * pkin(1) + t666;
t509 = t513 * t626 + t514 * t623;
t634 = -m(4) * t588 + mrSges(4,1) * t608 - mrSges(4,2) * t609 - t610 * t649 - t611 * t650 - t509;
t633 = -m(3) * t590 + (mrSges(3,2) * t630) + qJDD(1) * mrSges(3,3) - t634;
t506 = m(2) * t613 - (mrSges(2,1) * t630) - qJDD(1) * mrSges(2,2) + t633;
t655 = t501 * t628 + t506 * t625;
t654 = t578 * t668 - t579 * t660 - t615 * t657;
t653 = t578 * t657 - t579 * t658 + t615 * t667;
t652 = -t578 * t660 + t579 * t669 + t615 * t658;
t643 = -t501 * t625 + t506 * t628;
t642 = t508 * t627 - t624 * t523;
t596 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t627 - Ifges(4,4) * t624) * qJD(1);
t595 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t627 - Ifges(4,2) * t624) * qJD(1);
t594 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t627 - Ifges(4,6) * t624) * qJD(1);
t571 = Ifges(5,1) * t605 + Ifges(5,4) * t604 + Ifges(5,5) * t615;
t570 = Ifges(5,4) * t605 + Ifges(5,2) * t604 + Ifges(5,6) * t615;
t569 = Ifges(5,5) * t605 + Ifges(5,6) * t604 + Ifges(5,3) * t615;
t517 = mrSges(6,2) * t536 + mrSges(7,2) * t527 - mrSges(6,3) * t528 - mrSges(7,3) * t531 - qJ(6) * t524 - t548 * t660 + t549 * t669 + t578 * t653 + t603 * t658 + t615 * t654;
t516 = -mrSges(6,1) * t536 - mrSges(7,1) * t531 + mrSges(7,2) * t526 + mrSges(6,3) * t529 - pkin(5) * t524 - t548 * t668 + t549 * t660 + t579 * t653 + t603 * t657 + t615 * t652;
t502 = -m(3) * g(3) + t642;
t499 = mrSges(5,2) * t562 - mrSges(5,3) * t537 + Ifges(5,1) * t577 + Ifges(5,4) * t576 + Ifges(5,5) * t603 - qJ(5) * t515 - t516 * t622 + t517 * t656 + t569 * t604 - t570 * t615;
t498 = -mrSges(5,1) * t562 + mrSges(5,3) * t538 + Ifges(5,4) * t577 + Ifges(5,2) * t576 + Ifges(5,6) * t603 - pkin(4) * t632 + qJ(5) * t640 + t516 * t656 + t622 * t517 - t605 * t569 + t615 * t571;
t497 = Ifges(4,2) * t608 + Ifges(4,4) * t609 + t604 * t571 - t605 * t570 + mrSges(4,3) * t582 - mrSges(4,1) * t588 + qJD(3) * t596 - Ifges(5,6) * t576 - Ifges(5,5) * t577 - mrSges(5,1) * t537 + mrSges(5,2) * t538 + mrSges(6,2) * t529 - mrSges(7,3) * t526 + mrSges(7,1) * t527 - mrSges(6,1) * t528 - pkin(4) * t515 - pkin(3) * t509 + (-Ifges(5,3) + t667) * t603 + Ifges(4,6) * qJDD(3) + (qJ(6) * t554 - t652) * t578 + (pkin(5) * t554 + t654) * t579 + (mrSges(7,2) * qJ(6) + t657) * t548 + (mrSges(7,2) * pkin(5) - t658) * t549 - pkin(5) * t639 - t594 * t650 - qJ(6) * t646;
t496 = mrSges(4,2) * t588 - mrSges(4,3) * t581 + Ifges(4,1) * t609 + Ifges(4,4) * t608 + Ifges(4,5) * qJDD(3) - pkin(8) * t509 - qJD(3) * t595 - t498 * t623 + t499 * t626 - t594 * t649;
t495 = -qJ(2) * t502 - mrSges(2,3) * t612 + pkin(2) * t503 + mrSges(3,1) * t591 + t626 * t498 + pkin(3) * t631 + pkin(8) * t641 + t623 * t499 + mrSges(4,1) * t581 - mrSges(4,2) * t582 + Ifges(4,5) * t609 + Ifges(4,6) * t608 + Ifges(4,3) * qJDD(3) + (t659 * t630) + t661 * qJDD(1) + (t595 * t627 + t596 * t624) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t494 = -mrSges(3,1) * t590 + mrSges(2,3) * t613 - pkin(1) * t502 - pkin(2) * t634 - pkin(7) * t642 + g(3) * t663 - qJDD(1) * t659 - t624 * t496 - t627 * t497 + t630 * t661;
t1 = [-m(1) * g(1) + t643; -m(1) * g(2) + t655; (-m(1) - m(2) - m(3)) * g(3) + t642; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t655 - t494 * t625 + t495 * t628; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t643 + t628 * t494 + t625 * t495; pkin(1) * t635 + qJ(2) * t633 - t624 * t497 - pkin(7) * t503 + mrSges(2,1) * t612 - mrSges(2,2) * t613 + t627 * t496 + mrSges(3,2) * t591 - mrSges(3,3) * t590 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
