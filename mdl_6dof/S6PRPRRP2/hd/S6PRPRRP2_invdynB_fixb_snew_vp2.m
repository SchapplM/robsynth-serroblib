% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRP2
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
% Datum: 2019-05-04 23:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:35:20
% EndTime: 2019-05-04 23:35:30
% DurationCPUTime: 7.16s
% Computational Cost: add. (91874->274), mult. (167010->340), div. (0->0), fcn. (113197->12), ass. (0->120)
t627 = Ifges(6,1) + Ifges(7,1);
t621 = Ifges(6,4) - Ifges(7,5);
t626 = -Ifges(6,5) - Ifges(7,4);
t625 = Ifges(6,2) + Ifges(7,3);
t619 = Ifges(6,6) - Ifges(7,6);
t624 = -Ifges(6,3) - Ifges(7,2);
t623 = cos(qJ(5));
t622 = -mrSges(6,3) - mrSges(7,2);
t584 = sin(pkin(6));
t590 = sin(qJ(2));
t618 = t584 * t590;
t592 = cos(qJ(2));
t617 = t584 * t592;
t587 = cos(pkin(6));
t616 = t587 * t590;
t615 = t587 * t592;
t583 = sin(pkin(10));
t586 = cos(pkin(10));
t570 = g(1) * t583 - g(2) * t586;
t571 = -g(1) * t586 - g(2) * t583;
t581 = -g(3) + qJDD(1);
t525 = t570 * t615 - t571 * t590 + t581 * t617;
t523 = qJDD(2) * pkin(2) + t525;
t526 = t570 * t616 + t592 * t571 + t581 * t618;
t594 = qJD(2) ^ 2;
t524 = -pkin(2) * t594 + t526;
t582 = sin(pkin(11));
t585 = cos(pkin(11));
t519 = t582 * t523 + t585 * t524;
t517 = -pkin(3) * t594 + qJDD(2) * pkin(8) + t519;
t551 = -t570 * t584 + t587 * t581;
t550 = qJDD(3) + t551;
t589 = sin(qJ(4));
t591 = cos(qJ(4));
t513 = t591 * t517 + t589 * t550;
t566 = (-mrSges(5,1) * t591 + mrSges(5,2) * t589) * qJD(2);
t607 = qJD(2) * qJD(4);
t604 = t589 * t607;
t569 = qJDD(2) * t591 - t604;
t609 = qJD(2) * t589;
t572 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t609;
t567 = (-pkin(4) * t591 - pkin(9) * t589) * qJD(2);
t593 = qJD(4) ^ 2;
t608 = qJD(2) * t591;
t509 = -pkin(4) * t593 + qJDD(4) * pkin(9) + t567 * t608 + t513;
t518 = t585 * t523 - t582 * t524;
t516 = -qJDD(2) * pkin(3) - t594 * pkin(8) - t518;
t603 = t591 * t607;
t568 = qJDD(2) * t589 + t603;
t511 = (-t568 - t603) * pkin(9) + (-t569 + t604) * pkin(4) + t516;
t588 = sin(qJ(5));
t505 = t509 * t623 + t588 * t511;
t565 = t588 * qJD(4) + t609 * t623;
t537 = qJD(5) * t565 - qJDD(4) * t623 + t568 * t588;
t577 = qJD(5) - t608;
t547 = mrSges(6,1) * t577 - mrSges(6,3) * t565;
t562 = qJDD(5) - t569;
t564 = -qJD(4) * t623 + t588 * t609;
t541 = pkin(5) * t564 - qJ(6) * t565;
t576 = t577 ^ 2;
t502 = -pkin(5) * t576 + qJ(6) * t562 + 0.2e1 * qJD(6) * t577 - t541 * t564 + t505;
t548 = -mrSges(7,1) * t577 + mrSges(7,2) * t565;
t606 = m(7) * t502 + t562 * mrSges(7,3) + t577 * t548;
t542 = mrSges(7,1) * t564 - mrSges(7,3) * t565;
t610 = -mrSges(6,1) * t564 - mrSges(6,2) * t565 - t542;
t496 = m(6) * t505 - t562 * mrSges(6,2) + t537 * t622 - t577 * t547 + t564 * t610 + t606;
t504 = -t588 * t509 + t511 * t623;
t538 = -t564 * qJD(5) + t588 * qJDD(4) + t568 * t623;
t546 = -mrSges(6,2) * t577 - mrSges(6,3) * t564;
t503 = -t562 * pkin(5) - t576 * qJ(6) + t565 * t541 + qJDD(6) - t504;
t549 = -mrSges(7,2) * t564 + mrSges(7,3) * t577;
t598 = -m(7) * t503 + t562 * mrSges(7,1) + t577 * t549;
t497 = m(6) * t504 + t562 * mrSges(6,1) + t538 * t622 + t577 * t546 + t565 * t610 + t598;
t599 = t496 * t623 - t497 * t588;
t491 = m(5) * t513 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t569 - qJD(4) * t572 + t566 * t608 + t599;
t512 = -t589 * t517 + t591 * t550;
t573 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t608;
t508 = -qJDD(4) * pkin(4) - t593 * pkin(9) + t567 * t609 - t512;
t506 = -0.2e1 * qJD(6) * t565 + (t564 * t577 - t538) * qJ(6) + (t565 * t577 + t537) * pkin(5) + t508;
t500 = m(7) * t506 + mrSges(7,1) * t537 - t538 * mrSges(7,3) - t565 * t548 + t549 * t564;
t595 = -m(6) * t508 - t537 * mrSges(6,1) - mrSges(6,2) * t538 - t564 * t546 - t547 * t565 - t500;
t499 = m(5) * t512 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t568 + qJD(4) * t573 - t566 * t609 + t595;
t600 = t591 * t491 - t499 * t589;
t483 = m(4) * t519 - mrSges(4,1) * t594 - qJDD(2) * mrSges(4,2) + t600;
t494 = t588 * t496 + t497 * t623;
t596 = -m(5) * t516 + t569 * mrSges(5,1) - t568 * mrSges(5,2) - t572 * t609 + t573 * t608 - t494;
t488 = m(4) * t518 + qJDD(2) * mrSges(4,1) - t594 * mrSges(4,2) + t596;
t479 = t582 * t483 + t585 * t488;
t477 = m(3) * t525 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t594 + t479;
t601 = t585 * t483 - t488 * t582;
t478 = m(3) * t526 - mrSges(3,1) * t594 - qJDD(2) * mrSges(3,2) + t601;
t486 = t589 * t491 + t591 * t499;
t605 = m(4) * t550 + t486;
t485 = m(3) * t551 + t605;
t465 = t477 * t615 + t478 * t616 - t485 * t584;
t463 = m(2) * t570 + t465;
t469 = -t477 * t590 + t592 * t478;
t468 = m(2) * t571 + t469;
t614 = t586 * t463 + t583 * t468;
t613 = t625 * t564 - t621 * t565 - t619 * t577;
t612 = t619 * t564 + t626 * t565 + t624 * t577;
t611 = -t621 * t564 + t627 * t565 - t626 * t577;
t464 = t477 * t617 + t478 * t618 + t587 * t485;
t602 = -t463 * t583 + t586 * t468;
t492 = -mrSges(6,1) * t508 - mrSges(7,1) * t506 + mrSges(7,2) * t502 + mrSges(6,3) * t505 - pkin(5) * t500 - t625 * t537 + t621 * t538 + t619 * t562 + t612 * t565 + t611 * t577;
t493 = mrSges(6,2) * t508 + mrSges(7,2) * t503 - mrSges(6,3) * t504 - mrSges(7,3) * t506 - qJ(6) * t500 - t621 * t537 + t627 * t538 - t562 * t626 + t612 * t564 + t613 * t577;
t555 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t589 + Ifges(5,6) * t591) * qJD(2);
t556 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t589 + Ifges(5,2) * t591) * qJD(2);
t471 = mrSges(5,2) * t516 - mrSges(5,3) * t512 + Ifges(5,1) * t568 + Ifges(5,4) * t569 + Ifges(5,5) * qJDD(4) - pkin(9) * t494 - qJD(4) * t556 - t588 * t492 + t493 * t623 + t555 * t608;
t557 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t589 + Ifges(5,4) * t591) * qJD(2);
t480 = Ifges(5,4) * t568 + Ifges(5,2) * t569 + Ifges(5,6) * qJDD(4) - t555 * t609 + qJD(4) * t557 - mrSges(5,1) * t516 + mrSges(5,3) * t513 - mrSges(6,1) * t504 + mrSges(6,2) * t505 + mrSges(7,1) * t503 - mrSges(7,3) * t502 - pkin(5) * t598 - qJ(6) * t606 - pkin(4) * t494 + (pkin(5) * t542 + t613) * t565 + (qJ(6) * t542 - t611) * t564 + t624 * t562 + (mrSges(7,2) * pkin(5) + t626) * t538 + (mrSges(7,2) * qJ(6) + t619) * t537;
t461 = mrSges(4,2) * t550 - mrSges(4,3) * t518 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t594 - pkin(8) * t486 + t471 * t591 - t480 * t589;
t470 = Ifges(4,6) * qJDD(2) + t594 * Ifges(4,5) - mrSges(4,1) * t550 + mrSges(4,3) * t519 - Ifges(5,5) * t568 - Ifges(5,6) * t569 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t512 + mrSges(5,2) * t513 - t588 * t493 - t623 * t492 - pkin(4) * t595 - pkin(9) * t599 - pkin(3) * t486 + (-t556 * t589 + t557 * t591) * qJD(2);
t458 = -mrSges(3,1) * t551 + mrSges(3,3) * t526 + t594 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t605 + qJ(3) * t601 + t582 * t461 + t585 * t470;
t459 = mrSges(3,2) * t551 - mrSges(3,3) * t525 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t594 - qJ(3) * t479 + t461 * t585 - t470 * t582;
t597 = pkin(7) * t469 + t458 * t592 + t459 * t590;
t460 = mrSges(3,1) * t525 - mrSges(3,2) * t526 + mrSges(4,1) * t518 - mrSges(4,2) * t519 + t589 * t471 + t591 * t480 + pkin(3) * t596 + pkin(8) * t600 + pkin(2) * t479 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t457 = mrSges(2,2) * t581 - mrSges(2,3) * t570 - t590 * t458 + t592 * t459 + (-t464 * t584 - t465 * t587) * pkin(7);
t456 = -mrSges(2,1) * t581 + mrSges(2,3) * t571 - pkin(1) * t464 - t584 * t460 + t587 * t597;
t1 = [-m(1) * g(1) + t602; -m(1) * g(2) + t614; -m(1) * g(3) + m(2) * t581 + t464; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t614 - t583 * t456 + t586 * t457; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t602 + t586 * t456 + t583 * t457; -mrSges(1,1) * g(2) + mrSges(2,1) * t570 + mrSges(1,2) * g(1) - mrSges(2,2) * t571 + pkin(1) * t465 + t587 * t460 + t584 * t597;];
tauB  = t1;
