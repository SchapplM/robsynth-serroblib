% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-05-05 14:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:26:37
% EndTime: 2019-05-05 14:26:40
% DurationCPUTime: 1.61s
% Computational Cost: add. (10169->273), mult. (19316->306), div. (0->0), fcn. (8412->6), ass. (0->107)
t633 = 2 * qJD(1);
t632 = Ifges(5,1) + Ifges(6,2);
t621 = Ifges(5,4) + Ifges(6,6);
t620 = (Ifges(5,5) - Ifges(6,4));
t631 = Ifges(5,2) + Ifges(6,3);
t619 = (Ifges(5,6) - Ifges(6,5));
t630 = (-Ifges(5,3) - Ifges(6,1));
t583 = sin(qJ(1));
t586 = cos(qJ(1));
t559 = -t586 * g(1) - t583 * g(2);
t629 = qJDD(1) * qJ(2) + (qJD(2) * t633) + t559;
t558 = t583 * g(1) - t586 * g(2);
t588 = qJD(1) ^ 2;
t527 = -qJDD(1) * pkin(1) - t588 * qJ(2) + qJDD(2) - t558;
t598 = qJDD(1) * qJ(3) + (qJD(3) * t633) - t527;
t628 = -2 * qJD(5);
t627 = -m(3) - m(4);
t626 = pkin(4) + pkin(8);
t625 = pkin(8) * t588;
t582 = sin(qJ(4));
t624 = t582 * g(3);
t623 = t588 * pkin(7);
t622 = mrSges(2,1) - mrSges(3,2);
t523 = qJDD(3) + (-pkin(1) - qJ(3)) * t588 + t629;
t520 = -qJDD(1) * pkin(7) + t523;
t585 = cos(qJ(4));
t618 = t585 * t520;
t526 = t588 * pkin(1) - t629;
t514 = t618 + t624;
t611 = qJD(1) * qJD(4);
t564 = t582 * t611;
t551 = t585 * qJDD(1) - t564;
t613 = qJD(1) * t582;
t553 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t613;
t555 = mrSges(6,1) * t613 - (qJD(4) * mrSges(6,3));
t606 = t585 * t611;
t550 = t582 * qJDD(1) + t606;
t612 = t585 * qJD(1);
t557 = pkin(5) * t612 - (qJD(4) * pkin(8));
t576 = t582 ^ 2;
t590 = pkin(4) * t606 + t612 * t628 + t598 + (-t551 + t564) * qJ(5);
t499 = (-pkin(5) * t576 - pkin(7)) * t588 + t626 * t550 + t590 - t557 * t612;
t547 = (pkin(4) * t582 - qJ(5) * t585) * qJD(1);
t587 = qJD(4) ^ 2;
t597 = -t587 * qJ(5) + t547 * t612 + qJDD(5) - t618;
t502 = t551 * pkin(5) - t626 * qJDD(4) + (pkin(5) * t611 + t585 * t625 - g(3)) * t582 + t597;
t581 = sin(qJ(6));
t584 = cos(qJ(6));
t497 = -t581 * t499 + t584 * t502;
t545 = -t581 * qJD(4) + t584 * t613;
t513 = t545 * qJD(6) + t584 * qJDD(4) + t581 * t550;
t546 = t584 * qJD(4) + t581 * t613;
t517 = -t545 * mrSges(7,1) + t546 * mrSges(7,2);
t562 = qJD(6) + t612;
t524 = -t562 * mrSges(7,2) + t545 * mrSges(7,3);
t544 = qJDD(6) + t551;
t495 = m(7) * t497 + t544 * mrSges(7,1) - t513 * mrSges(7,3) - t546 * t517 + t562 * t524;
t498 = t584 * t499 + t581 * t502;
t512 = -t546 * qJD(6) - t581 * qJDD(4) + t584 * t550;
t525 = t562 * mrSges(7,1) - t546 * mrSges(7,3);
t496 = m(7) * t498 - t544 * mrSges(7,2) + t512 * mrSges(7,3) + t545 * t517 - t562 * t525;
t487 = t584 * t495 + t581 * t496;
t506 = -qJDD(4) * pkin(4) + t597 - t624;
t594 = -m(6) * t506 - t551 * mrSges(6,1) - t487;
t548 = (-mrSges(6,2) * t582 - mrSges(6,3) * t585) * qJD(1);
t602 = qJD(1) * (-t548 - (mrSges(5,1) * t582 + mrSges(5,2) * t585) * qJD(1));
t485 = m(5) * t514 - t551 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t553 - t555) * qJD(4) + t585 * t602 + t594;
t515 = -t585 * g(3) + t582 * t520;
t554 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t612;
t592 = -t587 * pkin(4) + qJDD(4) * qJ(5) - t547 * t613 + t515;
t505 = (qJD(4) * t628) - t592;
t556 = mrSges(6,1) * t612 + qJD(4) * mrSges(6,2);
t501 = -t576 * t625 - t550 * pkin(5) + ((2 * qJD(5)) + t557) * qJD(4) + t592;
t596 = -m(7) * t501 + t512 * mrSges(7,1) - t513 * mrSges(7,2) + t545 * t524 - t546 * t525;
t591 = -m(6) * t505 + qJDD(4) * mrSges(6,3) + qJD(4) * t556 - t596;
t492 = m(5) * t515 - qJDD(4) * mrSges(5,2) - qJD(4) * t554 + (-mrSges(5,3) - mrSges(6,1)) * t550 + t582 * t602 + t591;
t480 = t585 * t485 + t582 * t492;
t601 = -m(4) * t523 - qJDD(1) * mrSges(4,2) - t480;
t595 = -m(3) * t526 + (t588 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t601;
t478 = m(2) * t559 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t588) + t595;
t504 = t550 * pkin(4) + t590 - t623;
t616 = -t581 * t495 + t584 * t496;
t486 = m(6) * t504 - t550 * mrSges(6,2) - t551 * mrSges(6,3) - t555 * t613 - t556 * t612 + t616;
t519 = t598 - t623;
t593 = -m(5) * t519 - t550 * mrSges(5,1) - t551 * mrSges(5,2) - t553 * t613 - t554 * t612 - t486;
t483 = -m(4) * t598 - t588 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t593;
t589 = -m(3) * t527 + t588 * mrSges(3,3) - t483;
t482 = m(2) * t558 - t588 * mrSges(2,2) + t622 * qJDD(1) + t589;
t617 = t583 * t478 + t586 * t482;
t615 = -(t619 * qJD(4)) + (t582 * t631 - t585 * t621) * qJD(1);
t614 = (t620 * qJD(4)) + (-t582 * t621 + t585 * t632) * qJD(1);
t608 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t607 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t605 = t586 * t478 - t583 * t482;
t604 = -t582 * t485 + t585 * t492;
t603 = qJD(1) * ((t630 * qJD(4)) + (t582 * t619 - t585 * t620) * qJD(1));
t509 = Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t562;
t508 = Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t562;
t507 = Ifges(7,5) * t546 + Ifges(7,6) * t545 + Ifges(7,3) * t562;
t489 = mrSges(7,2) * t501 - mrSges(7,3) * t497 + Ifges(7,1) * t513 + Ifges(7,4) * t512 + Ifges(7,5) * t544 + t545 * t507 - t562 * t508;
t488 = -mrSges(7,1) * t501 + mrSges(7,3) * t498 + Ifges(7,4) * t513 + Ifges(7,2) * t512 + Ifges(7,6) * t544 - t546 * t507 + t562 * t509;
t479 = t627 * g(3) + t604;
t475 = mrSges(6,1) * t506 + mrSges(7,1) * t497 + mrSges(5,2) * t519 - mrSges(7,2) * t498 - mrSges(5,3) * t514 - mrSges(6,3) * t504 + Ifges(7,5) * t513 + Ifges(7,6) * t512 + Ifges(7,3) * t544 + pkin(5) * t487 - qJ(5) * t486 + t546 * t508 - t545 * t509 + t632 * t551 - t621 * t550 + t620 * qJDD(4) + t615 * qJD(4) + t582 * t603;
t474 = -mrSges(5,1) * t519 - mrSges(6,1) * t505 + mrSges(6,2) * t504 + mrSges(5,3) * t515 - pkin(4) * t486 - pkin(5) * t596 - pkin(8) * t616 + t614 * qJD(4) + t619 * qJDD(4) - t584 * t488 - t581 * t489 - t550 * t631 + t621 * t551 + t585 * t603;
t473 = -pkin(8) * t487 - pkin(1) * t479 + pkin(3) * t480 - pkin(2) * t601 + mrSges(2,3) * t559 - qJ(3) * t604 - t581 * t488 + t584 * t489 - mrSges(6,3) * t505 + mrSges(6,2) * t506 + mrSges(5,1) * t514 - mrSges(5,2) * t515 + pkin(4) * (-qJD(4) * t555 + t594) + qJ(5) * t591 + mrSges(4,1) * t523 - mrSges(3,1) * t526 + t620 * t551 + (-qJ(5) * mrSges(6,1) - t619) * t550 + (-pkin(4) * mrSges(6,2) - t630) * qJDD(4) + (-pkin(2) * mrSges(4,3) + t608) * t588 + t607 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t622) * g(3) + ((-pkin(4) * t548 - t615) * t585 + (-qJ(5) * t548 + t614) * t582) * qJD(1);
t472 = -qJ(2) * t479 - mrSges(2,3) * t558 + pkin(2) * t483 + mrSges(3,1) * t527 + t585 * t474 + pkin(3) * t593 + pkin(7) * t604 + t582 * t475 - mrSges(4,1) * t598 - t607 * t588 + t608 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t605; -m(1) * g(2) + t617; (-m(1) - m(2) + t627) * g(3) + t604; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t617 + t586 * t472 - t583 * t473; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t605 + t583 * t472 + t586 * t473; qJ(2) * (-(t588 * mrSges(4,3)) + t595) + pkin(1) * t589 + mrSges(2,1) * t558 - mrSges(2,2) * t559 - qJ(3) * t483 + mrSges(3,2) * t527 - mrSges(3,3) * t526 - t582 * t474 - pkin(7) * t480 + t585 * t475 + mrSges(4,2) * t523 + mrSges(4,3) * t598 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
