% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRRP2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:32:27
% EndTime: 2019-05-04 20:32:37
% DurationCPUTime: 9.99s
% Computational Cost: add. (176261->270), mult. (314728->344), div. (0->0), fcn. (240030->14), ass. (0->124)
t622 = Ifges(6,1) + Ifges(7,1);
t615 = Ifges(6,4) - Ifges(7,5);
t621 = -Ifges(6,5) - Ifges(7,4);
t620 = Ifges(6,2) + Ifges(7,3);
t613 = Ifges(6,6) - Ifges(7,6);
t619 = -Ifges(6,3) - Ifges(7,2);
t573 = sin(pkin(11));
t577 = cos(pkin(11));
t564 = -t577 * g(1) - t573 * g(2);
t572 = sin(pkin(12));
t576 = cos(pkin(12));
t563 = t573 * g(1) - t577 * g(2);
t571 = -g(3) + qJDD(1);
t575 = sin(pkin(6));
t579 = cos(pkin(6));
t591 = t563 * t579 + t571 * t575;
t523 = -t572 * t564 + t591 * t576;
t524 = t576 * t564 + t591 * t572;
t547 = -t575 * t563 + t579 * t571 + qJDD(2);
t584 = cos(qJ(3));
t578 = cos(pkin(7));
t582 = sin(qJ(3));
t608 = t578 * t582;
t574 = sin(pkin(7));
t609 = t574 * t582;
t517 = t523 * t608 + t584 * t524 + t547 * t609;
t586 = qJD(3) ^ 2;
t515 = -t586 * pkin(3) + qJDD(3) * pkin(9) + t517;
t519 = -t574 * t523 + t578 * t547;
t581 = sin(qJ(4));
t583 = cos(qJ(4));
t509 = t583 * t515 + t581 * t519;
t559 = (-mrSges(5,1) * t583 + mrSges(5,2) * t581) * qJD(3);
t600 = qJD(3) * qJD(4);
t598 = t581 * t600;
t562 = t583 * qJDD(3) - t598;
t602 = qJD(3) * t581;
t565 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t602;
t560 = (-pkin(4) * t583 - pkin(10) * t581) * qJD(3);
t585 = qJD(4) ^ 2;
t601 = t583 * qJD(3);
t507 = -t585 * pkin(4) + qJDD(4) * pkin(10) + t560 * t601 + t509;
t516 = -t582 * t524 + (t523 * t578 + t547 * t574) * t584;
t514 = -qJDD(3) * pkin(3) - t586 * pkin(9) - t516;
t597 = t583 * t600;
t561 = t581 * qJDD(3) + t597;
t511 = (-t561 - t597) * pkin(10) + (-t562 + t598) * pkin(4) + t514;
t580 = sin(qJ(5));
t617 = cos(qJ(5));
t503 = t617 * t507 + t580 * t511;
t558 = t580 * qJD(4) + t617 * t602;
t535 = t558 * qJD(5) - t617 * qJDD(4) + t580 * t561;
t568 = qJD(5) - t601;
t544 = t568 * mrSges(6,1) - t558 * mrSges(6,3);
t556 = qJDD(5) - t562;
t557 = -t617 * qJD(4) + t580 * t602;
t539 = t557 * pkin(5) - t558 * qJ(6);
t567 = t568 ^ 2;
t500 = -t567 * pkin(5) + t556 * qJ(6) + 0.2e1 * qJD(6) * t568 - t557 * t539 + t503;
t545 = -t568 * mrSges(7,1) + t558 * mrSges(7,2);
t599 = m(7) * t500 + t556 * mrSges(7,3) + t568 * t545;
t540 = t557 * mrSges(7,1) - t558 * mrSges(7,3);
t603 = -t557 * mrSges(6,1) - t558 * mrSges(6,2) - t540;
t616 = -mrSges(6,3) - mrSges(7,2);
t496 = m(6) * t503 - t556 * mrSges(6,2) + t616 * t535 - t568 * t544 + t603 * t557 + t599;
t502 = -t580 * t507 + t617 * t511;
t536 = -t557 * qJD(5) + t580 * qJDD(4) + t617 * t561;
t543 = -t568 * mrSges(6,2) - t557 * mrSges(6,3);
t501 = -t556 * pkin(5) - t567 * qJ(6) + t558 * t539 + qJDD(6) - t502;
t546 = -t557 * mrSges(7,2) + t568 * mrSges(7,3);
t593 = -m(7) * t501 + t556 * mrSges(7,1) + t568 * t546;
t497 = m(6) * t502 + t556 * mrSges(6,1) + t616 * t536 + t568 * t543 + t603 * t558 + t593;
t594 = t617 * t496 - t580 * t497;
t489 = m(5) * t509 - qJDD(4) * mrSges(5,2) + t562 * mrSges(5,3) - qJD(4) * t565 + t559 * t601 + t594;
t508 = -t581 * t515 + t583 * t519;
t566 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t601;
t506 = -qJDD(4) * pkin(4) - t585 * pkin(10) + t560 * t602 - t508;
t504 = -0.2e1 * qJD(6) * t558 + (t557 * t568 - t536) * qJ(6) + (t558 * t568 + t535) * pkin(5) + t506;
t498 = m(7) * t504 + t535 * mrSges(7,1) - t536 * mrSges(7,3) - t558 * t545 + t557 * t546;
t587 = -m(6) * t506 - t535 * mrSges(6,1) - t536 * mrSges(6,2) - t557 * t543 - t558 * t544 - t498;
t494 = m(5) * t508 + qJDD(4) * mrSges(5,1) - t561 * mrSges(5,3) + qJD(4) * t566 - t559 * t602 + t587;
t595 = t583 * t489 - t581 * t494;
t480 = m(4) * t517 - t586 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t595;
t483 = t581 * t489 + t583 * t494;
t482 = m(4) * t519 + t483;
t492 = t580 * t496 + t617 * t497;
t588 = -m(5) * t514 + t562 * mrSges(5,1) - t561 * mrSges(5,2) - t565 * t602 + t566 * t601 - t492;
t486 = m(4) * t516 + qJDD(3) * mrSges(4,1) - t586 * mrSges(4,2) + t588;
t610 = t486 * t584;
t469 = t480 * t608 - t574 * t482 + t578 * t610;
t465 = m(3) * t523 + t469;
t475 = t584 * t480 - t582 * t486;
t474 = m(3) * t524 + t475;
t618 = t465 * t576 + t474 * t572;
t468 = t480 * t609 + t578 * t482 + t574 * t610;
t467 = m(3) * t547 + t468;
t455 = -t575 * t467 + t618 * t579;
t453 = m(2) * t563 + t455;
t461 = -t572 * t465 + t576 * t474;
t460 = m(2) * t564 + t461;
t607 = t577 * t453 + t573 * t460;
t606 = t620 * t557 - t615 * t558 - t613 * t568;
t605 = t613 * t557 + t621 * t558 + t619 * t568;
t604 = -t615 * t557 + t622 * t558 - t621 * t568;
t454 = t579 * t467 + t618 * t575;
t596 = -t573 * t453 + t577 * t460;
t490 = -mrSges(6,1) * t506 - mrSges(7,1) * t504 + mrSges(7,2) * t500 + mrSges(6,3) * t503 - pkin(5) * t498 - t620 * t535 + t615 * t536 + t613 * t556 + t605 * t558 + t604 * t568;
t491 = mrSges(6,2) * t506 + mrSges(7,2) * t501 - mrSges(6,3) * t502 - mrSges(7,3) * t504 - qJ(6) * t498 - t615 * t535 + t622 * t536 - t556 * t621 + t605 * t557 + t606 * t568;
t549 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t581 + Ifges(5,6) * t583) * qJD(3);
t550 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t581 + Ifges(5,2) * t583) * qJD(3);
t470 = mrSges(5,2) * t514 - mrSges(5,3) * t508 + Ifges(5,1) * t561 + Ifges(5,4) * t562 + Ifges(5,5) * qJDD(4) - pkin(10) * t492 - qJD(4) * t550 - t580 * t490 + t617 * t491 + t549 * t601;
t551 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t581 + Ifges(5,4) * t583) * qJD(3);
t476 = Ifges(5,4) * t561 + Ifges(5,2) * t562 + Ifges(5,6) * qJDD(4) - t549 * t602 + qJD(4) * t551 - mrSges(5,1) * t514 + mrSges(5,3) * t509 - mrSges(6,1) * t502 + mrSges(6,2) * t503 + mrSges(7,1) * t501 - mrSges(7,3) * t500 - pkin(5) * t593 - qJ(6) * t599 - pkin(4) * t492 + (pkin(5) * t540 + t606) * t558 + (qJ(6) * t540 - t604) * t557 + t619 * t556 + (pkin(5) * mrSges(7,2) + t621) * t536 + (qJ(6) * mrSges(7,2) + t613) * t535;
t457 = mrSges(4,2) * t519 - mrSges(4,3) * t516 + Ifges(4,5) * qJDD(3) - t586 * Ifges(4,6) - pkin(9) * t483 + t583 * t470 - t581 * t476;
t462 = Ifges(4,6) * qJDD(3) + t586 * Ifges(4,5) - mrSges(4,1) * t519 + mrSges(4,3) * t517 - Ifges(5,5) * t561 - Ifges(5,6) * t562 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t508 + mrSges(5,2) * t509 - t580 * t491 - t617 * t490 - pkin(4) * t587 - pkin(10) * t594 - pkin(3) * t483 + (-t581 * t550 + t583 * t551) * qJD(3);
t590 = pkin(8) * t475 + t457 * t582 + t462 * t584;
t456 = mrSges(4,1) * t516 - mrSges(4,2) * t517 + Ifges(4,3) * qJDD(3) + pkin(3) * t588 + pkin(9) * t595 + t581 * t470 + t583 * t476;
t450 = -mrSges(3,1) * t547 + mrSges(3,3) * t524 - pkin(2) * t468 - t574 * t456 + t590 * t578;
t451 = mrSges(3,2) * t547 - mrSges(3,3) * t523 + t584 * t457 - t582 * t462 + (-t468 * t574 - t469 * t578) * pkin(8);
t589 = qJ(2) * t461 + t450 * t576 + t451 * t572;
t449 = mrSges(3,1) * t523 - mrSges(3,2) * t524 + pkin(2) * t469 + t578 * t456 + t590 * t574;
t448 = mrSges(2,2) * t571 - mrSges(2,3) * t563 - t572 * t450 + t576 * t451 + (-t454 * t575 - t455 * t579) * qJ(2);
t447 = -mrSges(2,1) * t571 + mrSges(2,3) * t564 - pkin(1) * t454 - t575 * t449 + t589 * t579;
t1 = [-m(1) * g(1) + t596; -m(1) * g(2) + t607; -m(1) * g(3) + m(2) * t571 + t454; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t607 - t573 * t447 + t577 * t448; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t596 + t577 * t447 + t573 * t448; -mrSges(1,1) * g(2) + mrSges(2,1) * t563 + mrSges(1,2) * g(1) - mrSges(2,2) * t564 + pkin(1) * t455 + t579 * t449 + t589 * t575;];
tauB  = t1;
