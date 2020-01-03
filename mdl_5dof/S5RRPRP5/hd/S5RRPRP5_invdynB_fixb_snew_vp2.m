% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:23
% EndTime: 2019-12-31 19:54:28
% DurationCPUTime: 4.67s
% Computational Cost: add. (45939->289), mult. (106890->358), div. (0->0), fcn. (71966->8), ass. (0->111)
t606 = Ifges(5,1) + Ifges(6,1);
t600 = Ifges(5,4) - Ifges(6,5);
t599 = Ifges(6,4) + Ifges(5,5);
t605 = Ifges(5,2) + Ifges(6,3);
t604 = -Ifges(6,2) - Ifges(5,3);
t598 = Ifges(5,6) - Ifges(6,6);
t603 = cos(qJ(4));
t578 = qJD(1) ^ 2;
t602 = pkin(2) * t578;
t601 = -mrSges(5,3) - mrSges(6,2);
t575 = sin(qJ(1));
t577 = cos(qJ(1));
t563 = -g(1) * t577 - g(2) * t575;
t552 = -pkin(1) * t578 + qJDD(1) * pkin(6) + t563;
t574 = sin(qJ(2));
t597 = t574 * t552;
t576 = cos(qJ(2));
t589 = qJD(1) * qJD(2);
t557 = qJDD(1) * t574 + t576 * t589;
t515 = qJDD(2) * pkin(2) - t557 * qJ(3) - t597 + (qJ(3) * t589 + t574 * t602 - g(3)) * t576;
t537 = -g(3) * t574 + t576 * t552;
t558 = qJDD(1) * t576 - t574 * t589;
t591 = qJD(1) * t574;
t559 = qJD(2) * pkin(2) - qJ(3) * t591;
t570 = t576 ^ 2;
t516 = qJ(3) * t558 - qJD(2) * t559 - t570 * t602 + t537;
t571 = sin(pkin(8));
t572 = cos(pkin(8));
t547 = (t571 * t576 + t572 * t574) * qJD(1);
t491 = -0.2e1 * qJD(3) * t547 + t572 * t515 - t571 * t516;
t535 = t557 * t572 + t558 * t571;
t546 = (-t571 * t574 + t572 * t576) * qJD(1);
t486 = (qJD(2) * t546 - t535) * pkin(7) + (t546 * t547 + qJDD(2)) * pkin(3) + t491;
t492 = 0.2e1 * qJD(3) * t546 + t571 * t515 + t572 * t516;
t534 = -t557 * t571 + t558 * t572;
t540 = qJD(2) * pkin(3) - pkin(7) * t547;
t545 = t546 ^ 2;
t488 = -pkin(3) * t545 + pkin(7) * t534 - qJD(2) * t540 + t492;
t573 = sin(qJ(4));
t484 = t573 * t486 + t603 * t488;
t528 = t546 * t573 + t547 * t603;
t497 = qJD(4) * t528 - t534 * t603 + t535 * t573;
t569 = qJD(2) + qJD(4);
t521 = mrSges(5,1) * t569 - mrSges(5,3) * t528;
t527 = -t546 * t603 + t547 * t573;
t568 = qJDD(2) + qJDD(4);
t509 = pkin(4) * t527 - qJ(5) * t528;
t567 = t569 ^ 2;
t479 = -pkin(4) * t567 + qJ(5) * t568 + 0.2e1 * qJD(5) * t569 - t509 * t527 + t484;
t522 = -mrSges(6,1) * t569 + mrSges(6,2) * t528;
t588 = m(6) * t479 + t568 * mrSges(6,3) + t569 * t522;
t510 = mrSges(6,1) * t527 - mrSges(6,3) * t528;
t592 = -mrSges(5,1) * t527 - mrSges(5,2) * t528 - t510;
t474 = m(5) * t484 - t568 * mrSges(5,2) + t497 * t601 - t569 * t521 + t527 * t592 + t588;
t483 = t486 * t603 - t488 * t573;
t498 = -qJD(4) * t527 + t534 * t573 + t535 * t603;
t520 = -mrSges(5,2) * t569 - mrSges(5,3) * t527;
t480 = -pkin(4) * t568 - qJ(5) * t567 + t509 * t528 + qJDD(5) - t483;
t523 = -mrSges(6,2) * t527 + mrSges(6,3) * t569;
t583 = -m(6) * t480 + t568 * mrSges(6,1) + t569 * t523;
t476 = m(5) * t483 + t568 * mrSges(5,1) + t498 * t601 + t569 * t520 + t528 * t592 + t583;
t469 = t474 * t573 + t476 * t603;
t531 = -mrSges(4,1) * t546 + mrSges(4,2) * t547;
t538 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t546;
t465 = m(4) * t491 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t535 + qJD(2) * t538 - t531 * t547 + t469;
t539 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t547;
t584 = t474 * t603 - t476 * t573;
t466 = m(4) * t492 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t534 - qJD(2) * t539 + t531 * t546 + t584;
t461 = t465 * t572 + t466 * t571;
t536 = -t576 * g(3) - t597;
t556 = (-mrSges(3,1) * t576 + mrSges(3,2) * t574) * qJD(1);
t590 = qJD(1) * t576;
t561 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t590;
t459 = m(3) * t536 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t557 + qJD(2) * t561 - t556 * t591 + t461;
t560 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t591;
t585 = -t465 * t571 + t466 * t572;
t460 = m(3) * t537 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t558 - qJD(2) * t560 + t556 * t590 + t585;
t586 = -t459 * t574 + t460 * t576;
t452 = m(2) * t563 - mrSges(2,1) * t578 - qJDD(1) * mrSges(2,2) + t586;
t562 = t575 * g(1) - t577 * g(2);
t582 = -qJDD(1) * pkin(1) - t562;
t551 = -t578 * pkin(6) + t582;
t519 = -t558 * pkin(2) + qJDD(3) + t559 * t591 + (-qJ(3) * t570 - pkin(6)) * t578 + t582;
t490 = -t534 * pkin(3) - t545 * pkin(7) + t547 * t540 + t519;
t482 = -0.2e1 * qJD(5) * t528 + (t527 * t569 - t498) * qJ(5) + (t528 * t569 + t497) * pkin(4) + t490;
t477 = m(6) * t482 + t497 * mrSges(6,1) - t498 * mrSges(6,3) - t528 * t522 + t527 * t523;
t581 = m(5) * t490 + t497 * mrSges(5,1) + t498 * mrSges(5,2) + t527 * t520 + t528 * t521 + t477;
t580 = m(4) * t519 - t534 * mrSges(4,1) + mrSges(4,2) * t535 - t546 * t538 + t539 * t547 + t581;
t579 = -m(3) * t551 + t558 * mrSges(3,1) - mrSges(3,2) * t557 - t560 * t591 + t561 * t590 - t580;
t471 = m(2) * t562 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t578 + t579;
t596 = t452 * t575 + t471 * t577;
t453 = t459 * t576 + t460 * t574;
t595 = t527 * t605 - t528 * t600 - t569 * t598;
t594 = t527 * t598 - t528 * t599 + t569 * t604;
t593 = -t527 * t600 + t528 * t606 + t569 * t599;
t587 = t452 * t577 - t471 * t575;
t550 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t574 + Ifges(3,4) * t576) * qJD(1);
t549 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t574 + Ifges(3,2) * t576) * qJD(1);
t548 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t574 + Ifges(3,6) * t576) * qJD(1);
t526 = Ifges(4,1) * t547 + Ifges(4,4) * t546 + Ifges(4,5) * qJD(2);
t525 = Ifges(4,4) * t547 + Ifges(4,2) * t546 + Ifges(4,6) * qJD(2);
t524 = Ifges(4,5) * t547 + Ifges(4,6) * t546 + Ifges(4,3) * qJD(2);
t468 = mrSges(5,2) * t490 + mrSges(6,2) * t480 - mrSges(5,3) * t483 - mrSges(6,3) * t482 - qJ(5) * t477 - t600 * t497 + t498 * t606 + t594 * t527 + t599 * t568 + t595 * t569;
t467 = -mrSges(5,1) * t490 - mrSges(6,1) * t482 + mrSges(6,2) * t479 + mrSges(5,3) * t484 - pkin(4) * t477 - t497 * t605 + t600 * t498 + t594 * t528 + t598 * t568 + t593 * t569;
t455 = mrSges(4,2) * t519 - mrSges(4,3) * t491 + Ifges(4,1) * t535 + Ifges(4,4) * t534 + Ifges(4,5) * qJDD(2) - pkin(7) * t469 - qJD(2) * t525 - t467 * t573 + t468 * t603 + t524 * t546;
t454 = -mrSges(4,1) * t519 + mrSges(4,3) * t492 + Ifges(4,4) * t535 + Ifges(4,2) * t534 + Ifges(4,6) * qJDD(2) - pkin(3) * t581 + pkin(7) * t584 + qJD(2) * t526 + t467 * t603 + t573 * t468 - t547 * t524;
t449 = mrSges(3,2) * t551 - mrSges(3,3) * t536 + Ifges(3,1) * t557 + Ifges(3,4) * t558 + Ifges(3,5) * qJDD(2) - qJ(3) * t461 - qJD(2) * t549 - t454 * t571 + t455 * t572 + t548 * t590;
t448 = t604 * t568 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t549 * t574 + t550 * t576) * qJD(1) + t578 * Ifges(2,5) - Ifges(3,5) * t557 - Ifges(3,6) * t558 + mrSges(2,3) * t563 + t546 * t526 - t547 * t525 - Ifges(4,6) * t534 - Ifges(4,5) * t535 - mrSges(3,1) * t536 + mrSges(3,2) * t537 - mrSges(4,1) * t491 + mrSges(4,2) * t492 - mrSges(5,1) * t483 + mrSges(5,2) * t484 - mrSges(6,3) * t479 + mrSges(6,1) * t480 - pkin(3) * t469 - pkin(2) * t461 + (mrSges(6,2) * pkin(4) - t599) * t498 + (mrSges(6,2) * qJ(5) + t598) * t497 + (qJ(5) * t510 - t593) * t527 + (pkin(4) * t510 + t595) * t528 - pkin(1) * t453 - qJ(5) * t588 - pkin(4) * t583;
t447 = -mrSges(3,1) * t551 + mrSges(3,3) * t537 + Ifges(3,4) * t557 + Ifges(3,2) * t558 + Ifges(3,6) * qJDD(2) - pkin(2) * t580 + qJ(3) * t585 + qJD(2) * t550 + t572 * t454 + t571 * t455 - t548 * t591;
t446 = -mrSges(2,2) * g(3) - mrSges(2,3) * t562 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t578 - pkin(6) * t453 - t447 * t574 + t449 * t576;
t1 = [-m(1) * g(1) + t587; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t453; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t596 + t446 * t577 - t448 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t587 + t575 * t446 + t577 * t448; -mrSges(1,1) * g(2) + mrSges(2,1) * t562 + mrSges(1,2) * g(1) - mrSges(2,2) * t563 + Ifges(2,3) * qJDD(1) + pkin(1) * t579 + pkin(6) * t586 + t576 * t447 + t574 * t449;];
tauB = t1;
