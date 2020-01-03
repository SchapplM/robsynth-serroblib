% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:38
% EndTime: 2019-12-31 17:59:40
% DurationCPUTime: 1.46s
% Computational Cost: add. (14812->221), mult. (26577->267), div. (0->0), fcn. (13786->8), ass. (0->96)
t560 = sin(qJ(1));
t563 = cos(qJ(1));
t539 = t560 * g(1) - t563 * g(2);
t530 = qJDD(1) * pkin(1) + t539;
t540 = -t563 * g(1) - t560 * g(2);
t565 = qJD(1) ^ 2;
t532 = -t565 * pkin(1) + t540;
t556 = sin(pkin(8));
t557 = cos(pkin(8));
t510 = t556 * t530 + t557 * t532;
t591 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t510;
t590 = -pkin(2) - pkin(6);
t589 = mrSges(3,1) - mrSges(4,2);
t588 = -Ifges(4,4) + Ifges(3,5);
t587 = Ifges(4,5) - Ifges(3,6);
t553 = -g(3) + qJDD(2);
t559 = sin(qJ(4));
t586 = t559 * t553;
t509 = t557 * t530 - t556 * t532;
t573 = -t565 * qJ(3) + qJDD(3) - t509;
t497 = t590 * qJDD(1) + t573;
t562 = cos(qJ(4));
t493 = t559 * t497 + t562 * t553;
t531 = (mrSges(5,1) * t559 + mrSges(5,2) * t562) * qJD(1);
t582 = qJD(1) * qJD(4);
t578 = t562 * t582;
t534 = -t559 * qJDD(1) - t578;
t584 = qJD(1) * t562;
t538 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t584;
t496 = t590 * t565 - t591;
t579 = t559 * t582;
t535 = t562 * qJDD(1) - t579;
t488 = (-t535 + t579) * pkin(7) + (-t534 + t578) * pkin(4) + t496;
t533 = (pkin(4) * t559 - pkin(7) * t562) * qJD(1);
t564 = qJD(4) ^ 2;
t583 = t559 * qJD(1);
t490 = -t564 * pkin(4) + qJDD(4) * pkin(7) - t533 * t583 + t493;
t558 = sin(qJ(5));
t561 = cos(qJ(5));
t486 = t561 * t488 - t558 * t490;
t528 = t561 * qJD(4) - t558 * t584;
t507 = t528 * qJD(5) + t558 * qJDD(4) + t561 * t535;
t529 = t558 * qJD(4) + t561 * t584;
t511 = -t528 * mrSges(6,1) + t529 * mrSges(6,2);
t541 = qJD(5) + t583;
t512 = -t541 * mrSges(6,2) + t528 * mrSges(6,3);
t527 = qJDD(5) - t534;
t483 = m(6) * t486 + t527 * mrSges(6,1) - t507 * mrSges(6,3) - t529 * t511 + t541 * t512;
t487 = t558 * t488 + t561 * t490;
t506 = -t529 * qJD(5) + t561 * qJDD(4) - t558 * t535;
t513 = t541 * mrSges(6,1) - t529 * mrSges(6,3);
t484 = m(6) * t487 - t527 * mrSges(6,2) + t506 * mrSges(6,3) + t528 * t511 - t541 * t513;
t574 = -t558 * t483 + t561 * t484;
t472 = m(5) * t493 - qJDD(4) * mrSges(5,2) + t534 * mrSges(5,3) - qJD(4) * t538 - t531 * t583 + t574;
t492 = t562 * t497 - t586;
t537 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t583;
t489 = -qJDD(4) * pkin(4) - t564 * pkin(7) + t586 + (qJD(1) * t533 - t497) * t562;
t570 = -m(6) * t489 + t506 * mrSges(6,1) - t507 * mrSges(6,2) + t528 * t512 - t529 * t513;
t479 = m(5) * t492 + qJDD(4) * mrSges(5,1) - t535 * mrSges(5,3) + qJD(4) * t537 - t531 * t584 + t570;
t466 = t559 * t472 + t562 * t479;
t500 = -qJDD(1) * pkin(2) + t573;
t572 = -m(4) * t500 + t565 * mrSges(4,3) - t466;
t461 = m(3) * t509 - t565 * mrSges(3,2) + t589 * qJDD(1) + t572;
t498 = t565 * pkin(2) + t591;
t474 = t561 * t483 + t558 * t484;
t571 = -m(5) * t496 + t534 * mrSges(5,1) - t535 * mrSges(5,2) - t537 * t583 - t538 * t584 - t474;
t568 = -m(4) * t498 + t565 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t571;
t469 = m(3) * t510 - t565 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t568;
t456 = t557 * t461 + t556 * t469;
t453 = m(2) * t539 + qJDD(1) * mrSges(2,1) - t565 * mrSges(2,2) + t456;
t576 = -t556 * t461 + t557 * t469;
t454 = m(2) * t540 - t565 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t576;
t585 = t563 * t453 + t560 * t454;
t577 = -t560 * t453 + t563 * t454;
t575 = t562 * t472 - t559 * t479;
t465 = m(4) * t553 + t575;
t464 = m(3) * t553 + t465;
t501 = Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * t541;
t503 = Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t541;
t477 = -mrSges(6,1) * t489 + mrSges(6,3) * t487 + Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t527 - t529 * t501 + t541 * t503;
t502 = Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t541;
t478 = mrSges(6,2) * t489 - mrSges(6,3) * t486 + Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t527 + t528 * t501 - t541 * t502;
t520 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t562 - Ifges(5,2) * t559) * qJD(1);
t521 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t562 - Ifges(5,4) * t559) * qJD(1);
t569 = mrSges(5,1) * t492 - mrSges(5,2) * t493 + Ifges(5,5) * t535 + Ifges(5,6) * t534 + Ifges(5,3) * qJDD(4) + pkin(4) * t570 + pkin(7) * t574 + t561 * t477 + t558 * t478 + t520 * t584 + t521 * t583;
t567 = mrSges(6,1) * t486 - mrSges(6,2) * t487 + Ifges(6,5) * t507 + Ifges(6,6) * t506 + Ifges(6,3) * t527 + t529 * t502 - t528 * t503;
t519 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t562 - Ifges(5,6) * t559) * qJD(1);
t458 = mrSges(5,2) * t496 - mrSges(5,3) * t492 + Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * qJDD(4) - pkin(7) * t474 - qJD(4) * t520 - t558 * t477 + t561 * t478 - t519 * t583;
t459 = -mrSges(5,1) * t496 + mrSges(5,3) * t493 + Ifges(5,4) * t535 + Ifges(5,2) * t534 + Ifges(5,6) * qJDD(4) - pkin(4) * t474 + qJD(4) * t521 - t519 * t584 - t567;
t463 = qJDD(1) * mrSges(4,2) - t572;
t566 = mrSges(2,1) * t539 + mrSges(3,1) * t509 - mrSges(2,2) * t540 - mrSges(3,2) * t510 + mrSges(4,2) * t500 - mrSges(4,3) * t498 + pkin(1) * t456 - pkin(2) * t463 - pkin(6) * t466 + qJ(3) * t568 + t562 * t458 - t559 * t459 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1);
t449 = t588 * qJDD(1) + (mrSges(3,2) - mrSges(4,3)) * t553 + t587 * t565 - mrSges(3,3) * t509 + mrSges(4,1) * t500 + pkin(3) * t466 - qJ(3) * t465 + t569;
t448 = -mrSges(4,1) * t498 + mrSges(3,3) * t510 - pkin(2) * t465 - pkin(3) * t571 - pkin(6) * t575 - t587 * qJDD(1) - t559 * t458 - t562 * t459 - t589 * t553 + t588 * t565;
t447 = -mrSges(2,2) * g(3) - mrSges(2,3) * t539 + Ifges(2,5) * qJDD(1) - t565 * Ifges(2,6) - qJ(2) * t456 - t556 * t448 + t557 * t449;
t446 = mrSges(2,1) * g(3) + mrSges(2,3) * t540 + t565 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t464 + qJ(2) * t576 + t557 * t448 + t556 * t449;
t1 = [-m(1) * g(1) + t577; -m(1) * g(2) + t585; (-m(1) - m(2)) * g(3) + t464; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t585 - t560 * t446 + t563 * t447; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t577 + t563 * t446 + t560 * t447; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t566; t566; t464; t463; t569; t567;];
tauJB = t1;
