% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:35
% EndTime: 2019-12-31 17:04:37
% DurationCPUTime: 2.60s
% Computational Cost: add. (26540->240), mult. (61392->310), div. (0->0), fcn. (39600->8), ass. (0->96)
t558 = sin(qJ(2));
t561 = cos(qJ(2));
t575 = qJD(1) * qJD(2);
t541 = t558 * qJDD(1) + t561 * t575;
t559 = sin(qJ(1));
t562 = cos(qJ(1));
t548 = -t562 * g(1) - t559 * g(2);
t563 = qJD(1) ^ 2;
t536 = -t563 * pkin(1) + qJDD(1) * pkin(5) + t548;
t579 = t558 * t536;
t581 = pkin(2) * t563;
t506 = qJDD(2) * pkin(2) - t541 * qJ(3) - t579 + (qJ(3) * t575 + t558 * t581 - g(3)) * t561;
t522 = -t558 * g(3) + t561 * t536;
t542 = t561 * qJDD(1) - t558 * t575;
t577 = qJD(1) * t558;
t544 = qJD(2) * pkin(2) - qJ(3) * t577;
t554 = t561 ^ 2;
t507 = t542 * qJ(3) - qJD(2) * t544 - t554 * t581 + t522;
t555 = sin(pkin(7));
t556 = cos(pkin(7));
t531 = (t555 * t561 + t556 * t558) * qJD(1);
t488 = -0.2e1 * qJD(3) * t531 + t556 * t506 - t555 * t507;
t520 = t556 * t541 + t555 * t542;
t530 = (-t555 * t558 + t556 * t561) * qJD(1);
t484 = (qJD(2) * t530 - t520) * pkin(6) + (t530 * t531 + qJDD(2)) * pkin(3) + t488;
t489 = 0.2e1 * qJD(3) * t530 + t555 * t506 + t556 * t507;
t519 = -t555 * t541 + t556 * t542;
t525 = qJD(2) * pkin(3) - t531 * pkin(6);
t529 = t530 ^ 2;
t485 = -t529 * pkin(3) + t519 * pkin(6) - qJD(2) * t525 + t489;
t557 = sin(qJ(4));
t560 = cos(qJ(4));
t482 = t560 * t484 - t557 * t485;
t514 = t560 * t530 - t557 * t531;
t496 = t514 * qJD(4) + t557 * t519 + t560 * t520;
t515 = t557 * t530 + t560 * t531;
t502 = -t514 * mrSges(5,1) + t515 * mrSges(5,2);
t552 = qJD(2) + qJD(4);
t509 = -t552 * mrSges(5,2) + t514 * mrSges(5,3);
t551 = qJDD(2) + qJDD(4);
t478 = m(5) * t482 + t551 * mrSges(5,1) - t496 * mrSges(5,3) - t515 * t502 + t552 * t509;
t483 = t557 * t484 + t560 * t485;
t495 = -t515 * qJD(4) + t560 * t519 - t557 * t520;
t510 = t552 * mrSges(5,1) - t515 * mrSges(5,3);
t479 = m(5) * t483 - t551 * mrSges(5,2) + t495 * mrSges(5,3) + t514 * t502 - t552 * t510;
t469 = t560 * t478 + t557 * t479;
t517 = -t530 * mrSges(4,1) + t531 * mrSges(4,2);
t523 = -qJD(2) * mrSges(4,2) + t530 * mrSges(4,3);
t467 = m(4) * t488 + qJDD(2) * mrSges(4,1) - t520 * mrSges(4,3) + qJD(2) * t523 - t531 * t517 + t469;
t524 = qJD(2) * mrSges(4,1) - t531 * mrSges(4,3);
t571 = -t557 * t478 + t560 * t479;
t468 = m(4) * t489 - qJDD(2) * mrSges(4,2) + t519 * mrSges(4,3) - qJD(2) * t524 + t530 * t517 + t571;
t463 = t556 * t467 + t555 * t468;
t512 = Ifges(4,4) * t531 + Ifges(4,2) * t530 + Ifges(4,6) * qJD(2);
t513 = Ifges(4,1) * t531 + Ifges(4,4) * t530 + Ifges(4,5) * qJD(2);
t521 = -t561 * g(3) - t579;
t533 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t558 + Ifges(3,2) * t561) * qJD(1);
t534 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t558 + Ifges(3,4) * t561) * qJD(1);
t498 = Ifges(5,4) * t515 + Ifges(5,2) * t514 + Ifges(5,6) * t552;
t499 = Ifges(5,1) * t515 + Ifges(5,4) * t514 + Ifges(5,5) * t552;
t566 = -mrSges(5,1) * t482 + mrSges(5,2) * t483 - Ifges(5,5) * t496 - Ifges(5,6) * t495 - Ifges(5,3) * t551 - t515 * t498 + t514 * t499;
t582 = mrSges(3,1) * t521 + mrSges(4,1) * t488 - mrSges(3,2) * t522 - mrSges(4,2) * t489 + Ifges(3,5) * t541 + Ifges(4,5) * t520 + Ifges(3,6) * t542 + Ifges(4,6) * t519 + pkin(2) * t463 + pkin(3) * t469 + (t558 * t533 - t561 * t534) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t531 * t512 - t530 * t513 - t566;
t540 = (-mrSges(3,1) * t561 + mrSges(3,2) * t558) * qJD(1);
t576 = qJD(1) * t561;
t546 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t576;
t461 = m(3) * t521 + qJDD(2) * mrSges(3,1) - t541 * mrSges(3,3) + qJD(2) * t546 - t540 * t577 + t463;
t545 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t577;
t572 = -t555 * t467 + t556 * t468;
t462 = m(3) * t522 - qJDD(2) * mrSges(3,2) + t542 * mrSges(3,3) - qJD(2) * t545 + t540 * t576 + t572;
t573 = -t558 * t461 + t561 * t462;
t453 = m(2) * t548 - t563 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t573;
t547 = t559 * g(1) - t562 * g(2);
t568 = -qJDD(1) * pkin(1) - t547;
t508 = -t542 * pkin(2) + qJDD(3) + t544 * t577 + (-qJ(3) * t554 - pkin(5)) * t563 + t568;
t487 = -t519 * pkin(3) - t529 * pkin(6) + t531 * t525 + t508;
t570 = m(5) * t487 - t495 * mrSges(5,1) + t496 * mrSges(5,2) - t514 * t509 + t515 * t510;
t480 = m(4) * t508 - t519 * mrSges(4,1) + t520 * mrSges(4,2) - t530 * t523 + t531 * t524 + t570;
t535 = -t563 * pkin(5) + t568;
t565 = -m(3) * t535 + t542 * mrSges(3,1) - t541 * mrSges(3,2) - t545 * t577 + t546 * t576 - t480;
t473 = m(2) * t547 + qJDD(1) * mrSges(2,1) - t563 * mrSges(2,2) + t565;
t578 = t559 * t453 + t562 * t473;
t455 = t561 * t461 + t558 * t462;
t574 = t562 * t453 - t559 * t473;
t497 = Ifges(5,5) * t515 + Ifges(5,6) * t514 + Ifges(5,3) * t552;
t470 = -mrSges(5,1) * t487 + mrSges(5,3) * t483 + Ifges(5,4) * t496 + Ifges(5,2) * t495 + Ifges(5,6) * t551 - t515 * t497 + t552 * t499;
t471 = mrSges(5,2) * t487 - mrSges(5,3) * t482 + Ifges(5,1) * t496 + Ifges(5,4) * t495 + Ifges(5,5) * t551 + t514 * t497 - t552 * t498;
t511 = Ifges(4,5) * t531 + Ifges(4,6) * t530 + Ifges(4,3) * qJD(2);
t456 = -mrSges(4,1) * t508 + mrSges(4,3) * t489 + Ifges(4,4) * t520 + Ifges(4,2) * t519 + Ifges(4,6) * qJDD(2) - pkin(3) * t570 + pkin(6) * t571 + qJD(2) * t513 + t560 * t470 + t557 * t471 - t531 * t511;
t457 = mrSges(4,2) * t508 - mrSges(4,3) * t488 + Ifges(4,1) * t520 + Ifges(4,4) * t519 + Ifges(4,5) * qJDD(2) - pkin(6) * t469 - qJD(2) * t512 - t557 * t470 + t560 * t471 + t530 * t511;
t532 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t558 + Ifges(3,6) * t561) * qJD(1);
t448 = -mrSges(3,1) * t535 + mrSges(3,3) * t522 + Ifges(3,4) * t541 + Ifges(3,2) * t542 + Ifges(3,6) * qJDD(2) - pkin(2) * t480 + qJ(3) * t572 + qJD(2) * t534 + t556 * t456 + t555 * t457 - t532 * t577;
t450 = mrSges(3,2) * t535 - mrSges(3,3) * t521 + Ifges(3,1) * t541 + Ifges(3,4) * t542 + Ifges(3,5) * qJDD(2) - qJ(3) * t463 - qJD(2) * t533 - t555 * t456 + t556 * t457 + t532 * t576;
t567 = mrSges(2,1) * t547 - mrSges(2,2) * t548 + Ifges(2,3) * qJDD(1) + pkin(1) * t565 + pkin(5) * t573 + t561 * t448 + t558 * t450;
t446 = mrSges(2,1) * g(3) + mrSges(2,3) * t548 + t563 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t455 - t582;
t445 = -mrSges(2,2) * g(3) - mrSges(2,3) * t547 + Ifges(2,5) * qJDD(1) - t563 * Ifges(2,6) - pkin(5) * t455 - t558 * t448 + t561 * t450;
t1 = [-m(1) * g(1) + t574; -m(1) * g(2) + t578; (-m(1) - m(2)) * g(3) + t455; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t578 + t562 * t445 - t559 * t446; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t574 + t559 * t445 + t562 * t446; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t567; t567; t582; t480; -t566;];
tauJB = t1;
