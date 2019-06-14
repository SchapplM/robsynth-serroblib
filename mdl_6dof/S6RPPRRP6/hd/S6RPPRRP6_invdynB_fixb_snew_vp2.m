% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:00:20
% EndTime: 2019-05-05 15:00:23
% DurationCPUTime: 1.82s
% Computational Cost: add. (13684->265), mult. (25066->300), div. (0->0), fcn. (12633->6), ass. (0->98)
t600 = Ifges(6,1) + Ifges(7,1);
t591 = Ifges(6,4) - Ifges(7,5);
t599 = -Ifges(6,5) - Ifges(7,4);
t598 = Ifges(6,2) + Ifges(7,3);
t589 = Ifges(6,6) - Ifges(7,6);
t597 = -Ifges(6,3) - Ifges(7,2);
t558 = sin(qJ(1));
t560 = cos(qJ(1));
t536 = t558 * g(1) - t560 * g(2);
t562 = qJD(1) ^ 2;
t516 = -qJDD(1) * pkin(1) - t562 * qJ(2) + qJDD(2) - t536;
t509 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t516;
t537 = -t560 * g(1) - t558 * g(2);
t596 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t537;
t595 = -m(3) - m(4);
t594 = cos(qJ(5));
t593 = mrSges(2,1) - mrSges(3,2);
t592 = -mrSges(6,3) - mrSges(7,2);
t515 = t562 * pkin(1) - t596;
t510 = qJDD(3) + (-pkin(1) - qJ(3)) * t562 + t596;
t507 = -qJDD(1) * pkin(7) + t510;
t557 = sin(qJ(4));
t559 = cos(qJ(4));
t498 = -t559 * g(3) + t557 * t507;
t530 = (mrSges(5,1) * t557 + mrSges(5,2) * t559) * qJD(1);
t581 = qJD(1) * qJD(4);
t574 = t559 * t581;
t532 = -t557 * qJDD(1) - t574;
t583 = qJD(1) * t559;
t535 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t583;
t506 = -t562 * pkin(7) - t509;
t575 = t557 * t581;
t533 = t559 * qJDD(1) - t575;
t481 = (-t533 + t575) * pkin(8) + (-t532 + t574) * pkin(4) + t506;
t531 = (pkin(4) * t557 - pkin(8) * t559) * qJD(1);
t561 = qJD(4) ^ 2;
t582 = t557 * qJD(1);
t484 = -t561 * pkin(4) + qJDD(4) * pkin(8) - t531 * t582 + t498;
t556 = sin(qJ(5));
t479 = t556 * t481 + t594 * t484;
t529 = t556 * qJD(4) + t594 * t583;
t495 = t529 * qJD(5) - t594 * qJDD(4) + t556 * t533;
t539 = qJD(5) + t582;
t512 = t539 * mrSges(6,1) - t529 * mrSges(6,3);
t527 = qJDD(5) - t532;
t528 = -t594 * qJD(4) + t556 * t583;
t502 = t528 * pkin(5) - t529 * qJ(6);
t538 = t539 ^ 2;
t475 = -t538 * pkin(5) + t527 * qJ(6) + 0.2e1 * qJD(6) * t539 - t528 * t502 + t479;
t513 = -t539 * mrSges(7,1) + t529 * mrSges(7,2);
t576 = m(7) * t475 + t527 * mrSges(7,3) + t539 * t513;
t503 = t528 * mrSges(7,1) - t529 * mrSges(7,3);
t584 = -t528 * mrSges(6,1) - t529 * mrSges(6,2) - t503;
t470 = m(6) * t479 - t527 * mrSges(6,2) + t592 * t495 - t539 * t512 + t584 * t528 + t576;
t478 = t594 * t481 - t556 * t484;
t496 = -t528 * qJD(5) + t556 * qJDD(4) + t594 * t533;
t511 = -t539 * mrSges(6,2) - t528 * mrSges(6,3);
t476 = -t527 * pkin(5) - t538 * qJ(6) + t529 * t502 + qJDD(6) - t478;
t514 = -t528 * mrSges(7,2) + t539 * mrSges(7,3);
t569 = -m(7) * t476 + t527 * mrSges(7,1) + t539 * t514;
t472 = m(6) * t478 + t527 * mrSges(6,1) + t592 * t496 + t539 * t511 + t584 * t529 + t569;
t571 = t594 * t470 - t556 * t472;
t464 = m(5) * t498 - qJDD(4) * mrSges(5,2) + t532 * mrSges(5,3) - qJD(4) * t535 - t530 * t582 + t571;
t497 = t557 * g(3) + t559 * t507;
t534 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t582;
t483 = -qJDD(4) * pkin(4) - t561 * pkin(8) + t531 * t583 - t497;
t477 = -0.2e1 * qJD(6) * t529 + (t528 * t539 - t496) * qJ(6) + (t529 * t539 + t495) * pkin(5) + t483;
t473 = m(7) * t477 + t495 * mrSges(7,1) - t496 * mrSges(7,3) - t529 * t513 + t528 * t514;
t563 = -m(6) * t483 - t495 * mrSges(6,1) - t496 * mrSges(6,2) - t528 * t511 - t529 * t512 - t473;
t467 = m(5) * t497 + qJDD(4) * mrSges(5,1) - t533 * mrSges(5,3) + qJD(4) * t534 - t530 * t583 + t563;
t456 = t557 * t464 + t559 * t467;
t570 = -m(4) * t510 - qJDD(1) * mrSges(4,2) - t456;
t565 = -m(3) * t515 + t562 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t570;
t454 = m(2) * t537 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t562 + t565;
t465 = t556 * t470 + t594 * t472;
t566 = -m(5) * t506 + t532 * mrSges(5,1) - t533 * mrSges(5,2) - t534 * t582 - t535 * t583 - t465;
t459 = m(4) * t509 - t562 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t566;
t564 = -m(3) * t516 + t562 * mrSges(3,3) - t459;
t458 = m(2) * t536 - t562 * mrSges(2,2) + t593 * qJDD(1) + t564;
t588 = t558 * t454 + t560 * t458;
t587 = t598 * t528 - t591 * t529 - t589 * t539;
t586 = t589 * t528 + t599 * t529 + t597 * t539;
t585 = -t591 * t528 + t600 * t529 - t599 * t539;
t578 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t577 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t573 = t560 * t454 - t558 * t458;
t572 = t559 * t464 - t557 * t467;
t520 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t559 - Ifges(5,4) * t557) * qJD(1);
t519 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t559 - Ifges(5,2) * t557) * qJD(1);
t518 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t559 - Ifges(5,6) * t557) * qJD(1);
t461 = mrSges(6,2) * t483 + mrSges(7,2) * t476 - mrSges(6,3) * t478 - mrSges(7,3) * t477 - qJ(6) * t473 - t591 * t495 + t600 * t496 - t527 * t599 + t586 * t528 + t587 * t539;
t460 = -mrSges(6,1) * t483 - mrSges(7,1) * t477 + mrSges(7,2) * t475 + mrSges(6,3) * t479 - pkin(5) * t473 - t598 * t495 + t591 * t496 + t589 * t527 + t586 * t529 + t585 * t539;
t455 = t595 * g(3) + t572;
t451 = Ifges(5,4) * t533 + Ifges(5,2) * t532 + Ifges(5,6) * qJDD(4) - t518 * t583 + qJD(4) * t520 - mrSges(5,1) * t506 + mrSges(5,3) * t498 - mrSges(6,1) * t478 + mrSges(6,2) * t479 + mrSges(7,1) * t476 - mrSges(7,3) * t475 - pkin(5) * t569 - qJ(6) * t576 - pkin(4) * t465 + (pkin(5) * t503 + t587) * t529 + (qJ(6) * t503 - t585) * t528 + t597 * t527 + (pkin(5) * mrSges(7,2) + t599) * t496 + (qJ(6) * mrSges(7,2) + t589) * t495;
t450 = mrSges(5,2) * t506 - mrSges(5,3) * t497 + Ifges(5,1) * t533 + Ifges(5,4) * t532 + Ifges(5,5) * qJDD(4) - pkin(8) * t465 - qJD(4) * t519 - t556 * t460 + t594 * t461 - t518 * t582;
t449 = -pkin(1) * t455 + pkin(3) * t456 + Ifges(5,3) * qJDD(4) + t594 * t460 + mrSges(5,1) * t497 - mrSges(5,2) * t498 + mrSges(4,1) * t510 - mrSges(3,1) * t515 + pkin(4) * t563 + Ifges(5,6) * t532 + Ifges(5,5) * t533 + mrSges(2,3) * t537 + t556 * t461 + pkin(8) * t571 - qJ(3) * t572 - pkin(2) * t570 + (t559 * t519 + t557 * t520) * qJD(1) + (-pkin(2) * mrSges(4,3) + t578) * t562 + t577 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t593) * g(3);
t448 = -qJ(2) * t455 - mrSges(2,3) * t536 + pkin(2) * t459 + mrSges(3,1) * t516 + t559 * t451 + pkin(3) * t566 + pkin(7) * t572 + t557 * t450 + mrSges(4,1) * t509 - t577 * t562 + t578 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t573; -m(1) * g(2) + t588; (-m(1) - m(2) + t595) * g(3) + t572; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t588 + t560 * t448 - t558 * t449; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t573 + t558 * t448 + t560 * t449; qJ(2) * (-t562 * mrSges(4,3) + t565) + pkin(1) * t564 + mrSges(2,1) * t536 - mrSges(2,2) * t537 - qJ(3) * t459 + mrSges(3,2) * t516 - mrSges(3,3) * t515 - t557 * t451 - pkin(7) * t456 + t559 * t450 + mrSges(4,2) * t510 - mrSges(4,3) * t509 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
