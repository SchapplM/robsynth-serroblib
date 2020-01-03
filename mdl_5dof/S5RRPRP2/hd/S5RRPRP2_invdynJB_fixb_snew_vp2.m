% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP2
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:35
% EndTime: 2019-12-31 19:49:36
% DurationCPUTime: 1.77s
% Computational Cost: add. (22855->205), mult. (29217->249), div. (0->0), fcn. (14185->8), ass. (0->88)
t597 = Ifges(5,1) + Ifges(6,1);
t590 = Ifges(5,4) - Ifges(6,5);
t589 = Ifges(5,5) + Ifges(6,4);
t596 = Ifges(5,2) + Ifges(6,3);
t588 = Ifges(5,6) - Ifges(6,6);
t595 = Ifges(5,3) + Ifges(6,2);
t554 = qJD(1) + qJD(2);
t561 = sin(qJ(4));
t564 = cos(qJ(4));
t528 = (-mrSges(6,1) * t564 - mrSges(6,3) * t561) * t554;
t553 = qJDD(1) + qJDD(2);
t579 = qJD(4) * t554;
t530 = t553 * t561 + t564 * t579;
t563 = sin(qJ(1));
t566 = cos(qJ(1));
t544 = t563 * g(1) - g(2) * t566;
t537 = qJDD(1) * pkin(1) + t544;
t545 = -g(1) * t566 - g(2) * t563;
t568 = qJD(1) ^ 2;
t538 = -pkin(1) * t568 + t545;
t562 = sin(qJ(2));
t565 = cos(qJ(2));
t513 = t565 * t537 - t538 * t562;
t510 = pkin(2) * t553 + t513;
t514 = t562 * t537 + t565 * t538;
t552 = t554 ^ 2;
t511 = -pkin(2) * t552 + t514;
t559 = sin(pkin(8));
t560 = cos(pkin(8));
t506 = t559 * t510 + t560 * t511;
t503 = -pkin(3) * t552 + pkin(7) * t553 + t506;
t527 = (-pkin(4) * t564 - qJ(5) * t561) * t554;
t567 = qJD(4) ^ 2;
t558 = -g(3) + qJDD(3);
t584 = t558 * t564;
t498 = -qJDD(4) * pkin(4) - qJ(5) * t567 - t584 + qJDD(5) + (t527 * t554 + t503) * t561;
t585 = t554 * t564;
t542 = mrSges(6,2) * t585 + qJD(4) * mrSges(6,3);
t573 = -m(6) * t498 + qJDD(4) * mrSges(6,1) + qJD(4) * t542;
t586 = t554 * t561;
t494 = t530 * mrSges(6,2) + t528 * t586 - t573;
t500 = t564 * t503 + t561 * t558;
t497 = -pkin(4) * t567 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t527 * t585 + t500;
t499 = -t503 * t561 + t584;
t531 = t553 * t564 - t561 * t579;
t540 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t586;
t574 = m(6) * t497 + qJDD(4) * mrSges(6,3) + qJD(4) * t540 + t528 * t585;
t580 = (t597 * t561 + t590 * t564) * t554 + t589 * qJD(4);
t582 = (-t590 * t561 - t596 * t564) * t554 - t588 * qJD(4);
t594 = -(t582 * t561 + t580 * t564) * t554 + t595 * qJDD(4) + t589 * t530 + t588 * t531 + mrSges(5,1) * t499 - mrSges(6,1) * t498 - mrSges(5,2) * t500 + mrSges(6,3) * t497 - pkin(4) * t494 + qJ(5) * (t531 * mrSges(6,2) + t574);
t591 = mrSges(5,3) + mrSges(6,2);
t529 = (-mrSges(5,1) * t564 + mrSges(5,2) * t561) * t554;
t539 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t586;
t490 = m(5) * t500 - qJDD(4) * mrSges(5,2) - qJD(4) * t539 + t529 * t585 + t591 * t531 + t574;
t541 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t585;
t491 = m(5) * t499 + qJDD(4) * mrSges(5,1) + qJD(4) * t541 + (-t528 - t529) * t586 - t591 * t530 + t573;
t575 = t564 * t490 - t491 * t561;
t480 = m(4) * t506 - mrSges(4,1) * t552 - mrSges(4,2) * t553 + t575;
t505 = t510 * t560 - t559 * t511;
t502 = -pkin(3) * t553 - pkin(7) * t552 - t505;
t495 = -t531 * pkin(4) - t530 * qJ(5) + (-0.2e1 * qJD(5) * t561 + (pkin(4) * t561 - qJ(5) * t564) * qJD(4)) * t554 + t502;
t492 = m(6) * t495 - t531 * mrSges(6,1) - t530 * mrSges(6,3) - t540 * t586 - t542 * t585;
t570 = -m(5) * t502 + t531 * mrSges(5,1) - t530 * mrSges(5,2) - t539 * t586 + t541 * t585 - t492;
t485 = m(4) * t505 + mrSges(4,1) * t553 - mrSges(4,2) * t552 + t570;
t473 = t559 * t480 + t560 * t485;
t470 = m(3) * t513 + mrSges(3,1) * t553 - mrSges(3,2) * t552 + t473;
t576 = t560 * t480 - t485 * t559;
t471 = m(3) * t514 - mrSges(3,1) * t552 - mrSges(3,2) * t553 + t576;
t464 = t565 * t470 + t562 * t471;
t461 = m(2) * t544 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t568 + t464;
t577 = -t562 * t470 + t565 * t471;
t462 = m(2) * t545 - mrSges(2,1) * t568 - qJDD(1) * mrSges(2,2) + t577;
t583 = t566 * t461 + t563 * t462;
t483 = t561 * t490 + t564 * t491;
t581 = (t589 * t561 + t588 * t564) * t554 + t595 * qJD(4);
t481 = m(4) * t558 + t483;
t578 = -t461 * t563 + t566 * t462;
t476 = -mrSges(5,1) * t502 - mrSges(6,1) * t495 + mrSges(6,2) * t497 + mrSges(5,3) * t500 - pkin(4) * t492 + t580 * qJD(4) + t588 * qJDD(4) + t590 * t530 + t596 * t531 - t581 * t586;
t477 = mrSges(5,2) * t502 + mrSges(6,2) * t498 - mrSges(5,3) * t499 - mrSges(6,3) * t495 - qJ(5) * t492 + t582 * qJD(4) + t589 * qJDD(4) + t597 * t530 + t590 * t531 + t581 * t585;
t572 = mrSges(3,1) * t513 + mrSges(4,1) * t505 - mrSges(3,2) * t514 - mrSges(4,2) * t506 + pkin(2) * t473 + pkin(3) * t570 + pkin(7) * t575 + t564 * t476 + t561 * t477 + (Ifges(4,3) + Ifges(3,3)) * t553;
t569 = mrSges(2,1) * t544 - mrSges(2,2) * t545 + Ifges(2,3) * qJDD(1) + pkin(1) * t464 + t572;
t466 = -mrSges(4,1) * t558 + mrSges(4,3) * t506 + t552 * Ifges(4,5) + Ifges(4,6) * t553 - pkin(3) * t483 - t594;
t465 = mrSges(4,2) * t558 - mrSges(4,3) * t505 + Ifges(4,5) * t553 - Ifges(4,6) * t552 - pkin(7) * t483 - t476 * t561 + t477 * t564;
t457 = -mrSges(3,2) * g(3) - mrSges(3,3) * t513 + Ifges(3,5) * t553 - Ifges(3,6) * t552 - qJ(3) * t473 + t465 * t560 - t466 * t559;
t456 = mrSges(3,1) * g(3) + mrSges(3,3) * t514 + t552 * Ifges(3,5) + Ifges(3,6) * t553 - pkin(2) * t481 + qJ(3) * t576 + t559 * t465 + t560 * t466;
t455 = -mrSges(2,2) * g(3) - mrSges(2,3) * t544 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t568 - pkin(6) * t464 - t456 * t562 + t457 * t565;
t454 = Ifges(2,6) * qJDD(1) + t568 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t545 + t562 * t457 + t565 * t456 - pkin(1) * (-m(3) * g(3) + t481) + pkin(6) * t577;
t1 = [-m(1) * g(1) + t578; -m(1) * g(2) + t583; (-m(1) - m(2) - m(3)) * g(3) + t481; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t583 - t563 * t454 + t566 * t455; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t578 + t566 * t454 + t563 * t455; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t569; t569; t572; t481; t594; t494;];
tauJB = t1;
