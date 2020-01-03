% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:17
% EndTime: 2019-12-31 21:49:18
% DurationCPUTime: 1.81s
% Computational Cost: add. (28623->207), mult. (29868->251), div. (0->0), fcn. (14500->8), ass. (0->91)
t599 = Ifges(5,1) + Ifges(6,1);
t590 = Ifges(5,4) - Ifges(6,5);
t589 = Ifges(5,5) + Ifges(6,4);
t598 = Ifges(5,2) + Ifges(6,3);
t588 = Ifges(5,6) - Ifges(6,6);
t597 = Ifges(5,3) + Ifges(6,2);
t555 = qJD(1) + qJD(2);
t551 = qJD(3) + t555;
t559 = sin(qJ(4));
t563 = cos(qJ(4));
t528 = (-mrSges(6,1) * t563 - mrSges(6,3) * t559) * t551;
t554 = qJDD(1) + qJDD(2);
t550 = qJDD(3) + t554;
t580 = qJD(4) * t551;
t530 = t559 * t550 + t563 * t580;
t562 = sin(qJ(1));
t566 = cos(qJ(1));
t545 = t562 * g(1) - t566 * g(2);
t541 = qJDD(1) * pkin(1) + t545;
t546 = -t566 * g(1) - t562 * g(2);
t568 = qJD(1) ^ 2;
t542 = -t568 * pkin(1) + t546;
t561 = sin(qJ(2));
t565 = cos(qJ(2));
t513 = t565 * t541 - t561 * t542;
t510 = t554 * pkin(2) + t513;
t514 = t561 * t541 + t565 * t542;
t553 = t555 ^ 2;
t511 = -t553 * pkin(2) + t514;
t560 = sin(qJ(3));
t564 = cos(qJ(3));
t506 = t560 * t510 + t564 * t511;
t549 = t551 ^ 2;
t503 = -t549 * pkin(3) + t550 * pkin(8) + t506;
t527 = (-pkin(4) * t563 - qJ(5) * t559) * t551;
t567 = qJD(4) ^ 2;
t592 = t563 * g(3);
t498 = -qJDD(4) * pkin(4) + t592 - t567 * qJ(5) + qJDD(5) + (t527 * t551 + t503) * t559;
t585 = t551 * t563;
t540 = mrSges(6,2) * t585 + qJD(4) * mrSges(6,3);
t574 = -m(6) * t498 + qJDD(4) * mrSges(6,1) + qJD(4) * t540;
t586 = t551 * t559;
t494 = t530 * mrSges(6,2) + t528 * t586 - t574;
t500 = -t559 * g(3) + t563 * t503;
t497 = -t567 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t527 * t585 + t500;
t499 = -t559 * t503 - t592;
t531 = t563 * t550 - t559 * t580;
t538 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t586;
t575 = m(6) * t497 + qJDD(4) * mrSges(6,3) + qJD(4) * t538 + t528 * t585;
t581 = (t599 * t559 + t590 * t563) * t551 + t589 * qJD(4);
t583 = (-t590 * t559 - t598 * t563) * t551 - t588 * qJD(4);
t596 = -(t583 * t559 + t581 * t563) * t551 + t597 * qJDD(4) + t589 * t530 + t588 * t531 + mrSges(5,1) * t499 - mrSges(6,1) * t498 - mrSges(5,2) * t500 + mrSges(6,3) * t497 - pkin(4) * t494 + qJ(5) * (t531 * mrSges(6,2) + t575);
t593 = -m(3) - m(4);
t591 = mrSges(5,3) + mrSges(6,2);
t529 = (-mrSges(5,1) * t563 + mrSges(5,2) * t559) * t551;
t537 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t586;
t490 = m(5) * t500 - qJDD(4) * mrSges(5,2) - qJD(4) * t537 + t529 * t585 + t591 * t531 + t575;
t539 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t585;
t491 = m(5) * t499 + qJDD(4) * mrSges(5,1) + qJD(4) * t539 + (-t528 - t529) * t586 - t591 * t530 + t574;
t576 = t563 * t490 - t559 * t491;
t481 = m(4) * t506 - t549 * mrSges(4,1) - t550 * mrSges(4,2) + t576;
t505 = t564 * t510 - t560 * t511;
t502 = -t550 * pkin(3) - t549 * pkin(8) - t505;
t495 = -t531 * pkin(4) - t530 * qJ(5) + (-0.2e1 * qJD(5) * t559 + (pkin(4) * t559 - qJ(5) * t563) * qJD(4)) * t551 + t502;
t492 = m(6) * t495 - t531 * mrSges(6,1) - t530 * mrSges(6,3) - t538 * t586 - t540 * t585;
t570 = -m(5) * t502 + t531 * mrSges(5,1) - t530 * mrSges(5,2) - t537 * t586 + t539 * t585 - t492;
t485 = m(4) * t505 + t550 * mrSges(4,1) - t549 * mrSges(4,2) + t570;
t474 = t560 * t481 + t564 * t485;
t471 = m(3) * t513 + t554 * mrSges(3,1) - t553 * mrSges(3,2) + t474;
t577 = t564 * t481 - t560 * t485;
t472 = m(3) * t514 - t553 * mrSges(3,1) - t554 * mrSges(3,2) + t577;
t465 = t565 * t471 + t561 * t472;
t462 = m(2) * t545 + qJDD(1) * mrSges(2,1) - t568 * mrSges(2,2) + t465;
t578 = -t561 * t471 + t565 * t472;
t463 = m(2) * t546 - t568 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t578;
t584 = t566 * t462 + t562 * t463;
t483 = t559 * t490 + t563 * t491;
t582 = (t589 * t559 + t588 * t563) * t551 + t597 * qJD(4);
t579 = -t562 * t462 + t566 * t463;
t477 = -mrSges(5,1) * t502 - mrSges(6,1) * t495 + mrSges(6,2) * t497 + mrSges(5,3) * t500 - pkin(4) * t492 + t581 * qJD(4) + t588 * qJDD(4) + t590 * t530 + t598 * t531 - t582 * t586;
t478 = mrSges(5,2) * t502 + mrSges(6,2) * t498 - mrSges(5,3) * t499 - mrSges(6,3) * t495 - qJ(5) * t492 + t583 * qJD(4) + t589 * qJDD(4) + t599 * t530 + t590 * t531 + t582 * t585;
t573 = mrSges(4,1) * t505 - mrSges(4,2) * t506 + Ifges(4,3) * t550 + pkin(3) * t570 + pkin(8) * t576 + t563 * t477 + t559 * t478;
t572 = mrSges(3,1) * t513 - mrSges(3,2) * t514 + Ifges(3,3) * t554 + pkin(2) * t474 + t573;
t569 = mrSges(2,1) * t545 - mrSges(2,2) * t546 + Ifges(2,3) * qJDD(1) + pkin(1) * t465 + t572;
t467 = mrSges(4,1) * g(3) + mrSges(4,3) * t506 + t549 * Ifges(4,5) + Ifges(4,6) * t550 - pkin(3) * t483 - t596;
t466 = -mrSges(4,2) * g(3) - mrSges(4,3) * t505 + Ifges(4,5) * t550 - t549 * Ifges(4,6) - pkin(8) * t483 - t559 * t477 + t563 * t478;
t458 = -mrSges(3,2) * g(3) - mrSges(3,3) * t513 + Ifges(3,5) * t554 - t553 * Ifges(3,6) - pkin(7) * t474 + t564 * t466 - t560 * t467;
t457 = Ifges(3,6) * t554 + t553 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t514 + t560 * t466 + t564 * t467 - pkin(2) * (-m(4) * g(3) + t483) + pkin(7) * t577;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t545 + Ifges(2,5) * qJDD(1) - t568 * Ifges(2,6) - pkin(6) * t465 - t561 * t457 + t565 * t458;
t455 = Ifges(2,6) * qJDD(1) + t568 * Ifges(2,5) + mrSges(2,3) * t546 + t561 * t458 + t565 * t457 - pkin(1) * t483 + pkin(6) * t578 + (-pkin(1) * t593 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t579; -m(1) * g(2) + t584; (-m(1) - m(2) + t593) * g(3) + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t584 - t562 * t455 + t566 * t456; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t579 + t566 * t455 + t562 * t456; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t569; t569; t572; t573; t596; t494;];
tauJB = t1;
