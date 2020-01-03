% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:59
% EndTime: 2019-12-31 16:59:00
% DurationCPUTime: 0.82s
% Computational Cost: add. (3605->190), mult. (7384->232), div. (0->0), fcn. (2964->4), ass. (0->75)
t596 = Ifges(3,1) + Ifges(4,1) + Ifges(5,1);
t582 = Ifges(3,4) - Ifges(4,5) - Ifges(5,4);
t581 = Ifges(3,5) + Ifges(4,4) - Ifges(5,5);
t595 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t580 = Ifges(3,6) - Ifges(4,6) + Ifges(5,6);
t594 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t558 = sin(qJ(2));
t560 = cos(qJ(2));
t527 = (-mrSges(4,1) * t560 - mrSges(4,3) * t558) * qJD(1);
t583 = qJD(1) * qJD(2);
t574 = t560 * t583;
t530 = t558 * qJDD(1) + t574;
t559 = sin(qJ(1));
t561 = cos(qJ(1));
t545 = -t561 * g(1) - t559 * g(2);
t563 = qJD(1) ^ 2;
t514 = -t563 * pkin(1) + qJDD(1) * pkin(5) + t545;
t496 = -t560 * g(3) - t558 * t514;
t526 = (-pkin(2) * t560 - qJ(3) * t558) * qJD(1);
t562 = qJD(2) ^ 2;
t585 = qJD(1) * t558;
t495 = -qJDD(2) * pkin(2) - t562 * qJ(3) + t526 * t585 + qJDD(3) - t496;
t578 = -0.2e1 * qJD(1) * qJD(4);
t491 = t558 * t578 + (-t530 + t574) * qJ(4) + (-t558 * t560 * t563 - qJDD(2)) * pkin(3) + t495;
t528 = (mrSges(5,1) * t560 + mrSges(5,2) * t558) * qJD(1);
t584 = qJD(1) * t560;
t541 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t584;
t486 = m(5) * t491 - qJDD(2) * mrSges(5,1) - t530 * mrSges(5,3) - qJD(2) * t541 - t528 * t585;
t543 = mrSges(4,2) * t584 + qJD(2) * mrSges(4,3);
t566 = -m(4) * t495 + qJDD(2) * mrSges(4,1) + qJD(2) * t543 - t486;
t484 = t530 * mrSges(4,2) + t527 * t585 - t566;
t497 = -t558 * g(3) + t560 * t514;
t590 = 2 * qJD(3);
t494 = -t562 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t590 + t526 * t584 + t497;
t531 = t560 * qJDD(1) - t558 * t583;
t537 = -qJD(2) * pkin(3) - qJ(4) * t585;
t557 = t560 ^ 2;
t490 = -t557 * t563 * pkin(3) - t531 * qJ(4) + qJD(2) * t537 + t560 * t578 + t494;
t540 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t585;
t538 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t585;
t571 = m(5) * t490 + qJDD(2) * mrSges(5,2) - t531 * mrSges(5,3) + qJD(2) * t538;
t568 = m(4) * t494 + qJDD(2) * mrSges(4,3) + qJD(2) * t540 + t527 * t584 + t571;
t575 = t581 * qJD(2) + (t596 * t558 + t582 * t560) * qJD(1);
t576 = -t580 * qJD(2) + (-t582 * t558 - t595 * t560) * qJD(1);
t593 = -(t576 * t558 + t575 * t560) * qJD(1) + t594 * qJDD(2) + t581 * t530 + t580 * t531 + mrSges(3,1) * t496 - mrSges(4,1) * t495 - mrSges(5,1) * t491 - mrSges(3,2) * t497 + mrSges(5,2) * t490 + mrSges(4,3) * t494 - pkin(2) * t484 - pkin(3) * t486 + qJ(3) * (t531 * mrSges(4,2) - t528 * t584 + t568);
t589 = t563 * pkin(5);
t588 = mrSges(3,3) + mrSges(4,2);
t587 = qJ(3) * t560;
t529 = (-mrSges(3,1) * t560 + mrSges(3,2) * t558) * qJD(1);
t539 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t585;
t480 = m(3) * t497 - qJDD(2) * mrSges(3,2) - qJD(2) * t539 + t588 * t531 + (-t528 + t529) * t584 + t568;
t542 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t584;
t481 = m(3) * t496 + qJDD(2) * mrSges(3,1) + qJD(2) * t542 - t588 * t530 + (-t527 - t529) * t585 + t566;
t572 = t560 * t480 - t558 * t481;
t471 = m(2) * t545 - t563 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t572;
t544 = t559 * g(1) - t561 * g(2);
t570 = qJDD(1) * pkin(1) + t544;
t569 = -t530 * qJ(3) - t570;
t488 = qJDD(4) + (-qJ(4) * t557 + pkin(5)) * t563 + (pkin(2) + pkin(3)) * t531 + (qJD(2) * t587 + (-pkin(2) * qJD(2) + t537 + t590) * t558) * qJD(1) - t569;
t485 = m(5) * t488 + t531 * mrSges(5,1) + t530 * mrSges(5,2) + t538 * t585 + t541 * t584;
t492 = -t531 * pkin(2) - t589 + (-0.2e1 * qJD(3) * t558 + (pkin(2) * t558 - t587) * qJD(2)) * qJD(1) + t569;
t482 = m(4) * t492 - t531 * mrSges(4,1) - t530 * mrSges(4,3) - t540 * t585 - t543 * t584 - t485;
t513 = -t570 - t589;
t564 = -m(3) * t513 + t531 * mrSges(3,1) - t530 * mrSges(3,2) - t539 * t585 + t542 * t584 - t482;
t475 = m(2) * t544 + qJDD(1) * mrSges(2,1) - t563 * mrSges(2,2) + t564;
t586 = t559 * t471 + t561 * t475;
t473 = t558 * t480 + t560 * t481;
t577 = -t594 * qJD(2) + (-t581 * t558 - t580 * t560) * qJD(1);
t573 = t561 * t471 - t559 * t475;
t466 = -mrSges(3,1) * t513 + mrSges(3,3) * t497 - mrSges(4,1) * t492 + mrSges(4,2) * t494 + mrSges(5,1) * t488 - mrSges(5,3) * t490 + pkin(3) * t485 - qJ(4) * t571 - pkin(2) * t482 + t595 * t531 + t582 * t530 + t580 * qJDD(2) + t575 * qJD(2) + (qJ(4) * t528 * t560 + t577 * t558) * qJD(1);
t468 = mrSges(3,2) * t513 + mrSges(4,2) * t495 + mrSges(5,2) * t488 - mrSges(3,3) * t496 - mrSges(4,3) * t492 - mrSges(5,3) * t491 - qJ(3) * t482 - qJ(4) * t486 + t576 * qJD(2) + t581 * qJDD(2) + t596 * t530 + t582 * t531 - t577 * t584;
t567 = mrSges(2,1) * t544 - mrSges(2,2) * t545 + Ifges(2,3) * qJDD(1) + pkin(1) * t564 + pkin(5) * t572 + t560 * t466 + t558 * t468;
t464 = mrSges(2,1) * g(3) + mrSges(2,3) * t545 + t563 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t473 - t593;
t463 = -mrSges(2,2) * g(3) - mrSges(2,3) * t544 + Ifges(2,5) * qJDD(1) - t563 * Ifges(2,6) - pkin(5) * t473 - t558 * t466 + t560 * t468;
t1 = [-m(1) * g(1) + t573; -m(1) * g(2) + t586; (-m(1) - m(2)) * g(3) + t473; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t586 + t561 * t463 - t559 * t464; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t573 + t559 * t463 + t561 * t464; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t567; t567; t593; t484; t485;];
tauJB = t1;
