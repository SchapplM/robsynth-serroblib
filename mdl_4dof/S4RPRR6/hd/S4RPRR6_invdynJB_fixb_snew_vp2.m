% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR6
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:38
% EndTime: 2019-12-31 16:52:40
% DurationCPUTime: 2.24s
% Computational Cost: add. (22804->218), mult. (55057->276), div. (0->0), fcn. (37980->8), ass. (0->98)
t564 = qJD(1) ^ 2;
t557 = cos(pkin(7));
t590 = pkin(2) * t557;
t556 = sin(pkin(7));
t589 = mrSges(3,2) * t556;
t553 = t557 ^ 2;
t588 = t553 * t564;
t560 = sin(qJ(1));
t563 = cos(qJ(1));
t543 = -t563 * g(1) - t560 * g(2);
t538 = -t564 * pkin(1) + qJDD(1) * qJ(2) + t543;
t584 = qJD(1) * qJD(2);
t582 = -t557 * g(3) - 0.2e1 * t556 * t584;
t515 = (-pkin(5) * qJDD(1) + t564 * t590 - t538) * t556 + t582;
t529 = -t556 * g(3) + (t538 + 0.2e1 * t584) * t557;
t583 = qJDD(1) * t557;
t516 = -pkin(2) * t588 + pkin(5) * t583 + t529;
t559 = sin(qJ(3));
t562 = cos(qJ(3));
t502 = t562 * t515 - t559 * t516;
t572 = t556 * t562 + t557 * t559;
t571 = -t556 * t559 + t557 * t562;
t536 = t571 * qJD(1);
t585 = t536 * qJD(3);
t527 = t572 * qJDD(1) + t585;
t537 = t572 * qJD(1);
t491 = (-t527 + t585) * pkin(6) + (t536 * t537 + qJDD(3)) * pkin(3) + t502;
t503 = t559 * t515 + t562 * t516;
t526 = -t537 * qJD(3) + t571 * qJDD(1);
t532 = qJD(3) * pkin(3) - t537 * pkin(6);
t535 = t536 ^ 2;
t492 = -t535 * pkin(3) + t526 * pkin(6) - qJD(3) * t532 + t503;
t558 = sin(qJ(4));
t561 = cos(qJ(4));
t489 = t561 * t491 - t558 * t492;
t520 = t561 * t536 - t558 * t537;
t501 = t520 * qJD(4) + t558 * t526 + t561 * t527;
t521 = t558 * t536 + t561 * t537;
t509 = -t520 * mrSges(5,1) + t521 * mrSges(5,2);
t554 = qJD(3) + qJD(4);
t512 = -t554 * mrSges(5,2) + t520 * mrSges(5,3);
t551 = qJDD(3) + qJDD(4);
t486 = m(5) * t489 + t551 * mrSges(5,1) - t501 * mrSges(5,3) - t521 * t509 + t554 * t512;
t490 = t558 * t491 + t561 * t492;
t500 = -t521 * qJD(4) + t561 * t526 - t558 * t527;
t513 = t554 * mrSges(5,1) - t521 * mrSges(5,3);
t487 = m(5) * t490 - t551 * mrSges(5,2) + t500 * mrSges(5,3) + t520 * t509 - t554 * t513;
t476 = t561 * t486 + t558 * t487;
t523 = -t536 * mrSges(4,1) + t537 * mrSges(4,2);
t530 = -qJD(3) * mrSges(4,2) + t536 * mrSges(4,3);
t474 = m(4) * t502 + qJDD(3) * mrSges(4,1) - t527 * mrSges(4,3) + qJD(3) * t530 - t537 * t523 + t476;
t531 = qJD(3) * mrSges(4,1) - t537 * mrSges(4,3);
t578 = -t558 * t486 + t561 * t487;
t475 = m(4) * t503 - qJDD(3) * mrSges(4,2) + t526 * mrSges(4,3) - qJD(3) * t531 + t536 * t523 + t578;
t470 = t562 * t474 + t559 * t475;
t528 = -t556 * t538 + t582;
t570 = mrSges(3,3) * qJDD(1) + t564 * (-mrSges(3,1) * t557 + t589);
t468 = m(3) * t528 - t570 * t556 + t470;
t579 = -t559 * t474 + t562 * t475;
t469 = m(3) * t529 + t570 * t557 + t579;
t580 = -t556 * t468 + t557 * t469;
t460 = m(2) * t543 - t564 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t580;
t542 = t560 * g(1) - t563 * g(2);
t577 = qJDD(2) - t542;
t534 = -qJDD(1) * pkin(1) - t564 * qJ(2) + t577;
t552 = t556 ^ 2;
t525 = (-pkin(1) - t590) * qJDD(1) + (-qJ(2) + (-t552 - t553) * pkin(5)) * t564 + t577;
t495 = -t526 * pkin(3) - t535 * pkin(6) + t537 * t532 + t525;
t573 = m(5) * t495 - t500 * mrSges(5,1) + t501 * mrSges(5,2) - t520 * t512 + t521 * t513;
t567 = m(4) * t525 - t526 * mrSges(4,1) + t527 * mrSges(4,2) - t536 * t530 + t537 * t531 + t573;
t566 = -m(3) * t534 + mrSges(3,1) * t583 - t567 + (t552 * t564 + t588) * mrSges(3,3);
t480 = -t564 * mrSges(2,2) + m(2) * t542 + t566 + (mrSges(2,1) - t589) * qJDD(1);
t587 = t560 * t460 + t563 * t480;
t462 = t557 * t468 + t556 * t469;
t574 = Ifges(3,5) * t556 + Ifges(3,6) * t557;
t586 = t564 * t574;
t581 = t563 * t460 - t560 * t480;
t576 = Ifges(3,1) * t556 + Ifges(3,4) * t557;
t575 = Ifges(3,4) * t556 + Ifges(3,2) * t557;
t504 = Ifges(5,5) * t521 + Ifges(5,6) * t520 + Ifges(5,3) * t554;
t506 = Ifges(5,1) * t521 + Ifges(5,4) * t520 + Ifges(5,5) * t554;
t477 = -mrSges(5,1) * t495 + mrSges(5,3) * t490 + Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t551 - t521 * t504 + t554 * t506;
t505 = Ifges(5,4) * t521 + Ifges(5,2) * t520 + Ifges(5,6) * t554;
t478 = mrSges(5,2) * t495 - mrSges(5,3) * t489 + Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t551 + t520 * t504 - t554 * t505;
t517 = Ifges(4,5) * t537 + Ifges(4,6) * t536 + Ifges(4,3) * qJD(3);
t519 = Ifges(4,1) * t537 + Ifges(4,4) * t536 + Ifges(4,5) * qJD(3);
t463 = -mrSges(4,1) * t525 + mrSges(4,3) * t503 + Ifges(4,4) * t527 + Ifges(4,2) * t526 + Ifges(4,6) * qJDD(3) - pkin(3) * t573 + pkin(6) * t578 + qJD(3) * t519 + t561 * t477 + t558 * t478 - t537 * t517;
t518 = Ifges(4,4) * t537 + Ifges(4,2) * t536 + Ifges(4,6) * qJD(3);
t464 = mrSges(4,2) * t525 - mrSges(4,3) * t502 + Ifges(4,1) * t527 + Ifges(4,4) * t526 + Ifges(4,5) * qJDD(3) - pkin(6) * t476 - qJD(3) * t518 - t558 * t477 + t561 * t478 + t536 * t517;
t455 = -mrSges(3,1) * t534 + mrSges(3,3) * t529 - pkin(2) * t567 + pkin(5) * t579 + t575 * qJDD(1) + t562 * t463 + t559 * t464 - t556 * t586;
t457 = mrSges(3,2) * t534 - mrSges(3,3) * t528 - pkin(5) * t470 + t576 * qJDD(1) - t559 * t463 + t562 * t464 + t557 * t586;
t482 = qJDD(1) * t589 - t566;
t569 = mrSges(2,1) * t542 - mrSges(2,2) * t543 + Ifges(2,3) * qJDD(1) - pkin(1) * t482 + qJ(2) * t580 + t557 * t455 + t556 * t457;
t568 = -mrSges(5,1) * t489 + mrSges(5,2) * t490 - Ifges(5,5) * t501 - Ifges(5,6) * t500 - Ifges(5,3) * t551 - t521 * t505 + t520 * t506;
t565 = mrSges(4,1) * t502 - mrSges(4,2) * t503 + Ifges(4,5) * t527 + Ifges(4,6) * t526 + Ifges(4,3) * qJDD(3) + pkin(3) * t476 + t537 * t518 - t536 * t519 - t568;
t453 = -mrSges(3,1) * t528 + mrSges(3,2) * t529 - pkin(1) * t462 + (Ifges(2,6) - t574) * qJDD(1) + mrSges(2,3) * t543 - pkin(2) * t470 + mrSges(2,1) * g(3) - t565 + (-t556 * t575 + t557 * t576 + Ifges(2,5)) * t564;
t452 = -mrSges(2,2) * g(3) - mrSges(2,3) * t542 + Ifges(2,5) * qJDD(1) - t564 * Ifges(2,6) - qJ(2) * t462 - t556 * t455 + t557 * t457;
t1 = [-m(1) * g(1) + t581; -m(1) * g(2) + t587; (-m(1) - m(2)) * g(3) + t462; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t587 + t563 * t452 - t560 * t453; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t581 + t560 * t452 + t563 * t453; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t569; t569; t482; t565; -t568;];
tauJB = t1;
