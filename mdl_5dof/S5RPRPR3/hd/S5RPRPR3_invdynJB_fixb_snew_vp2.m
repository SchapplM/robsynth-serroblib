% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:47
% EndTime: 2022-01-23 09:20:50
% DurationCPUTime: 2.73s
% Computational Cost: add. (36512->207), mult. (51081->271), div. (0->0), fcn. (28386->10), ass. (0->104)
t601 = 2 * qJD(4);
t565 = sin(qJ(1));
t568 = cos(qJ(1));
t543 = t565 * g(1) - g(2) * t568;
t537 = qJDD(1) * pkin(1) + t543;
t544 = -g(1) * t568 - g(2) * t565;
t569 = qJD(1) ^ 2;
t538 = -pkin(1) * t569 + t544;
t560 = sin(pkin(8));
t562 = cos(pkin(8));
t523 = t562 * t537 - t538 * t560;
t517 = qJDD(1) * pkin(2) + t523;
t524 = t560 * t537 + t562 * t538;
t518 = -pkin(2) * t569 + t524;
t564 = sin(qJ(3));
t567 = cos(qJ(3));
t513 = t564 * t517 + t567 * t518;
t555 = (qJD(1) + qJD(3));
t553 = t555 ^ 2;
t554 = qJDD(1) + qJDD(3);
t510 = -pkin(3) * t553 + qJ(4) * t554 + t513;
t600 = (t555 * t601) + t510;
t559 = sin(pkin(9));
t599 = mrSges(5,2) * t559;
t597 = mrSges(5,3) * t554;
t596 = t559 * t555;
t563 = sin(qJ(5));
t595 = t559 * t563;
t566 = cos(qJ(5));
t594 = t559 * t566;
t561 = cos(pkin(9));
t593 = t561 * t554;
t592 = t561 * t555;
t558 = -g(3) + qJDD(2);
t591 = t561 * t558;
t506 = t559 * t558 + t561 * t600;
t531 = (-mrSges(5,1) * t561 + t599) * t555;
t579 = -pkin(4) * t561 - pkin(7) * t559;
t533 = t579 * t555;
t504 = t533 * t592 + t506;
t512 = t567 * t517 - t564 * t518;
t574 = -t553 * qJ(4) + qJDD(4) - t512;
t507 = (-pkin(3) + t579) * t554 + t574;
t501 = -t504 * t563 + t507 * t566;
t541 = qJD(5) - t592;
t587 = t555 * t595;
t526 = -mrSges(6,2) * t541 - mrSges(6,3) * t587;
t528 = (mrSges(6,1) * t563 + mrSges(6,2) * t566) * t596;
t588 = qJD(5) * t555;
t530 = (t554 * t566 - t563 * t588) * t559;
t540 = qJDD(5) - t593;
t586 = t555 * t594;
t499 = m(6) * t501 + mrSges(6,1) * t540 - mrSges(6,3) * t530 + t526 * t541 - t528 * t586;
t502 = t504 * t566 + t507 * t563;
t527 = mrSges(6,1) * t541 - mrSges(6,3) * t586;
t529 = (-t554 * t563 - t566 * t588) * t559;
t500 = m(6) * t502 - mrSges(6,2) * t540 + mrSges(6,3) * t529 - t527 * t541 - t528 * t587;
t580 = -t563 * t499 + t500 * t566;
t490 = m(5) * t506 + (t531 * t555 + t597) * t561 + t580;
t505 = -t559 * t600 + t591;
t503 = -t591 + (t510 + (t601 + t533) * t555) * t559;
t575 = -m(6) * t503 + t529 * mrSges(6,1) - t530 * mrSges(6,2);
t497 = m(5) * t505 + (-t597 + (-t526 * t563 - t527 * t566 - t531) * t555) * t559 + t575;
t581 = t490 * t561 - t497 * t559;
t482 = m(4) * t513 - mrSges(4,1) * t553 - mrSges(4,2) * t554 + t581;
t493 = t566 * t499 + t563 * t500;
t509 = -pkin(3) * t554 + t574;
t572 = -m(5) * t509 + mrSges(5,1) * t593 - t493 + (t559 ^ 2 + t561 ^ 2) * mrSges(5,3) * t553;
t487 = m(4) * t512 - t553 * mrSges(4,2) + (mrSges(4,1) - t599) * t554 + t572;
t475 = t482 * t564 + t487 * t567;
t472 = m(3) * t523 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t569 + t475;
t582 = t482 * t567 - t487 * t564;
t473 = m(3) * t524 - mrSges(3,1) * t569 - qJDD(1) * mrSges(3,2) + t582;
t467 = t472 * t562 + t473 * t560;
t464 = m(2) * t543 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t569 + t467;
t583 = -t472 * t560 + t473 * t562;
t465 = m(2) * t544 - mrSges(2,1) * t569 - qJDD(1) * mrSges(2,2) + t583;
t590 = t464 * t568 + t465 * t565;
t485 = t490 * t559 + t497 * t561;
t585 = m(4) * t558 + t485;
t584 = -t464 * t565 + t465 * t568;
t483 = m(3) * t558 + t585;
t578 = Ifges(5,1) * t559 + Ifges(5,4) * t561;
t577 = Ifges(5,5) * t559 + Ifges(5,6) * t561;
t521 = Ifges(6,6) * t541 + (Ifges(6,4) * t566 - Ifges(6,2) * t563) * t596;
t522 = Ifges(6,5) * t541 + (Ifges(6,1) * t566 - Ifges(6,4) * t563) * t596;
t576 = t521 * t566 + t522 * t563;
t520 = Ifges(6,3) * t541 + (Ifges(6,5) * t566 - Ifges(6,6) * t563) * t596;
t494 = -mrSges(6,1) * t503 + mrSges(6,3) * t502 + Ifges(6,4) * t530 + Ifges(6,2) * t529 + Ifges(6,6) * t540 - t520 * t586 + t522 * t541;
t495 = mrSges(6,2) * t503 - mrSges(6,3) * t501 + Ifges(6,1) * t530 + Ifges(6,4) * t529 + Ifges(6,5) * t540 - t520 * t587 - t521 * t541;
t532 = t577 * t555;
t477 = mrSges(5,2) * t509 - mrSges(5,3) * t505 - pkin(7) * t493 - t563 * t494 + t566 * t495 + t532 * t592 + t554 * t578;
t571 = mrSges(6,1) * t501 - mrSges(6,2) * t502 + Ifges(6,5) * t530 + Ifges(6,6) * t529 + Ifges(6,3) * t540;
t479 = Ifges(5,2) * t593 - mrSges(5,1) * t509 + mrSges(5,3) * t506 - pkin(4) * t493 + (Ifges(5,4) * t554 + (-t532 - t576) * t555) * t559 - t571;
t492 = t554 * t599 - t572;
t573 = mrSges(4,1) * t512 - mrSges(4,2) * t513 + Ifges(4,3) * t554 - pkin(3) * t492 + qJ(4) * t581 + t477 * t559 + t479 * t561;
t570 = mrSges(2,1) * t543 + mrSges(3,1) * t523 - mrSges(2,2) * t544 - mrSges(3,2) * t524 + pkin(1) * t467 + pkin(2) * t475 + t573 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t468 = t553 * Ifges(4,5) - mrSges(4,1) * t558 + mrSges(4,3) * t513 - mrSges(5,1) * t505 + mrSges(5,2) * t506 - t563 * t495 - t566 * t494 - pkin(4) * t575 - pkin(7) * t580 - pkin(3) * t485 + (Ifges(4,6) - t577) * t554 + (-pkin(4) * (-t526 * t595 - t527 * t594) + (-t559 * (Ifges(5,4) * t559 + Ifges(5,2) * t561) + t561 * t578) * t555) * t555;
t460 = mrSges(4,2) * t558 - mrSges(4,3) * t512 + Ifges(4,5) * t554 - Ifges(4,6) * t553 - qJ(4) * t485 + t477 * t561 - t479 * t559;
t459 = mrSges(3,2) * t558 - mrSges(3,3) * t523 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t569 - pkin(6) * t475 + t460 * t567 - t468 * t564;
t458 = -mrSges(3,1) * t558 + mrSges(3,3) * t524 + t569 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t585 + pkin(6) * t582 + t564 * t460 + t567 * t468;
t457 = -mrSges(2,2) * g(3) - mrSges(2,3) * t543 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t569 - qJ(2) * t467 - t458 * t560 + t459 * t562;
t456 = mrSges(2,1) * g(3) + mrSges(2,3) * t544 + t569 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t483 + qJ(2) * t583 + t562 * t458 + t560 * t459;
t1 = [-m(1) * g(1) + t584; -m(1) * g(2) + t590; (-m(1) - m(2)) * g(3) + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t590 - t456 * t565 + t457 * t568; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t584 + t568 * t456 + t565 * t457; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t570; t570; t483; t573; t492; t576 * t596 + t571;];
tauJB = t1;
