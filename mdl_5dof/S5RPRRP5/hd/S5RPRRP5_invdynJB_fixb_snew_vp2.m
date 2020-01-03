% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:50
% EndTime: 2019-12-31 18:40:52
% DurationCPUTime: 1.73s
% Computational Cost: add. (21288->204), mult. (28835->249), div. (0->0), fcn. (13996->8), ass. (0->89)
t600 = Ifges(5,1) + Ifges(6,1);
t593 = Ifges(5,4) - Ifges(6,5);
t592 = Ifges(5,5) + Ifges(6,4);
t599 = Ifges(5,2) + Ifges(6,3);
t591 = Ifges(5,6) - Ifges(6,6);
t598 = Ifges(5,3) + Ifges(6,2);
t555 = qJD(1) + qJD(3);
t563 = sin(qJ(4));
t566 = cos(qJ(4));
t529 = (-mrSges(6,1) * t566 - mrSges(6,3) * t563) * t555;
t554 = qJDD(1) + qJDD(3);
t582 = qJD(4) * t555;
t531 = t563 * t554 + t566 * t582;
t565 = sin(qJ(1));
t568 = cos(qJ(1));
t545 = t565 * g(1) - t568 * g(2);
t538 = qJDD(1) * pkin(1) + t545;
t546 = -t568 * g(1) - t565 * g(2);
t570 = qJD(1) ^ 2;
t539 = -t570 * pkin(1) + t546;
t561 = sin(pkin(8));
t562 = cos(pkin(8));
t514 = t562 * t538 - t561 * t539;
t511 = qJDD(1) * pkin(2) + t514;
t515 = t561 * t538 + t562 * t539;
t512 = -t570 * pkin(2) + t515;
t564 = sin(qJ(3));
t567 = cos(qJ(3));
t507 = t564 * t511 + t567 * t512;
t553 = t555 ^ 2;
t504 = -t553 * pkin(3) + t554 * pkin(7) + t507;
t528 = (-pkin(4) * t566 - qJ(5) * t563) * t555;
t569 = qJD(4) ^ 2;
t560 = -g(3) + qJDD(2);
t587 = t566 * t560;
t499 = -qJDD(4) * pkin(4) - t569 * qJ(5) - t587 + qJDD(5) + (t528 * t555 + t504) * t563;
t588 = t555 * t566;
t543 = mrSges(6,2) * t588 + qJD(4) * mrSges(6,3);
t575 = -m(6) * t499 + qJDD(4) * mrSges(6,1) + qJD(4) * t543;
t589 = t555 * t563;
t495 = t531 * mrSges(6,2) + t529 * t589 - t575;
t501 = t566 * t504 + t563 * t560;
t498 = -t569 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t528 * t588 + t501;
t500 = -t563 * t504 + t587;
t532 = t566 * t554 - t563 * t582;
t541 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t589;
t576 = m(6) * t498 + qJDD(4) * mrSges(6,3) + qJD(4) * t541 + t529 * t588;
t583 = (t600 * t563 + t593 * t566) * t555 + t592 * qJD(4);
t585 = (-t593 * t563 - t599 * t566) * t555 - t591 * qJD(4);
t597 = -(t585 * t563 + t583 * t566) * t555 + t598 * qJDD(4) + t592 * t531 + t591 * t532 + mrSges(5,1) * t500 - mrSges(6,1) * t499 - mrSges(5,2) * t501 + mrSges(6,3) * t498 - pkin(4) * t495 + qJ(5) * (t532 * mrSges(6,2) + t576);
t594 = mrSges(5,3) + mrSges(6,2);
t530 = (-mrSges(5,1) * t566 + mrSges(5,2) * t563) * t555;
t540 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t589;
t491 = m(5) * t501 - qJDD(4) * mrSges(5,2) - qJD(4) * t540 + t530 * t588 + t594 * t532 + t576;
t542 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t588;
t492 = m(5) * t500 + qJDD(4) * mrSges(5,1) + qJD(4) * t542 + (-t529 - t530) * t589 - t594 * t531 + t575;
t577 = t566 * t491 - t563 * t492;
t481 = m(4) * t507 - t553 * mrSges(4,1) - t554 * mrSges(4,2) + t577;
t506 = t567 * t511 - t564 * t512;
t503 = -t554 * pkin(3) - t553 * pkin(7) - t506;
t496 = -t532 * pkin(4) - t531 * qJ(5) + (-0.2e1 * qJD(5) * t563 + (pkin(4) * t563 - qJ(5) * t566) * qJD(4)) * t555 + t503;
t493 = m(6) * t496 - t532 * mrSges(6,1) - t531 * mrSges(6,3) - t541 * t589 - t543 * t588;
t572 = -m(5) * t503 + t532 * mrSges(5,1) - t531 * mrSges(5,2) - t540 * t589 + t542 * t588 - t493;
t486 = m(4) * t506 + t554 * mrSges(4,1) - t553 * mrSges(4,2) + t572;
t474 = t564 * t481 + t567 * t486;
t471 = m(3) * t514 + qJDD(1) * mrSges(3,1) - t570 * mrSges(3,2) + t474;
t578 = t567 * t481 - t564 * t486;
t472 = m(3) * t515 - t570 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t578;
t465 = t562 * t471 + t561 * t472;
t462 = m(2) * t545 + qJDD(1) * mrSges(2,1) - t570 * mrSges(2,2) + t465;
t579 = -t561 * t471 + t562 * t472;
t463 = m(2) * t546 - t570 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t579;
t586 = t568 * t462 + t565 * t463;
t484 = t563 * t491 + t566 * t492;
t584 = (t592 * t563 + t591 * t566) * t555 + t598 * qJD(4);
t581 = m(4) * t560 + t484;
t580 = -t565 * t462 + t568 * t463;
t482 = m(3) * t560 + t581;
t477 = -mrSges(5,1) * t503 - mrSges(6,1) * t496 + mrSges(6,2) * t498 + mrSges(5,3) * t501 - pkin(4) * t493 + t583 * qJD(4) + t591 * qJDD(4) + t593 * t531 + t599 * t532 - t584 * t589;
t478 = mrSges(5,2) * t503 + mrSges(6,2) * t499 - mrSges(5,3) * t500 - mrSges(6,3) * t496 - qJ(5) * t493 + t585 * qJD(4) + t592 * qJDD(4) + t600 * t531 + t593 * t532 + t584 * t588;
t574 = mrSges(4,1) * t506 - mrSges(4,2) * t507 + Ifges(4,3) * t554 + pkin(3) * t572 + pkin(7) * t577 + t566 * t477 + t563 * t478;
t571 = mrSges(2,1) * t545 + mrSges(3,1) * t514 - mrSges(2,2) * t546 - mrSges(3,2) * t515 + pkin(1) * t465 + pkin(2) * t474 + t574 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t467 = -mrSges(4,1) * t560 + mrSges(4,3) * t507 + t553 * Ifges(4,5) + Ifges(4,6) * t554 - pkin(3) * t484 - t597;
t466 = mrSges(4,2) * t560 - mrSges(4,3) * t506 + Ifges(4,5) * t554 - t553 * Ifges(4,6) - pkin(7) * t484 - t563 * t477 + t566 * t478;
t458 = mrSges(3,2) * t560 - mrSges(3,3) * t514 + Ifges(3,5) * qJDD(1) - t570 * Ifges(3,6) - pkin(6) * t474 + t567 * t466 - t564 * t467;
t457 = -mrSges(3,1) * t560 + mrSges(3,3) * t515 + t570 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t581 + pkin(6) * t578 + t564 * t466 + t567 * t467;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t545 + Ifges(2,5) * qJDD(1) - t570 * Ifges(2,6) - qJ(2) * t465 - t561 * t457 + t562 * t458;
t455 = mrSges(2,1) * g(3) + mrSges(2,3) * t546 + t570 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t482 + qJ(2) * t579 + t562 * t457 + t561 * t458;
t1 = [-m(1) * g(1) + t580; -m(1) * g(2) + t586; (-m(1) - m(2)) * g(3) + t482; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t586 - t565 * t455 + t568 * t456; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t580 + t568 * t455 + t565 * t456; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t571; t571; t482; t574; t597; t495;];
tauJB = t1;
