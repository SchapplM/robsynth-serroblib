% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:47
% EndTime: 2019-12-05 16:45:49
% DurationCPUTime: 1.38s
% Computational Cost: add. (15998->195), mult. (19883->240), div. (0->0), fcn. (10464->8), ass. (0->87)
t580 = Ifges(5,1) + Ifges(6,1);
t572 = Ifges(5,4) - Ifges(6,5);
t571 = Ifges(5,5) + Ifges(6,4);
t579 = Ifges(5,2) + Ifges(6,3);
t570 = Ifges(5,6) - Ifges(6,6);
t578 = Ifges(5,3) + Ifges(6,2);
t536 = qJD(2) + qJD(3);
t541 = sin(qJ(4));
t544 = cos(qJ(4));
t514 = (-mrSges(6,1) * t544 - mrSges(6,3) * t541) * t536;
t535 = qJDD(2) + qJDD(3);
t560 = qJD(4) * t536;
t516 = t541 * t535 + t544 * t560;
t540 = sin(pkin(8));
t568 = cos(pkin(8));
t527 = -t568 * g(1) - t540 * g(2);
t539 = -g(3) + qJDD(1);
t543 = sin(qJ(2));
t546 = cos(qJ(2));
t500 = -t543 * t527 + t546 * t539;
t497 = qJDD(2) * pkin(2) + t500;
t501 = t546 * t527 + t543 * t539;
t548 = qJD(2) ^ 2;
t498 = -t548 * pkin(2) + t501;
t542 = sin(qJ(3));
t545 = cos(qJ(3));
t493 = t542 * t497 + t545 * t498;
t534 = t536 ^ 2;
t490 = -t534 * pkin(3) + t535 * pkin(7) + t493;
t513 = (-pkin(4) * t544 - qJ(5) * t541) * t536;
t547 = qJD(4) ^ 2;
t526 = t540 * g(1) - t568 * g(2);
t565 = t544 * t526;
t485 = -qJDD(4) * pkin(4) - t547 * qJ(5) + t565 + qJDD(5) + (t513 * t536 + t490) * t541;
t566 = t536 * t544;
t525 = mrSges(6,2) * t566 + qJD(4) * mrSges(6,3);
t553 = -m(6) * t485 + qJDD(4) * mrSges(6,1) + qJD(4) * t525;
t567 = t536 * t541;
t481 = t516 * mrSges(6,2) + t514 * t567 - t553;
t487 = t544 * t490 - t541 * t526;
t484 = -t547 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t513 * t566 + t487;
t486 = -t541 * t490 - t565;
t517 = t544 * t535 - t541 * t560;
t523 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t567;
t554 = m(6) * t484 + qJDD(4) * mrSges(6,3) + qJD(4) * t523 + t514 * t566;
t561 = (t580 * t541 + t572 * t544) * t536 + t571 * qJD(4);
t563 = (-t572 * t541 - t579 * t544) * t536 - t570 * qJD(4);
t577 = -(t563 * t541 + t561 * t544) * t536 + t578 * qJDD(4) + t571 * t516 + t570 * t517 + mrSges(5,1) * t486 - mrSges(6,1) * t485 - mrSges(5,2) * t487 + mrSges(6,3) * t484 - pkin(4) * t481 + qJ(5) * (t517 * mrSges(6,2) + t554);
t574 = m(3) + m(4);
t573 = mrSges(5,3) + mrSges(6,2);
t515 = (-mrSges(5,1) * t544 + mrSges(5,2) * t541) * t536;
t522 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t567;
t477 = m(5) * t487 - qJDD(4) * mrSges(5,2) - qJD(4) * t522 + t515 * t566 + t573 * t517 + t554;
t524 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t566;
t478 = m(5) * t486 + qJDD(4) * mrSges(5,1) + qJD(4) * t524 + (-t514 - t515) * t567 - t573 * t516 + t553;
t555 = t544 * t477 - t541 * t478;
t464 = m(4) * t493 - t534 * mrSges(4,1) - t535 * mrSges(4,2) + t555;
t492 = t545 * t497 - t542 * t498;
t489 = -t535 * pkin(3) - t534 * pkin(7) - t492;
t482 = -t517 * pkin(4) - t516 * qJ(5) + (-0.2e1 * qJD(5) * t541 + (pkin(4) * t541 - qJ(5) * t544) * qJD(4)) * t536 + t489;
t479 = m(6) * t482 - t517 * mrSges(6,1) - t516 * mrSges(6,3) - t523 * t567 - t525 * t566;
t550 = -m(5) * t489 + t517 * mrSges(5,1) - t516 * mrSges(5,2) - t522 * t567 + t524 * t566 - t479;
t471 = m(4) * t492 + t535 * mrSges(4,1) - t534 * mrSges(4,2) + t550;
t457 = t542 * t464 + t545 * t471;
t455 = m(3) * t500 + qJDD(2) * mrSges(3,1) - t548 * mrSges(3,2) + t457;
t556 = t545 * t464 - t542 * t471;
t456 = m(3) * t501 - t548 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t556;
t557 = -t543 * t455 + t546 * t456;
t448 = m(2) * t527 + t557;
t468 = t541 * t477 + t544 * t478;
t466 = (m(2) + t574) * t526 - t468;
t564 = t540 * t448 + t568 * t466;
t449 = t546 * t455 + t543 * t456;
t562 = (t571 * t541 + t570 * t544) * t536 + t578 * qJD(4);
t559 = m(2) * t539 + t449;
t558 = t568 * t448 - t540 * t466;
t460 = -mrSges(5,1) * t489 - mrSges(6,1) * t482 + mrSges(6,2) * t484 + mrSges(5,3) * t487 - pkin(4) * t479 + t561 * qJD(4) + t570 * qJDD(4) + t572 * t516 + t579 * t517 - t562 * t567;
t461 = mrSges(5,2) * t489 + mrSges(6,2) * t485 - mrSges(5,3) * t486 - mrSges(6,3) * t482 - qJ(5) * t479 + t563 * qJD(4) + t571 * qJDD(4) + t580 * t516 + t572 * t517 + t562 * t566;
t552 = mrSges(4,1) * t492 - mrSges(4,2) * t493 + Ifges(4,3) * t535 + pkin(3) * t550 + pkin(7) * t555 + t544 * t460 + t541 * t461;
t549 = mrSges(3,1) * t500 - mrSges(3,2) * t501 + Ifges(3,3) * qJDD(2) + pkin(2) * t457 + t552;
t451 = mrSges(4,1) * t526 + mrSges(4,3) * t493 + t534 * Ifges(4,5) + Ifges(4,6) * t535 - pkin(3) * t468 - t577;
t450 = -mrSges(4,2) * t526 - mrSges(4,3) * t492 + Ifges(4,5) * t535 - t534 * Ifges(4,6) - pkin(7) * t468 - t541 * t460 + t544 * t461;
t445 = -mrSges(3,2) * t526 - mrSges(3,3) * t500 + Ifges(3,5) * qJDD(2) - t548 * Ifges(3,6) - pkin(6) * t457 + t545 * t450 - t542 * t451;
t444 = Ifges(3,6) * qJDD(2) + t548 * Ifges(3,5) + mrSges(3,1) * t526 + mrSges(3,3) * t501 + t542 * t450 + t545 * t451 - pkin(2) * (-m(4) * t526 + t468) + pkin(6) * t556;
t443 = -mrSges(2,1) * t539 + mrSges(2,3) * t527 - pkin(1) * t449 - t549;
t442 = mrSges(2,2) * t539 - mrSges(2,3) * t526 - pkin(5) * t449 - t543 * t444 + t546 * t445;
t1 = [-m(1) * g(1) + t558; -m(1) * g(2) + t564; -m(1) * g(3) + t559; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t564 + t568 * t442 - t540 * t443; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t558 + t540 * t442 + t568 * t443; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t527 + t543 * t445 + t546 * t444 - pkin(1) * t468 + pkin(5) * t557 + (pkin(1) * t574 + mrSges(2,1)) * t526; t559; t549; t552; t577; t481;];
tauJB = t1;
