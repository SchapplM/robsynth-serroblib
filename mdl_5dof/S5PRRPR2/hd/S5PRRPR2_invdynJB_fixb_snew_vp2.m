% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:31
% EndTime: 2019-12-05 16:17:34
% DurationCPUTime: 2.45s
% Computational Cost: add. (31169->196), mult. (43139->262), div. (0->0), fcn. (27676->10), ass. (0->103)
t580 = 2 * qJD(4);
t537 = sin(pkin(8));
t539 = cos(pkin(8));
t518 = t537 * g(1) - t539 * g(2);
t520 = -t539 * g(1) - t537 * g(2);
t542 = sin(qJ(2));
t545 = cos(qJ(2));
t504 = t545 * t518 - t542 * t520;
t501 = qJDD(2) * pkin(2) + t504;
t505 = t542 * t518 + t545 * t520;
t546 = qJD(2) ^ 2;
t502 = -t546 * pkin(2) + t505;
t541 = sin(qJ(3));
t544 = cos(qJ(3));
t494 = t541 * t501 + t544 * t502;
t533 = (qJD(2) + qJD(3));
t531 = t533 ^ 2;
t532 = qJDD(2) + qJDD(3);
t491 = -t531 * pkin(3) + t532 * qJ(4) + t494;
t579 = (t533 * t580) + t491;
t536 = sin(pkin(9));
t578 = mrSges(5,2) * t536;
t576 = mrSges(5,3) * t532;
t575 = t536 * t533;
t540 = sin(qJ(5));
t574 = t536 * t540;
t543 = cos(qJ(5));
t573 = t536 * t543;
t538 = cos(pkin(9));
t572 = t538 * t532;
t571 = t538 * t533;
t535 = -g(3) + qJDD(1);
t570 = t538 * t535;
t487 = t536 * t535 + t579 * t538;
t512 = (-mrSges(5,1) * t538 + t578) * t533;
t557 = -pkin(4) * t538 - pkin(7) * t536;
t514 = t557 * t533;
t485 = t514 * t571 + t487;
t493 = t544 * t501 - t541 * t502;
t551 = -t531 * qJ(4) + qJDD(4) - t493;
t488 = (-pkin(3) + t557) * t532 + t551;
t482 = -t540 * t485 + t543 * t488;
t521 = qJD(5) - t571;
t566 = t533 * t574;
t507 = -t521 * mrSges(6,2) - mrSges(6,3) * t566;
t509 = (mrSges(6,1) * t540 + mrSges(6,2) * t543) * t575;
t567 = qJD(5) * t533;
t511 = (t532 * t543 - t540 * t567) * t536;
t519 = qJDD(5) - t572;
t565 = t533 * t573;
t480 = m(6) * t482 + t519 * mrSges(6,1) - t511 * mrSges(6,3) + t521 * t507 - t509 * t565;
t483 = t543 * t485 + t540 * t488;
t508 = t521 * mrSges(6,1) - mrSges(6,3) * t565;
t510 = (-t532 * t540 - t543 * t567) * t536;
t481 = m(6) * t483 - t519 * mrSges(6,2) + t510 * mrSges(6,3) - t521 * t508 - t509 * t566;
t559 = -t540 * t480 + t543 * t481;
t471 = m(5) * t487 + (t512 * t533 + t576) * t538 + t559;
t486 = -t579 * t536 + t570;
t484 = -t570 + (t491 + (t580 + t514) * t533) * t536;
t552 = -m(6) * t484 + t510 * mrSges(6,1) - t511 * mrSges(6,2);
t478 = m(5) * t486 + (-t576 + (-t507 * t540 - t508 * t543 - t512) * t533) * t536 + t552;
t560 = t538 * t471 - t536 * t478;
t464 = m(4) * t494 - t531 * mrSges(4,1) - t532 * mrSges(4,2) + t560;
t474 = t543 * t480 + t540 * t481;
t490 = -t532 * pkin(3) + t551;
t549 = -m(5) * t490 + mrSges(5,1) * t572 - t474 + (t536 ^ 2 + t538 ^ 2) * mrSges(5,3) * t531;
t468 = m(4) * t493 - t531 * mrSges(4,2) + (mrSges(4,1) - t578) * t532 + t549;
t457 = t541 * t464 + t544 * t468;
t454 = m(3) * t504 + qJDD(2) * mrSges(3,1) - t546 * mrSges(3,2) + t457;
t561 = t544 * t464 - t541 * t468;
t455 = m(3) * t505 - t546 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t561;
t449 = t545 * t454 + t542 * t455;
t447 = m(2) * t518 + t449;
t562 = -t542 * t454 + t545 * t455;
t448 = m(2) * t520 + t562;
t569 = t539 * t447 + t537 * t448;
t466 = t536 * t471 + t538 * t478;
t564 = m(4) * t535 + t466;
t563 = -t537 * t447 + t539 * t448;
t558 = m(3) * t535 + t564;
t556 = Ifges(5,1) * t536 + Ifges(5,4) * t538;
t555 = Ifges(5,5) * t536 + Ifges(5,6) * t538;
t554 = m(2) * t535 + t558;
t496 = Ifges(6,6) * t521 + (Ifges(6,4) * t543 - Ifges(6,2) * t540) * t575;
t497 = Ifges(6,5) * t521 + (Ifges(6,1) * t543 - Ifges(6,4) * t540) * t575;
t553 = t496 * t543 + t497 * t540;
t495 = Ifges(6,3) * t521 + (Ifges(6,5) * t543 - Ifges(6,6) * t540) * t575;
t475 = -mrSges(6,1) * t484 + mrSges(6,3) * t483 + Ifges(6,4) * t511 + Ifges(6,2) * t510 + Ifges(6,6) * t519 - t495 * t565 + t521 * t497;
t476 = mrSges(6,2) * t484 - mrSges(6,3) * t482 + Ifges(6,1) * t511 + Ifges(6,4) * t510 + Ifges(6,5) * t519 - t495 * t566 - t521 * t496;
t513 = t555 * t533;
t459 = mrSges(5,2) * t490 - mrSges(5,3) * t486 - pkin(7) * t474 - t540 * t475 + t543 * t476 + t513 * t571 + t532 * t556;
t548 = mrSges(6,1) * t482 - mrSges(6,2) * t483 + Ifges(6,5) * t511 + Ifges(6,6) * t510 + Ifges(6,3) * t519;
t461 = Ifges(5,2) * t572 - mrSges(5,1) * t490 + mrSges(5,3) * t487 - pkin(4) * t474 + (Ifges(5,4) * t532 + (-t513 - t553) * t533) * t536 - t548;
t473 = t532 * t578 - t549;
t550 = mrSges(4,1) * t493 - mrSges(4,2) * t494 + Ifges(4,3) * t532 - pkin(3) * t473 + qJ(4) * t560 + t536 * t459 + t538 * t461;
t547 = mrSges(3,1) * t504 - mrSges(3,2) * t505 + Ifges(3,3) * qJDD(2) + pkin(2) * t457 + t550;
t450 = t531 * Ifges(4,5) - mrSges(4,1) * t535 + mrSges(4,3) * t494 - mrSges(5,1) * t486 + mrSges(5,2) * t487 - t540 * t476 - t543 * t475 - pkin(4) * t552 - pkin(7) * t559 - pkin(3) * t466 + (Ifges(4,6) - t555) * t532 + (-pkin(4) * (-t507 * t574 - t508 * t573) + (-t536 * (Ifges(5,4) * t536 + Ifges(5,2) * t538) + t538 * t556) * t533) * t533;
t443 = mrSges(4,2) * t535 - mrSges(4,3) * t493 + Ifges(4,5) * t532 - t531 * Ifges(4,6) - qJ(4) * t466 + t538 * t459 - t536 * t461;
t442 = mrSges(3,2) * t535 - mrSges(3,3) * t504 + Ifges(3,5) * qJDD(2) - t546 * Ifges(3,6) - pkin(6) * t457 + t544 * t443 - t541 * t450;
t441 = -mrSges(3,1) * t535 + mrSges(3,3) * t505 + t546 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t564 + pkin(6) * t561 + t541 * t443 + t544 * t450;
t440 = mrSges(2,2) * t535 - mrSges(2,3) * t518 - pkin(5) * t449 - t542 * t441 + t545 * t442;
t439 = -mrSges(2,1) * t535 + mrSges(2,3) * t520 - pkin(1) * t558 + pkin(5) * t562 + t545 * t441 + t542 * t442;
t1 = [-m(1) * g(1) + t563; -m(1) * g(2) + t569; -m(1) * g(3) + t554; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t569 - t537 * t439 + t539 * t440; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t563 + t539 * t439 + t537 * t440; -mrSges(1,1) * g(2) + mrSges(2,1) * t518 + mrSges(1,2) * g(1) - mrSges(2,2) * t520 + pkin(1) * t449 + t547; t554; t547; t550; t473; t553 * t575 + t548;];
tauJB = t1;
