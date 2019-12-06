% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:20
% EndTime: 2019-12-05 15:19:30
% DurationCPUTime: 7.60s
% Computational Cost: add. (93701->221), mult. (166157->297), div. (0->0), fcn. (129443->14), ass. (0->110)
t570 = sin(pkin(10));
t574 = cos(pkin(10));
t563 = -t574 * g(1) - t570 * g(2);
t569 = sin(pkin(11));
t573 = cos(pkin(11));
t562 = t570 * g(1) - t574 * g(2);
t568 = -g(3) + qJDD(1);
t572 = sin(pkin(5));
t576 = cos(pkin(5));
t590 = t562 * t576 + t568 * t572;
t536 = -t569 * t563 + t590 * t573;
t537 = t573 * t563 + t590 * t569;
t549 = -t572 * t562 + t576 * t568 + qJDD(2);
t582 = cos(qJ(3));
t575 = cos(pkin(6));
t579 = sin(qJ(3));
t603 = t575 * t579;
t571 = sin(pkin(6));
t604 = t571 * t579;
t530 = t536 * t603 + t582 * t537 + t549 * t604;
t584 = qJD(3) ^ 2;
t528 = -t584 * pkin(3) + qJDD(3) * pkin(8) + t530;
t532 = -t571 * t536 + t575 * t549;
t578 = sin(qJ(4));
t581 = cos(qJ(4));
t524 = t581 * t528 + t578 * t532;
t559 = (-pkin(4) * t581 - pkin(9) * t578) * qJD(3);
t583 = qJD(4) ^ 2;
t599 = t581 * qJD(3);
t522 = -t583 * pkin(4) + qJDD(4) * pkin(9) + t559 * t599 + t524;
t529 = -t579 * t537 + (t536 * t575 + t549 * t571) * t582;
t527 = -qJDD(3) * pkin(3) - t584 * pkin(8) - t529;
t598 = qJD(3) * qJD(4);
t596 = t581 * t598;
t560 = t578 * qJDD(3) + t596;
t597 = t578 * t598;
t561 = t581 * qJDD(3) - t597;
t525 = (-t560 - t596) * pkin(9) + (-t561 + t597) * pkin(4) + t527;
t577 = sin(qJ(5));
t580 = cos(qJ(5));
t518 = -t577 * t522 + t580 * t525;
t600 = qJD(3) * t578;
t556 = t580 * qJD(4) - t577 * t600;
t544 = t556 * qJD(5) + t577 * qJDD(4) + t580 * t560;
t557 = t577 * qJD(4) + t580 * t600;
t545 = -t556 * mrSges(6,1) + t557 * mrSges(6,2);
t566 = qJD(5) - t599;
t547 = -t566 * mrSges(6,2) + t556 * mrSges(6,3);
t555 = qJDD(5) - t561;
t516 = m(6) * t518 + t555 * mrSges(6,1) - t544 * mrSges(6,3) - t557 * t545 + t566 * t547;
t519 = t580 * t522 + t577 * t525;
t543 = -t557 * qJD(5) + t580 * qJDD(4) - t577 * t560;
t548 = t566 * mrSges(6,1) - t557 * mrSges(6,3);
t517 = m(6) * t519 - t555 * mrSges(6,2) + t543 * mrSges(6,3) + t556 * t545 - t566 * t548;
t510 = -t577 * t516 + t580 * t517;
t558 = (-mrSges(5,1) * t581 + mrSges(5,2) * t578) * qJD(3);
t564 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t600;
t508 = m(5) * t524 - qJDD(4) * mrSges(5,2) + t561 * mrSges(5,3) - qJD(4) * t564 + t558 * t599 + t510;
t602 = t581 * t532;
t521 = -qJDD(4) * pkin(4) - t583 * pkin(9) - t602 + (qJD(3) * t559 + t528) * t578;
t520 = -m(6) * t521 + t543 * mrSges(6,1) - t544 * mrSges(6,2) + t556 * t547 - t557 * t548;
t523 = -t578 * t528 + t602;
t565 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t599;
t514 = m(5) * t523 + qJDD(4) * mrSges(5,1) - t560 * mrSges(5,3) + qJD(4) * t565 - t558 * t600 + t520;
t594 = t581 * t508 - t578 * t514;
t499 = m(4) * t530 - t584 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t594;
t502 = t578 * t508 + t581 * t514;
t501 = m(4) * t532 + t502;
t509 = t580 * t516 + t577 * t517;
t587 = -m(5) * t527 + t561 * mrSges(5,1) - t560 * mrSges(5,2) - t564 * t600 + t565 * t599 - t509;
t505 = m(4) * t529 + qJDD(3) * mrSges(4,1) - t584 * mrSges(4,2) + t587;
t605 = t505 * t582;
t488 = t499 * t603 - t571 * t501 + t575 * t605;
t484 = m(3) * t536 + t488;
t493 = t582 * t499 - t579 * t505;
t492 = m(3) * t537 + t493;
t609 = t484 * t573 + t492 * t569;
t538 = Ifges(6,5) * t557 + Ifges(6,6) * t556 + Ifges(6,3) * t566;
t540 = Ifges(6,1) * t557 + Ifges(6,4) * t556 + Ifges(6,5) * t566;
t511 = -mrSges(6,1) * t521 + mrSges(6,3) * t519 + Ifges(6,4) * t544 + Ifges(6,2) * t543 + Ifges(6,6) * t555 - t557 * t538 + t566 * t540;
t539 = Ifges(6,4) * t557 + Ifges(6,2) * t556 + Ifges(6,6) * t566;
t512 = mrSges(6,2) * t521 - mrSges(6,3) * t518 + Ifges(6,1) * t544 + Ifges(6,4) * t543 + Ifges(6,5) * t555 + t556 * t538 - t566 * t539;
t551 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t578 + Ifges(5,2) * t581) * qJD(3);
t552 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t578 + Ifges(5,4) * t581) * qJD(3);
t608 = mrSges(5,1) * t523 - mrSges(5,2) * t524 + Ifges(5,5) * t560 + Ifges(5,6) * t561 + Ifges(5,3) * qJDD(4) + pkin(4) * t520 + pkin(9) * t510 + t580 * t511 + t577 * t512 + (t578 * t551 - t581 * t552) * qJD(3);
t487 = t499 * t604 + t575 * t501 + t571 * t605;
t486 = m(3) * t549 + t487;
t474 = -t572 * t486 + t609 * t576;
t472 = m(2) * t562 + t474;
t478 = -t569 * t484 + t573 * t492;
t477 = m(2) * t563 + t478;
t601 = t574 * t472 + t570 * t477;
t473 = t576 * t486 + t609 * t572;
t595 = -t570 * t472 + t574 * t477;
t593 = m(2) * t568 + t473;
t550 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t578 + Ifges(5,6) * t581) * qJD(3);
t494 = mrSges(5,2) * t527 - mrSges(5,3) * t523 + Ifges(5,1) * t560 + Ifges(5,4) * t561 + Ifges(5,5) * qJDD(4) - pkin(9) * t509 - qJD(4) * t551 - t577 * t511 + t580 * t512 + t550 * t599;
t586 = mrSges(6,1) * t518 - mrSges(6,2) * t519 + Ifges(6,5) * t544 + Ifges(6,6) * t543 + Ifges(6,3) * t555 + t557 * t539 - t556 * t540;
t495 = -mrSges(5,1) * t527 + mrSges(5,3) * t524 + Ifges(5,4) * t560 + Ifges(5,2) * t561 + Ifges(5,6) * qJDD(4) - pkin(4) * t509 + qJD(4) * t552 - t550 * t600 - t586;
t480 = mrSges(4,2) * t532 - mrSges(4,3) * t529 + Ifges(4,5) * qJDD(3) - t584 * Ifges(4,6) - pkin(8) * t502 + t581 * t494 - t578 * t495;
t481 = -mrSges(4,1) * t532 + mrSges(4,3) * t530 + t584 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t502 - t608;
t589 = pkin(7) * t493 + t480 * t579 + t481 * t582;
t479 = mrSges(4,1) * t529 - mrSges(4,2) * t530 + Ifges(4,3) * qJDD(3) + pkin(3) * t587 + pkin(8) * t594 + t578 * t494 + t581 * t495;
t469 = -mrSges(3,1) * t549 + mrSges(3,3) * t537 - pkin(2) * t487 - t571 * t479 + t589 * t575;
t470 = mrSges(3,2) * t549 - mrSges(3,3) * t536 + t582 * t480 - t579 * t481 + (-t487 * t571 - t488 * t575) * pkin(7);
t588 = qJ(2) * t478 + t469 * t573 + t470 * t569;
t468 = mrSges(3,1) * t536 - mrSges(3,2) * t537 + pkin(2) * t488 + t575 * t479 + t589 * t571;
t467 = mrSges(2,2) * t568 - mrSges(2,3) * t562 - t569 * t469 + t573 * t470 + (-t473 * t572 - t474 * t576) * qJ(2);
t466 = -mrSges(2,1) * t568 + mrSges(2,3) * t563 - pkin(1) * t473 - t572 * t468 + t588 * t576;
t1 = [-m(1) * g(1) + t595; -m(1) * g(2) + t601; -m(1) * g(3) + t593; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t601 - t570 * t466 + t574 * t467; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t595 + t574 * t466 + t570 * t467; -mrSges(1,1) * g(2) + mrSges(2,1) * t562 + mrSges(1,2) * g(1) - mrSges(2,2) * t563 + pkin(1) * t474 + t576 * t468 + t588 * t572; t593; t486; t479; t608; t586;];
tauJB = t1;
