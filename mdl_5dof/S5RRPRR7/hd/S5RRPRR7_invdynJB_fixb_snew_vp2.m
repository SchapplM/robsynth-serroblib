% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:37
% EndTime: 2019-12-31 20:15:39
% DurationCPUTime: 1.79s
% Computational Cost: add. (26457->224), mult. (32105->273), div. (0->0), fcn. (17394->8), ass. (0->97)
t580 = sin(qJ(1));
t584 = cos(qJ(1));
t557 = t580 * g(1) - t584 * g(2);
t551 = qJDD(1) * pkin(1) + t557;
t558 = -t584 * g(1) - t580 * g(2);
t585 = qJD(1) ^ 2;
t552 = -t585 * pkin(1) + t558;
t579 = sin(qJ(2));
t583 = cos(qJ(2));
t530 = t579 * t551 + t583 * t552;
t571 = qJDD(1) + qJDD(2);
t573 = (qJD(1) + qJD(2));
t595 = t571 * qJ(3) + (2 * qJD(3) * t573) + t530;
t611 = -m(3) - m(4);
t610 = -pkin(2) - pkin(7);
t609 = mrSges(3,1) - mrSges(4,2);
t608 = Ifges(3,5) - Ifges(4,4);
t607 = (-Ifges(3,6) + Ifges(4,5));
t578 = sin(qJ(4));
t606 = t573 * t578;
t582 = cos(qJ(4));
t605 = t573 * t582;
t529 = t583 * t551 - t579 * t552;
t569 = t573 ^ 2;
t594 = -t569 * qJ(3) + qJDD(3) - t529;
t519 = t610 * t571 + t594;
t513 = t578 * g(3) + t582 * t519;
t602 = qJD(4) * t573;
t600 = t578 * t602;
t547 = t582 * t571 - t600;
t503 = (-t547 - t600) * pkin(8) + (-t569 * t578 * t582 + qJDD(4)) * pkin(4) + t513;
t514 = -t582 * g(3) + t578 * t519;
t546 = -t578 * t571 - t582 * t602;
t555 = qJD(4) * pkin(4) - pkin(8) * t605;
t576 = t578 ^ 2;
t505 = -t576 * t569 * pkin(4) + t546 * pkin(8) - qJD(4) * t555 + t514;
t577 = sin(qJ(5));
t581 = cos(qJ(5));
t500 = t581 * t503 - t577 * t505;
t538 = (-t577 * t582 - t578 * t581) * t573;
t511 = t538 * qJD(5) + t577 * t546 + t581 * t547;
t539 = (-t577 * t578 + t581 * t582) * t573;
t527 = -t538 * mrSges(6,1) + t539 * mrSges(6,2);
t572 = qJD(4) + qJD(5);
t531 = -t572 * mrSges(6,2) + t538 * mrSges(6,3);
t570 = qJDD(4) + qJDD(5);
t497 = m(6) * t500 + t570 * mrSges(6,1) - t511 * mrSges(6,3) - t539 * t527 + t572 * t531;
t501 = t577 * t503 + t581 * t505;
t510 = -t539 * qJD(5) + t581 * t546 - t577 * t547;
t532 = t572 * mrSges(6,1) - t539 * mrSges(6,3);
t498 = m(6) * t501 - t570 * mrSges(6,2) + t510 * mrSges(6,3) + t538 * t527 - t572 * t532;
t487 = t581 * t497 + t577 * t498;
t545 = (mrSges(5,1) * t578 + mrSges(5,2) * t582) * t573;
t553 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t606;
t484 = m(5) * t513 + qJDD(4) * mrSges(5,1) - t547 * mrSges(5,3) + qJD(4) * t553 - t545 * t605 + t487;
t554 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t605;
t596 = -t577 * t497 + t581 * t498;
t485 = m(5) * t514 - qJDD(4) * mrSges(5,2) + t546 * mrSges(5,3) - qJD(4) * t554 - t545 * t606 + t596;
t482 = t582 * t484 + t578 * t485;
t526 = -t571 * pkin(2) + t594;
t592 = -m(4) * t526 + (t569 * mrSges(4,3)) - t482;
t478 = m(3) * t529 - (t569 * mrSges(3,2)) + t609 * t571 + t592;
t523 = t569 * pkin(2) - t595;
t516 = t610 * t569 + t595;
t504 = t555 * t605 - t546 * pkin(4) + (-pkin(8) * t576 + t610) * t569 + t595;
t593 = m(6) * t504 - t510 * mrSges(6,1) + t511 * mrSges(6,2) - t538 * t531 + t539 * t532;
t590 = -m(5) * t516 + t546 * mrSges(5,1) - t547 * mrSges(5,2) - t553 * t606 - t554 * t605 - t593;
t588 = -m(4) * t523 + (t569 * mrSges(4,2)) + t571 * mrSges(4,3) - t590;
t492 = m(3) * t530 - (t569 * mrSges(3,1)) - t571 * mrSges(3,2) + t588;
t473 = t583 * t478 + t579 * t492;
t470 = m(2) * t557 + qJDD(1) * mrSges(2,1) - t585 * mrSges(2,2) + t473;
t598 = -t579 * t478 + t583 * t492;
t471 = m(2) * t558 - t585 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t598;
t604 = t584 * t470 + t580 * t471;
t599 = -t580 * t470 + t584 * t471;
t597 = t578 * t484 - t582 * t485;
t521 = Ifges(6,4) * t539 + Ifges(6,2) * t538 + Ifges(6,6) * t572;
t522 = Ifges(6,1) * t539 + Ifges(6,4) * t538 + Ifges(6,5) * t572;
t591 = mrSges(6,1) * t500 - mrSges(6,2) * t501 + Ifges(6,5) * t511 + Ifges(6,6) * t510 + Ifges(6,3) * t570 + t539 * t521 - t538 * t522;
t520 = Ifges(6,5) * t539 + Ifges(6,6) * t538 + Ifges(6,3) * t572;
t488 = -mrSges(6,1) * t504 + mrSges(6,3) * t501 + Ifges(6,4) * t511 + Ifges(6,2) * t510 + Ifges(6,6) * t570 - t539 * t520 + t572 * t522;
t489 = mrSges(6,2) * t504 - mrSges(6,3) * t500 + Ifges(6,1) * t511 + Ifges(6,4) * t510 + Ifges(6,5) * t570 + t538 * t520 - t572 * t521;
t535 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t582 - Ifges(5,6) * t578) * t573;
t537 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t582 - Ifges(5,4) * t578) * t573;
t474 = -mrSges(5,1) * t516 + mrSges(5,3) * t514 + Ifges(5,4) * t547 + Ifges(5,2) * t546 + Ifges(5,6) * qJDD(4) - pkin(4) * t593 + pkin(8) * t596 + qJD(4) * t537 + t581 * t488 + t577 * t489 - t535 * t605;
t536 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t582 - Ifges(5,2) * t578) * t573;
t476 = mrSges(5,2) * t516 - mrSges(5,3) * t513 + Ifges(5,1) * t547 + Ifges(5,4) * t546 + Ifges(5,5) * qJDD(4) - pkin(8) * t487 - qJD(4) * t536 - t577 * t488 + t581 * t489 - t535 * t606;
t480 = t571 * mrSges(4,2) - t592;
t589 = mrSges(3,1) * t529 - mrSges(3,2) * t530 + mrSges(4,2) * t526 - mrSges(4,3) * t523 - pkin(2) * t480 - pkin(7) * t482 + qJ(3) * t588 - t578 * t474 + t582 * t476 + (Ifges(3,3) + Ifges(4,1)) * t571;
t587 = mrSges(5,1) * t513 - mrSges(5,2) * t514 + Ifges(5,5) * t547 + Ifges(5,6) * t546 + Ifges(5,3) * qJDD(4) + pkin(4) * t487 + t536 * t605 + t537 * t606 + t591;
t586 = mrSges(2,1) * t557 - mrSges(2,2) * t558 + Ifges(2,3) * qJDD(1) + pkin(1) * t473 + t589;
t481 = -m(4) * g(3) - t597;
t466 = t587 + (-mrSges(3,2) + mrSges(4,3)) * g(3) + t608 * t571 + (t607 * t569) + mrSges(4,1) * t526 - mrSges(3,3) * t529 - qJ(3) * t481 + pkin(3) * t482;
t465 = -mrSges(4,1) * t523 + mrSges(3,3) * t530 - pkin(2) * t481 - pkin(3) * t590 + pkin(7) * t597 + t609 * g(3) - t582 * t474 - t578 * t476 + t608 * t569 - t607 * t571;
t464 = -mrSges(2,2) * g(3) - mrSges(2,3) * t557 + Ifges(2,5) * qJDD(1) - t585 * Ifges(2,6) - pkin(6) * t473 - t579 * t465 + t583 * t466;
t463 = Ifges(2,6) * qJDD(1) + t585 * Ifges(2,5) + mrSges(2,3) * t558 + t579 * t466 + t583 * t465 + pkin(1) * t597 + pkin(6) * t598 + (-pkin(1) * t611 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t599; -m(1) * g(2) + t604; (-m(1) - m(2) + t611) * g(3) - t597; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t604 - t580 * t463 + t584 * t464; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t599 + t584 * t463 + t580 * t464; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t586; t586; t589; t480; t587; t591;];
tauJB = t1;
