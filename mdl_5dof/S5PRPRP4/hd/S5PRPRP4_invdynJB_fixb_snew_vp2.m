% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:17
% EndTime: 2019-12-05 15:35:21
% DurationCPUTime: 1.36s
% Computational Cost: add. (11376->193), mult. (19366->238), div. (0->0), fcn. (10199->8), ass. (0->84)
t561 = Ifges(5,1) + Ifges(6,1);
t553 = Ifges(5,4) - Ifges(6,5);
t552 = -Ifges(5,5) - Ifges(6,4);
t560 = Ifges(5,2) + Ifges(6,3);
t551 = Ifges(6,6) - Ifges(5,6);
t559 = Ifges(5,3) + Ifges(6,2);
t527 = sin(qJ(4));
t529 = cos(qJ(4));
t502 = (-mrSges(6,1) * t529 - mrSges(6,3) * t527) * qJD(2);
t541 = qJD(2) * qJD(4);
t504 = t527 * qJDD(2) + t529 * t541;
t524 = sin(pkin(7));
t526 = cos(pkin(7));
t510 = -t526 * g(1) - t524 * g(2);
t522 = -g(3) + qJDD(1);
t528 = sin(qJ(2));
t530 = cos(qJ(2));
t485 = -t528 * t510 + t530 * t522;
t483 = qJDD(2) * pkin(2) + t485;
t486 = t530 * t510 + t528 * t522;
t532 = qJD(2) ^ 2;
t484 = -t532 * pkin(2) + t486;
t523 = sin(pkin(8));
t525 = cos(pkin(8));
t479 = t523 * t483 + t525 * t484;
t477 = -t532 * pkin(3) + qJDD(2) * pkin(6) + t479;
t501 = (-pkin(4) * t529 - qJ(5) * t527) * qJD(2);
t531 = qJD(4) ^ 2;
t509 = t524 * g(1) - t526 * g(2);
t508 = qJDD(3) - t509;
t548 = t529 * t508;
t472 = -qJDD(4) * pkin(4) - t531 * qJ(5) - t548 + qJDD(5) + (qJD(2) * t501 + t477) * t527;
t542 = qJD(2) * t529;
t514 = mrSges(6,2) * t542 + qJD(4) * mrSges(6,3);
t535 = -m(6) * t472 + qJDD(4) * mrSges(6,1) + qJD(4) * t514;
t543 = qJD(2) * t527;
t468 = t504 * mrSges(6,2) + t502 * t543 - t535;
t474 = t529 * t477 + t527 * t508;
t471 = -t531 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t501 * t542 + t474;
t473 = -t527 * t477 + t548;
t505 = t529 * qJDD(2) - t527 * t541;
t512 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t543;
t536 = m(6) * t471 + qJDD(4) * mrSges(6,3) + qJD(4) * t512 + t502 * t542;
t544 = -t552 * qJD(4) + (t561 * t527 + t553 * t529) * qJD(2);
t546 = t551 * qJD(4) + (-t553 * t527 - t560 * t529) * qJD(2);
t558 = -(t546 * t527 + t544 * t529) * qJD(2) + t559 * qJDD(4) - t552 * t504 - t551 * t505 + mrSges(5,1) * t473 - mrSges(6,1) * t472 - mrSges(5,2) * t474 + mrSges(6,3) * t471 - pkin(4) * t468 + qJ(5) * (t505 * mrSges(6,2) + t536);
t503 = (-mrSges(5,1) * t529 + mrSges(5,2) * t527) * qJD(2);
t511 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t543;
t554 = mrSges(5,3) + mrSges(6,2);
t464 = m(5) * t474 - qJDD(4) * mrSges(5,2) - qJD(4) * t511 + t503 * t542 + t554 * t505 + t536;
t513 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t542;
t465 = m(5) * t473 + qJDD(4) * mrSges(5,1) + qJD(4) * t513 - t554 * t504 + (-t502 - t503) * t543 + t535;
t457 = t529 * t464 - t527 * t465;
t452 = m(4) * t479 - t532 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t457;
t478 = t525 * t483 - t523 * t484;
t476 = -qJDD(2) * pkin(3) - t532 * pkin(6) - t478;
t469 = -t505 * pkin(4) - t504 * qJ(5) + (-0.2e1 * qJD(5) * t527 + (pkin(4) * t527 - qJ(5) * t529) * qJD(4)) * qJD(2) + t476;
t466 = m(6) * t469 - t505 * mrSges(6,1) - t504 * mrSges(6,3) - t512 * t543 - t514 * t542;
t460 = -m(5) * t476 + t505 * mrSges(5,1) - t504 * mrSges(5,2) - t511 * t543 + t513 * t542 - t466;
t459 = m(4) * t478 + qJDD(2) * mrSges(4,1) - t532 * mrSges(4,2) + t460;
t447 = t523 * t452 + t525 * t459;
t545 = t559 * qJD(4) + (-t552 * t527 - t551 * t529) * qJD(2);
t448 = -mrSges(5,1) * t476 - mrSges(6,1) * t469 + mrSges(6,2) * t471 + mrSges(5,3) * t474 - pkin(4) * t466 + t544 * qJD(4) - t551 * qJDD(4) + t553 * t504 + t560 * t505 - t545 * t543;
t449 = mrSges(5,2) * t476 + mrSges(6,2) * t472 - mrSges(5,3) * t473 - mrSges(6,3) * t469 - qJ(5) * t466 + t546 * qJD(4) - t552 * qJDD(4) + t561 * t504 + t553 * t505 + t545 * t542;
t557 = mrSges(3,1) * t485 + mrSges(4,1) * t478 - mrSges(3,2) * t486 - mrSges(4,2) * t479 + pkin(2) * t447 + pkin(3) * t460 + pkin(6) * t457 + t529 * t448 + t527 * t449 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t445 = m(3) * t485 + qJDD(2) * mrSges(3,1) - t532 * mrSges(3,2) + t447;
t537 = t525 * t452 - t523 * t459;
t446 = m(3) * t486 - t532 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t537;
t538 = -t528 * t445 + t530 * t446;
t438 = m(2) * t510 + t538;
t456 = t527 * t464 + t529 * t465;
t455 = m(4) * t508 + t456;
t454 = (m(2) + m(3)) * t509 - t455;
t547 = t524 * t438 + t526 * t454;
t439 = t530 * t445 + t528 * t446;
t540 = m(2) * t522 + t439;
t539 = t526 * t438 - t524 * t454;
t441 = -mrSges(4,1) * t508 + mrSges(4,3) * t479 + t532 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t456 - t558;
t440 = mrSges(4,2) * t508 - mrSges(4,3) * t478 + Ifges(4,5) * qJDD(2) - t532 * Ifges(4,6) - pkin(6) * t456 - t527 * t448 + t529 * t449;
t435 = -mrSges(3,2) * t509 - mrSges(3,3) * t485 + Ifges(3,5) * qJDD(2) - t532 * Ifges(3,6) - qJ(3) * t447 + t525 * t440 - t523 * t441;
t434 = mrSges(3,1) * t509 + mrSges(3,3) * t486 + t532 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t455 + qJ(3) * t537 + t523 * t440 + t525 * t441;
t433 = -mrSges(2,1) * t522 + mrSges(2,3) * t510 - pkin(1) * t439 - t557;
t432 = mrSges(2,2) * t522 - mrSges(2,3) * t509 - pkin(5) * t439 - t528 * t434 + t530 * t435;
t1 = [-m(1) * g(1) + t539; -m(1) * g(2) + t547; -m(1) * g(3) + t540; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t547 + t526 * t432 - t524 * t433; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t539 + t524 * t432 + t526 * t433; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t509 - mrSges(2,2) * t510 + t528 * t435 + t530 * t434 + pkin(1) * (m(3) * t509 - t455) + pkin(5) * t538; t540; t557; t455; t558; t468;];
tauJB = t1;
