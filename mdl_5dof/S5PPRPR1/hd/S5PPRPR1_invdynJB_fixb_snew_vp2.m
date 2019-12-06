% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:04
% EndTime: 2019-12-05 15:01:06
% DurationCPUTime: 1.94s
% Computational Cost: add. (20249->184), mult. (38994->238), div. (0->0), fcn. (26769->10), ass. (0->92)
t512 = qJD(3) ^ 2;
t504 = sin(pkin(7));
t507 = cos(pkin(7));
t493 = -t507 * g(1) - t504 * g(2);
t501 = -g(3) + qJDD(1);
t503 = sin(pkin(8));
t506 = cos(pkin(8));
t483 = -t503 * t493 + t506 * t501;
t484 = t506 * t493 + t503 * t501;
t509 = sin(qJ(3));
t511 = cos(qJ(3));
t469 = t509 * t483 + t511 * t484;
t467 = -t512 * pkin(3) + qJDD(3) * qJ(4) + t469;
t502 = sin(pkin(9));
t492 = t504 * g(1) - t507 * g(2);
t491 = qJDD(2) - t492;
t505 = cos(pkin(9));
t530 = qJD(3) * qJD(4);
t532 = t505 * t491 - 0.2e1 * t502 * t530;
t536 = pkin(4) * t505;
t460 = (-pkin(6) * qJDD(3) + t512 * t536 - t467) * t502 + t532;
t463 = t502 * t491 + (t467 + 0.2e1 * t530) * t505;
t529 = qJDD(3) * t505;
t500 = t505 ^ 2;
t534 = t500 * t512;
t461 = -pkin(4) * t534 + pkin(6) * t529 + t463;
t508 = sin(qJ(5));
t510 = cos(qJ(5));
t458 = t510 * t460 - t508 * t461;
t518 = -t502 * t508 + t505 * t510;
t485 = t518 * qJD(3);
t519 = t502 * t510 + t505 * t508;
t486 = t519 * qJD(3);
t474 = -t485 * mrSges(6,1) + t486 * mrSges(6,2);
t477 = t485 * qJD(5) + t519 * qJDD(3);
t481 = -qJD(5) * mrSges(6,2) + t485 * mrSges(6,3);
t455 = m(6) * t458 + qJDD(5) * mrSges(6,1) - t477 * mrSges(6,3) + qJD(5) * t481 - t486 * t474;
t459 = t508 * t460 + t510 * t461;
t476 = -t486 * qJD(5) + t518 * qJDD(3);
t482 = qJD(5) * mrSges(6,1) - t486 * mrSges(6,3);
t456 = m(6) * t459 - qJDD(5) * mrSges(6,2) + t476 * mrSges(6,3) - qJD(5) * t482 + t485 * t474;
t447 = t510 * t455 + t508 * t456;
t462 = -t502 * t467 + t532;
t535 = mrSges(5,2) * t502;
t517 = mrSges(5,3) * qJDD(3) + t512 * (-mrSges(5,1) * t505 + t535);
t445 = m(5) * t462 - t517 * t502 + t447;
t524 = -t508 * t455 + t510 * t456;
t446 = m(5) * t463 + t517 * t505 + t524;
t440 = t505 * t445 + t502 * t446;
t439 = (m(3) + m(4)) * t491 + t440;
t538 = t502 ^ 2;
t441 = -t502 * t445 + t505 * t446;
t436 = m(4) * t469 - t512 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t441;
t468 = t511 * t483 - t509 * t484;
t520 = qJDD(4) - t468;
t466 = -qJDD(3) * pkin(3) - t512 * qJ(4) + t520;
t464 = (-pkin(3) - t536) * qJDD(3) + (-qJ(4) + (-t500 - t538) * pkin(6)) * t512 + t520;
t516 = m(6) * t464 - t476 * mrSges(6,1) + t477 * mrSges(6,2) - t485 * t481 + t486 * t482;
t515 = -m(5) * t466 + mrSges(5,1) * t529 - t516 + (t512 * t538 + t534) * mrSges(5,3);
t451 = m(4) * t468 - t512 * mrSges(4,2) + (mrSges(4,1) - t535) * qJDD(3) + t515;
t432 = t509 * t436 + t511 * t451;
t429 = m(3) * t483 + t432;
t525 = t511 * t436 - t509 * t451;
t430 = m(3) * t484 + t525;
t526 = -t503 * t429 + t506 * t430;
t423 = m(2) * t493 + t526;
t438 = m(2) * t492 - t439;
t533 = t504 * t423 + t507 * t438;
t424 = t506 * t429 + t503 * t430;
t521 = Ifges(5,5) * t502 + Ifges(5,6) * t505;
t531 = t512 * t521;
t528 = m(2) * t501 + t424;
t527 = t507 * t423 - t504 * t438;
t523 = Ifges(5,1) * t502 + Ifges(5,4) * t505;
t522 = Ifges(5,4) * t502 + Ifges(5,2) * t505;
t470 = Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * qJD(5);
t472 = Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * qJD(5);
t448 = -mrSges(6,1) * t464 + mrSges(6,3) * t459 + Ifges(6,4) * t477 + Ifges(6,2) * t476 + Ifges(6,6) * qJDD(5) + qJD(5) * t472 - t486 * t470;
t471 = Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * qJD(5);
t449 = mrSges(6,2) * t464 - mrSges(6,3) * t458 + Ifges(6,1) * t477 + Ifges(6,4) * t476 + Ifges(6,5) * qJDD(5) - qJD(5) * t471 + t485 * t470;
t431 = -mrSges(5,1) * t466 + mrSges(5,3) * t463 - pkin(4) * t516 + pkin(6) * t524 + t522 * qJDD(3) + t510 * t448 + t508 * t449 - t502 * t531;
t433 = mrSges(5,2) * t466 - mrSges(5,3) * t462 - pkin(6) * t447 + t523 * qJDD(3) - t508 * t448 + t510 * t449 + t505 * t531;
t457 = qJDD(3) * t535 - t515;
t514 = mrSges(4,1) * t468 - mrSges(4,2) * t469 + Ifges(4,3) * qJDD(3) - pkin(3) * t457 + qJ(4) * t441 + t505 * t431 + t502 * t433;
t513 = mrSges(6,1) * t458 - mrSges(6,2) * t459 + Ifges(6,5) * t477 + Ifges(6,6) * t476 + Ifges(6,3) * qJDD(5) + t486 * t471 - t485 * t472;
t425 = -mrSges(4,1) * t491 - mrSges(5,1) * t462 + mrSges(5,2) * t463 + mrSges(4,3) * t469 - pkin(3) * t440 - pkin(4) * t447 + (Ifges(4,6) - t521) * qJDD(3) - t513 + (-t502 * t522 + t505 * t523 + Ifges(4,5)) * t512;
t420 = mrSges(4,2) * t491 - mrSges(4,3) * t468 + Ifges(4,5) * qJDD(3) - t512 * Ifges(4,6) - qJ(4) * t440 - t502 * t431 + t505 * t433;
t419 = mrSges(3,2) * t491 - mrSges(3,3) * t483 - pkin(5) * t432 + t511 * t420 - t509 * t425;
t418 = -mrSges(2,1) * t501 - mrSges(3,1) * t483 + mrSges(3,2) * t484 + mrSges(2,3) * t493 - pkin(1) * t424 - pkin(2) * t432 - t514;
t417 = -mrSges(3,1) * t491 + mrSges(3,3) * t484 + t509 * t420 + t511 * t425 - pkin(2) * (m(4) * t491 + t440) + pkin(5) * t525;
t416 = mrSges(2,2) * t501 - mrSges(2,3) * t492 - qJ(2) * t424 - t503 * t417 + t506 * t419;
t1 = [-m(1) * g(1) + t527; -m(1) * g(2) + t533; -m(1) * g(3) + t528; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t533 + t507 * t416 - t504 * t418; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t527 + t504 * t416 + t507 * t418; -mrSges(1,1) * g(2) + mrSges(2,1) * t492 + mrSges(1,2) * g(1) - mrSges(2,2) * t493 - pkin(1) * t439 + qJ(2) * t526 + t506 * t417 + t503 * t419; t528; t439; t514; t457; t513;];
tauJB = t1;
