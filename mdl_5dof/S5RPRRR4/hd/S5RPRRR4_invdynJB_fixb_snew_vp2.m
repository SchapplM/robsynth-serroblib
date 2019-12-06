% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:31
% EndTime: 2019-12-05 18:14:34
% DurationCPUTime: 2.47s
% Computational Cost: add. (40644->185), mult. (52826->233), div. (0->0), fcn. (28272->10), ass. (0->86)
t505 = sin(qJ(1));
t509 = cos(qJ(1));
t482 = t509 * g(2) + t505 * g(3);
t478 = qJDD(1) * pkin(1) + t482;
t481 = t505 * g(2) - t509 * g(3);
t510 = qJD(1) ^ 2;
t479 = -t510 * pkin(1) + t481;
t500 = sin(pkin(9));
t501 = cos(pkin(9));
t463 = t501 * t478 - t500 * t479;
t460 = qJDD(1) * pkin(2) + t463;
t464 = t500 * t478 + t501 * t479;
t461 = -t510 * pkin(2) + t464;
t504 = sin(qJ(3));
t508 = cos(qJ(3));
t455 = t508 * t460 - t504 * t461;
t495 = qJDD(1) + qJDD(3);
t452 = t495 * pkin(3) + t455;
t456 = t504 * t460 + t508 * t461;
t496 = qJD(1) + qJD(3);
t494 = t496 ^ 2;
t453 = -t494 * pkin(3) + t456;
t503 = sin(qJ(4));
t507 = cos(qJ(4));
t449 = t503 * t452 + t507 * t453;
t487 = qJD(4) + t496;
t485 = t487 ^ 2;
t486 = qJDD(4) + t495;
t446 = -t485 * pkin(4) + t486 * pkin(8) + t449;
t499 = -g(1) + qJDD(2);
t502 = sin(qJ(5));
t506 = cos(qJ(5));
t443 = -t502 * t446 + t506 * t499;
t444 = t506 * t446 + t502 * t499;
t466 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t502 + Ifges(6,2) * t506) * t487;
t467 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t502 + Ifges(6,4) * t506) * t487;
t525 = qJD(5) * t487;
t471 = t502 * t486 + t506 * t525;
t472 = t506 * t486 - t502 * t525;
t528 = mrSges(6,1) * t443 - mrSges(6,2) * t444 + Ifges(6,5) * t471 + Ifges(6,6) * t472 + Ifges(6,3) * qJDD(5) + (t466 * t502 - t467 * t506) * t487;
t527 = t487 * t502;
t526 = t487 * t506;
t470 = (-mrSges(6,1) * t506 + mrSges(6,2) * t502) * t487;
t477 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t526;
t441 = m(6) * t443 + qJDD(5) * mrSges(6,1) - t471 * mrSges(6,3) + qJD(5) * t477 - t470 * t527;
t476 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t527;
t442 = m(6) * t444 - qJDD(5) * mrSges(6,2) + t472 * mrSges(6,3) - qJD(5) * t476 + t470 * t526;
t519 = -t502 * t441 + t506 * t442;
t427 = m(5) * t449 - t485 * mrSges(5,1) - t486 * mrSges(5,2) + t519;
t448 = t507 * t452 - t503 * t453;
t445 = -t486 * pkin(4) - t485 * pkin(8) - t448;
t514 = -m(6) * t445 + t472 * mrSges(6,1) - t471 * mrSges(6,2) - t476 * t527 + t477 * t526;
t436 = m(5) * t448 + t486 * mrSges(5,1) - t485 * mrSges(5,2) + t514;
t424 = t503 * t427 + t507 * t436;
t420 = m(4) * t455 + t495 * mrSges(4,1) - t494 * mrSges(4,2) + t424;
t520 = t507 * t427 - t503 * t436;
t421 = m(4) * t456 - t494 * mrSges(4,1) - t495 * mrSges(4,2) + t520;
t415 = t508 * t420 + t504 * t421;
t412 = m(3) * t463 + qJDD(1) * mrSges(3,1) - t510 * mrSges(3,2) + t415;
t521 = -t504 * t420 + t508 * t421;
t413 = m(3) * t464 - t510 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t521;
t406 = t501 * t412 + t500 * t413;
t430 = t506 * t441 + t502 * t442;
t524 = m(5) * t499 + t430;
t522 = -t500 * t412 + t501 * t413;
t403 = m(2) * t481 - t510 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t522;
t404 = m(2) * t482 + qJDD(1) * mrSges(2,1) - t510 * mrSges(2,2) + t406;
t523 = t509 * t403 - t505 * t404;
t518 = m(4) * t499 + t524;
t428 = m(3) * t499 + t518;
t517 = -t505 * t403 - t509 * t404;
t465 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t502 + Ifges(6,6) * t506) * t487;
t433 = -mrSges(6,1) * t445 + mrSges(6,3) * t444 + Ifges(6,4) * t471 + Ifges(6,2) * t472 + Ifges(6,6) * qJDD(5) + qJD(5) * t467 - t465 * t527;
t434 = mrSges(6,2) * t445 - mrSges(6,3) * t443 + Ifges(6,1) * t471 + Ifges(6,4) * t472 + Ifges(6,5) * qJDD(5) - qJD(5) * t466 + t465 * t526;
t515 = mrSges(5,1) * t448 - mrSges(5,2) * t449 + Ifges(5,3) * t486 + pkin(4) * t514 + pkin(8) * t519 + t506 * t433 + t502 * t434;
t512 = mrSges(4,1) * t455 - mrSges(4,2) * t456 + Ifges(4,3) * t495 + pkin(3) * t424 + t515;
t511 = mrSges(2,1) * t482 + mrSges(3,1) * t463 - mrSges(2,2) * t481 - mrSges(3,2) * t464 + pkin(1) * t406 + pkin(2) * t415 + t512 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t422 = -mrSges(5,1) * t499 + mrSges(5,3) * t449 + t485 * Ifges(5,5) + Ifges(5,6) * t486 - pkin(4) * t430 - t528;
t416 = mrSges(5,2) * t499 - mrSges(5,3) * t448 + Ifges(5,5) * t486 - t485 * Ifges(5,6) - pkin(8) * t430 - t502 * t433 + t506 * t434;
t408 = mrSges(4,2) * t499 - mrSges(4,3) * t455 + Ifges(4,5) * t495 - t494 * Ifges(4,6) - pkin(7) * t424 + t507 * t416 - t503 * t422;
t407 = -mrSges(4,1) * t499 + mrSges(4,3) * t456 + t494 * Ifges(4,5) + Ifges(4,6) * t495 - pkin(3) * t524 + pkin(7) * t520 + t503 * t416 + t507 * t422;
t401 = mrSges(3,2) * t499 - mrSges(3,3) * t463 + Ifges(3,5) * qJDD(1) - t510 * Ifges(3,6) - pkin(6) * t415 - t504 * t407 + t508 * t408;
t400 = -mrSges(3,1) * t499 + mrSges(3,3) * t464 + t510 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t518 + pkin(6) * t521 + t508 * t407 + t504 * t408;
t399 = -mrSges(2,2) * g(1) - mrSges(2,3) * t482 + Ifges(2,5) * qJDD(1) - t510 * Ifges(2,6) - qJ(2) * t406 - t500 * t400 + t501 * t401;
t398 = mrSges(2,1) * g(1) + mrSges(2,3) * t481 + t510 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t428 + qJ(2) * t522 + t501 * t400 + t500 * t401;
t1 = [(-m(1) - m(2)) * g(1) + t428; -m(1) * g(2) + t517; -m(1) * g(3) + t523; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t511; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t523 - t509 * t398 - t505 * t399; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t517 - t505 * t398 + t509 * t399; t511; t428; t512; t515; t528;];
tauJB = t1;
