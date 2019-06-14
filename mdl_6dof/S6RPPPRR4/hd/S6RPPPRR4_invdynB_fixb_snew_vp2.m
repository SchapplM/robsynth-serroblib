% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:45:57
% EndTime: 2019-05-05 13:45:59
% DurationCPUTime: 1.64s
% Computational Cost: add. (17914->253), mult. (30357->292), div. (0->0), fcn. (13102->8), ass. (0->102)
t474 = qJD(1) ^ 2;
t469 = sin(qJ(1));
t472 = cos(qJ(1));
t450 = -t472 * g(1) - t469 * g(2);
t480 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t450;
t504 = -pkin(1) - pkin(2);
t431 = t504 * t474 + t480;
t449 = t469 * g(1) - t472 * g(2);
t479 = -t474 * qJ(2) + qJDD(2) - t449;
t433 = t504 * qJDD(1) + t479;
t465 = sin(pkin(9));
t466 = cos(pkin(9));
t417 = t466 * t431 + t465 * t433;
t506 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t417;
t468 = sin(qJ(5));
t493 = t468 * qJD(1);
t447 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t493;
t471 = cos(qJ(5));
t494 = qJD(1) * t471;
t448 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t494;
t505 = (t447 * t468 + t448 * t471) * qJD(1);
t503 = mrSges(2,1) + mrSges(3,1);
t502 = mrSges(4,2) - mrSges(5,3);
t501 = Ifges(3,4) + Ifges(2,5);
t500 = Ifges(5,4) - Ifges(4,5);
t499 = Ifges(5,5) - Ifges(4,6);
t498 = Ifges(2,6) - Ifges(3,6);
t460 = g(3) + qJDD(3);
t497 = t468 * t460;
t434 = -t474 * pkin(1) + t480;
t416 = -t465 * t431 + t466 * t433;
t415 = qJDD(1) * pkin(3) - t474 * qJ(4) + qJDD(4) - t416;
t413 = qJDD(1) * pkin(7) + t415;
t409 = t468 * t413 + t471 * t460;
t443 = (-mrSges(6,1) * t468 - mrSges(6,2) * t471) * qJD(1);
t492 = qJD(1) * qJD(5);
t489 = t471 * t492;
t445 = t468 * qJDD(1) + t489;
t412 = (-pkin(3) - pkin(7)) * t474 + t506;
t490 = t468 * t492;
t446 = -t471 * qJDD(1) + t490;
t405 = (-t446 - t490) * pkin(8) + (-t445 - t489) * pkin(5) + t412;
t444 = (-pkin(5) * t468 + pkin(8) * t471) * qJD(1);
t473 = qJD(5) ^ 2;
t407 = -t473 * pkin(5) + qJDD(5) * pkin(8) + t444 * t493 + t409;
t467 = sin(qJ(6));
t470 = cos(qJ(6));
t403 = t470 * t405 - t467 * t407;
t441 = t470 * qJD(5) + t467 * t494;
t424 = t441 * qJD(6) + t467 * qJDD(5) + t470 * t446;
t442 = t467 * qJD(5) - t470 * t494;
t425 = -t441 * mrSges(7,1) + t442 * mrSges(7,2);
t451 = qJD(6) - t493;
t429 = -t451 * mrSges(7,2) + t441 * mrSges(7,3);
t440 = qJDD(6) - t445;
t401 = m(7) * t403 + t440 * mrSges(7,1) - t424 * mrSges(7,3) - t442 * t425 + t451 * t429;
t404 = t467 * t405 + t470 * t407;
t423 = -t442 * qJD(6) + t470 * qJDD(5) - t467 * t446;
t430 = t451 * mrSges(7,1) - t442 * mrSges(7,3);
t402 = m(7) * t404 - t440 * mrSges(7,2) + t423 * mrSges(7,3) + t441 * t425 - t451 * t430;
t486 = -t467 * t401 + t470 * t402;
t392 = m(6) * t409 - qJDD(5) * mrSges(6,2) + t445 * mrSges(6,3) - qJD(5) * t448 + t443 * t493 + t486;
t408 = t471 * t413 - t497;
t406 = -qJDD(5) * pkin(5) - t473 * pkin(8) + t497 + (-qJD(1) * t444 - t413) * t471;
t478 = -m(7) * t406 + t423 * mrSges(7,1) - t424 * mrSges(7,2) + t441 * t429 - t442 * t430;
t397 = m(6) * t408 + qJDD(5) * mrSges(6,1) - t446 * mrSges(6,3) + qJD(5) * t447 + t443 * t494 + t478;
t388 = t468 * t392 + t471 * t397;
t476 = -m(5) * t415 + qJDD(1) * mrSges(5,2) + t474 * mrSges(5,3) - t388;
t385 = m(4) * t416 - qJDD(1) * mrSges(4,1) - t474 * mrSges(4,2) + t476;
t414 = t474 * pkin(3) - t506;
t393 = t470 * t401 + t467 * t402;
t481 = -m(6) * t412 + t445 * mrSges(6,1) - t446 * mrSges(6,2) - t393;
t475 = -m(5) * t414 + t474 * mrSges(5,2) - t481;
t390 = m(4) * t417 - t474 * mrSges(4,1) + t502 * qJDD(1) + t475 - t505;
t487 = -t465 * t385 + t466 * t390;
t482 = m(3) * t434 + qJDD(1) * mrSges(3,3) + t487;
t380 = m(2) * t450 - qJDD(1) * mrSges(2,2) - t503 * t474 + t482;
t382 = t466 * t385 + t465 * t390;
t435 = -qJDD(1) * pkin(1) + t479;
t477 = -m(3) * t435 + qJDD(1) * mrSges(3,1) + t474 * mrSges(3,3) - t382;
t381 = m(2) * t449 + qJDD(1) * mrSges(2,1) - t474 * mrSges(2,2) + t477;
t496 = t469 * t380 + t472 * t381;
t495 = t471 * t392 - t468 * t397;
t488 = t472 * t380 - t469 * t381;
t483 = (-m(4) - m(5)) * t460 - t495;
t438 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t471 + Ifges(6,4) * t468) * qJD(1);
t437 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t471 + Ifges(6,2) * t468) * qJD(1);
t436 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t471 + Ifges(6,6) * t468) * qJD(1);
t420 = Ifges(7,1) * t442 + Ifges(7,4) * t441 + Ifges(7,5) * t451;
t419 = Ifges(7,4) * t442 + Ifges(7,2) * t441 + Ifges(7,6) * t451;
t418 = Ifges(7,5) * t442 + Ifges(7,6) * t441 + Ifges(7,3) * t451;
t395 = mrSges(7,2) * t406 - mrSges(7,3) * t403 + Ifges(7,1) * t424 + Ifges(7,4) * t423 + Ifges(7,5) * t440 + t441 * t418 - t451 * t419;
t394 = -mrSges(7,1) * t406 + mrSges(7,3) * t404 + Ifges(7,4) * t424 + Ifges(7,2) * t423 + Ifges(7,6) * t440 - t442 * t418 + t451 * t420;
t387 = m(5) * t460 + t495;
t386 = -m(3) * g(3) + t483;
t384 = -mrSges(6,1) * t412 - mrSges(7,1) * t403 + mrSges(7,2) * t404 + mrSges(6,3) * t409 + Ifges(6,4) * t446 - Ifges(7,5) * t424 + Ifges(6,2) * t445 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t423 - Ifges(7,3) * t440 - pkin(5) * t393 + qJD(5) * t438 - t442 * t419 + t441 * t420 + t436 * t494;
t383 = mrSges(6,2) * t412 - mrSges(6,3) * t408 + Ifges(6,1) * t446 + Ifges(6,4) * t445 + Ifges(6,5) * qJDD(5) - pkin(8) * t393 - qJD(5) * t437 - t467 * t394 + t470 * t395 + t436 * t493;
t376 = -qJ(4) * t387 - mrSges(4,3) * t416 + pkin(4) * t388 + mrSges(5,1) * t415 - mrSges(6,2) * t409 + t467 * t395 + t470 * t394 + pkin(5) * t478 + pkin(8) * t486 + Ifges(6,5) * t446 + Ifges(6,6) * t445 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t408 + t499 * t474 + t502 * t460 + t500 * qJDD(1) + (-t471 * t437 - t468 * t438) * qJD(1);
t375 = mrSges(4,3) * t417 - mrSges(5,1) * t414 - t468 * t383 - t471 * t384 - pkin(4) * (t481 + t505) - pkin(7) * t495 - pkin(3) * t387 - t500 * t474 + (-mrSges(4,1) + mrSges(5,2)) * t460 + t499 * qJDD(1);
t374 = mrSges(3,2) * t435 - mrSges(2,3) * t449 - qJ(2) * t386 - qJ(3) * t382 - t465 * t375 + t466 * t376 - t498 * t474 + t501 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t373 = mrSges(3,2) * t434 + mrSges(2,3) * t450 - pkin(1) * t386 - pkin(2) * t483 + t503 * g(3) - qJ(3) * t487 + t498 * qJDD(1) - t466 * t375 - t465 * t376 + t501 * t474;
t1 = [-m(1) * g(1) + t488; -m(1) * g(2) + t496; (-m(1) - m(2) - m(3)) * g(3) + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t496 - t469 * t373 + t472 * t374; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t488 + t472 * t373 + t469 * t374; qJ(2) * (-t474 * mrSges(3,1) + t482) + pkin(1) * t477 - mrSges(2,2) * t450 + mrSges(2,1) * t449 - pkin(2) * t382 + mrSges(3,3) * t434 - mrSges(3,1) * t435 - pkin(3) * t476 - qJ(4) * (-t447 * t493 - t448 * t494 + t475) + pkin(7) * t388 - mrSges(4,1) * t416 + mrSges(4,2) * t417 - t471 * t383 + t468 * t384 - mrSges(5,2) * t415 + mrSges(5,3) * t414 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (qJ(4) * mrSges(5,3) + Ifges(5,1) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB  = t1;
