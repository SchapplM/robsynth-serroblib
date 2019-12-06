% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:27
% EndTime: 2019-12-05 16:06:30
% DurationCPUTime: 2.24s
% Computational Cost: add. (19182->235), mult. (41262->292), div. (0->0), fcn. (25110->8), ass. (0->96)
t511 = Ifges(5,1) + Ifges(6,1);
t506 = Ifges(5,4) - Ifges(6,5);
t505 = Ifges(5,5) + Ifges(6,4);
t510 = Ifges(5,2) + Ifges(6,3);
t509 = -Ifges(6,2) - Ifges(5,3);
t504 = Ifges(5,6) - Ifges(6,6);
t508 = -2 * qJD(4);
t507 = -mrSges(5,3) - mrSges(6,2);
t503 = cos(pkin(8));
t475 = sin(pkin(7));
t476 = cos(pkin(7));
t461 = g(1) * t475 - g(2) * t476;
t462 = -g(1) * t476 - g(2) * t475;
t478 = sin(qJ(2));
t480 = cos(qJ(2));
t438 = t478 * t461 + t480 * t462;
t482 = qJD(2) ^ 2;
t434 = -pkin(2) * t482 + qJDD(2) * pkin(6) + t438;
t473 = -g(3) + qJDD(1);
t477 = sin(qJ(3));
t479 = cos(qJ(3));
t415 = -t434 * t477 + t479 * t473;
t495 = qJD(2) * qJD(3);
t492 = t479 * t495;
t459 = qJDD(2) * t477 + t492;
t412 = (-t459 + t492) * qJ(4) + (t477 * t479 * t482 + qJDD(3)) * pkin(3) + t415;
t416 = t479 * t434 + t477 * t473;
t460 = qJDD(2) * t479 - t477 * t495;
t497 = qJD(2) * t477;
t463 = qJD(3) * pkin(3) - qJ(4) * t497;
t472 = t479 ^ 2;
t413 = -pkin(3) * t472 * t482 + qJ(4) * t460 - qJD(3) * t463 + t416;
t474 = sin(pkin(8));
t496 = qJD(2) * t479;
t446 = t474 * t497 - t503 * t496;
t409 = t474 * t412 + t503 * t413 + t446 * t508;
t435 = t459 * t474 - t503 * t460;
t447 = (t474 * t479 + t503 * t477) * qJD(2);
t442 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t447;
t426 = pkin(4) * t446 - qJ(5) * t447;
t481 = qJD(3) ^ 2;
t404 = -pkin(4) * t481 + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t426 * t446 + t409;
t443 = -qJD(3) * mrSges(6,1) + mrSges(6,2) * t447;
t494 = m(6) * t404 + qJDD(3) * mrSges(6,3) + qJD(3) * t443;
t427 = mrSges(6,1) * t446 - mrSges(6,3) * t447;
t498 = -mrSges(5,1) * t446 - mrSges(5,2) * t447 - t427;
t400 = m(5) * t409 - qJDD(3) * mrSges(5,2) - qJD(3) * t442 + t507 * t435 + t498 * t446 + t494;
t485 = t503 * t412 - t474 * t413;
t408 = t447 * t508 + t485;
t436 = t503 * t459 + t474 * t460;
t441 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t446;
t405 = -qJDD(3) * pkin(4) - t481 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t426) * t447 - t485;
t444 = -mrSges(6,2) * t446 + qJD(3) * mrSges(6,3);
t487 = -m(6) * t405 + qJDD(3) * mrSges(6,1) + qJD(3) * t444;
t401 = m(5) * t408 + qJDD(3) * mrSges(5,1) + qJD(3) * t441 + t507 * t436 + t498 * t447 + t487;
t394 = t474 * t400 + t503 * t401;
t458 = (-mrSges(4,1) * t479 + mrSges(4,2) * t477) * qJD(2);
t465 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t496;
t390 = m(4) * t415 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t459 + qJD(3) * t465 - t458 * t497 + t394;
t464 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t497;
t488 = t503 * t400 - t401 * t474;
t391 = m(4) * t416 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t460 - qJD(3) * t464 + t458 * t496 + t488;
t489 = -t390 * t477 + t479 * t391;
t385 = m(3) * t438 - mrSges(3,1) * t482 - qJDD(2) * mrSges(3,2) + t489;
t437 = t461 * t480 - t478 * t462;
t486 = -qJDD(2) * pkin(2) - t437;
t433 = -pkin(6) * t482 + t486;
t414 = -pkin(3) * t460 + qJDD(4) + t463 * t497 + (-qJ(4) * t472 - pkin(6)) * t482 + t486;
t407 = -0.2e1 * qJD(5) * t447 + (qJD(3) * t446 - t436) * qJ(5) + (qJD(3) * t447 + t435) * pkin(4) + t414;
t402 = m(6) * t407 + t435 * mrSges(6,1) - t436 * mrSges(6,3) - t447 * t443 + t446 * t444;
t484 = m(5) * t414 + t435 * mrSges(5,1) + mrSges(5,2) * t436 + t446 * t441 + t442 * t447 + t402;
t483 = -m(4) * t433 + t460 * mrSges(4,1) - mrSges(4,2) * t459 - t464 * t497 + t465 * t496 - t484;
t396 = m(3) * t437 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t482 + t483;
t382 = t478 * t385 + t480 * t396;
t380 = m(2) * t461 + t382;
t490 = t480 * t385 - t478 * t396;
t381 = m(2) * t462 + t490;
t502 = t476 * t380 + t475 * t381;
t386 = t479 * t390 + t477 * t391;
t501 = -t504 * qJD(3) + t510 * t446 - t506 * t447;
t500 = t509 * qJD(3) + t504 * t446 - t505 * t447;
t499 = t505 * qJD(3) - t506 * t446 + t511 * t447;
t493 = m(3) * t473 + t386;
t491 = -t380 * t475 + t476 * t381;
t450 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t477 + Ifges(4,4) * t479) * qJD(2);
t449 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t477 + Ifges(4,2) * t479) * qJD(2);
t448 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t477 + Ifges(4,6) * t479) * qJD(2);
t393 = mrSges(5,2) * t414 + mrSges(6,2) * t405 - mrSges(5,3) * t408 - mrSges(6,3) * t407 - qJ(5) * t402 + t501 * qJD(3) + t505 * qJDD(3) - t506 * t435 + t511 * t436 + t500 * t446;
t392 = -mrSges(5,1) * t414 - mrSges(6,1) * t407 + mrSges(6,2) * t404 + mrSges(5,3) * t409 - pkin(4) * t402 + t499 * qJD(3) + t504 * qJDD(3) - t510 * t435 + t506 * t436 + t500 * t447;
t376 = mrSges(4,2) * t433 - mrSges(4,3) * t415 + Ifges(4,1) * t459 + Ifges(4,4) * t460 + Ifges(4,5) * qJDD(3) - qJ(4) * t394 - qJD(3) * t449 - t474 * t392 + t503 * t393 + t448 * t496;
t375 = -mrSges(4,1) * t433 + mrSges(4,3) * t416 + Ifges(4,4) * t459 + Ifges(4,2) * t460 + Ifges(4,6) * qJDD(3) - pkin(3) * t484 + qJ(4) * t488 + qJD(3) * t450 + t503 * t392 + t474 * t393 - t448 * t497;
t374 = Ifges(3,6) * qJDD(2) + t482 * Ifges(3,5) - mrSges(3,1) * t473 - Ifges(4,6) * t460 - qJ(5) * t494 - pkin(4) * t487 - Ifges(4,5) * t459 + mrSges(3,3) * t438 - mrSges(5,1) * t408 + mrSges(5,2) * t409 - mrSges(4,1) * t415 + mrSges(4,2) * t416 - mrSges(6,3) * t404 + mrSges(6,1) * t405 - pkin(3) * t394 - pkin(2) * t386 + (pkin(4) * t427 + t501) * t447 + (qJ(5) * t427 - t499) * t446 + (pkin(4) * mrSges(6,2) - t505) * t436 + (qJ(5) * mrSges(6,2) + t504) * t435 + (-t449 * t477 + t450 * t479) * qJD(2) + (-Ifges(4,3) + t509) * qJDD(3);
t373 = mrSges(3,2) * t473 - mrSges(3,3) * t437 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t482 - pkin(6) * t386 - t375 * t477 + t376 * t479;
t372 = mrSges(2,2) * t473 - mrSges(2,3) * t461 - pkin(5) * t382 + t373 * t480 - t374 * t478;
t371 = -mrSges(2,1) * t473 + mrSges(2,3) * t462 - pkin(1) * t493 + pkin(5) * t490 + t478 * t373 + t480 * t374;
t1 = [-m(1) * g(1) + t491; -m(1) * g(2) + t502; -m(1) * g(3) + m(2) * t473 + t493; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t502 - t475 * t371 + t476 * t372; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t491 + t476 * t371 + t475 * t372; -mrSges(1,1) * g(2) + mrSges(2,1) * t461 + mrSges(3,1) * t437 + mrSges(1,2) * g(1) - mrSges(2,2) * t462 - mrSges(3,2) * t438 + Ifges(3,3) * qJDD(2) + pkin(1) * t382 + pkin(2) * t483 + pkin(6) * t489 + t479 * t375 + t477 * t376;];
tauB = t1;
