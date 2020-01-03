% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR8
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:29
% EndTime: 2019-12-31 17:42:30
% DurationCPUTime: 1.73s
% Computational Cost: add. (24639->175), mult. (30899->222), div. (0->0), fcn. (18662->10), ass. (0->81)
t464 = sin(pkin(8));
t488 = cos(pkin(8));
t452 = -t488 * g(1) - t464 * g(2);
t462 = -g(3) + qJDD(1);
t468 = sin(qJ(2));
t471 = cos(qJ(2));
t436 = -t468 * t452 + t471 * t462;
t434 = qJDD(2) * pkin(2) + t436;
t437 = t471 * t452 + t468 * t462;
t472 = qJD(2) ^ 2;
t435 = -t472 * pkin(2) + t437;
t467 = sin(qJ(3));
t470 = cos(qJ(3));
t429 = t470 * t434 - t467 * t435;
t460 = qJDD(2) + qJDD(3);
t426 = t460 * pkin(3) + t429;
t430 = t467 * t434 + t470 * t435;
t461 = qJD(2) + qJD(3);
t459 = t461 ^ 2;
t427 = -t459 * pkin(3) + t430;
t463 = sin(pkin(9));
t465 = cos(pkin(9));
t423 = t463 * t426 + t465 * t427;
t420 = -t459 * pkin(4) + t460 * pkin(7) + t423;
t451 = t464 * g(1) - t488 * g(2);
t450 = qJDD(4) - t451;
t466 = sin(qJ(5));
t469 = cos(qJ(5));
t417 = -t466 * t420 + t469 * t450;
t418 = t469 * t420 + t466 * t450;
t439 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t466 + Ifges(6,2) * t469) * t461;
t440 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t466 + Ifges(6,4) * t469) * t461;
t484 = qJD(5) * t461;
t444 = t466 * t460 + t469 * t484;
t445 = t469 * t460 - t466 * t484;
t490 = mrSges(6,1) * t417 - mrSges(6,2) * t418 + Ifges(6,5) * t444 + Ifges(6,6) * t445 + Ifges(6,3) * qJDD(5) + (t439 * t466 - t440 * t469) * t461;
t489 = m(3) + m(4);
t487 = t461 * t466;
t486 = t461 * t469;
t443 = (-mrSges(6,1) * t469 + mrSges(6,2) * t466) * t461;
t449 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t486;
t413 = m(6) * t417 + qJDD(5) * mrSges(6,1) - t444 * mrSges(6,3) + qJD(5) * t449 - t443 * t487;
t448 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t487;
t414 = m(6) * t418 - qJDD(5) * mrSges(6,2) + t445 * mrSges(6,3) - qJD(5) * t448 + t443 * t486;
t478 = -t466 * t413 + t469 * t414;
t397 = m(5) * t423 - t459 * mrSges(5,1) - t460 * mrSges(5,2) + t478;
t422 = t465 * t426 - t463 * t427;
t419 = -t460 * pkin(4) - t459 * pkin(7) - t422;
t476 = -m(6) * t419 + t445 * mrSges(6,1) - t444 * mrSges(6,2) - t448 * t487 + t449 * t486;
t409 = m(5) * t422 + t460 * mrSges(5,1) - t459 * mrSges(5,2) + t476;
t394 = t463 * t397 + t465 * t409;
t390 = m(4) * t429 + t460 * mrSges(4,1) - t459 * mrSges(4,2) + t394;
t479 = t465 * t397 - t463 * t409;
t391 = m(4) * t430 - t459 * mrSges(4,1) - t460 * mrSges(4,2) + t479;
t385 = t470 * t390 + t467 * t391;
t383 = m(3) * t436 + qJDD(2) * mrSges(3,1) - t472 * mrSges(3,2) + t385;
t480 = -t467 * t390 + t470 * t391;
t384 = m(3) * t437 - t472 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t480;
t481 = -t468 * t383 + t471 * t384;
t376 = m(2) * t452 + t481;
t402 = t469 * t413 + t466 * t414;
t400 = m(5) * t450 + t402;
t399 = (m(2) + t489) * t451 - t400;
t485 = t464 * t376 + t488 * t399;
t377 = t471 * t383 + t468 * t384;
t483 = m(2) * t462 + t377;
t482 = t488 * t376 - t464 * t399;
t438 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t466 + Ifges(6,6) * t469) * t461;
t406 = -mrSges(6,1) * t419 + mrSges(6,3) * t418 + Ifges(6,4) * t444 + Ifges(6,2) * t445 + Ifges(6,6) * qJDD(5) + qJD(5) * t440 - t438 * t487;
t407 = mrSges(6,2) * t419 - mrSges(6,3) * t417 + Ifges(6,1) * t444 + Ifges(6,4) * t445 + Ifges(6,5) * qJDD(5) - qJD(5) * t439 + t438 * t486;
t474 = mrSges(4,1) * t429 + mrSges(5,1) * t422 - mrSges(4,2) * t430 - mrSges(5,2) * t423 + pkin(3) * t394 + pkin(4) * t476 + pkin(7) * t478 + t469 * t406 + t466 * t407 + (Ifges(5,3) + Ifges(4,3)) * t460;
t473 = mrSges(3,1) * t436 - mrSges(3,2) * t437 + Ifges(3,3) * qJDD(2) + pkin(2) * t385 + t474;
t392 = -mrSges(5,1) * t450 + mrSges(5,3) * t423 + t459 * Ifges(5,5) + Ifges(5,6) * t460 - pkin(4) * t402 - t490;
t386 = mrSges(5,2) * t450 - mrSges(5,3) * t422 + Ifges(5,5) * t460 - t459 * Ifges(5,6) - pkin(7) * t402 - t466 * t406 + t469 * t407;
t379 = -mrSges(4,2) * t451 - mrSges(4,3) * t429 + Ifges(4,5) * t460 - t459 * Ifges(4,6) - qJ(4) * t394 + t465 * t386 - t463 * t392;
t378 = mrSges(4,1) * t451 + mrSges(4,3) * t430 + t459 * Ifges(4,5) + Ifges(4,6) * t460 - pkin(3) * t400 + qJ(4) * t479 + t463 * t386 + t465 * t392;
t373 = -mrSges(2,1) * t462 + mrSges(2,3) * t452 - pkin(1) * t377 - t473;
t372 = -mrSges(3,2) * t451 - mrSges(3,3) * t436 + Ifges(3,5) * qJDD(2) - t472 * Ifges(3,6) - pkin(6) * t385 - t467 * t378 + t470 * t379;
t371 = Ifges(3,6) * qJDD(2) + t472 * Ifges(3,5) + mrSges(3,1) * t451 + mrSges(3,3) * t437 + t467 * t379 + t470 * t378 - pkin(2) * (-m(4) * t451 + t400) + pkin(6) * t480;
t370 = mrSges(2,2) * t462 - mrSges(2,3) * t451 - pkin(5) * t377 - t468 * t371 + t471 * t372;
t1 = [-m(1) * g(1) + t482; -m(1) * g(2) + t485; -m(1) * g(3) + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t485 + t488 * t370 - t464 * t373; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t482 + t464 * t370 + t488 * t373; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t452 + t468 * t372 + t471 * t371 - pkin(1) * t400 + pkin(5) * t481 + (pkin(1) * t489 + mrSges(2,1)) * t451; t483; t473; t474; t400; t490;];
tauJB = t1;
