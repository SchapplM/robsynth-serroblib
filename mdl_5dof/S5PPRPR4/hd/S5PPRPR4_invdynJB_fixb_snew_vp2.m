% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:22
% EndTime: 2019-12-31 17:32:23
% DurationCPUTime: 1.10s
% Computational Cost: add. (10471->177), mult. (21224->222), div. (0->0), fcn. (13455->8), ass. (0->85)
t462 = qJD(3) ^ 2;
t454 = sin(pkin(8));
t488 = t454 ^ 2;
t456 = cos(pkin(8));
t487 = pkin(4) * t456;
t486 = -mrSges(2,2) + mrSges(3,3);
t485 = mrSges(5,2) * t454;
t452 = t456 ^ 2;
t484 = t452 * t462;
t455 = sin(pkin(7));
t457 = cos(pkin(7));
t441 = g(1) * t455 - g(2) * t457;
t439 = qJDD(2) - t441;
t442 = -g(1) * t457 - g(2) * t455;
t459 = sin(qJ(3));
t461 = cos(qJ(3));
t429 = t459 * t439 + t461 * t442;
t425 = -pkin(3) * t462 + qJDD(3) * qJ(4) + t429;
t453 = g(3) - qJDD(1);
t480 = qJD(3) * qJD(4);
t482 = t456 * t453 - 0.2e1 * t454 * t480;
t412 = (-pkin(6) * qJDD(3) + t462 * t487 - t425) * t454 + t482;
t416 = t454 * t453 + (t425 + 0.2e1 * t480) * t456;
t479 = qJDD(3) * t456;
t413 = -pkin(4) * t484 + pkin(6) * t479 + t416;
t458 = sin(qJ(5));
t460 = cos(qJ(5));
t410 = t412 * t460 - t413 * t458;
t469 = -t454 * t458 + t456 * t460;
t432 = t469 * qJD(3);
t470 = t454 * t460 + t456 * t458;
t433 = t470 * qJD(3);
t423 = -mrSges(6,1) * t432 + mrSges(6,2) * t433;
t427 = t432 * qJD(5) + t470 * qJDD(3);
t430 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t432;
t408 = m(6) * t410 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t427 + qJD(5) * t430 - t423 * t433;
t411 = t412 * t458 + t413 * t460;
t426 = -qJD(5) * t433 + t469 * qJDD(3);
t431 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t433;
t409 = m(6) * t411 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t426 - qJD(5) * t431 + t423 * t432;
t400 = t460 * t408 + t458 * t409;
t415 = -t425 * t454 + t482;
t468 = mrSges(5,3) * qJDD(3) + t462 * (-mrSges(5,1) * t456 + t485);
t398 = m(5) * t415 - t468 * t454 + t400;
t476 = -t408 * t458 + t460 * t409;
t399 = m(5) * t416 + t468 * t456 + t476;
t396 = -t398 * t454 + t456 * t399;
t393 = m(4) * t429 - mrSges(4,1) * t462 - qJDD(3) * mrSges(4,2) + t396;
t428 = t439 * t461 - t459 * t442;
t471 = qJDD(4) - t428;
t422 = -qJDD(3) * pkin(3) - qJ(4) * t462 + t471;
t414 = (-pkin(3) - t487) * qJDD(3) + (-qJ(4) + (-t452 - t488) * pkin(6)) * t462 + t471;
t467 = m(6) * t414 - t426 * mrSges(6,1) + mrSges(6,2) * t427 - t432 * t430 + t431 * t433;
t465 = -m(5) * t422 + mrSges(5,1) * t479 - t467 + (t462 * t488 + t484) * mrSges(5,3);
t403 = m(4) * t428 - mrSges(4,2) * t462 + (mrSges(4,1) - t485) * qJDD(3) + t465;
t390 = t393 * t459 + t403 * t461;
t389 = m(3) * t439 + t390;
t387 = m(2) * t441 - t389;
t477 = t461 * t393 - t459 * t403;
t475 = m(3) * t442 + t477;
t388 = m(2) * t442 + t475;
t483 = t457 * t387 + t455 * t388;
t472 = Ifges(5,5) * t454 + Ifges(5,6) * t456;
t481 = t462 * t472;
t478 = -t387 * t455 + t457 * t388;
t474 = Ifges(5,1) * t454 + Ifges(5,4) * t456;
t473 = Ifges(5,4) * t454 + Ifges(5,2) * t456;
t395 = t398 * t456 + t399 * t454;
t394 = -t395 + (-m(3) - m(4)) * t453;
t466 = -m(2) * t453 + t394;
t417 = Ifges(6,5) * t433 + Ifges(6,6) * t432 + Ifges(6,3) * qJD(5);
t419 = Ifges(6,1) * t433 + Ifges(6,4) * t432 + Ifges(6,5) * qJD(5);
t401 = -mrSges(6,1) * t414 + mrSges(6,3) * t411 + Ifges(6,4) * t427 + Ifges(6,2) * t426 + Ifges(6,6) * qJDD(5) + qJD(5) * t419 - t417 * t433;
t418 = Ifges(6,4) * t433 + Ifges(6,2) * t432 + Ifges(6,6) * qJD(5);
t402 = mrSges(6,2) * t414 - mrSges(6,3) * t410 + Ifges(6,1) * t427 + Ifges(6,4) * t426 + Ifges(6,5) * qJDD(5) - qJD(5) * t418 + t417 * t432;
t383 = -mrSges(5,1) * t422 + mrSges(5,3) * t416 - pkin(4) * t467 + pkin(6) * t476 + t473 * qJDD(3) + t460 * t401 + t458 * t402 - t454 * t481;
t391 = mrSges(5,2) * t422 - mrSges(5,3) * t415 - pkin(6) * t400 + t474 * qJDD(3) - t401 * t458 + t402 * t460 + t456 * t481;
t404 = qJDD(3) * t485 - t465;
t464 = mrSges(4,1) * t428 - mrSges(4,2) * t429 + Ifges(4,3) * qJDD(3) - pkin(3) * t404 + qJ(4) * t396 + t456 * t383 + t454 * t391;
t463 = mrSges(6,1) * t410 - mrSges(6,2) * t411 + Ifges(6,5) * t427 + Ifges(6,6) * t426 + Ifges(6,3) * qJDD(5) + t433 * t418 - t432 * t419;
t382 = -mrSges(4,1) * t453 - mrSges(5,1) * t415 + mrSges(5,2) * t416 + mrSges(4,3) * t429 - pkin(3) * t395 - pkin(4) * t400 + (Ifges(4,6) - t472) * qJDD(3) - t463 + (-t454 * t473 + t456 * t474 + Ifges(4,5)) * t462;
t381 = mrSges(4,2) * t453 - mrSges(4,3) * t428 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t462 - qJ(4) * t395 - t383 * t454 + t391 * t456;
t380 = mrSges(3,2) * t439 - mrSges(2,3) * t441 - pkin(5) * t390 - qJ(2) * t394 + t381 * t461 - t382 * t459 + t486 * t453;
t379 = -t459 * t381 - t461 * t382 + pkin(2) * t395 - pkin(5) * t477 - pkin(1) * t394 + (m(4) * pkin(2) + mrSges(2,1) + mrSges(3,1)) * t453 + (mrSges(3,2) + mrSges(2,3)) * t442;
t1 = [-m(1) * g(1) + t478; -m(1) * g(2) + t483; -m(1) * g(3) + t466; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t483 - t455 * t379 + t457 * t380; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t478 + t457 * t379 + t455 * t380; -mrSges(1,1) * g(2) + mrSges(2,1) * t441 - mrSges(3,1) * t439 + mrSges(1,2) * g(1) - pkin(1) * t389 - pkin(2) * t390 + qJ(2) * t475 + t486 * t442 - t464; t466; t389; t464; t404; t463;];
tauJB = t1;
