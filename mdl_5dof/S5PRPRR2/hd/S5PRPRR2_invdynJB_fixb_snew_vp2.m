% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:58
% EndTime: 2019-12-05 15:45:01
% DurationCPUTime: 1.73s
% Computational Cost: add. (23107->175), mult. (30637->222), div. (0->0), fcn. (18502->10), ass. (0->81)
t458 = sin(pkin(8));
t460 = cos(pkin(8));
t447 = -t460 * g(1) - t458 * g(2);
t456 = -g(3) + qJDD(1);
t463 = sin(qJ(2));
t466 = cos(qJ(2));
t431 = -t463 * t447 + t466 * t456;
t429 = qJDD(2) * pkin(2) + t431;
t432 = t466 * t447 + t463 * t456;
t467 = qJD(2) ^ 2;
t430 = -t467 * pkin(2) + t432;
t457 = sin(pkin(9));
t459 = cos(pkin(9));
t424 = t459 * t429 - t457 * t430;
t422 = qJDD(2) * pkin(3) + t424;
t425 = t457 * t429 + t459 * t430;
t423 = -t467 * pkin(3) + t425;
t462 = sin(qJ(4));
t465 = cos(qJ(4));
t419 = t462 * t422 + t465 * t423;
t455 = qJD(2) + qJD(4);
t453 = t455 ^ 2;
t454 = qJDD(2) + qJDD(4);
t416 = -t453 * pkin(4) + t454 * pkin(7) + t419;
t446 = t458 * g(1) - t460 * g(2);
t445 = qJDD(3) - t446;
t461 = sin(qJ(5));
t464 = cos(qJ(5));
t413 = -t461 * t416 + t464 * t445;
t414 = t464 * t416 + t461 * t445;
t434 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t461 + Ifges(6,2) * t464) * t455;
t435 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t461 + Ifges(6,4) * t464) * t455;
t480 = qJD(5) * t455;
t439 = t461 * t454 + t464 * t480;
t440 = t464 * t454 - t461 * t480;
t486 = mrSges(6,1) * t413 - mrSges(6,2) * t414 + Ifges(6,5) * t439 + Ifges(6,6) * t440 + Ifges(6,3) * qJDD(5) + (t434 * t461 - t435 * t464) * t455;
t438 = (-mrSges(6,1) * t464 + mrSges(6,2) * t461) * t455;
t482 = t455 * t464;
t444 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t482;
t483 = t455 * t461;
t409 = m(6) * t413 + qJDD(5) * mrSges(6,1) - t439 * mrSges(6,3) + qJD(5) * t444 - t438 * t483;
t443 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t483;
t410 = m(6) * t414 - qJDD(5) * mrSges(6,2) + t440 * mrSges(6,3) - qJD(5) * t443 + t438 * t482;
t473 = -t461 * t409 + t464 * t410;
t393 = m(5) * t419 - t453 * mrSges(5,1) - t454 * mrSges(5,2) + t473;
t418 = t465 * t422 - t462 * t423;
t415 = -t454 * pkin(4) - t453 * pkin(7) - t418;
t470 = -m(6) * t415 + t440 * mrSges(6,1) - t439 * mrSges(6,2) - t443 * t483 + t444 * t482;
t405 = m(5) * t418 + t454 * mrSges(5,1) - t453 * mrSges(5,2) + t470;
t390 = t462 * t393 + t465 * t405;
t388 = m(4) * t424 + qJDD(2) * mrSges(4,1) - t467 * mrSges(4,2) + t390;
t474 = t465 * t393 - t462 * t405;
t389 = m(4) * t425 - t467 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t474;
t382 = t459 * t388 + t457 * t389;
t433 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t461 + Ifges(6,6) * t464) * t455;
t402 = -mrSges(6,1) * t415 + mrSges(6,3) * t414 + Ifges(6,4) * t439 + Ifges(6,2) * t440 + Ifges(6,6) * qJDD(5) + qJD(5) * t435 - t433 * t483;
t403 = mrSges(6,2) * t415 - mrSges(6,3) * t413 + Ifges(6,1) * t439 + Ifges(6,4) * t440 + Ifges(6,5) * qJDD(5) - qJD(5) * t434 + t433 * t482;
t471 = -mrSges(5,1) * t418 + mrSges(5,2) * t419 - Ifges(5,3) * t454 - pkin(4) * t470 - pkin(7) * t473 - t464 * t402 - t461 * t403;
t485 = mrSges(3,1) * t431 + mrSges(4,1) * t424 - mrSges(3,2) * t432 - mrSges(4,2) * t425 + pkin(2) * t382 + pkin(3) * t390 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) - t471;
t380 = m(3) * t431 + qJDD(2) * mrSges(3,1) - t467 * mrSges(3,2) + t382;
t475 = -t457 * t388 + t459 * t389;
t381 = m(3) * t432 - t467 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t475;
t476 = -t463 * t380 + t466 * t381;
t373 = m(2) * t447 + t476;
t398 = t464 * t409 + t461 * t410;
t479 = m(5) * t445 + t398;
t396 = m(4) * t445 + t479;
t395 = (m(2) + m(3)) * t446 - t396;
t481 = t458 * t373 + t460 * t395;
t374 = t466 * t380 + t463 * t381;
t478 = m(2) * t456 + t374;
t477 = t460 * t373 - t458 * t395;
t387 = -mrSges(5,1) * t445 + mrSges(5,3) * t419 + t453 * Ifges(5,5) + Ifges(5,6) * t454 - pkin(4) * t398 - t486;
t383 = mrSges(5,2) * t445 - mrSges(5,3) * t418 + Ifges(5,5) * t454 - t453 * Ifges(5,6) - pkin(7) * t398 - t461 * t402 + t464 * t403;
t376 = mrSges(4,2) * t445 - mrSges(4,3) * t424 + Ifges(4,5) * qJDD(2) - t467 * Ifges(4,6) - pkin(6) * t390 + t465 * t383 - t462 * t387;
t375 = -mrSges(4,1) * t445 + mrSges(4,3) * t425 + t467 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t479 + pkin(6) * t474 + t462 * t383 + t465 * t387;
t370 = -mrSges(2,1) * t456 + mrSges(2,3) * t447 - pkin(1) * t374 - t485;
t369 = -mrSges(3,2) * t446 - mrSges(3,3) * t431 + Ifges(3,5) * qJDD(2) - t467 * Ifges(3,6) - qJ(3) * t382 - t457 * t375 + t459 * t376;
t368 = mrSges(3,1) * t446 + mrSges(3,3) * t432 + t467 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t396 + qJ(3) * t475 + t459 * t375 + t457 * t376;
t367 = mrSges(2,2) * t456 - mrSges(2,3) * t446 - pkin(5) * t374 - t463 * t368 + t466 * t369;
t1 = [-m(1) * g(1) + t477; -m(1) * g(2) + t481; -m(1) * g(3) + t478; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t481 + t460 * t367 - t458 * t370; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t477 + t458 * t367 + t460 * t370; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t446 - mrSges(2,2) * t447 + t463 * t369 + t466 * t368 + pkin(1) * (m(3) * t446 - t396) + pkin(5) * t476; t478; t485; t396; -t471; t486;];
tauJB = t1;
