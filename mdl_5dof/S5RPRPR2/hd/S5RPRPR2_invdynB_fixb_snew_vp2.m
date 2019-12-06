% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:31
% EndTime: 2019-12-05 17:49:33
% DurationCPUTime: 2.60s
% Computational Cost: add. (37791->206), mult. (54393->259), div. (0->0), fcn. (31463->10), ass. (0->95)
t438 = qJD(1) + qJD(3);
t434 = t438 ^ 2;
t442 = cos(pkin(9));
t474 = pkin(4) * t442;
t440 = sin(pkin(9));
t473 = mrSges(5,2) * t440;
t437 = t442 ^ 2;
t472 = t434 * t437;
t435 = qJDD(1) + qJDD(3);
t471 = t435 * t442;
t458 = Ifges(5,5) * t440 + Ifges(5,6) * t442;
t470 = t434 * t458;
t446 = sin(qJ(1));
t449 = cos(qJ(1));
t423 = t449 * g(2) + t446 * g(3);
t418 = qJDD(1) * pkin(1) + t423;
t422 = t446 * g(2) - t449 * g(3);
t450 = qJD(1) ^ 2;
t419 = -t450 * pkin(1) + t422;
t441 = sin(pkin(8));
t443 = cos(pkin(8));
t407 = t443 * t418 - t441 * t419;
t405 = qJDD(1) * pkin(2) + t407;
t408 = t441 * t418 + t443 * t419;
t406 = -t450 * pkin(2) + t408;
t445 = sin(qJ(3));
t448 = cos(qJ(3));
t393 = t445 * t405 + t448 * t406;
t391 = -t434 * pkin(3) + t435 * qJ(4) + t393;
t439 = -g(1) + qJDD(2);
t468 = qJD(4) * t438;
t469 = t442 * t439 - 0.2e1 * t440 * t468;
t384 = (-pkin(7) * t435 + t434 * t474 - t391) * t440 + t469;
t388 = t440 * t439 + (t391 + 0.2e1 * t468) * t442;
t385 = -pkin(4) * t472 + pkin(7) * t471 + t388;
t444 = sin(qJ(5));
t447 = cos(qJ(5));
t382 = t447 * t384 - t444 * t385;
t453 = -t440 * t444 + t442 * t447;
t411 = t453 * t438;
t454 = t440 * t447 + t442 * t444;
t412 = t454 * t438;
t399 = -t411 * mrSges(6,1) + t412 * mrSges(6,2);
t401 = t411 * qJD(5) + t454 * t435;
t409 = -qJD(5) * mrSges(6,2) + t411 * mrSges(6,3);
t380 = m(6) * t382 + qJDD(5) * mrSges(6,1) - t401 * mrSges(6,3) + qJD(5) * t409 - t412 * t399;
t383 = t444 * t384 + t447 * t385;
t400 = -t412 * qJD(5) + t453 * t435;
t410 = qJD(5) * mrSges(6,1) - t412 * mrSges(6,3);
t381 = m(6) * t383 - qJDD(5) * mrSges(6,2) + t400 * mrSges(6,3) - qJD(5) * t410 + t411 * t399;
t372 = t447 * t380 + t444 * t381;
t387 = -t440 * t391 + t469;
t456 = mrSges(5,3) * t435 + (-mrSges(5,1) * t442 + t473) * t434;
t370 = m(5) * t387 - t456 * t440 + t372;
t462 = -t444 * t380 + t447 * t381;
t371 = m(5) * t388 + t456 * t442 + t462;
t463 = -t440 * t370 + t442 * t371;
t365 = m(4) * t393 - t434 * mrSges(4,1) - t435 * mrSges(4,2) + t463;
t392 = t448 * t405 - t445 * t406;
t457 = qJDD(4) - t392;
t390 = -t435 * pkin(3) - t434 * qJ(4) + t457;
t436 = t440 ^ 2;
t386 = (-pkin(3) - t474) * t435 + (-qJ(4) + (-t436 - t437) * pkin(7)) * t434 + t457;
t452 = m(6) * t386 - t400 * mrSges(6,1) + t401 * mrSges(6,2) - t411 * t409 + t412 * t410;
t451 = -m(5) * t390 + mrSges(5,1) * t471 - t452 + (t434 * t436 + t472) * mrSges(5,3);
t376 = m(4) * t392 - t434 * mrSges(4,2) + (mrSges(4,1) - t473) * t435 + t451;
t361 = t445 * t365 + t448 * t376;
t358 = m(3) * t407 + qJDD(1) * mrSges(3,1) - t450 * mrSges(3,2) + t361;
t464 = t448 * t365 - t445 * t376;
t359 = m(3) * t408 - t450 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t464;
t353 = t443 * t358 + t441 * t359;
t366 = t442 * t370 + t440 * t371;
t467 = m(4) * t439 + t366;
t465 = -t441 * t358 + t443 * t359;
t351 = m(2) * t422 - t450 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t465;
t352 = m(2) * t423 + qJDD(1) * mrSges(2,1) - t450 * mrSges(2,2) + t353;
t466 = t449 * t351 - t446 * t352;
t461 = m(3) * t439 + t467;
t460 = Ifges(5,1) * t440 + Ifges(5,4) * t442;
t459 = Ifges(5,4) * t440 + Ifges(5,2) * t442;
t455 = -t446 * t351 - t449 * t352;
t396 = Ifges(6,1) * t412 + Ifges(6,4) * t411 + Ifges(6,5) * qJD(5);
t395 = Ifges(6,4) * t412 + Ifges(6,2) * t411 + Ifges(6,6) * qJD(5);
t394 = Ifges(6,5) * t412 + Ifges(6,6) * t411 + Ifges(6,3) * qJD(5);
t374 = mrSges(6,2) * t386 - mrSges(6,3) * t382 + Ifges(6,1) * t401 + Ifges(6,4) * t400 + Ifges(6,5) * qJDD(5) - qJD(5) * t395 + t411 * t394;
t373 = -mrSges(6,1) * t386 + mrSges(6,3) * t383 + Ifges(6,4) * t401 + Ifges(6,2) * t400 + Ifges(6,6) * qJDD(5) + qJD(5) * t396 - t412 * t394;
t362 = mrSges(5,2) * t390 - mrSges(5,3) * t387 - pkin(7) * t372 - t444 * t373 + t447 * t374 + t460 * t435 + t442 * t470;
t360 = -mrSges(5,1) * t390 + mrSges(5,3) * t388 - pkin(4) * t452 + pkin(7) * t462 + t447 * t373 + t444 * t374 + t459 * t435 - t440 * t470;
t354 = -mrSges(4,1) * t439 - mrSges(5,1) * t387 - mrSges(6,1) * t382 + mrSges(5,2) * t388 + mrSges(6,2) * t383 + mrSges(4,3) * t393 - Ifges(6,5) * t401 - Ifges(6,6) * t400 - Ifges(6,3) * qJDD(5) - pkin(3) * t366 - pkin(4) * t372 - t412 * t395 + t411 * t396 + (Ifges(4,6) - t458) * t435 + (-t440 * t459 + t442 * t460 + Ifges(4,5)) * t434;
t349 = mrSges(4,2) * t439 - mrSges(4,3) * t392 + Ifges(4,5) * t435 - t434 * Ifges(4,6) - qJ(4) * t366 - t440 * t360 + t442 * t362;
t348 = mrSges(3,2) * t439 - mrSges(3,3) * t407 + Ifges(3,5) * qJDD(1) - t450 * Ifges(3,6) - pkin(6) * t361 + t448 * t349 - t445 * t354;
t347 = -mrSges(3,1) * t439 + mrSges(3,3) * t408 + t450 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t467 + pkin(6) * t464 + t445 * t349 + t448 * t354;
t346 = -mrSges(2,2) * g(1) - mrSges(2,3) * t423 + Ifges(2,5) * qJDD(1) - t450 * Ifges(2,6) - qJ(2) * t353 - t441 * t347 + t443 * t348;
t345 = mrSges(2,1) * g(1) + mrSges(2,3) * t422 + t450 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t461 + qJ(2) * t465 + t443 * t347 + t441 * t348;
t1 = [(-m(1) - m(2)) * g(1) + t461; -m(1) * g(2) + t455; -m(1) * g(3) + t466; pkin(1) * t353 + pkin(2) * t361 + mrSges(3,1) * t407 - mrSges(3,2) * t408 + t440 * t362 + t442 * t360 + pkin(3) * (-t435 * t473 + t451) + qJ(4) * t463 + mrSges(4,1) * t392 - mrSges(4,2) * t393 + mrSges(2,1) * t423 - mrSges(2,2) * t422 + Ifges(4,3) * t435 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t466 - t449 * t345 - t446 * t346; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t455 - t446 * t345 + t449 * t346;];
tauB = t1;
