% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:20
% EndTime: 2019-12-05 18:18:23
% DurationCPUTime: 2.53s
% Computational Cost: add. (39959->207), mult. (54393->259), div. (0->0), fcn. (31463->10), ass. (0->94)
t438 = qJD(1) + qJD(2);
t434 = t438 ^ 2;
t442 = cos(pkin(9));
t473 = pkin(4) * t442;
t440 = sin(pkin(9));
t472 = mrSges(5,2) * t440;
t437 = t442 ^ 2;
t471 = t434 * t437;
t435 = qJDD(1) + qJDD(2);
t470 = t435 * t442;
t458 = Ifges(5,5) * t440 + Ifges(5,6) * t442;
t469 = t434 * t458;
t446 = sin(qJ(1));
t449 = cos(qJ(1));
t424 = t449 * g(2) + t446 * g(3);
t419 = qJDD(1) * pkin(1) + t424;
t423 = t446 * g(2) - t449 * g(3);
t450 = qJD(1) ^ 2;
t420 = -t450 * pkin(1) + t423;
t445 = sin(qJ(2));
t448 = cos(qJ(2));
t408 = t448 * t419 - t445 * t420;
t406 = t435 * pkin(2) + t408;
t409 = t445 * t419 + t448 * t420;
t407 = -t434 * pkin(2) + t409;
t441 = sin(pkin(8));
t443 = cos(pkin(8));
t394 = t441 * t406 + t443 * t407;
t392 = -t434 * pkin(3) + t435 * qJ(4) + t394;
t439 = -g(1) + qJDD(3);
t467 = qJD(4) * t438;
t468 = t442 * t439 - 0.2e1 * t440 * t467;
t385 = (-pkin(7) * t435 + t434 * t473 - t392) * t440 + t468;
t389 = t440 * t439 + (t392 + 0.2e1 * t467) * t442;
t386 = -pkin(4) * t471 + pkin(7) * t470 + t389;
t444 = sin(qJ(5));
t447 = cos(qJ(5));
t383 = t447 * t385 - t444 * t386;
t453 = -t440 * t444 + t442 * t447;
t412 = t453 * t438;
t454 = t440 * t447 + t442 * t444;
t413 = t454 * t438;
t400 = -t412 * mrSges(6,1) + t413 * mrSges(6,2);
t402 = t412 * qJD(5) + t454 * t435;
t410 = -qJD(5) * mrSges(6,2) + t412 * mrSges(6,3);
t381 = m(6) * t383 + qJDD(5) * mrSges(6,1) - t402 * mrSges(6,3) + qJD(5) * t410 - t413 * t400;
t384 = t444 * t385 + t447 * t386;
t401 = -t413 * qJD(5) + t453 * t435;
t411 = qJD(5) * mrSges(6,1) - t413 * mrSges(6,3);
t382 = m(6) * t384 - qJDD(5) * mrSges(6,2) + t401 * mrSges(6,3) - qJD(5) * t411 + t412 * t400;
t373 = t447 * t381 + t444 * t382;
t388 = -t440 * t392 + t468;
t456 = mrSges(5,3) * t435 + (-mrSges(5,1) * t442 + t472) * t434;
t371 = m(5) * t388 - t456 * t440 + t373;
t461 = -t444 * t381 + t447 * t382;
t372 = m(5) * t389 + t456 * t442 + t461;
t462 = -t440 * t371 + t442 * t372;
t366 = m(4) * t394 - t434 * mrSges(4,1) - t435 * mrSges(4,2) + t462;
t393 = t443 * t406 - t441 * t407;
t457 = qJDD(4) - t393;
t391 = -t435 * pkin(3) - t434 * qJ(4) + t457;
t436 = t440 ^ 2;
t387 = (-pkin(3) - t473) * t435 + (-qJ(4) + (-t436 - t437) * pkin(7)) * t434 + t457;
t452 = m(6) * t387 - t401 * mrSges(6,1) + t402 * mrSges(6,2) - t412 * t410 + t413 * t411;
t451 = -m(5) * t391 + mrSges(5,1) * t470 - t452 + (t434 * t436 + t471) * mrSges(5,3);
t377 = m(4) * t393 - t434 * mrSges(4,2) + (mrSges(4,1) - t472) * t435 + t451;
t362 = t441 * t366 + t443 * t377;
t359 = m(3) * t408 + t435 * mrSges(3,1) - t434 * mrSges(3,2) + t362;
t463 = t443 * t366 - t441 * t377;
t360 = m(3) * t409 - t434 * mrSges(3,1) - t435 * mrSges(3,2) + t463;
t354 = t448 * t359 + t445 * t360;
t367 = t442 * t371 + t440 * t372;
t466 = m(4) * t439 + t367;
t464 = -t445 * t359 + t448 * t360;
t352 = m(2) * t423 - t450 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t464;
t353 = m(2) * t424 + qJDD(1) * mrSges(2,1) - t450 * mrSges(2,2) + t354;
t465 = t449 * t352 - t446 * t353;
t460 = Ifges(5,1) * t440 + Ifges(5,4) * t442;
t459 = Ifges(5,4) * t440 + Ifges(5,2) * t442;
t455 = -t446 * t352 - t449 * t353;
t397 = Ifges(6,1) * t413 + Ifges(6,4) * t412 + Ifges(6,5) * qJD(5);
t396 = Ifges(6,4) * t413 + Ifges(6,2) * t412 + Ifges(6,6) * qJD(5);
t395 = Ifges(6,5) * t413 + Ifges(6,6) * t412 + Ifges(6,3) * qJD(5);
t375 = mrSges(6,2) * t387 - mrSges(6,3) * t383 + Ifges(6,1) * t402 + Ifges(6,4) * t401 + Ifges(6,5) * qJDD(5) - qJD(5) * t396 + t412 * t395;
t374 = -mrSges(6,1) * t387 + mrSges(6,3) * t384 + Ifges(6,4) * t402 + Ifges(6,2) * t401 + Ifges(6,6) * qJDD(5) + qJD(5) * t397 - t413 * t395;
t363 = mrSges(5,2) * t391 - mrSges(5,3) * t388 - pkin(7) * t373 - t444 * t374 + t447 * t375 + t460 * t435 + t442 * t469;
t361 = -mrSges(5,1) * t391 + mrSges(5,3) * t389 - pkin(4) * t452 + pkin(7) * t461 + t447 * t374 + t444 * t375 + t459 * t435 - t440 * t469;
t355 = -mrSges(4,1) * t439 - mrSges(5,1) * t388 - mrSges(6,1) * t383 + mrSges(5,2) * t389 + mrSges(6,2) * t384 + mrSges(4,3) * t394 - Ifges(6,5) * t402 - Ifges(6,6) * t401 - Ifges(6,3) * qJDD(5) - pkin(3) * t367 - pkin(4) * t373 - t413 * t396 + t412 * t397 + (Ifges(4,6) - t458) * t435 + (-t440 * t459 + t442 * t460 + Ifges(4,5)) * t434;
t350 = mrSges(4,2) * t439 - mrSges(4,3) * t393 + Ifges(4,5) * t435 - t434 * Ifges(4,6) - qJ(4) * t367 - t440 * t361 + t442 * t363;
t349 = -mrSges(3,2) * g(1) - mrSges(3,3) * t408 + Ifges(3,5) * t435 - t434 * Ifges(3,6) - qJ(3) * t362 + t443 * t350 - t441 * t355;
t348 = mrSges(3,1) * g(1) + mrSges(3,3) * t409 + t434 * Ifges(3,5) + Ifges(3,6) * t435 - pkin(2) * t466 + qJ(3) * t463 + t441 * t350 + t443 * t355;
t347 = -mrSges(2,2) * g(1) - mrSges(2,3) * t424 + Ifges(2,5) * qJDD(1) - t450 * Ifges(2,6) - pkin(6) * t354 - t445 * t348 + t448 * t349;
t346 = Ifges(2,6) * qJDD(1) + t450 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t423 + t445 * t349 + t448 * t348 - pkin(1) * (-m(3) * g(1) + t466) + pkin(6) * t464;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t466; -m(1) * g(2) + t455; -m(1) * g(3) + t465; pkin(1) * t354 + pkin(2) * t362 + mrSges(3,1) * t408 - mrSges(3,2) * t409 + t440 * t363 + t442 * t361 + pkin(3) * t451 + qJ(4) * t462 - mrSges(4,2) * t394 + mrSges(4,1) * t393 + mrSges(2,1) * t424 - mrSges(2,2) * t423 + Ifges(2,3) * qJDD(1) - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(3) * t472 + Ifges(3,3) + Ifges(4,3)) * t435; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t465 - t449 * t346 - t446 * t347; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t455 - t446 * t346 + t449 * t347;];
tauB = t1;
