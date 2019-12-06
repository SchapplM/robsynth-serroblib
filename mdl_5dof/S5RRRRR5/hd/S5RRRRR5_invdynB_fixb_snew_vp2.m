% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:23
% EndTime: 2019-12-05 18:58:26
% DurationCPUTime: 2.77s
% Computational Cost: add. (57651->232), mult. (59398->293), div. (0->0), fcn. (33553->10), ass. (0->98)
t468 = -m(3) - m(4);
t441 = qJD(1) + qJD(2);
t433 = qJD(3) + t441;
t431 = t433 ^ 2;
t467 = pkin(4) * t431;
t444 = sin(qJ(4));
t466 = t433 * t444;
t449 = cos(qJ(4));
t465 = t433 * t449;
t447 = sin(qJ(1));
t452 = cos(qJ(1));
t430 = t452 * g(2) + t447 * g(3);
t427 = qJDD(1) * pkin(1) + t430;
t429 = t447 * g(2) - t452 * g(3);
t453 = qJD(1) ^ 2;
t428 = -t453 * pkin(1) + t429;
t446 = sin(qJ(2));
t451 = cos(qJ(2));
t407 = t451 * t427 - t446 * t428;
t439 = qJDD(1) + qJDD(2);
t405 = t439 * pkin(2) + t407;
t408 = t446 * t427 + t451 * t428;
t437 = t441 ^ 2;
t406 = -t437 * pkin(2) + t408;
t445 = sin(qJ(3));
t450 = cos(qJ(3));
t393 = t445 * t405 + t450 * t406;
t432 = qJDD(3) + t439;
t391 = -t431 * pkin(3) + t432 * pkin(8) + t393;
t464 = t444 * t391;
t463 = qJD(4) * t433;
t419 = t444 * t432 + t449 * t463;
t384 = qJDD(4) * pkin(4) - t419 * pkin(9) - t464 + (pkin(9) * t463 + t444 * t467 - g(1)) * t449;
t388 = -t444 * g(1) + t449 * t391;
t420 = t449 * t432 - t444 * t463;
t426 = qJD(4) * pkin(4) - pkin(9) * t466;
t442 = t449 ^ 2;
t385 = t420 * pkin(9) - qJD(4) * t426 - t442 * t467 + t388;
t443 = sin(qJ(5));
t448 = cos(qJ(5));
t382 = t448 * t384 - t443 * t385;
t414 = (-t443 * t444 + t448 * t449) * t433;
t396 = t414 * qJD(5) + t448 * t419 + t443 * t420;
t415 = (t443 * t449 + t444 * t448) * t433;
t401 = -t414 * mrSges(6,1) + t415 * mrSges(6,2);
t440 = qJD(4) + qJD(5);
t409 = -t440 * mrSges(6,2) + t414 * mrSges(6,3);
t438 = qJDD(4) + qJDD(5);
t380 = m(6) * t382 + t438 * mrSges(6,1) - t396 * mrSges(6,3) - t415 * t401 + t440 * t409;
t383 = t443 * t384 + t448 * t385;
t395 = -t415 * qJD(5) - t443 * t419 + t448 * t420;
t410 = t440 * mrSges(6,1) - t415 * mrSges(6,3);
t381 = m(6) * t383 - t438 * mrSges(6,2) + t395 * mrSges(6,3) + t414 * t401 - t440 * t410;
t372 = t448 * t380 + t443 * t381;
t387 = -t449 * g(1) - t464;
t418 = (-mrSges(5,1) * t449 + mrSges(5,2) * t444) * t433;
t425 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t465;
t370 = m(5) * t387 + qJDD(4) * mrSges(5,1) - t419 * mrSges(5,3) + qJD(4) * t425 - t418 * t466 + t372;
t424 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t466;
t458 = -t443 * t380 + t448 * t381;
t371 = m(5) * t388 - qJDD(4) * mrSges(5,2) + t420 * mrSges(5,3) - qJD(4) * t424 + t418 * t465 + t458;
t459 = -t444 * t370 + t449 * t371;
t365 = m(4) * t393 - t431 * mrSges(4,1) - t432 * mrSges(4,2) + t459;
t392 = t450 * t405 - t445 * t406;
t456 = -t432 * pkin(3) - t392;
t390 = -t431 * pkin(8) + t456;
t386 = t426 * t466 - t420 * pkin(4) + (-pkin(9) * t442 - pkin(8)) * t431 + t456;
t455 = m(6) * t386 - t395 * mrSges(6,1) + t396 * mrSges(6,2) - t414 * t409 + t415 * t410;
t454 = -m(5) * t390 + t420 * mrSges(5,1) - t419 * mrSges(5,2) - t424 * t466 + t425 * t465 - t455;
t376 = m(4) * t392 + t432 * mrSges(4,1) - t431 * mrSges(4,2) + t454;
t361 = t445 * t365 + t450 * t376;
t359 = m(3) * t407 + t439 * mrSges(3,1) - t437 * mrSges(3,2) + t361;
t460 = t450 * t365 - t445 * t376;
t360 = m(3) * t408 - t437 * mrSges(3,1) - t439 * mrSges(3,2) + t460;
t353 = t451 * t359 + t446 * t360;
t366 = t449 * t370 + t444 * t371;
t461 = -t446 * t359 + t451 * t360;
t351 = m(2) * t429 - t453 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t461;
t352 = m(2) * t430 + qJDD(1) * mrSges(2,1) - t453 * mrSges(2,2) + t353;
t462 = t452 * t351 - t447 * t352;
t457 = -t447 * t351 - t452 * t352;
t413 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t444 + Ifges(5,4) * t449) * t433;
t412 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t444 + Ifges(5,2) * t449) * t433;
t411 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t444 + Ifges(5,6) * t449) * t433;
t399 = Ifges(6,1) * t415 + Ifges(6,4) * t414 + Ifges(6,5) * t440;
t398 = Ifges(6,4) * t415 + Ifges(6,2) * t414 + Ifges(6,6) * t440;
t397 = Ifges(6,5) * t415 + Ifges(6,6) * t414 + Ifges(6,3) * t440;
t374 = mrSges(6,2) * t386 - mrSges(6,3) * t382 + Ifges(6,1) * t396 + Ifges(6,4) * t395 + Ifges(6,5) * t438 + t414 * t397 - t440 * t398;
t373 = -mrSges(6,1) * t386 + mrSges(6,3) * t383 + Ifges(6,4) * t396 + Ifges(6,2) * t395 + Ifges(6,6) * t438 - t415 * t397 + t440 * t399;
t362 = mrSges(5,2) * t390 - mrSges(5,3) * t387 + Ifges(5,1) * t419 + Ifges(5,4) * t420 + Ifges(5,5) * qJDD(4) - pkin(9) * t372 - qJD(4) * t412 - t443 * t373 + t448 * t374 + t411 * t465;
t355 = -mrSges(5,1) * t390 + mrSges(5,3) * t388 + Ifges(5,4) * t419 + Ifges(5,2) * t420 + Ifges(5,6) * qJDD(4) - pkin(4) * t455 + pkin(9) * t458 + qJD(4) * t413 + t448 * t373 + t443 * t374 - t411 * t466;
t354 = Ifges(4,6) * t432 + t431 * Ifges(4,5) + mrSges(4,1) * g(1) + mrSges(4,3) * t393 - Ifges(5,5) * t419 - Ifges(5,6) * t420 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t387 + mrSges(5,2) * t388 - Ifges(6,5) * t396 - Ifges(6,6) * t395 - Ifges(6,3) * t438 - t415 * t398 + t414 * t399 - mrSges(6,1) * t382 + mrSges(6,2) * t383 - pkin(4) * t372 - pkin(3) * t366 + (-t444 * t412 + t449 * t413) * t433;
t349 = -mrSges(4,2) * g(1) - mrSges(4,3) * t392 + Ifges(4,5) * t432 - t431 * Ifges(4,6) - pkin(8) * t366 - t444 * t355 + t449 * t362;
t348 = -mrSges(3,2) * g(1) - mrSges(3,3) * t407 + Ifges(3,5) * t439 - t437 * Ifges(3,6) - pkin(7) * t361 + t450 * t349 - t445 * t354;
t347 = Ifges(3,6) * t439 + t437 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t408 + t445 * t349 + t450 * t354 - pkin(2) * (-m(4) * g(1) + t366) + pkin(7) * t460;
t346 = -mrSges(2,2) * g(1) - mrSges(2,3) * t430 + Ifges(2,5) * qJDD(1) - t453 * Ifges(2,6) - pkin(6) * t353 - t446 * t347 + t451 * t348;
t345 = Ifges(2,6) * qJDD(1) + t453 * Ifges(2,5) + mrSges(2,3) * t429 + t446 * t348 + t451 * t347 - pkin(1) * t366 + pkin(6) * t461 + (-pkin(1) * t468 + mrSges(2,1)) * g(1);
t1 = [(-m(1) - m(2) + t468) * g(1) + t366; -m(1) * g(2) + t457; -m(1) * g(3) + t462; pkin(1) * t353 + pkin(2) * t361 + mrSges(3,1) * t407 - mrSges(3,2) * t408 + t444 * t362 + t449 * t355 + pkin(3) * t454 + pkin(8) * t459 + mrSges(4,1) * t392 - mrSges(4,2) * t393 + mrSges(2,1) * t430 - mrSges(2,2) * t429 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t439 + Ifges(4,3) * t432 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t462 - t452 * t345 - t447 * t346; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t457 - t447 * t345 + t452 * t346;];
tauB = t1;
