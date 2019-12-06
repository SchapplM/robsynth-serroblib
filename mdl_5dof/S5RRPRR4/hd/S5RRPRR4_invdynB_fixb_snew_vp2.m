% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:32:01
% EndTime: 2019-12-05 18:32:04
% DurationCPUTime: 2.68s
% Computational Cost: add. (46120->229), mult. (59398->290), div. (0->0), fcn. (33553->10), ass. (0->95)
t441 = qJD(1) + qJD(2);
t447 = sin(qJ(4));
t468 = t441 * t447;
t451 = cos(qJ(4));
t467 = t441 * t451;
t449 = sin(qJ(1));
t453 = cos(qJ(1));
t430 = t453 * g(2) + t449 * g(3);
t424 = qJDD(1) * pkin(1) + t430;
t429 = t449 * g(2) - t453 * g(3);
t454 = qJD(1) ^ 2;
t425 = -t454 * pkin(1) + t429;
t448 = sin(qJ(2));
t452 = cos(qJ(2));
t407 = t452 * t424 - t448 * t425;
t439 = qJDD(1) + qJDD(2);
t405 = t439 * pkin(2) + t407;
t408 = t448 * t424 + t452 * t425;
t437 = t441 ^ 2;
t406 = -t437 * pkin(2) + t408;
t444 = sin(pkin(9));
t445 = cos(pkin(9));
t393 = t444 * t405 + t445 * t406;
t391 = -t437 * pkin(3) + t439 * pkin(7) + t393;
t443 = -g(1) + qJDD(3);
t387 = -t447 * t391 + t451 * t443;
t466 = qJD(4) * t441;
t464 = t451 * t466;
t419 = t447 * t439 + t464;
t384 = (-t419 + t464) * pkin(8) + (t437 * t447 * t451 + qJDD(4)) * pkin(4) + t387;
t388 = t451 * t391 + t447 * t443;
t420 = t451 * t439 - t447 * t466;
t428 = qJD(4) * pkin(4) - pkin(8) * t468;
t442 = t451 ^ 2;
t385 = -t442 * t437 * pkin(4) + t420 * pkin(8) - qJD(4) * t428 + t388;
t446 = sin(qJ(5));
t450 = cos(qJ(5));
t382 = t450 * t384 - t446 * t385;
t414 = (-t446 * t447 + t450 * t451) * t441;
t396 = t414 * qJD(5) + t450 * t419 + t446 * t420;
t415 = (t446 * t451 + t447 * t450) * t441;
t401 = -t414 * mrSges(6,1) + t415 * mrSges(6,2);
t440 = qJD(4) + qJD(5);
t409 = -t440 * mrSges(6,2) + t414 * mrSges(6,3);
t438 = qJDD(4) + qJDD(5);
t380 = m(6) * t382 + t438 * mrSges(6,1) - t396 * mrSges(6,3) - t415 * t401 + t440 * t409;
t383 = t446 * t384 + t450 * t385;
t395 = -t415 * qJD(5) - t446 * t419 + t450 * t420;
t410 = t440 * mrSges(6,1) - t415 * mrSges(6,3);
t381 = m(6) * t383 - t438 * mrSges(6,2) + t395 * mrSges(6,3) + t414 * t401 - t440 * t410;
t372 = t450 * t380 + t446 * t381;
t418 = (-mrSges(5,1) * t451 + mrSges(5,2) * t447) * t441;
t427 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t467;
t370 = m(5) * t387 + qJDD(4) * mrSges(5,1) - t419 * mrSges(5,3) + qJD(4) * t427 - t418 * t468 + t372;
t426 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t468;
t459 = -t446 * t380 + t450 * t381;
t371 = m(5) * t388 - qJDD(4) * mrSges(5,2) + t420 * mrSges(5,3) - qJD(4) * t426 + t418 * t467 + t459;
t460 = -t447 * t370 + t451 * t371;
t365 = m(4) * t393 - t437 * mrSges(4,1) - t439 * mrSges(4,2) + t460;
t392 = t445 * t405 - t444 * t406;
t457 = -t439 * pkin(3) - t392;
t390 = -t437 * pkin(7) + t457;
t386 = t428 * t468 - t420 * pkin(4) + (-pkin(8) * t442 - pkin(7)) * t437 + t457;
t456 = m(6) * t386 - t395 * mrSges(6,1) + t396 * mrSges(6,2) - t414 * t409 + t415 * t410;
t455 = -m(5) * t390 + t420 * mrSges(5,1) - t419 * mrSges(5,2) - t426 * t468 + t427 * t467 - t456;
t376 = m(4) * t392 + t439 * mrSges(4,1) - t437 * mrSges(4,2) + t455;
t361 = t444 * t365 + t445 * t376;
t359 = m(3) * t407 + t439 * mrSges(3,1) - t437 * mrSges(3,2) + t361;
t461 = t445 * t365 - t444 * t376;
t360 = m(3) * t408 - t437 * mrSges(3,1) - t439 * mrSges(3,2) + t461;
t353 = t452 * t359 + t448 * t360;
t366 = t451 * t370 + t447 * t371;
t465 = m(4) * t443 + t366;
t462 = -t448 * t359 + t452 * t360;
t351 = m(2) * t429 - t454 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t462;
t352 = m(2) * t430 + qJDD(1) * mrSges(2,1) - t454 * mrSges(2,2) + t353;
t463 = t453 * t351 - t449 * t352;
t458 = -t449 * t351 - t453 * t352;
t413 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t447 + Ifges(5,4) * t451) * t441;
t412 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t447 + Ifges(5,2) * t451) * t441;
t411 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t447 + Ifges(5,6) * t451) * t441;
t399 = Ifges(6,1) * t415 + Ifges(6,4) * t414 + Ifges(6,5) * t440;
t398 = Ifges(6,4) * t415 + Ifges(6,2) * t414 + Ifges(6,6) * t440;
t397 = Ifges(6,5) * t415 + Ifges(6,6) * t414 + Ifges(6,3) * t440;
t374 = mrSges(6,2) * t386 - mrSges(6,3) * t382 + Ifges(6,1) * t396 + Ifges(6,4) * t395 + Ifges(6,5) * t438 + t414 * t397 - t440 * t398;
t373 = -mrSges(6,1) * t386 + mrSges(6,3) * t383 + Ifges(6,4) * t396 + Ifges(6,2) * t395 + Ifges(6,6) * t438 - t415 * t397 + t440 * t399;
t362 = mrSges(5,2) * t390 - mrSges(5,3) * t387 + Ifges(5,1) * t419 + Ifges(5,4) * t420 + Ifges(5,5) * qJDD(4) - pkin(8) * t372 - qJD(4) * t412 - t446 * t373 + t450 * t374 + t411 * t467;
t355 = -mrSges(5,1) * t390 + mrSges(5,3) * t388 + Ifges(5,4) * t419 + Ifges(5,2) * t420 + Ifges(5,6) * qJDD(4) - pkin(4) * t456 + pkin(8) * t459 + qJD(4) * t413 + t450 * t373 + t446 * t374 - t411 * t468;
t354 = Ifges(4,6) * t439 + t437 * Ifges(4,5) - mrSges(4,1) * t443 + mrSges(4,3) * t393 - Ifges(5,5) * t419 - Ifges(5,6) * t420 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t387 + mrSges(5,2) * t388 - Ifges(6,5) * t396 - Ifges(6,6) * t395 - Ifges(6,3) * t438 - t415 * t398 + t414 * t399 - mrSges(6,1) * t382 + mrSges(6,2) * t383 - pkin(4) * t372 - pkin(3) * t366 + (-t447 * t412 + t451 * t413) * t441;
t349 = mrSges(4,2) * t443 - mrSges(4,3) * t392 + Ifges(4,5) * t439 - t437 * Ifges(4,6) - pkin(7) * t366 - t447 * t355 + t451 * t362;
t348 = -mrSges(3,2) * g(1) - mrSges(3,3) * t407 + Ifges(3,5) * t439 - t437 * Ifges(3,6) - qJ(3) * t361 + t445 * t349 - t444 * t354;
t347 = mrSges(3,1) * g(1) + mrSges(3,3) * t408 + t437 * Ifges(3,5) + Ifges(3,6) * t439 - pkin(2) * t465 + qJ(3) * t461 + t444 * t349 + t445 * t354;
t346 = -mrSges(2,2) * g(1) - mrSges(2,3) * t430 + Ifges(2,5) * qJDD(1) - t454 * Ifges(2,6) - pkin(6) * t353 - t448 * t347 + t452 * t348;
t345 = Ifges(2,6) * qJDD(1) + t454 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t429 + t448 * t348 + t452 * t347 - pkin(1) * (-m(3) * g(1) + t465) + pkin(6) * t462;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t465; -m(1) * g(2) + t458; -m(1) * g(3) + t463; pkin(1) * t353 + pkin(2) * t361 + mrSges(3,1) * t407 - mrSges(3,2) * t408 + t451 * t355 + pkin(3) * t455 + pkin(7) * t460 + t447 * t362 + mrSges(4,1) * t392 - mrSges(4,2) * t393 + mrSges(2,1) * t430 - mrSges(2,2) * t429 + Ifges(2,3) * qJDD(1) - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(4,3) + Ifges(3,3)) * t439; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t463 - t453 * t345 - t449 * t346; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t458 - t449 * t345 + t453 * t346;];
tauB = t1;
