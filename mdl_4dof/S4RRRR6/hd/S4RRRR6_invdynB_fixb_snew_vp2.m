% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:40
% EndTime: 2019-12-31 17:29:44
% DurationCPUTime: 3.79s
% Computational Cost: add. (42283->251), mult. (91557->336), div. (0->0), fcn. (67836->10), ass. (0->109)
t440 = sin(pkin(4));
t444 = sin(qJ(2));
t448 = cos(qJ(2));
t460 = qJD(1) * qJD(2);
t430 = (-qJDD(1) * t448 + t444 * t460) * t440;
t470 = pkin(6) * t440;
t441 = cos(pkin(4));
t469 = t441 * g(3);
t468 = t440 * t444;
t467 = t440 * t448;
t466 = t441 * t444;
t465 = t441 * t448;
t445 = sin(qJ(1));
t449 = cos(qJ(1));
t433 = t445 * g(1) - t449 * g(2);
t450 = qJD(1) ^ 2;
t425 = qJDD(1) * pkin(1) + t450 * t470 + t433;
t434 = -t449 * g(1) - t445 * g(2);
t426 = -t450 * pkin(1) + qJDD(1) * t470 + t434;
t463 = t425 * t466 + t448 * t426;
t404 = -g(3) * t468 + t463;
t437 = t441 * qJD(1) + qJD(2);
t462 = qJD(1) * t440;
t459 = t444 * t462;
t423 = t437 * mrSges(3,1) - mrSges(3,3) * t459;
t427 = (-mrSges(3,1) * t448 + mrSges(3,2) * t444) * t462;
t436 = t441 * qJDD(1) + qJDD(2);
t428 = (-pkin(2) * t448 - pkin(7) * t444) * t462;
t435 = t437 ^ 2;
t461 = qJD(1) * t448;
t390 = -t435 * pkin(2) + t436 * pkin(7) + (-g(3) * t444 + t428 * t461) * t440 + t463;
t429 = (qJDD(1) * t444 + t448 * t460) * t440;
t391 = t430 * pkin(2) - t429 * pkin(7) - t469 + (-t425 + (pkin(2) * t444 - pkin(7) * t448) * t437 * qJD(1)) * t440;
t443 = sin(qJ(3));
t447 = cos(qJ(3));
t379 = t447 * t390 + t443 * t391;
t419 = t443 * t437 + t447 * t459;
t401 = -t419 * qJD(3) - t443 * t429 + t447 * t436;
t418 = t447 * t437 - t443 * t459;
t405 = -t418 * mrSges(4,1) + t419 * mrSges(4,2);
t458 = t440 * t461;
t432 = qJD(3) - t458;
t410 = t432 * mrSges(4,1) - t419 * mrSges(4,3);
t422 = qJDD(3) + t430;
t406 = -t418 * pkin(3) - t419 * pkin(8);
t431 = t432 ^ 2;
t376 = -t431 * pkin(3) + t422 * pkin(8) + t418 * t406 + t379;
t403 = -g(3) * t467 + t425 * t465 - t444 * t426;
t389 = -t436 * pkin(2) - t435 * pkin(7) + t428 * t459 - t403;
t402 = t418 * qJD(3) + t447 * t429 + t443 * t436;
t377 = (-t418 * t432 - t402) * pkin(8) + (t419 * t432 - t401) * pkin(3) + t389;
t442 = sin(qJ(4));
t446 = cos(qJ(4));
t373 = -t442 * t376 + t446 * t377;
t407 = -t442 * t419 + t446 * t432;
t382 = t407 * qJD(4) + t446 * t402 + t442 * t422;
t408 = t446 * t419 + t442 * t432;
t392 = -t407 * mrSges(5,1) + t408 * mrSges(5,2);
t417 = qJD(4) - t418;
t393 = -t417 * mrSges(5,2) + t407 * mrSges(5,3);
t399 = qJDD(4) - t401;
t371 = m(5) * t373 + t399 * mrSges(5,1) - t382 * mrSges(5,3) - t408 * t392 + t417 * t393;
t374 = t446 * t376 + t442 * t377;
t381 = -t408 * qJD(4) - t442 * t402 + t446 * t422;
t394 = t417 * mrSges(5,1) - t408 * mrSges(5,3);
t372 = m(5) * t374 - t399 * mrSges(5,2) + t381 * mrSges(5,3) + t407 * t392 - t417 * t394;
t455 = -t442 * t371 + t446 * t372;
t364 = m(4) * t379 - t422 * mrSges(4,2) + t401 * mrSges(4,3) + t418 * t405 - t432 * t410 + t455;
t378 = -t443 * t390 + t447 * t391;
t409 = -t432 * mrSges(4,2) + t418 * mrSges(4,3);
t375 = -t422 * pkin(3) - t431 * pkin(8) + t419 * t406 - t378;
t452 = -m(5) * t375 + t381 * mrSges(5,1) - t382 * mrSges(5,2) + t407 * t393 - t408 * t394;
t369 = m(4) * t378 + t422 * mrSges(4,1) - t402 * mrSges(4,3) - t419 * t405 + t432 * t409 + t452;
t456 = t447 * t364 - t443 * t369;
t355 = m(3) * t404 - t436 * mrSges(3,2) - t430 * mrSges(3,3) - t437 * t423 + t427 * t458 + t456;
t358 = t443 * t364 + t447 * t369;
t414 = -t440 * t425 - t469;
t424 = -t437 * mrSges(3,2) + mrSges(3,3) * t458;
t357 = m(3) * t414 + t430 * mrSges(3,1) + t429 * mrSges(3,2) + (t423 * t444 - t424 * t448) * t462 + t358;
t365 = t446 * t371 + t442 * t372;
t451 = -m(4) * t389 + t401 * mrSges(4,1) - t402 * mrSges(4,2) + t418 * t409 - t419 * t410 - t365;
t361 = m(3) * t403 + t436 * mrSges(3,1) - t429 * mrSges(3,3) + t437 * t424 - t427 * t459 + t451;
t345 = t355 * t466 - t440 * t357 + t361 * t465;
t343 = m(2) * t433 + qJDD(1) * mrSges(2,1) - t450 * mrSges(2,2) + t345;
t349 = t448 * t355 - t444 * t361;
t348 = m(2) * t434 - t450 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t349;
t464 = t449 * t343 + t445 * t348;
t344 = t355 * t468 + t441 * t357 + t361 * t467;
t457 = -t445 * t343 + t449 * t348;
t383 = Ifges(5,5) * t408 + Ifges(5,6) * t407 + Ifges(5,3) * t417;
t385 = Ifges(5,1) * t408 + Ifges(5,4) * t407 + Ifges(5,5) * t417;
t366 = -mrSges(5,1) * t375 + mrSges(5,3) * t374 + Ifges(5,4) * t382 + Ifges(5,2) * t381 + Ifges(5,6) * t399 - t408 * t383 + t417 * t385;
t384 = Ifges(5,4) * t408 + Ifges(5,2) * t407 + Ifges(5,6) * t417;
t367 = mrSges(5,2) * t375 - mrSges(5,3) * t373 + Ifges(5,1) * t382 + Ifges(5,4) * t381 + Ifges(5,5) * t399 + t407 * t383 - t417 * t384;
t395 = Ifges(4,5) * t419 + Ifges(4,6) * t418 + Ifges(4,3) * t432;
t396 = Ifges(4,4) * t419 + Ifges(4,2) * t418 + Ifges(4,6) * t432;
t350 = mrSges(4,2) * t389 - mrSges(4,3) * t378 + Ifges(4,1) * t402 + Ifges(4,4) * t401 + Ifges(4,5) * t422 - pkin(8) * t365 - t442 * t366 + t446 * t367 + t418 * t395 - t432 * t396;
t397 = Ifges(4,1) * t419 + Ifges(4,4) * t418 + Ifges(4,5) * t432;
t351 = -mrSges(4,1) * t389 - mrSges(5,1) * t373 + mrSges(5,2) * t374 + mrSges(4,3) * t379 + Ifges(4,4) * t402 - Ifges(5,5) * t382 + Ifges(4,2) * t401 + Ifges(4,6) * t422 - Ifges(5,6) * t381 - Ifges(5,3) * t399 - pkin(3) * t365 - t408 * t384 + t407 * t385 - t419 * t395 + t432 * t397;
t411 = Ifges(3,3) * t437 + (Ifges(3,5) * t444 + Ifges(3,6) * t448) * t462;
t412 = Ifges(3,6) * t437 + (Ifges(3,4) * t444 + Ifges(3,2) * t448) * t462;
t340 = mrSges(3,2) * t414 - mrSges(3,3) * t403 + Ifges(3,1) * t429 - Ifges(3,4) * t430 + Ifges(3,5) * t436 - pkin(7) * t358 + t447 * t350 - t443 * t351 + t411 * t458 - t437 * t412;
t413 = Ifges(3,5) * t437 + (Ifges(3,1) * t444 + Ifges(3,4) * t448) * t462;
t341 = Ifges(3,4) * t429 - Ifges(3,2) * t430 + Ifges(3,6) * t436 - t411 * t459 + t437 * t413 - mrSges(3,1) * t414 + mrSges(3,3) * t404 - Ifges(4,5) * t402 - Ifges(4,6) * t401 - Ifges(4,3) * t422 - t419 * t396 + t418 * t397 - mrSges(4,1) * t378 + mrSges(4,2) * t379 - t442 * t367 - t446 * t366 - pkin(3) * t452 - pkin(8) * t455 - pkin(2) * t358;
t453 = pkin(6) * t349 + t340 * t444 + t341 * t448;
t339 = Ifges(3,5) * t429 - Ifges(3,6) * t430 + Ifges(3,3) * t436 + mrSges(3,1) * t403 - mrSges(3,2) * t404 + t443 * t350 + t447 * t351 + pkin(2) * t451 + pkin(7) * t456 + (t412 * t444 - t413 * t448) * t462;
t338 = -mrSges(2,2) * g(3) - mrSges(2,3) * t433 + Ifges(2,5) * qJDD(1) - t450 * Ifges(2,6) + t448 * t340 - t444 * t341 + (-t344 * t440 - t345 * t441) * pkin(6);
t337 = mrSges(2,1) * g(3) + mrSges(2,3) * t434 + t450 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t344 - t440 * t339 + t453 * t441;
t1 = [-m(1) * g(1) + t457; -m(1) * g(2) + t464; (-m(1) - m(2)) * g(3) + t344; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t464 - t445 * t337 + t449 * t338; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t457 + t449 * t337 + t445 * t338; -mrSges(1,1) * g(2) + mrSges(2,1) * t433 + mrSges(1,2) * g(1) - mrSges(2,2) * t434 + Ifges(2,3) * qJDD(1) + pkin(1) * t345 + t441 * t339 + t453 * t440;];
tauB = t1;
