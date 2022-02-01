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
% m [6x1]
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:01:58
% EndTime: 2022-01-20 12:02:01
% DurationCPUTime: 3.00s
% Computational Cost: add. (57651->232), mult. (59398->293), div. (0->0), fcn. (33553->10), ass. (0->98)
t462 = -m(3) - m(4);
t435 = qJD(1) + qJD(2);
t429 = qJD(3) + t435;
t427 = t429 ^ 2;
t461 = pkin(4) * t427;
t438 = sin(qJ(4));
t460 = t429 * t438;
t443 = cos(qJ(4));
t459 = t429 * t443;
t441 = sin(qJ(1));
t446 = cos(qJ(1));
t425 = t441 * g(1) - g(2) * t446;
t423 = qJDD(1) * pkin(1) + t425;
t426 = -g(1) * t446 - g(2) * t441;
t447 = qJD(1) ^ 2;
t424 = -pkin(1) * t447 + t426;
t440 = sin(qJ(2));
t445 = cos(qJ(2));
t403 = t445 * t423 - t424 * t440;
t433 = qJDD(1) + qJDD(2);
t401 = pkin(2) * t433 + t403;
t404 = t440 * t423 + t445 * t424;
t431 = t435 ^ 2;
t402 = -pkin(2) * t431 + t404;
t439 = sin(qJ(3));
t444 = cos(qJ(3));
t389 = t439 * t401 + t444 * t402;
t428 = qJDD(3) + t433;
t387 = -pkin(3) * t427 + pkin(8) * t428 + t389;
t458 = t438 * t387;
t456 = qJD(4) * t429;
t415 = t428 * t438 + t443 * t456;
t380 = qJDD(4) * pkin(4) - t415 * pkin(9) - t458 + (pkin(9) * t456 + t438 * t461 - g(3)) * t443;
t384 = -g(3) * t438 + t443 * t387;
t416 = t428 * t443 - t438 * t456;
t422 = qJD(4) * pkin(4) - pkin(9) * t460;
t436 = t443 ^ 2;
t381 = pkin(9) * t416 - qJD(4) * t422 - t436 * t461 + t384;
t437 = sin(qJ(5));
t442 = cos(qJ(5));
t378 = t380 * t442 - t381 * t437;
t410 = (-t437 * t438 + t442 * t443) * t429;
t392 = qJD(5) * t410 + t415 * t442 + t416 * t437;
t411 = (t437 * t443 + t438 * t442) * t429;
t397 = -mrSges(6,1) * t410 + mrSges(6,2) * t411;
t434 = qJD(4) + qJD(5);
t405 = -mrSges(6,2) * t434 + mrSges(6,3) * t410;
t432 = qJDD(4) + qJDD(5);
t376 = m(6) * t378 + mrSges(6,1) * t432 - mrSges(6,3) * t392 - t397 * t411 + t405 * t434;
t379 = t380 * t437 + t381 * t442;
t391 = -qJD(5) * t411 - t415 * t437 + t416 * t442;
t406 = mrSges(6,1) * t434 - mrSges(6,3) * t411;
t377 = m(6) * t379 - mrSges(6,2) * t432 + mrSges(6,3) * t391 + t397 * t410 - t406 * t434;
t368 = t442 * t376 + t437 * t377;
t383 = -t443 * g(3) - t458;
t414 = (-mrSges(5,1) * t443 + mrSges(5,2) * t438) * t429;
t421 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t459;
t366 = m(5) * t383 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t415 + qJD(4) * t421 - t414 * t460 + t368;
t420 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t460;
t451 = -t376 * t437 + t442 * t377;
t367 = m(5) * t384 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t416 - qJD(4) * t420 + t414 * t459 + t451;
t452 = -t366 * t438 + t443 * t367;
t361 = m(4) * t389 - mrSges(4,1) * t427 - mrSges(4,2) * t428 + t452;
t388 = t401 * t444 - t439 * t402;
t450 = -pkin(3) * t428 - t388;
t386 = -pkin(8) * t427 + t450;
t382 = t422 * t460 - t416 * pkin(4) + (-pkin(9) * t436 - pkin(8)) * t427 + t450;
t449 = m(6) * t382 - t391 * mrSges(6,1) + mrSges(6,2) * t392 - t410 * t405 + t406 * t411;
t448 = -m(5) * t386 + t416 * mrSges(5,1) - mrSges(5,2) * t415 - t420 * t460 + t421 * t459 - t449;
t372 = m(4) * t388 + mrSges(4,1) * t428 - mrSges(4,2) * t427 + t448;
t357 = t439 * t361 + t444 * t372;
t355 = m(3) * t403 + mrSges(3,1) * t433 - mrSges(3,2) * t431 + t357;
t453 = t444 * t361 - t372 * t439;
t356 = m(3) * t404 - mrSges(3,1) * t431 - mrSges(3,2) * t433 + t453;
t349 = t445 * t355 + t440 * t356;
t347 = m(2) * t425 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t447 + t349;
t454 = -t355 * t440 + t445 * t356;
t348 = m(2) * t426 - mrSges(2,1) * t447 - qJDD(1) * mrSges(2,2) + t454;
t457 = t446 * t347 + t441 * t348;
t362 = t443 * t366 + t438 * t367;
t455 = -t347 * t441 + t446 * t348;
t409 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t438 + Ifges(5,4) * t443) * t429;
t408 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t438 + Ifges(5,2) * t443) * t429;
t407 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t438 + Ifges(5,6) * t443) * t429;
t395 = Ifges(6,1) * t411 + Ifges(6,4) * t410 + Ifges(6,5) * t434;
t394 = Ifges(6,4) * t411 + Ifges(6,2) * t410 + Ifges(6,6) * t434;
t393 = Ifges(6,5) * t411 + Ifges(6,6) * t410 + Ifges(6,3) * t434;
t370 = mrSges(6,2) * t382 - mrSges(6,3) * t378 + Ifges(6,1) * t392 + Ifges(6,4) * t391 + Ifges(6,5) * t432 + t393 * t410 - t394 * t434;
t369 = -mrSges(6,1) * t382 + mrSges(6,3) * t379 + Ifges(6,4) * t392 + Ifges(6,2) * t391 + Ifges(6,6) * t432 - t393 * t411 + t395 * t434;
t358 = mrSges(5,2) * t386 - mrSges(5,3) * t383 + Ifges(5,1) * t415 + Ifges(5,4) * t416 + Ifges(5,5) * qJDD(4) - pkin(9) * t368 - qJD(4) * t408 - t369 * t437 + t370 * t442 + t407 * t459;
t351 = -mrSges(5,1) * t386 + mrSges(5,3) * t384 + Ifges(5,4) * t415 + Ifges(5,2) * t416 + Ifges(5,6) * qJDD(4) - pkin(4) * t449 + pkin(9) * t451 + qJD(4) * t409 + t442 * t369 + t437 * t370 - t407 * t460;
t350 = Ifges(4,6) * t428 + t427 * Ifges(4,5) + mrSges(4,1) * g(3) + mrSges(4,3) * t389 - Ifges(5,5) * t415 - Ifges(5,6) * t416 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t383 + mrSges(5,2) * t384 - Ifges(6,5) * t392 - Ifges(6,6) * t391 - Ifges(6,3) * t432 - t411 * t394 + t410 * t395 - mrSges(6,1) * t378 + mrSges(6,2) * t379 - pkin(4) * t368 - pkin(3) * t362 + (-t408 * t438 + t409 * t443) * t429;
t343 = -mrSges(4,2) * g(3) - mrSges(4,3) * t388 + Ifges(4,5) * t428 - Ifges(4,6) * t427 - pkin(8) * t362 - t351 * t438 + t358 * t443;
t342 = -mrSges(3,2) * g(3) - mrSges(3,3) * t403 + Ifges(3,5) * t433 - Ifges(3,6) * t431 - pkin(7) * t357 + t343 * t444 - t350 * t439;
t341 = Ifges(3,6) * t433 + t431 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t404 + t439 * t343 + t444 * t350 - pkin(2) * (-m(4) * g(3) + t362) + pkin(7) * t453;
t340 = -mrSges(2,2) * g(3) - mrSges(2,3) * t425 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t447 - pkin(6) * t349 - t341 * t440 + t342 * t445;
t339 = Ifges(2,6) * qJDD(1) + t447 * Ifges(2,5) + mrSges(2,3) * t426 + t440 * t342 + t445 * t341 - pkin(1) * t362 + pkin(6) * t454 + (-pkin(1) * t462 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t455; -m(1) * g(2) + t457; (-m(1) - m(2) + t462) * g(3) + t362; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t457 - t441 * t339 + t446 * t340; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t455 + t446 * t339 + t441 * t340; pkin(1) * t349 + mrSges(2,1) * t425 - mrSges(2,2) * t426 + pkin(2) * t357 + mrSges(3,1) * t403 - mrSges(3,2) * t404 + pkin(8) * t452 + t438 * t358 + t443 * t351 + pkin(3) * t448 + mrSges(4,1) * t388 - mrSges(4,2) * t389 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t428 + Ifges(3,3) * t433 + Ifges(2,3) * qJDD(1);];
tauB = t1;
