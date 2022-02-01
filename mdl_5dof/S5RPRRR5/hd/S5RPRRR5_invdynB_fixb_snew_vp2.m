% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:50
% EndTime: 2022-01-20 09:48:53
% DurationCPUTime: 3.01s
% Computational Cost: add. (43952->228), mult. (59398->290), div. (0->0), fcn. (33553->10), ass. (0->96)
t435 = qJD(1) + qJD(3);
t441 = sin(qJ(4));
t463 = t435 * t441;
t445 = cos(qJ(4));
t462 = t435 * t445;
t443 = sin(qJ(1));
t447 = cos(qJ(1));
t424 = t443 * g(1) - t447 * g(2);
t419 = qJDD(1) * pkin(1) + t424;
t425 = -t447 * g(1) - t443 * g(2);
t448 = qJD(1) ^ 2;
t420 = -t448 * pkin(1) + t425;
t438 = sin(pkin(9));
t439 = cos(pkin(9));
t402 = t439 * t419 - t438 * t420;
t400 = qJDD(1) * pkin(2) + t402;
t403 = t438 * t419 + t439 * t420;
t401 = -t448 * pkin(2) + t403;
t442 = sin(qJ(3));
t446 = cos(qJ(3));
t388 = t442 * t400 + t446 * t401;
t431 = t435 ^ 2;
t433 = qJDD(1) + qJDD(3);
t386 = -t431 * pkin(3) + t433 * pkin(7) + t388;
t437 = -g(3) + qJDD(2);
t382 = -t441 * t386 + t445 * t437;
t460 = qJD(4) * t435;
t458 = t445 * t460;
t414 = t441 * t433 + t458;
t379 = (-t414 + t458) * pkin(8) + (t431 * t441 * t445 + qJDD(4)) * pkin(4) + t382;
t383 = t445 * t386 + t441 * t437;
t415 = t445 * t433 - t441 * t460;
t423 = qJD(4) * pkin(4) - pkin(8) * t463;
t436 = t445 ^ 2;
t380 = -t436 * t431 * pkin(4) + t415 * pkin(8) - qJD(4) * t423 + t383;
t440 = sin(qJ(5));
t444 = cos(qJ(5));
t377 = t444 * t379 - t440 * t380;
t409 = (-t440 * t441 + t444 * t445) * t435;
t391 = t409 * qJD(5) + t444 * t414 + t440 * t415;
t410 = (t440 * t445 + t441 * t444) * t435;
t396 = -t409 * mrSges(6,1) + t410 * mrSges(6,2);
t434 = qJD(4) + qJD(5);
t404 = -t434 * mrSges(6,2) + t409 * mrSges(6,3);
t432 = qJDD(4) + qJDD(5);
t375 = m(6) * t377 + t432 * mrSges(6,1) - t391 * mrSges(6,3) - t410 * t396 + t434 * t404;
t378 = t440 * t379 + t444 * t380;
t390 = -t410 * qJD(5) - t440 * t414 + t444 * t415;
t405 = t434 * mrSges(6,1) - t410 * mrSges(6,3);
t376 = m(6) * t378 - t432 * mrSges(6,2) + t390 * mrSges(6,3) + t409 * t396 - t434 * t405;
t367 = t444 * t375 + t440 * t376;
t413 = (-mrSges(5,1) * t445 + mrSges(5,2) * t441) * t435;
t422 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t462;
t365 = m(5) * t382 + qJDD(4) * mrSges(5,1) - t414 * mrSges(5,3) + qJD(4) * t422 - t413 * t463 + t367;
t421 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t463;
t453 = -t440 * t375 + t444 * t376;
t366 = m(5) * t383 - qJDD(4) * mrSges(5,2) + t415 * mrSges(5,3) - qJD(4) * t421 + t413 * t462 + t453;
t454 = -t441 * t365 + t445 * t366;
t360 = m(4) * t388 - t431 * mrSges(4,1) - t433 * mrSges(4,2) + t454;
t387 = t446 * t400 - t442 * t401;
t451 = -t433 * pkin(3) - t387;
t385 = -t431 * pkin(7) + t451;
t381 = t423 * t463 - t415 * pkin(4) + (-pkin(8) * t436 - pkin(7)) * t431 + t451;
t450 = m(6) * t381 - t390 * mrSges(6,1) + t391 * mrSges(6,2) - t409 * t404 + t410 * t405;
t449 = -m(5) * t385 + t415 * mrSges(5,1) - t414 * mrSges(5,2) - t421 * t463 + t422 * t462 - t450;
t371 = m(4) * t387 + t433 * mrSges(4,1) - t431 * mrSges(4,2) + t449;
t356 = t442 * t360 + t446 * t371;
t354 = m(3) * t402 + qJDD(1) * mrSges(3,1) - t448 * mrSges(3,2) + t356;
t455 = t446 * t360 - t442 * t371;
t355 = m(3) * t403 - t448 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t455;
t348 = t439 * t354 + t438 * t355;
t346 = m(2) * t424 + qJDD(1) * mrSges(2,1) - t448 * mrSges(2,2) + t348;
t456 = -t438 * t354 + t439 * t355;
t347 = m(2) * t425 - t448 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t456;
t461 = t447 * t346 + t443 * t347;
t361 = t445 * t365 + t441 * t366;
t459 = m(4) * t437 + t361;
t457 = -t443 * t346 + t447 * t347;
t452 = m(3) * t437 + t459;
t408 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t441 + Ifges(5,4) * t445) * t435;
t407 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t441 + Ifges(5,2) * t445) * t435;
t406 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t441 + Ifges(5,6) * t445) * t435;
t394 = Ifges(6,1) * t410 + Ifges(6,4) * t409 + Ifges(6,5) * t434;
t393 = Ifges(6,4) * t410 + Ifges(6,2) * t409 + Ifges(6,6) * t434;
t392 = Ifges(6,5) * t410 + Ifges(6,6) * t409 + Ifges(6,3) * t434;
t369 = mrSges(6,2) * t381 - mrSges(6,3) * t377 + Ifges(6,1) * t391 + Ifges(6,4) * t390 + Ifges(6,5) * t432 + t409 * t392 - t434 * t393;
t368 = -mrSges(6,1) * t381 + mrSges(6,3) * t378 + Ifges(6,4) * t391 + Ifges(6,2) * t390 + Ifges(6,6) * t432 - t410 * t392 + t434 * t394;
t357 = mrSges(5,2) * t385 - mrSges(5,3) * t382 + Ifges(5,1) * t414 + Ifges(5,4) * t415 + Ifges(5,5) * qJDD(4) - pkin(8) * t367 - qJD(4) * t407 - t440 * t368 + t444 * t369 + t406 * t462;
t350 = -mrSges(5,1) * t385 + mrSges(5,3) * t383 + Ifges(5,4) * t414 + Ifges(5,2) * t415 + Ifges(5,6) * qJDD(4) - pkin(4) * t450 + pkin(8) * t453 + qJD(4) * t408 + t444 * t368 + t440 * t369 - t406 * t463;
t349 = Ifges(4,6) * t433 + t431 * Ifges(4,5) - mrSges(4,1) * t437 + mrSges(4,3) * t388 - Ifges(5,5) * t414 - Ifges(5,6) * t415 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t382 + mrSges(5,2) * t383 - Ifges(6,5) * t391 - Ifges(6,6) * t390 - Ifges(6,3) * t432 - t410 * t393 + t409 * t394 - mrSges(6,1) * t377 + mrSges(6,2) * t378 - pkin(4) * t367 - pkin(3) * t361 + (-t441 * t407 + t445 * t408) * t435;
t342 = mrSges(4,2) * t437 - mrSges(4,3) * t387 + Ifges(4,5) * t433 - t431 * Ifges(4,6) - pkin(7) * t361 - t441 * t350 + t445 * t357;
t341 = mrSges(3,2) * t437 - mrSges(3,3) * t402 + Ifges(3,5) * qJDD(1) - t448 * Ifges(3,6) - pkin(6) * t356 + t446 * t342 - t442 * t349;
t340 = -mrSges(3,1) * t437 + mrSges(3,3) * t403 + t448 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t459 + pkin(6) * t455 + t442 * t342 + t446 * t349;
t339 = -mrSges(2,2) * g(3) - mrSges(2,3) * t424 + Ifges(2,5) * qJDD(1) - t448 * Ifges(2,6) - qJ(2) * t348 - t438 * t340 + t439 * t341;
t338 = mrSges(2,1) * g(3) + mrSges(2,3) * t425 + t448 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t452 + qJ(2) * t456 + t439 * t340 + t438 * t341;
t1 = [-m(1) * g(1) + t457; -m(1) * g(2) + t461; (-m(1) - m(2)) * g(3) + t452; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t461 - t443 * t338 + t447 * t339; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t457 + t447 * t338 + t443 * t339; pkin(1) * t348 + mrSges(2,1) * t424 - mrSges(2,2) * t425 + pkin(2) * t356 + mrSges(3,1) * t402 - mrSges(3,2) * t403 + t441 * t357 + t445 * t350 + pkin(3) * t449 + pkin(7) * t454 + mrSges(4,1) * t387 - mrSges(4,2) * t388 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t433 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
