% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR3
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:45
% EndTime: 2022-01-23 09:20:47
% DurationCPUTime: 2.45s
% Computational Cost: add. (32310->206), mult. (45201->270), div. (0->0), fcn. (25113->10), ass. (0->99)
t455 = 2 * qJD(4);
t422 = sin(qJ(1));
t425 = cos(qJ(1));
t403 = t422 * g(1) - t425 * g(2);
t398 = qJDD(1) * pkin(1) + t403;
t404 = -t425 * g(1) - t422 * g(2);
t426 = qJD(1) ^ 2;
t399 = -t426 * pkin(1) + t404;
t417 = sin(pkin(8));
t419 = cos(pkin(8));
t384 = t419 * t398 - t417 * t399;
t379 = qJDD(1) * pkin(2) + t384;
t385 = t417 * t398 + t419 * t399;
t380 = -t426 * pkin(2) + t385;
t421 = sin(qJ(3));
t424 = cos(qJ(3));
t375 = t421 * t379 + t424 * t380;
t414 = (qJD(1) + qJD(3));
t412 = t414 ^ 2;
t413 = qJDD(1) + qJDD(3);
t373 = -t412 * pkin(3) + t413 * qJ(4) + t375;
t454 = (t414 * t455) + t373;
t416 = sin(pkin(9));
t453 = mrSges(5,2) * t416;
t451 = mrSges(5,3) * t413;
t450 = t416 * t414;
t420 = sin(qJ(5));
t449 = t416 * t420;
t423 = cos(qJ(5));
t448 = t416 * t423;
t418 = cos(pkin(9));
t447 = t418 * t413;
t446 = t418 * t414;
t415 = -g(3) + qJDD(2);
t445 = t418 * t415;
t369 = t416 * t415 + t454 * t418;
t392 = (-mrSges(5,1) * t418 + t453) * t414;
t432 = -pkin(4) * t418 - pkin(7) * t416;
t394 = t432 * t414;
t367 = t394 * t446 + t369;
t374 = t424 * t379 - t421 * t380;
t428 = -t412 * qJ(4) + qJDD(4) - t374;
t370 = (-pkin(3) + t432) * t413 + t428;
t364 = -t420 * t367 + t423 * t370;
t401 = qJD(5) - t446;
t441 = t414 * t449;
t387 = -t401 * mrSges(6,2) - mrSges(6,3) * t441;
t389 = (mrSges(6,1) * t420 + mrSges(6,2) * t423) * t450;
t442 = qJD(5) * t414;
t391 = (t413 * t423 - t420 * t442) * t416;
t400 = qJDD(5) - t447;
t440 = t414 * t448;
t362 = m(6) * t364 + t400 * mrSges(6,1) - t391 * mrSges(6,3) + t401 * t387 - t389 * t440;
t365 = t423 * t367 + t420 * t370;
t388 = t401 * mrSges(6,1) - mrSges(6,3) * t440;
t390 = (-t413 * t420 - t423 * t442) * t416;
t363 = m(6) * t365 - t400 * mrSges(6,2) + t390 * mrSges(6,3) - t401 * t388 - t389 * t441;
t434 = -t420 * t362 + t423 * t363;
t355 = m(5) * t369 + (t392 * t414 + t451) * t418 + t434;
t368 = -t454 * t416 + t445;
t366 = -t445 + (t373 + (t455 + t394) * t414) * t416;
t429 = -m(6) * t366 + t390 * mrSges(6,1) - t391 * mrSges(6,2);
t360 = m(5) * t368 + (-t451 + (-t387 * t420 - t388 * t423 - t392) * t414) * t416 + t429;
t435 = t418 * t355 - t416 * t360;
t349 = m(4) * t375 - t412 * mrSges(4,1) - t413 * mrSges(4,2) + t435;
t356 = t423 * t362 + t420 * t363;
t372 = -t413 * pkin(3) + t428;
t427 = -m(5) * t372 + mrSges(5,1) * t447 - t356 + (t416 ^ 2 + t418 ^ 2) * mrSges(5,3) * t412;
t352 = m(4) * t374 - t412 * mrSges(4,2) + (mrSges(4,1) - t453) * t413 + t427;
t344 = t421 * t349 + t424 * t352;
t342 = m(3) * t384 + qJDD(1) * mrSges(3,1) - t426 * mrSges(3,2) + t344;
t436 = t424 * t349 - t421 * t352;
t343 = m(3) * t385 - t426 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t436;
t337 = t419 * t342 + t417 * t343;
t335 = m(2) * t403 + qJDD(1) * mrSges(2,1) - t426 * mrSges(2,2) + t337;
t437 = -t417 * t342 + t419 * t343;
t336 = m(2) * t404 - t426 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t437;
t444 = t425 * t335 + t422 * t336;
t350 = t416 * t355 + t418 * t360;
t439 = m(4) * t415 + t350;
t438 = -t422 * t335 + t425 * t336;
t433 = m(3) * t415 + t439;
t431 = Ifges(5,1) * t416 + Ifges(5,4) * t418;
t430 = Ifges(5,5) * t416 + Ifges(5,6) * t418;
t393 = t430 * t414;
t383 = Ifges(6,5) * t401 + (Ifges(6,1) * t423 - Ifges(6,4) * t420) * t450;
t382 = Ifges(6,6) * t401 + (Ifges(6,4) * t423 - Ifges(6,2) * t420) * t450;
t381 = Ifges(6,3) * t401 + (Ifges(6,5) * t423 - Ifges(6,6) * t420) * t450;
t358 = mrSges(6,2) * t366 - mrSges(6,3) * t364 + Ifges(6,1) * t391 + Ifges(6,4) * t390 + Ifges(6,5) * t400 - t381 * t441 - t401 * t382;
t357 = -mrSges(6,1) * t366 + mrSges(6,3) * t365 + Ifges(6,4) * t391 + Ifges(6,2) * t390 + Ifges(6,6) * t400 - t381 * t440 + t401 * t383;
t346 = Ifges(5,2) * t447 - mrSges(5,1) * t372 - mrSges(6,1) * t364 + mrSges(6,2) * t365 + mrSges(5,3) * t369 - Ifges(6,5) * t391 - Ifges(6,6) * t390 - Ifges(6,3) * t400 - pkin(4) * t356 + (Ifges(5,4) * t413 + (-t382 * t423 - t383 * t420 - t393) * t414) * t416;
t345 = mrSges(5,2) * t372 - mrSges(5,3) * t368 - pkin(7) * t356 - t420 * t357 + t423 * t358 + t393 * t446 + t431 * t413;
t338 = t412 * Ifges(4,5) - mrSges(4,1) * t415 + mrSges(4,3) * t375 - mrSges(5,1) * t368 + mrSges(5,2) * t369 - t420 * t358 - t423 * t357 - pkin(4) * t429 - pkin(7) * t434 - pkin(3) * t350 + (Ifges(4,6) - t430) * t413 + (-pkin(4) * (-t387 * t449 - t388 * t448) + (-t416 * (Ifges(5,4) * t416 + Ifges(5,2) * t418) + t418 * t431) * t414) * t414;
t331 = mrSges(4,2) * t415 - mrSges(4,3) * t374 + Ifges(4,5) * t413 - t412 * Ifges(4,6) - qJ(4) * t350 + t418 * t345 - t416 * t346;
t330 = mrSges(3,2) * t415 - mrSges(3,3) * t384 + Ifges(3,5) * qJDD(1) - t426 * Ifges(3,6) - pkin(6) * t344 + t424 * t331 - t421 * t338;
t329 = -mrSges(3,1) * t415 + mrSges(3,3) * t385 + t426 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t439 + pkin(6) * t436 + t421 * t331 + t424 * t338;
t328 = -mrSges(2,2) * g(3) - mrSges(2,3) * t403 + Ifges(2,5) * qJDD(1) - t426 * Ifges(2,6) - qJ(2) * t337 - t417 * t329 + t419 * t330;
t327 = mrSges(2,1) * g(3) + mrSges(2,3) * t404 + t426 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t433 + qJ(2) * t437 + t419 * t329 + t417 * t330;
t1 = [-m(1) * g(1) + t438; -m(1) * g(2) + t444; (-m(1) - m(2)) * g(3) + t433; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t444 - t422 * t327 + t425 * t328; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t438 + t425 * t327 + t422 * t328; pkin(1) * t337 + mrSges(2,1) * t403 - mrSges(2,2) * t404 + pkin(2) * t344 + mrSges(3,1) * t384 - mrSges(3,2) * t385 + pkin(3) * (-t413 * t453 + t427) + qJ(4) * t435 + t416 * t345 + t418 * t346 + mrSges(4,1) * t374 - mrSges(4,2) * t375 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t413 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
