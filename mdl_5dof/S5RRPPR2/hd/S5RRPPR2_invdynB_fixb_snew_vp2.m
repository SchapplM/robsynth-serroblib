% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:40
% EndTime: 2022-01-20 10:05:43
% DurationCPUTime: 2.64s
% Computational Cost: add. (34334->207), mult. (45201->270), div. (0->0), fcn. (25113->10), ass. (0->98)
t454 = 2 * qJD(4);
t422 = sin(qJ(1));
t425 = cos(qJ(1));
t404 = t422 * g(1) - g(2) * t425;
t399 = qJDD(1) * pkin(1) + t404;
t405 = -g(1) * t425 - g(2) * t422;
t426 = qJD(1) ^ 2;
t400 = -pkin(1) * t426 + t405;
t421 = sin(qJ(2));
t424 = cos(qJ(2));
t385 = t424 * t399 - t400 * t421;
t413 = qJDD(1) + qJDD(2);
t380 = pkin(2) * t413 + t385;
t386 = t421 * t399 + t424 * t400;
t414 = (qJD(1) + qJD(2));
t412 = t414 ^ 2;
t381 = -pkin(2) * t412 + t386;
t417 = sin(pkin(8));
t419 = cos(pkin(8));
t376 = t417 * t380 + t419 * t381;
t374 = -pkin(3) * t412 + qJ(4) * t413 + t376;
t453 = (t414 * t454) + t374;
t416 = sin(pkin(9));
t452 = mrSges(5,2) * t416;
t450 = mrSges(5,3) * t413;
t449 = t416 * t414;
t420 = sin(qJ(5));
t448 = t416 * t420;
t423 = cos(qJ(5));
t447 = t416 * t423;
t418 = cos(pkin(9));
t446 = t418 * t413;
t445 = t418 * t414;
t415 = -g(3) + qJDD(3);
t444 = t418 * t415;
t370 = t416 * t415 + t418 * t453;
t393 = (-mrSges(5,1) * t418 + t452) * t414;
t432 = -pkin(4) * t418 - pkin(7) * t416;
t395 = t432 * t414;
t368 = t395 * t445 + t370;
t375 = t419 * t380 - t417 * t381;
t428 = -t412 * qJ(4) + qJDD(4) - t375;
t371 = (-pkin(3) + t432) * t413 + t428;
t365 = -t368 * t420 + t371 * t423;
t402 = qJD(5) - t445;
t440 = t414 * t448;
t388 = -mrSges(6,2) * t402 - mrSges(6,3) * t440;
t390 = (mrSges(6,1) * t420 + mrSges(6,2) * t423) * t449;
t441 = qJD(5) * t414;
t392 = (t413 * t423 - t420 * t441) * t416;
t401 = qJDD(5) - t446;
t439 = t414 * t447;
t363 = m(6) * t365 + mrSges(6,1) * t401 - mrSges(6,3) * t392 + t388 * t402 - t390 * t439;
t366 = t368 * t423 + t371 * t420;
t389 = mrSges(6,1) * t402 - mrSges(6,3) * t439;
t391 = (-t413 * t420 - t423 * t441) * t416;
t364 = m(6) * t366 - mrSges(6,2) * t401 + mrSges(6,3) * t391 - t389 * t402 - t390 * t440;
t433 = -t420 * t363 + t423 * t364;
t356 = m(5) * t370 + (t393 * t414 + t450) * t418 + t433;
t369 = -t416 * t453 + t444;
t367 = -t444 + (t374 + (t454 + t395) * t414) * t416;
t429 = -m(6) * t367 + t391 * mrSges(6,1) - t392 * mrSges(6,2);
t361 = m(5) * t369 + (-t450 + (-t388 * t420 - t389 * t423 - t393) * t414) * t416 + t429;
t434 = t418 * t356 - t361 * t416;
t350 = m(4) * t376 - mrSges(4,1) * t412 - mrSges(4,2) * t413 + t434;
t357 = t423 * t363 + t420 * t364;
t373 = -pkin(3) * t413 + t428;
t427 = -m(5) * t373 + mrSges(5,1) * t446 - t357 + (t416 ^ 2 + t418 ^ 2) * mrSges(5,3) * t412;
t353 = m(4) * t375 - t412 * mrSges(4,2) + (mrSges(4,1) - t452) * t413 + t427;
t345 = t350 * t417 + t419 * t353;
t343 = m(3) * t385 + mrSges(3,1) * t413 - mrSges(3,2) * t412 + t345;
t435 = t350 * t419 - t353 * t417;
t344 = m(3) * t386 - mrSges(3,1) * t412 - mrSges(3,2) * t413 + t435;
t338 = t343 * t424 + t344 * t421;
t336 = m(2) * t404 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t426 + t338;
t436 = -t343 * t421 + t344 * t424;
t337 = m(2) * t405 - mrSges(2,1) * t426 - qJDD(1) * mrSges(2,2) + t436;
t443 = t336 * t425 + t337 * t422;
t351 = t416 * t356 + t418 * t361;
t438 = m(4) * t415 + t351;
t437 = -t336 * t422 + t337 * t425;
t431 = Ifges(5,1) * t416 + Ifges(5,4) * t418;
t430 = Ifges(5,5) * t416 + Ifges(5,6) * t418;
t394 = t430 * t414;
t384 = Ifges(6,5) * t402 + (Ifges(6,1) * t423 - Ifges(6,4) * t420) * t449;
t383 = Ifges(6,6) * t402 + (Ifges(6,4) * t423 - Ifges(6,2) * t420) * t449;
t382 = Ifges(6,3) * t402 + (Ifges(6,5) * t423 - Ifges(6,6) * t420) * t449;
t359 = mrSges(6,2) * t367 - mrSges(6,3) * t365 + Ifges(6,1) * t392 + Ifges(6,4) * t391 + Ifges(6,5) * t401 - t382 * t440 - t383 * t402;
t358 = -mrSges(6,1) * t367 + mrSges(6,3) * t366 + Ifges(6,4) * t392 + Ifges(6,2) * t391 + Ifges(6,6) * t401 - t382 * t439 + t384 * t402;
t347 = Ifges(5,2) * t446 - mrSges(5,1) * t373 - mrSges(6,1) * t365 + mrSges(6,2) * t366 + mrSges(5,3) * t370 - Ifges(6,5) * t392 - Ifges(6,6) * t391 - Ifges(6,3) * t401 - pkin(4) * t357 + (Ifges(5,4) * t413 + (-t383 * t423 - t384 * t420 - t394) * t414) * t416;
t346 = mrSges(5,2) * t373 - mrSges(5,3) * t369 - pkin(7) * t357 - t420 * t358 + t423 * t359 + t394 * t445 + t413 * t431;
t339 = t412 * Ifges(4,5) - mrSges(4,1) * t415 + mrSges(4,3) * t376 - mrSges(5,1) * t369 + mrSges(5,2) * t370 - t420 * t359 - t423 * t358 - pkin(4) * t429 - pkin(7) * t433 - pkin(3) * t351 + (Ifges(4,6) - t430) * t413 + (-pkin(4) * (-t388 * t448 - t389 * t447) + (-t416 * (Ifges(5,4) * t416 + Ifges(5,2) * t418) + t418 * t431) * t414) * t414;
t332 = mrSges(4,2) * t415 - mrSges(4,3) * t375 + Ifges(4,5) * t413 - Ifges(4,6) * t412 - qJ(4) * t351 + t346 * t418 - t347 * t416;
t331 = -mrSges(3,2) * g(3) - mrSges(3,3) * t385 + Ifges(3,5) * t413 - Ifges(3,6) * t412 - qJ(3) * t345 + t332 * t419 - t339 * t417;
t330 = mrSges(3,1) * g(3) + mrSges(3,3) * t386 + t412 * Ifges(3,5) + Ifges(3,6) * t413 - pkin(2) * t438 + qJ(3) * t435 + t417 * t332 + t419 * t339;
t329 = -mrSges(2,2) * g(3) - mrSges(2,3) * t404 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t426 - pkin(6) * t338 - t330 * t421 + t331 * t424;
t328 = Ifges(2,6) * qJDD(1) + t426 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t405 + t421 * t331 + t424 * t330 - pkin(1) * (-m(3) * g(3) + t438) + pkin(6) * t436;
t1 = [-m(1) * g(1) + t437; -m(1) * g(2) + t443; (-m(1) - m(2) - m(3)) * g(3) + t438; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t443 - t328 * t422 + t329 * t425; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t437 + t425 * t328 + t422 * t329; pkin(1) * t338 + mrSges(2,1) * t404 - mrSges(2,2) * t405 + pkin(2) * t345 + mrSges(3,1) * t385 - mrSges(3,2) * t386 + qJ(4) * t434 + t416 * t346 + t418 * t347 + pkin(3) * t427 + mrSges(4,1) * t375 - mrSges(4,2) * t376 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (-pkin(3) * t452 + Ifges(3,3) + Ifges(4,3)) * t413;];
tauB = t1;
