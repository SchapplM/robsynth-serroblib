% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:51
% EndTime: 2019-12-05 15:59:53
% DurationCPUTime: 1.28s
% Computational Cost: add. (11832->213), mult. (21747->264), div. (0->0), fcn. (12620->8), ass. (0->89)
t409 = sin(pkin(8));
t433 = cos(pkin(8));
t394 = -g(1) * t433 - g(2) * t409;
t406 = -g(3) + qJDD(1);
t412 = sin(qJ(2));
t415 = cos(qJ(2));
t378 = t415 * t394 + t412 * t406;
t422 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t378;
t439 = m(3) + m(4);
t438 = -pkin(2) - pkin(6);
t416 = qJD(2) ^ 2;
t437 = pkin(4) * t416;
t436 = mrSges(3,1) - mrSges(4,2);
t435 = Ifges(3,5) - Ifges(4,4);
t434 = (-Ifges(3,6) + Ifges(4,5));
t377 = -t412 * t394 + t415 * t406;
t421 = -t416 * qJ(3) + qJDD(3) - t377;
t372 = qJDD(2) * t438 + t421;
t414 = cos(qJ(4));
t369 = t414 * t372;
t411 = sin(qJ(4));
t428 = qJD(2) * qJD(4);
t392 = qJDD(2) * t414 - t411 * t428;
t393 = g(1) * t409 - g(2) * t433;
t354 = (qJDD(4) * pkin(4)) - t392 * pkin(7) + t369 + (-pkin(7) * t428 - t414 * t437 + t393) * t411;
t360 = t411 * t372 - t414 * t393;
t391 = -qJDD(2) * t411 - t414 * t428;
t429 = qJD(2) * t414;
t397 = (qJD(4) * pkin(4)) - pkin(7) * t429;
t405 = t411 ^ 2;
t355 = pkin(7) * t391 - qJD(4) * t397 - t405 * t437 + t360;
t410 = sin(qJ(5));
t413 = cos(qJ(5));
t352 = t354 * t413 - t355 * t410;
t382 = (-t410 * t414 - t411 * t413) * qJD(2);
t362 = qJD(5) * t382 + t391 * t410 + t392 * t413;
t383 = (-t410 * t411 + t413 * t414) * qJD(2);
t367 = -mrSges(6,1) * t382 + mrSges(6,2) * t383;
t402 = qJD(4) + qJD(5);
t375 = -mrSges(6,2) * t402 + mrSges(6,3) * t382;
t401 = qJDD(4) + qJDD(5);
t350 = m(6) * t352 + mrSges(6,1) * t401 - mrSges(6,3) * t362 - t367 * t383 + t375 * t402;
t353 = t354 * t410 + t355 * t413;
t361 = -qJD(5) * t383 + t391 * t413 - t392 * t410;
t376 = mrSges(6,1) * t402 - mrSges(6,3) * t383;
t351 = m(6) * t353 - mrSges(6,2) * t401 + mrSges(6,3) * t361 + t367 * t382 - t376 * t402;
t341 = t413 * t350 + t410 * t351;
t359 = t411 * t393 + t369;
t390 = (mrSges(5,1) * t411 + mrSges(5,2) * t414) * qJD(2);
t430 = qJD(2) * t411;
t395 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t430;
t339 = m(5) * t359 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t392 + qJD(4) * t395 - t390 * t429 + t341;
t396 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t429;
t423 = -t350 * t410 + t413 * t351;
t340 = m(5) * t360 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t391 - qJD(4) * t396 - t390 * t430 + t423;
t336 = t414 * t339 + t411 * t340;
t374 = -qJDD(2) * pkin(2) + t421;
t419 = -m(4) * t374 + (t416 * mrSges(4,3)) - t336;
t332 = m(3) * t377 - (t416 * mrSges(3,2)) + qJDD(2) * t436 + t419;
t373 = pkin(2) * t416 - t422;
t371 = t416 * t438 + t422;
t357 = t397 * t429 - t391 * pkin(4) + (-pkin(7) * t405 + t438) * t416 + t422;
t420 = m(6) * t357 - mrSges(6,1) * t361 + t362 * mrSges(6,2) - t375 * t382 + t383 * t376;
t418 = -m(5) * t371 + mrSges(5,1) * t391 - t392 * mrSges(5,2) - t395 * t430 - t396 * t429 - t420;
t417 = -m(4) * t373 + (t416 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t418;
t346 = m(3) * t378 - (mrSges(3,1) * t416) - qJDD(2) * mrSges(3,2) + t417;
t424 = -t332 * t412 + t415 * t346;
t327 = m(2) * t394 + t424;
t431 = -t411 * t339 + t414 * t340;
t334 = (m(2) + t439) * t393 - t431;
t432 = t327 * t409 + t334 * t433;
t329 = t332 * t415 + t412 * t346;
t425 = t327 * t433 - t334 * t409;
t381 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t414 - Ifges(5,4) * t411) * qJD(2);
t380 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t414 - Ifges(5,2) * t411) * qJD(2);
t379 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t414 - Ifges(5,6) * t411) * qJD(2);
t365 = Ifges(6,1) * t383 + Ifges(6,4) * t382 + Ifges(6,5) * t402;
t364 = Ifges(6,4) * t383 + Ifges(6,2) * t382 + Ifges(6,6) * t402;
t363 = Ifges(6,5) * t383 + Ifges(6,6) * t382 + Ifges(6,3) * t402;
t343 = mrSges(6,2) * t357 - mrSges(6,3) * t352 + Ifges(6,1) * t362 + Ifges(6,4) * t361 + Ifges(6,5) * t401 + t363 * t382 - t364 * t402;
t342 = -mrSges(6,1) * t357 + mrSges(6,3) * t353 + Ifges(6,4) * t362 + Ifges(6,2) * t361 + Ifges(6,6) * t401 - t363 * t383 + t365 * t402;
t335 = -m(4) * t393 + t431;
t330 = mrSges(5,2) * t371 - mrSges(5,3) * t359 + Ifges(5,1) * t392 + Ifges(5,4) * t391 + (Ifges(5,5) * qJDD(4)) - pkin(7) * t341 - qJD(4) * t380 - t342 * t410 + t343 * t413 - t379 * t430;
t328 = -mrSges(5,1) * t371 + mrSges(5,3) * t360 + Ifges(5,4) * t392 + Ifges(5,2) * t391 + (Ifges(5,6) * qJDD(4)) - pkin(4) * t420 + pkin(7) * t423 + qJD(4) * t381 + t413 * t342 + t410 * t343 - t379 * t429;
t324 = (Ifges(5,3) * qJDD(4)) + Ifges(6,3) * t401 + Ifges(5,5) * t392 - t382 * t365 + t383 * t364 + Ifges(5,6) * t391 + mrSges(4,1) * t374 - mrSges(3,3) * t377 + mrSges(5,1) * t359 - mrSges(5,2) * t360 + Ifges(6,6) * t361 + Ifges(6,5) * t362 + mrSges(6,1) * t352 - mrSges(6,2) * t353 + pkin(4) * t341 - qJ(3) * t335 + pkin(3) * t336 + (t434 * t416) + (-mrSges(3,2) + mrSges(4,3)) * t393 + t435 * qJDD(2) + (t380 * t414 + t381 * t411) * qJD(2);
t323 = -mrSges(4,1) * t373 + mrSges(3,3) * t378 - pkin(2) * t335 - pkin(3) * t418 - pkin(6) * t431 - qJDD(2) * t434 - t414 * t328 - t411 * t330 + t393 * t436 + t416 * t435;
t322 = -pkin(1) * t329 + mrSges(2,3) * t394 - pkin(2) * t419 - qJ(3) * t417 + t411 * t328 + pkin(6) * t336 - mrSges(3,1) * t377 + mrSges(3,2) * t378 - mrSges(4,2) * t374 + mrSges(4,3) * t373 - t414 * t330 - mrSges(2,1) * t406 + (mrSges(4,2) * pkin(2) - Ifges(4,1) - Ifges(3,3)) * qJDD(2);
t321 = mrSges(2,2) * t406 - mrSges(2,3) * t393 - pkin(5) * t329 - t323 * t412 + t324 * t415;
t1 = [-m(1) * g(1) + t425; -m(1) * g(2) + t432; -m(1) * g(3) + m(2) * t406 + t329; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t432 + t321 * t433 - t322 * t409; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t425 + t409 * t321 + t322 * t433; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t394 + t412 * t324 + t415 * t323 - pkin(1) * t431 + pkin(5) * t424 + (pkin(1) * t439 + mrSges(2,1)) * t393;];
tauB = t1;
