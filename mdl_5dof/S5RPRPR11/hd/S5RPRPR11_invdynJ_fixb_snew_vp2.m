% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:23
% EndTime: 2019-12-31 18:27:25
% DurationCPUTime: 1.45s
% Computational Cost: add. (7003->230), mult. (17008->280), div. (0->0), fcn. (11632->8), ass. (0->103)
t420 = cos(pkin(8));
t427 = qJD(1) ^ 2;
t423 = sin(qJ(1));
t425 = cos(qJ(1));
t450 = t423 * g(1) - t425 * g(2);
t442 = -qJDD(2) + t450;
t414 = t420 ^ 2;
t419 = sin(pkin(8));
t449 = t419 ^ 2 + t414;
t385 = -(t449 * pkin(6) + qJ(2)) * t427 - (pkin(2) * t420 + pkin(1)) * qJDD(1) - t442;
t422 = sin(qJ(3));
t462 = cos(qJ(3));
t434 = t462 * t419 + t420 * t422;
t443 = t420 * t462;
t448 = qJD(1) * t419;
t400 = -qJD(1) * t443 + t422 * t448;
t446 = t400 * qJD(3);
t387 = t434 * qJDD(1) - t446;
t470 = t385 + (-t387 + t446) * qJ(4);
t469 = Ifges(4,1) + Ifges(5,1);
t459 = Ifges(4,4) - Ifges(5,5);
t458 = Ifges(4,5) + Ifges(5,4);
t468 = -Ifges(4,2) - Ifges(5,3);
t457 = Ifges(4,6) - Ifges(5,6);
t467 = Ifges(4,3) + Ifges(5,2);
t463 = 2 * qJD(4);
t461 = pkin(2) * t427;
t460 = -mrSges(4,3) - mrSges(5,2);
t456 = pkin(6) * qJDD(1);
t437 = -g(1) * t425 - g(2) * t423;
t402 = -pkin(1) * t427 + qJDD(1) * qJ(2) + t437;
t445 = qJD(1) * qJD(2);
t441 = -g(3) * t420 - 0.2e1 * t419 * t445;
t370 = (t420 * t461 - t402 - t456) * t419 + t441;
t389 = -g(3) * t419 + (t402 + 0.2e1 * t445) * t420;
t371 = -t414 * t461 + t420 * t456 + t389;
t357 = t422 * t370 + t462 * t371;
t401 = t434 * qJD(1);
t447 = qJD(3) * t401;
t386 = t447 + (t419 * t422 - t443) * qJDD(1);
t393 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t401;
t381 = pkin(3) * t400 - qJ(4) * t401;
t426 = qJD(3) ^ 2;
t351 = -pkin(3) * t426 + qJDD(3) * qJ(4) + qJD(3) * t463 - t400 * t381 + t357;
t394 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t401;
t356 = t462 * t370 - t422 * t371;
t352 = -qJDD(3) * pkin(3) - t426 * qJ(4) + t401 * t381 + qJDD(4) - t356;
t346 = (-t387 - t446) * pkin(7) + (t400 * t401 - qJDD(3)) * pkin(4) + t352;
t396 = -qJD(3) * pkin(4) - pkin(7) * t401;
t399 = t400 ^ 2;
t347 = -pkin(4) * t399 + pkin(7) * t386 + qJD(3) * t396 + t351;
t421 = sin(qJ(5));
t424 = cos(qJ(5));
t343 = t346 * t424 - t347 * t421;
t378 = t400 * t424 - t401 * t421;
t355 = qJD(5) * t378 + t386 * t421 + t387 * t424;
t379 = t400 * t421 + t401 * t424;
t362 = -mrSges(6,1) * t378 + mrSges(6,2) * t379;
t415 = -qJD(3) + qJD(5);
t366 = -mrSges(6,2) * t415 + mrSges(6,3) * t378;
t412 = -qJDD(3) + qJDD(5);
t341 = m(6) * t343 + mrSges(6,1) * t412 - mrSges(6,3) * t355 - t362 * t379 + t366 * t415;
t344 = t346 * t421 + t347 * t424;
t354 = -qJD(5) * t379 + t386 * t424 - t387 * t421;
t367 = mrSges(6,1) * t415 - mrSges(6,3) * t379;
t342 = m(6) * t344 - mrSges(6,2) * t412 + mrSges(6,3) * t354 + t362 * t378 - t367 * t415;
t439 = -t341 * t421 + t424 * t342;
t433 = m(5) * t351 + qJDD(3) * mrSges(5,3) + qJD(3) * t394 + t439;
t382 = mrSges(5,1) * t400 - mrSges(5,3) * t401;
t451 = -mrSges(4,1) * t400 - mrSges(4,2) * t401 - t382;
t332 = m(4) * t357 - qJDD(3) * mrSges(4,2) - qJD(3) * t393 + t460 * t386 + t451 * t400 + t433;
t392 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t400;
t335 = t341 * t424 + t342 * t421;
t395 = -mrSges(5,2) * t400 + qJD(3) * mrSges(5,3);
t431 = -m(5) * t352 + qJDD(3) * mrSges(5,1) + qJD(3) * t395 - t335;
t333 = m(4) * t356 + qJDD(3) * mrSges(4,1) + qJD(3) * t392 + t460 * t387 + t451 * t401 + t431;
t455 = t422 * t332 + t462 * t333;
t454 = -t467 * qJD(3) + t457 * t400 - t458 * t401;
t453 = t457 * qJD(3) + t468 * t400 + t459 * t401;
t452 = t458 * qJD(3) - t459 * t400 + t469 * t401;
t440 = t462 * t332 - t422 * t333;
t436 = -mrSges(3,1) * t420 + mrSges(3,2) * t419;
t435 = mrSges(3,3) * qJDD(1) + t427 * t436;
t345 = -pkin(7) * t399 + (-pkin(3) - pkin(4)) * t386 + (-pkin(3) * qJD(3) + t396 + t463) * t401 - t470;
t432 = -m(6) * t345 + t354 * mrSges(6,1) - t355 * mrSges(6,2) + t378 * t366 - t379 * t367;
t349 = -0.2e1 * qJD(4) * t401 + (t386 + t447) * pkin(3) + t470;
t430 = m(5) * t349 + t386 * mrSges(5,1) + t400 * t395 + t432;
t359 = Ifges(6,4) * t379 + Ifges(6,2) * t378 + Ifges(6,6) * t415;
t360 = Ifges(6,1) * t379 + Ifges(6,4) * t378 + Ifges(6,5) * t415;
t429 = mrSges(6,1) * t343 - mrSges(6,2) * t344 + Ifges(6,5) * t355 + Ifges(6,6) * t354 + Ifges(6,3) * t412 + t379 * t359 - t378 * t360;
t428 = m(4) * t385 + t386 * mrSges(4,1) + t400 * t392 + (t393 - t394) * t401 + (mrSges(4,2) - mrSges(5,3)) * t387 + t430;
t404 = (Ifges(3,5) * t419 + Ifges(3,6) * t420) * qJD(1);
t398 = -qJDD(1) * pkin(1) - qJ(2) * t427 - t442;
t388 = -t402 * t419 + t441;
t358 = Ifges(6,5) * t379 + Ifges(6,6) * t378 + Ifges(6,3) * t415;
t339 = -t387 * mrSges(5,3) - t401 * t394 + t430;
t338 = -t449 * t427 * mrSges(3,3) + m(3) * t398 + t436 * qJDD(1) + t428;
t337 = mrSges(6,2) * t345 - mrSges(6,3) * t343 + Ifges(6,1) * t355 + Ifges(6,4) * t354 + Ifges(6,5) * t412 + t358 * t378 - t359 * t415;
t336 = -mrSges(6,1) * t345 + mrSges(6,3) * t344 + Ifges(6,4) * t355 + Ifges(6,2) * t354 + Ifges(6,6) * t412 - t358 * t379 + t360 * t415;
t334 = mrSges(5,2) * t387 + t382 * t401 - t431;
t328 = mrSges(4,2) * t385 + mrSges(5,2) * t352 - mrSges(4,3) * t356 - mrSges(5,3) * t349 - pkin(7) * t335 - qJ(4) * t339 - t453 * qJD(3) + t458 * qJDD(3) - t336 * t421 + t337 * t424 - t459 * t386 + t469 * t387 + t454 * t400;
t327 = -mrSges(4,1) * t385 - mrSges(5,1) * t349 + mrSges(5,2) * t351 + mrSges(4,3) * t357 - pkin(3) * t339 - pkin(4) * t432 - pkin(7) * t439 + t452 * qJD(3) + t457 * qJDD(3) - t424 * t336 - t421 * t337 + t468 * t386 + t459 * t387 + t454 * t401;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t450 - mrSges(2,2) * t437 + t419 * (t420 * qJD(1) * t404 + mrSges(3,2) * t398 - mrSges(3,3) * t388 + t462 * t328 - t422 * t327 - pkin(6) * t455 + (Ifges(3,1) * t419 + Ifges(3,4) * t420) * qJDD(1)) + t420 * (-t404 * t448 - mrSges(3,1) * t398 + mrSges(3,3) * t389 + t422 * t328 + t462 * t327 - pkin(2) * t428 + pkin(6) * t440 + (Ifges(3,4) * t419 + Ifges(3,2) * t420) * qJDD(1)) - pkin(1) * t338 + qJ(2) * ((m(3) * t389 + t435 * t420 + t440) * t420 + (-m(3) * t388 + t435 * t419 - t455) * t419); t338; -t429 - mrSges(4,2) * t357 + mrSges(4,1) * t356 + mrSges(5,3) * t351 - mrSges(5,1) * t352 - pkin(4) * t335 - pkin(3) * t334 + qJ(4) * t433 + t467 * qJDD(3) + (-qJ(4) * mrSges(5,2) - t457) * t386 + t458 * t387 + (-qJ(4) * t382 + t452) * t400 + t453 * t401; t334; t429;];
tauJ = t1;
