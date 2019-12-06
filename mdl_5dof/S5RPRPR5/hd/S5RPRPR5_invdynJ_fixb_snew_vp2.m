% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR5
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:56:23
% EndTime: 2019-12-05 17:56:28
% DurationCPUTime: 2.43s
% Computational Cost: add. (15870->240), mult. (40702->322), div. (0->0), fcn. (27700->10), ass. (0->105)
t417 = qJD(1) ^ 2;
t413 = sin(qJ(1));
t416 = cos(qJ(1));
t428 = t413 * g(2) - g(3) * t416;
t443 = -pkin(1) * t417 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t428;
t408 = sin(pkin(8));
t410 = cos(pkin(8));
t370 = -t410 * g(1) - t443 * t408;
t371 = -g(1) * t408 + t443 * t410;
t425 = -pkin(2) * t410 - pkin(6) * t408;
t395 = t425 * qJD(1);
t436 = t410 * qJD(1);
t361 = t395 * t436 + t371;
t424 = g(2) * t416 + g(3) * t413;
t420 = -qJ(2) * t417 + qJDD(2) - t424;
t372 = (-pkin(1) + t425) * qJDD(1) + t420;
t415 = cos(qJ(3));
t369 = t415 * t372;
t412 = sin(qJ(3));
t434 = qJD(1) * qJD(3);
t389 = (qJDD(1) * t415 - t412 * t434) * t408;
t433 = t410 * qJDD(1);
t399 = qJDD(3) - t433;
t400 = qJD(3) - t436;
t438 = qJD(1) * t408;
t405 = t408 ^ 2;
t440 = t405 * t417;
t342 = pkin(3) * t399 - qJ(4) * t389 + t369 + (-pkin(3) * t415 * t440 - qJ(4) * t400 * t438 - t361) * t412;
t351 = t415 * t361 + t412 * t372;
t430 = t415 * t438;
t385 = pkin(3) * t400 - qJ(4) * t430;
t388 = (-qJDD(1) * t412 - t415 * t434) * t408;
t432 = t412 ^ 2 * t440;
t343 = -pkin(3) * t432 + qJ(4) * t388 - t385 * t400 + t351;
t407 = sin(pkin(9));
t409 = cos(pkin(9));
t381 = (-t412 * t407 + t415 * t409) * t438;
t329 = -0.2e1 * qJD(4) * t381 + t409 * t342 - t343 * t407;
t364 = t388 * t407 + t389 * t409;
t380 = (-t415 * t407 - t412 * t409) * t438;
t327 = (t380 * t400 - t364) * pkin(7) + (t380 * t381 + t399) * pkin(4) + t329;
t330 = 0.2e1 * qJD(4) * t380 + t407 * t342 + t409 * t343;
t363 = t388 * t409 - t389 * t407;
t367 = pkin(4) * t400 - pkin(7) * t381;
t379 = t380 ^ 2;
t328 = -pkin(4) * t379 + pkin(7) * t363 - t367 * t400 + t330;
t411 = sin(qJ(5));
t414 = cos(qJ(5));
t325 = t327 * t414 - t328 * t411;
t358 = t380 * t414 - t381 * t411;
t338 = qJD(5) * t358 + t363 * t411 + t364 * t414;
t359 = t380 * t411 + t381 * t414;
t348 = -mrSges(6,1) * t358 + mrSges(6,2) * t359;
t398 = qJD(5) + t400;
t352 = -mrSges(6,2) * t398 + mrSges(6,3) * t358;
t397 = qJDD(5) + t399;
t321 = m(6) * t325 + mrSges(6,1) * t397 - mrSges(6,3) * t338 - t348 * t359 + t352 * t398;
t326 = t327 * t411 + t328 * t414;
t337 = -qJD(5) * t359 + t363 * t414 - t364 * t411;
t353 = mrSges(6,1) * t398 - mrSges(6,3) * t359;
t322 = m(6) * t326 - mrSges(6,2) * t397 + mrSges(6,3) * t337 + t348 * t358 - t353 * t398;
t315 = t414 * t321 + t411 * t322;
t362 = -mrSges(5,1) * t380 + mrSges(5,2) * t381;
t365 = -mrSges(5,2) * t400 + mrSges(5,3) * t380;
t313 = m(5) * t329 + mrSges(5,1) * t399 - mrSges(5,3) * t364 - t362 * t381 + t365 * t400 + t315;
t366 = mrSges(5,1) * t400 - mrSges(5,3) * t381;
t426 = -t321 * t411 + t414 * t322;
t314 = m(5) * t330 - mrSges(5,2) * t399 + mrSges(5,3) * t363 + t362 * t380 - t366 * t400 + t426;
t309 = t409 * t313 + t407 * t314;
t350 = -t361 * t412 + t369;
t355 = Ifges(5,4) * t381 + Ifges(5,2) * t380 + Ifges(5,6) * t400;
t356 = Ifges(5,1) * t381 + Ifges(5,4) * t380 + Ifges(5,5) * t400;
t345 = Ifges(6,4) * t359 + Ifges(6,2) * t358 + Ifges(6,6) * t398;
t346 = Ifges(6,1) * t359 + Ifges(6,4) * t358 + Ifges(6,5) * t398;
t419 = -mrSges(6,1) * t325 + mrSges(6,2) * t326 - Ifges(6,5) * t338 - Ifges(6,6) * t337 - Ifges(6,3) * t397 - t359 * t345 + t358 * t346;
t442 = -mrSges(4,1) * t350 - mrSges(5,1) * t329 + mrSges(4,2) * t351 + mrSges(5,2) * t330 - Ifges(4,5) * t389 - Ifges(5,5) * t364 - Ifges(4,6) * t388 - Ifges(5,6) * t363 - pkin(3) * t309 - pkin(4) * t315 - t381 * t355 + t380 * t356 - (Ifges(4,3) + Ifges(5,3)) * t399 + t419;
t437 = qJDD(1) * mrSges(3,3);
t431 = t412 * t438;
t427 = -t313 * t407 + t409 * t314;
t360 = t395 * t438 - t370;
t423 = -mrSges(3,1) * t410 + mrSges(3,2) * t408;
t384 = -mrSges(4,2) * t400 - mrSges(4,3) * t431;
t387 = (t412 * mrSges(4,1) + t415 * mrSges(4,2)) * t438;
t307 = m(4) * t350 + mrSges(4,1) * t399 - mrSges(4,3) * t389 + t384 * t400 - t387 * t430 + t309;
t386 = mrSges(4,1) * t400 - mrSges(4,3) * t430;
t308 = m(4) * t351 - mrSges(4,2) * t399 + mrSges(4,3) * t388 - t386 * t400 - t387 * t431 + t427;
t304 = t307 * t415 + t308 * t412;
t374 = Ifges(4,6) * t400 + (t415 * Ifges(4,4) - t412 * Ifges(4,2)) * t438;
t375 = Ifges(4,5) * t400 + (t415 * Ifges(4,1) - t412 * Ifges(4,4)) * t438;
t422 = t374 * t415 + t375 * t412;
t349 = -pkin(3) * t388 - qJ(4) * t432 + t385 * t430 + qJDD(4) + t360;
t332 = -pkin(4) * t363 - pkin(7) * t379 + t367 * t381 + t349;
t421 = m(6) * t332 - t337 * mrSges(6,1) + t338 * mrSges(6,2) - t358 * t352 + t359 * t353;
t323 = m(5) * t349 - t363 * mrSges(5,1) + t364 * mrSges(5,2) - t380 * t365 + t381 * t366 + t421;
t394 = (Ifges(3,5) * t408 + Ifges(3,6) * t410) * qJD(1);
t393 = t423 * qJD(1);
t391 = -qJDD(1) * pkin(1) + t420;
t354 = Ifges(5,5) * t381 + Ifges(5,6) * t380 + Ifges(5,3) * t400;
t344 = Ifges(6,5) * t359 + Ifges(6,6) * t358 + Ifges(6,3) * t398;
t317 = mrSges(6,2) * t332 - mrSges(6,3) * t325 + Ifges(6,1) * t338 + Ifges(6,4) * t337 + Ifges(6,5) * t397 + t344 * t358 - t345 * t398;
t316 = -mrSges(6,1) * t332 + mrSges(6,3) * t326 + Ifges(6,4) * t338 + Ifges(6,2) * t337 + Ifges(6,6) * t397 - t344 * t359 + t346 * t398;
t306 = mrSges(5,2) * t349 - mrSges(5,3) * t329 + Ifges(5,1) * t364 + Ifges(5,4) * t363 + Ifges(5,5) * t399 - pkin(7) * t315 - t316 * t411 + t317 * t414 + t354 * t380 - t355 * t400;
t305 = -mrSges(5,1) * t349 + mrSges(5,3) * t330 + Ifges(5,4) * t364 + Ifges(5,2) * t363 + Ifges(5,6) * t399 - pkin(4) * t421 + pkin(7) * t426 + t414 * t316 + t411 * t317 - t381 * t354 + t400 * t356;
t303 = m(3) * t391 + t423 * qJDD(1) + (-t410 ^ 2 - t405) * t417 * mrSges(3,3) + t304;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t424 - mrSges(2,2) * t428 + t408 * (mrSges(3,2) * t391 - mrSges(3,3) * t370 + t415 * (mrSges(4,2) * t360 - mrSges(4,3) * t350 + Ifges(4,1) * t389 + Ifges(4,4) * t388 + Ifges(4,5) * t399 - qJ(4) * t309 - t305 * t407 + t306 * t409 - t374 * t400) - t412 * (-mrSges(4,1) * t360 + mrSges(4,3) * t351 + Ifges(4,4) * t389 + Ifges(4,2) * t388 + Ifges(4,6) * t399 - pkin(3) * t323 + qJ(4) * t427 + t409 * t305 + t407 * t306 + t400 * t375) - pkin(6) * t304 + (Ifges(3,1) * t408 + t410 * Ifges(3,4)) * qJDD(1) + t394 * t436) + t410 * ((Ifges(3,4) * qJDD(1) + (-t394 - t422) * qJD(1)) * t408 + Ifges(3,2) * t433 - mrSges(3,1) * t391 + mrSges(3,3) * t371 - pkin(2) * t304 + t442) - pkin(1) * t303 + qJ(2) * ((m(3) * t371 - t412 * t307 + t415 * t308 + (qJD(1) * t393 + t437) * t410) * t410 + (t323 + (t384 * t412 + t386 * t415 + t393) * t438 + t408 * t437 - t388 * mrSges(4,1) + t389 * mrSges(4,2) - m(3) * t370 + m(4) * t360) * t408); t303; t422 * t438 - t442; t323; -t419;];
tauJ = t1;
