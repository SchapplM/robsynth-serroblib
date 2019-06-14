% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:31:05
% EndTime: 2019-05-05 15:31:08
% DurationCPUTime: 1.38s
% Computational Cost: add. (11531->236), mult. (21588->292), div. (0->0), fcn. (12998->10), ass. (0->100)
t416 = sin(qJ(1));
t420 = cos(qJ(1));
t434 = t416 * g(1) - t420 * g(2);
t389 = qJDD(1) * pkin(1) + t434;
t422 = qJD(1) ^ 2;
t431 = -t420 * g(1) - t416 * g(2);
t391 = -t422 * pkin(1) + t431;
t411 = sin(pkin(10));
t412 = cos(pkin(10));
t368 = t411 * t389 + t412 * t391;
t441 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t368;
t440 = -pkin(2) - pkin(7);
t352 = t440 * t422 - t441;
t419 = cos(qJ(4));
t438 = qJD(1) * qJD(4);
t401 = t419 * t438;
t415 = sin(qJ(4));
t393 = -t415 * qJDD(1) - t401;
t435 = t415 * t438;
t394 = t419 * qJDD(1) - t435;
t339 = (-t394 + t435) * pkin(8) + (-t393 + t401) * pkin(4) + t352;
t367 = t412 * t389 - t411 * t391;
t429 = -t422 * qJ(3) + qJDD(3) - t367;
t355 = t440 * qJDD(1) + t429;
t408 = -g(3) + qJDD(2);
t349 = t415 * t355 + t419 * t408;
t392 = (t415 * pkin(4) - t419 * pkin(8)) * qJD(1);
t403 = t415 * qJD(1);
t421 = qJD(4) ^ 2;
t346 = -t421 * pkin(4) + qJDD(4) * pkin(8) - t392 * t403 + t349;
t414 = sin(qJ(5));
t418 = cos(qJ(5));
t328 = t418 * t339 - t414 * t346;
t439 = qJD(1) * t419;
t387 = t418 * qJD(4) - t414 * t439;
t364 = t387 * qJD(5) + t414 * qJDD(4) + t418 * t394;
t386 = qJDD(5) - t393;
t388 = t414 * qJD(4) + t418 * t439;
t398 = t403 + qJD(5);
t326 = (t387 * t398 - t364) * pkin(9) + (t387 * t388 + t386) * pkin(5) + t328;
t329 = t414 * t339 + t418 * t346;
t363 = -t388 * qJD(5) + t418 * qJDD(4) - t414 * t394;
t372 = t398 * pkin(5) - t388 * pkin(9);
t385 = t387 ^ 2;
t327 = -t385 * pkin(5) + t363 * pkin(9) - t398 * t372 + t329;
t413 = sin(qJ(6));
t417 = cos(qJ(6));
t324 = t417 * t326 - t413 * t327;
t365 = t417 * t387 - t413 * t388;
t336 = t365 * qJD(6) + t413 * t363 + t417 * t364;
t366 = t413 * t387 + t417 * t388;
t347 = -t365 * mrSges(7,1) + t366 * mrSges(7,2);
t397 = qJD(6) + t398;
t353 = -t397 * mrSges(7,2) + t365 * mrSges(7,3);
t381 = qJDD(6) + t386;
t321 = m(7) * t324 + t381 * mrSges(7,1) - t336 * mrSges(7,3) - t366 * t347 + t397 * t353;
t325 = t413 * t326 + t417 * t327;
t335 = -t366 * qJD(6) + t417 * t363 - t413 * t364;
t354 = t397 * mrSges(7,1) - t366 * mrSges(7,3);
t322 = m(7) * t325 - t381 * mrSges(7,2) + t335 * mrSges(7,3) + t365 * t347 - t397 * t354;
t314 = t417 * t321 + t413 * t322;
t369 = -t387 * mrSges(6,1) + t388 * mrSges(6,2);
t370 = -t398 * mrSges(6,2) + t387 * mrSges(6,3);
t312 = m(6) * t328 + t386 * mrSges(6,1) - t364 * mrSges(6,3) - t388 * t369 + t398 * t370 + t314;
t371 = t398 * mrSges(6,1) - t388 * mrSges(6,3);
t432 = -t413 * t321 + t417 * t322;
t313 = m(6) * t329 - t386 * mrSges(6,2) + t363 * mrSges(6,3) + t387 * t369 - t398 * t371 + t432;
t308 = t418 * t312 + t414 * t313;
t433 = -t414 * t312 + t418 * t313;
t348 = t419 * t355 - t415 * t408;
t390 = (t415 * mrSges(5,1) + t419 * mrSges(5,2)) * qJD(1);
t396 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t439;
t307 = m(5) * t349 - qJDD(4) * mrSges(5,2) + t393 * mrSges(5,3) - qJD(4) * t396 - t390 * t403 + t433;
t395 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t403;
t345 = -qJDD(4) * pkin(4) - t421 * pkin(8) + t392 * t439 - t348;
t330 = -t363 * pkin(5) - t385 * pkin(9) + t388 * t372 + t345;
t427 = m(7) * t330 - t335 * mrSges(7,1) + t336 * mrSges(7,2) - t365 * t353 + t366 * t354;
t424 = -m(6) * t345 + t363 * mrSges(6,1) - t364 * mrSges(6,2) + t387 * t370 - t388 * t371 - t427;
t317 = m(5) * t348 + qJDD(4) * mrSges(5,1) - t394 * mrSges(5,3) + qJD(4) * t395 - t390 * t439 + t424;
t430 = t415 * t307 + t419 * t317;
t357 = -qJDD(1) * pkin(2) + t429;
t428 = m(4) * t357 - t422 * mrSges(4,3) + t430;
t342 = Ifges(7,4) * t366 + Ifges(7,2) * t365 + Ifges(7,6) * t397;
t343 = Ifges(7,1) * t366 + Ifges(7,4) * t365 + Ifges(7,5) * t397;
t426 = -mrSges(7,1) * t324 + mrSges(7,2) * t325 - Ifges(7,5) * t336 - Ifges(7,6) * t335 - Ifges(7,3) * t381 - t366 * t342 + t365 * t343;
t356 = t422 * pkin(2) + t441;
t425 = -m(4) * t356 + m(5) * t352 - t393 * mrSges(5,1) + t422 * mrSges(4,2) + t394 * mrSges(5,2) + qJDD(1) * mrSges(4,3) + t395 * t403 + t396 * t439 + t308;
t359 = Ifges(6,4) * t388 + Ifges(6,2) * t387 + Ifges(6,6) * t398;
t360 = Ifges(6,1) * t388 + Ifges(6,4) * t387 + Ifges(6,5) * t398;
t423 = mrSges(6,1) * t328 - mrSges(6,2) * t329 + Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t386 + pkin(5) * t314 + t388 * t359 - t387 * t360 - t426;
t380 = (Ifges(5,5) * qJD(4)) + (t419 * Ifges(5,1) - t415 * Ifges(5,4)) * qJD(1);
t379 = (Ifges(5,6) * qJD(4)) + (t419 * Ifges(5,4) - t415 * Ifges(5,2)) * qJD(1);
t358 = Ifges(6,5) * t388 + Ifges(6,6) * t387 + Ifges(6,3) * t398;
t341 = Ifges(7,5) * t366 + Ifges(7,6) * t365 + Ifges(7,3) * t397;
t316 = mrSges(7,2) * t330 - mrSges(7,3) * t324 + Ifges(7,1) * t336 + Ifges(7,4) * t335 + Ifges(7,5) * t381 + t365 * t341 - t397 * t342;
t315 = -mrSges(7,1) * t330 + mrSges(7,3) * t325 + Ifges(7,4) * t336 + Ifges(7,2) * t335 + Ifges(7,6) * t381 - t366 * t341 + t397 * t343;
t306 = mrSges(6,2) * t345 - mrSges(6,3) * t328 + Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t386 - pkin(9) * t314 - t413 * t315 + t417 * t316 + t387 * t358 - t398 * t359;
t305 = qJDD(1) * mrSges(4,2) + t428;
t304 = -mrSges(6,1) * t345 + mrSges(6,3) * t329 + Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t386 - pkin(5) * t427 + pkin(9) * t432 + t417 * t315 + t413 * t316 - t388 * t358 + t398 * t360;
t1 = [pkin(1) * (t411 * (m(3) * t368 - t422 * mrSges(3,1) + t425) + t412 * (m(3) * t367 - t422 * mrSges(3,2) - t428)) + mrSges(2,1) * t434 - mrSges(2,2) * t431 - pkin(2) * t305 + qJ(3) * t425 + t419 * (mrSges(5,2) * t352 - mrSges(5,3) * t348 + Ifges(5,1) * t394 + Ifges(5,4) * t393 + Ifges(5,5) * qJDD(4) - pkin(8) * t308 - qJD(4) * t379 - t414 * t304 + t418 * t306) - t415 * (-mrSges(5,1) * t352 + mrSges(5,3) * t349 + Ifges(5,4) * t394 + Ifges(5,2) * t393 + Ifges(5,6) * qJDD(4) - pkin(4) * t308 + qJD(4) * t380 - t423) - pkin(7) * t430 + mrSges(3,1) * t367 - mrSges(3,2) * t368 + mrSges(4,2) * t357 - mrSges(4,3) * t356 + (pkin(1) * (-t411 * mrSges(3,2) + t412 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1); t419 * t307 - t415 * t317 + (m(3) + m(4)) * t408; t305; Ifges(5,5) * t394 + Ifges(5,6) * t393 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t348 - mrSges(5,2) * t349 + t414 * t306 + t418 * t304 + pkin(4) * t424 + pkin(8) * t433 + (t419 * t379 + t415 * t380) * qJD(1); t423; -t426;];
tauJ  = t1;
