% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-05-04 20:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:04:06
% EndTime: 2019-05-04 20:04:09
% DurationCPUTime: 2.18s
% Computational Cost: add. (20760->227), mult. (39339->300), div. (0->0), fcn. (30702->16), ass. (0->105)
t405 = sin(pkin(11));
t410 = cos(pkin(11));
t397 = -g(1) * t410 - g(2) * t405;
t404 = sin(pkin(12));
t409 = cos(pkin(12));
t396 = g(1) * t405 - g(2) * t410;
t402 = -g(3) + qJDD(1);
t407 = sin(pkin(6));
t412 = cos(pkin(6));
t424 = t396 * t412 + t402 * t407;
t363 = -t397 * t404 + t424 * t409;
t374 = -t396 * t407 + t402 * t412 + qJDD(2);
t406 = sin(pkin(7));
t411 = cos(pkin(7));
t436 = t363 * t411 + t374 * t406;
t364 = t397 * t409 + t424 * t404;
t415 = sin(qJ(3));
t418 = cos(qJ(3));
t346 = -t415 * t364 + t436 * t418;
t347 = t418 * t364 + t436 * t415;
t420 = qJD(3) ^ 2;
t345 = -pkin(3) * t420 + qJDD(3) * pkin(9) + t347;
t357 = -t363 * t406 + t374 * t411;
t414 = sin(qJ(4));
t417 = cos(qJ(4));
t338 = t417 * t345 + t414 * t357;
t392 = (-t417 * pkin(4) - t414 * qJ(5)) * qJD(3);
t419 = qJD(4) ^ 2;
t433 = qJD(3) * t417;
t336 = -pkin(4) * t419 + qJDD(4) * qJ(5) + t392 * t433 + t338;
t344 = -qJDD(3) * pkin(3) - pkin(9) * t420 - t346;
t431 = qJD(3) * qJD(4);
t430 = t417 * t431;
t394 = qJDD(3) * t414 + t430;
t401 = t414 * t431;
t395 = qJDD(3) * t417 - t401;
t341 = (-t394 - t430) * qJ(5) + (-t395 + t401) * pkin(4) + t344;
t403 = sin(pkin(13));
t408 = cos(pkin(13));
t432 = t414 * qJD(3);
t389 = qJD(4) * t403 + t408 * t432;
t331 = -0.2e1 * qJD(5) * t389 - t336 * t403 + t408 * t341;
t378 = qJDD(4) * t403 + t394 * t408;
t388 = qJD(4) * t408 - t403 * t432;
t329 = (-t388 * t433 - t378) * pkin(10) + (t388 * t389 - t395) * pkin(5) + t331;
t332 = 0.2e1 * qJD(5) * t388 + t408 * t336 + t403 * t341;
t377 = qJDD(4) * t408 - t394 * t403;
t379 = -pkin(5) * t433 - pkin(10) * t389;
t387 = t388 ^ 2;
t330 = -pkin(5) * t387 + pkin(10) * t377 + t379 * t433 + t332;
t413 = sin(qJ(6));
t416 = cos(qJ(6));
t326 = t329 * t416 - t330 * t413;
t369 = t388 * t416 - t389 * t413;
t350 = qJD(6) * t369 + t377 * t413 + t378 * t416;
t370 = t388 * t413 + t389 * t416;
t356 = -mrSges(7,1) * t369 + mrSges(7,2) * t370;
t400 = qJD(6) - t433;
t359 = -mrSges(7,2) * t400 + mrSges(7,3) * t369;
t391 = qJDD(6) - t395;
t323 = m(7) * t326 + mrSges(7,1) * t391 - mrSges(7,3) * t350 - t356 * t370 + t359 * t400;
t327 = t329 * t413 + t330 * t416;
t349 = -qJD(6) * t370 + t377 * t416 - t378 * t413;
t360 = mrSges(7,1) * t400 - mrSges(7,3) * t370;
t324 = m(7) * t327 - mrSges(7,2) * t391 + mrSges(7,3) * t349 + t356 * t369 - t360 * t400;
t317 = t416 * t323 + t413 * t324;
t393 = (-t417 * mrSges(5,1) + t414 * mrSges(5,2)) * qJD(3);
t398 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t432;
t371 = -mrSges(6,1) * t388 + mrSges(6,2) * t389;
t375 = mrSges(6,2) * t433 + mrSges(6,3) * t388;
t315 = m(6) * t331 - mrSges(6,1) * t395 - mrSges(6,3) * t378 - t371 * t389 - t375 * t433 + t317;
t376 = -mrSges(6,1) * t433 - mrSges(6,3) * t389;
t427 = -t323 * t413 + t416 * t324;
t316 = m(6) * t332 + mrSges(6,2) * t395 + mrSges(6,3) * t377 + t371 * t388 + t376 * t433 + t427;
t428 = -t315 * t403 + t408 * t316;
t312 = m(5) * t338 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t395 - qJD(4) * t398 + t393 * t433 + t428;
t337 = -t414 * t345 + t357 * t417;
t335 = -qJDD(4) * pkin(4) - qJ(5) * t419 + t392 * t432 + qJDD(5) - t337;
t333 = -pkin(5) * t377 - pkin(10) * t387 + t379 * t389 + t335;
t423 = m(7) * t333 - t349 * mrSges(7,1) + mrSges(7,2) * t350 - t369 * t359 + t360 * t370;
t328 = m(6) * t335 - t377 * mrSges(6,1) + mrSges(6,2) * t378 - t388 * t375 + t376 * t389 + t423;
t399 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t433;
t325 = m(5) * t337 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t394 + qJD(4) * t399 - t393 * t432 - t328;
t429 = t417 * t312 - t414 * t325;
t306 = m(4) * t347 - mrSges(4,1) * t420 - qJDD(3) * mrSges(4,2) + t429;
t313 = t408 * t315 + t403 * t316;
t422 = -m(5) * t344 + t395 * mrSges(5,1) - t394 * mrSges(5,2) - t398 * t432 + t399 * t433 - t313;
t310 = m(4) * t346 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t420 + t422;
t426 = t306 * t415 + t310 * t418;
t352 = Ifges(7,4) * t370 + Ifges(7,2) * t369 + Ifges(7,6) * t400;
t353 = Ifges(7,1) * t370 + Ifges(7,4) * t369 + Ifges(7,5) * t400;
t421 = mrSges(7,1) * t326 - mrSges(7,2) * t327 + Ifges(7,5) * t350 + Ifges(7,6) * t349 + Ifges(7,3) * t391 + t370 * t352 - t369 * t353;
t385 = Ifges(5,5) * qJD(4) + (t414 * Ifges(5,1) + t417 * Ifges(5,4)) * qJD(3);
t384 = Ifges(5,6) * qJD(4) + (t414 * Ifges(5,4) + Ifges(5,2) * t417) * qJD(3);
t367 = Ifges(6,1) * t389 + Ifges(6,4) * t388 - Ifges(6,5) * t433;
t366 = Ifges(6,4) * t389 + Ifges(6,2) * t388 - Ifges(6,6) * t433;
t365 = Ifges(6,5) * t389 + Ifges(6,6) * t388 - Ifges(6,3) * t433;
t351 = Ifges(7,5) * t370 + Ifges(7,6) * t369 + Ifges(7,3) * t400;
t319 = mrSges(7,2) * t333 - mrSges(7,3) * t326 + Ifges(7,1) * t350 + Ifges(7,4) * t349 + Ifges(7,5) * t391 + t351 * t369 - t352 * t400;
t318 = -mrSges(7,1) * t333 + mrSges(7,3) * t327 + Ifges(7,4) * t350 + Ifges(7,2) * t349 + Ifges(7,6) * t391 - t351 * t370 + t353 * t400;
t309 = mrSges(6,2) * t335 - mrSges(6,3) * t331 + Ifges(6,1) * t378 + Ifges(6,4) * t377 - Ifges(6,5) * t395 - pkin(10) * t317 - t318 * t413 + t319 * t416 + t365 * t388 + t366 * t433;
t308 = -mrSges(6,1) * t335 + mrSges(6,3) * t332 + Ifges(6,4) * t378 + Ifges(6,2) * t377 - Ifges(6,6) * t395 - pkin(5) * t423 + pkin(10) * t427 + t416 * t318 + t413 * t319 - t389 * t365 - t367 * t433;
t307 = m(4) * t357 + t312 * t414 + t325 * t417;
t305 = m(3) * t374 + t307 * t411 + t426 * t406;
t1 = [m(2) * t402 + t412 * t305 + (t404 * (m(3) * t364 + t306 * t418 - t310 * t415) + t409 * (m(3) * t363 - t406 * t307 + t426 * t411)) * t407; t305; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t346 - mrSges(4,2) * t347 + t414 * (mrSges(5,2) * t344 - mrSges(5,3) * t337 + Ifges(5,1) * t394 + Ifges(5,4) * t395 + Ifges(5,5) * qJDD(4) - qJ(5) * t313 - qJD(4) * t384 - t403 * t308 + t408 * t309) + t417 * (-mrSges(5,1) * t344 - mrSges(6,1) * t331 + mrSges(6,2) * t332 + mrSges(5,3) * t338 + Ifges(5,4) * t394 - Ifges(6,5) * t378 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t377 - pkin(4) * t313 - pkin(5) * t317 + qJD(4) * t385 - t389 * t366 + t388 * t367 - t421 + (Ifges(5,2) + Ifges(6,3)) * t395) + pkin(3) * t422 + pkin(9) * t429; Ifges(5,5) * t394 + Ifges(5,6) * t395 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t337 - mrSges(5,2) * t338 + t403 * t309 + t408 * t308 - pkin(4) * t328 + qJ(5) * t428 + (t414 * t384 - t417 * t385) * qJD(3); t328; t421;];
tauJ  = t1;
