% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR13
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:53:34
% EndTime: 2019-05-05 20:54:45
% DurationCPUTime: 58.60s
% Computational Cost: add. (914239->405), mult. (2896973->532), div. (0->0), fcn. (2435665->14), ass. (0->179)
t371 = cos(pkin(12));
t369 = sin(pkin(7));
t432 = cos(pkin(6));
t410 = t432 * t369;
t370 = sin(pkin(6));
t431 = cos(pkin(7));
t413 = t370 * t431;
t347 = (t371 * t413 + t410) * qJD(1) * pkin(9);
t375 = sin(qJ(1));
t378 = cos(qJ(1));
t364 = -t378 * g(1) - t375 * g(2);
t379 = qJD(1) ^ 2;
t430 = qJ(2) * t370;
t352 = -t379 * pkin(1) + qJDD(1) * t430 + t364;
t368 = sin(pkin(12));
t435 = pkin(9) * t369;
t399 = -pkin(2) * t371 - t368 * t435;
t416 = t431 * pkin(9);
t422 = qJD(1) * t370;
t393 = qJD(1) * t399 * t422 + qJDD(1) * t416;
t363 = t375 * g(1) - t378 * g(2);
t351 = qJDD(1) * pkin(1) + t379 * t430 + t363;
t412 = t371 * t432;
t415 = qJD(2) * t422;
t426 = t370 * t371;
t401 = -g(3) * t426 + t351 * t412 - 0.2e1 * t368 * t415;
t407 = qJDD(1) * t432;
t409 = qJD(1) * t432;
t281 = pkin(2) * t407 + t347 * t409 + (-t370 * t393 - t352) * t368 + t401;
t353 = (-pkin(9) * t368 * t413 + pkin(2) * t432) * qJD(1);
t414 = t368 * t432;
t418 = t351 * t414 + (t352 + 0.2e1 * t415) * t371;
t282 = t407 * t435 - t353 * t409 + (-g(3) * t368 + t371 * t393) * t370 + t418;
t402 = -g(3) * t432 + qJDD(2);
t293 = (-t351 + t399 * qJDD(1) + (-t347 * t371 + t353 * t368) * qJD(1)) * t370 + t402;
t374 = sin(qJ(3));
t411 = t374 * t431;
t427 = t369 * t374;
t436 = cos(qJ(3));
t255 = t281 * t411 + t282 * t436 + t293 * t427;
t406 = t431 * t436;
t428 = t368 * t370;
t437 = t374 * t428 - t406 * t426 - t436 * t410;
t331 = t437 * qJD(1);
t382 = t374 * t410 + (t368 * t436 + t371 * t411) * t370;
t332 = t382 * qJD(1);
t311 = t331 * pkin(3) - t332 * qJ(4);
t438 = -t369 * t426 + t432 * t431;
t348 = -t438 * qJD(1) - qJD(3);
t344 = t348 ^ 2;
t345 = t438 * qJDD(1) + qJDD(3);
t249 = t344 * pkin(3) - t345 * qJ(4) + 0.2e1 * qJD(4) * t348 + t331 * t311 - t255;
t434 = mrSges(4,1) - mrSges(5,2);
t433 = -Ifges(5,6) - Ifges(4,4);
t429 = t331 * t348;
t417 = t369 * t436;
t254 = t281 * t406 - t374 * t282 + t293 * t417;
t252 = -t345 * pkin(3) - t344 * qJ(4) + t332 * t311 + qJDD(4) - t254;
t315 = -t331 * qJD(3) + qJDD(1) * t382;
t244 = (t331 * t332 - t345) * pkin(10) + (t315 - t429) * pkin(4) + t252;
t314 = t332 * qJD(3) + t437 * qJDD(1);
t325 = t332 * pkin(4) + t348 * pkin(10);
t330 = t331 ^ 2;
t263 = -t369 * t281 + t431 * t293;
t387 = (-t315 - t429) * qJ(4) + t263 + (-t348 * pkin(3) - 0.2e1 * qJD(4)) * t332;
t248 = -t330 * pkin(4) - t332 * t325 + (pkin(3) + pkin(10)) * t314 + t387;
t373 = sin(qJ(5));
t377 = cos(qJ(5));
t241 = t373 * t244 + t377 * t248;
t319 = t377 * t331 + t373 * t348;
t320 = t373 * t331 - t377 * t348;
t285 = -t319 * pkin(5) - t320 * pkin(11);
t310 = qJDD(5) + t315;
t329 = qJD(5) + t332;
t328 = t329 ^ 2;
t238 = -t328 * pkin(5) + t310 * pkin(11) + t319 * t285 + t241;
t246 = -t314 * pkin(4) - t330 * pkin(10) - t348 * t325 - t249;
t276 = -t320 * qJD(5) + t377 * t314 - t373 * t345;
t277 = t319 * qJD(5) + t373 * t314 + t377 * t345;
t242 = (-t319 * t329 - t277) * pkin(11) + (t320 * t329 - t276) * pkin(5) + t246;
t372 = sin(qJ(6));
t376 = cos(qJ(6));
t233 = -t372 * t238 + t376 * t242;
t291 = -t372 * t320 + t376 * t329;
t258 = t291 * qJD(6) + t376 * t277 + t372 * t310;
t292 = t376 * t320 + t372 * t329;
t265 = -t291 * mrSges(7,1) + t292 * mrSges(7,2);
t316 = qJD(6) - t319;
t266 = -t316 * mrSges(7,2) + t291 * mrSges(7,3);
t274 = qJDD(6) - t276;
t231 = m(7) * t233 + t274 * mrSges(7,1) - t258 * mrSges(7,3) - t292 * t265 + t316 * t266;
t234 = t376 * t238 + t372 * t242;
t257 = -t292 * qJD(6) - t372 * t277 + t376 * t310;
t267 = t316 * mrSges(7,1) - t292 * mrSges(7,3);
t232 = m(7) * t234 - t274 * mrSges(7,2) + t257 * mrSges(7,3) + t291 * t265 - t316 * t267;
t221 = t376 * t231 + t372 * t232;
t300 = -Ifges(5,1) * t348 - Ifges(5,4) * t332 + Ifges(5,5) * t331;
t425 = -Ifges(4,5) * t332 + Ifges(4,6) * t331 + Ifges(4,3) * t348 - t300;
t298 = -Ifges(5,4) * t348 - Ifges(5,2) * t332 + Ifges(5,6) * t331;
t424 = -Ifges(4,1) * t332 + Ifges(4,4) * t331 + Ifges(4,5) * t348 + t298;
t322 = t331 * mrSges(5,1) + t348 * mrSges(5,3);
t423 = t348 * mrSges(4,2) - t331 * mrSges(4,3) - t322;
t312 = t331 * mrSges(4,1) + t332 * mrSges(4,2);
t284 = -t319 * mrSges(6,1) + t320 * mrSges(6,2);
t295 = t329 * mrSges(6,1) - t320 * mrSges(6,3);
t408 = -t372 * t231 + t376 * t232;
t218 = m(6) * t241 - t310 * mrSges(6,2) + t276 * mrSges(6,3) + t319 * t284 - t329 * t295 + t408;
t240 = t377 * t244 - t373 * t248;
t294 = -t329 * mrSges(6,2) + t319 * mrSges(6,3);
t237 = -t310 * pkin(5) - t328 * pkin(11) + t320 * t285 - t240;
t392 = -m(7) * t237 + t257 * mrSges(7,1) - t258 * mrSges(7,2) + t291 * t266 - t292 * t267;
t227 = m(6) * t240 + t310 * mrSges(6,1) - t277 * mrSges(6,3) - t320 * t284 + t329 * t294 + t392;
t211 = t373 * t218 + t377 * t227;
t313 = -t331 * mrSges(5,2) - t332 * mrSges(5,3);
t391 = -m(5) * t252 - t315 * mrSges(5,1) - t332 * t313 - t211;
t207 = m(4) * t254 - t315 * mrSges(4,3) - t332 * t312 + t345 * t434 - t348 * t423 + t391;
t324 = -t348 * mrSges(4,1) - t332 * mrSges(4,3);
t212 = t377 * t218 - t373 * t227;
t253 = t314 * pkin(3) + t387;
t323 = t332 * mrSges(5,1) - t348 * mrSges(5,2);
t397 = m(5) * t253 - t315 * mrSges(5,3) - t332 * t323 + t212;
t209 = m(4) * t263 + t315 * mrSges(4,2) + t314 * t434 + t332 * t324 + t331 * t423 + t397;
t219 = -m(6) * t246 + t276 * mrSges(6,1) - t277 * mrSges(6,2) + t319 * t294 - t320 * t295 - t221;
t384 = -m(5) * t249 + t345 * mrSges(5,3) - t348 * t323 - t219;
t216 = t384 + t348 * t324 - t345 * mrSges(4,2) + m(4) * t255 + (-t312 - t313) * t331 + (-mrSges(4,3) - mrSges(5,1)) * t314;
t196 = t207 * t417 + t431 * t209 + t216 * t427;
t197 = t207 * t406 - t369 * t209 + t216 * t411;
t317 = -t368 * t352 + t401;
t404 = -mrSges(3,1) * t371 + mrSges(3,2) * t368;
t350 = t404 * t422;
t395 = -mrSges(3,2) * t432 + mrSges(3,3) * t426;
t355 = t395 * qJD(1);
t396 = mrSges(3,1) * t432 - mrSges(3,3) * t428;
t194 = m(3) * t317 + t396 * qJDD(1) + (-t350 * t428 + t355 * t432) * qJD(1) + t197;
t201 = -t374 * t207 + t216 * t436;
t318 = -g(3) * t428 + t418;
t354 = t396 * qJD(1);
t200 = m(3) * t318 + t395 * qJDD(1) + (t350 * t426 - t354 * t432) * qJD(1) + t201;
t191 = -t368 * t194 + t371 * t200;
t333 = -t370 * t351 + t402;
t195 = m(3) * t333 + (t404 * qJDD(1) + (t354 * t368 - t355 * t371) * qJD(1)) * t370 + t196;
t186 = t194 * t412 - t370 * t195 + t200 * t414;
t403 = Ifges(3,5) * t368 + Ifges(3,6) * t371;
t299 = Ifges(4,4) * t332 - Ifges(4,2) * t331 - Ifges(4,6) * t348;
t259 = Ifges(7,5) * t292 + Ifges(7,6) * t291 + Ifges(7,3) * t316;
t261 = Ifges(7,1) * t292 + Ifges(7,4) * t291 + Ifges(7,5) * t316;
t225 = -mrSges(7,1) * t237 + mrSges(7,3) * t234 + Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t274 - t292 * t259 + t316 * t261;
t260 = Ifges(7,4) * t292 + Ifges(7,2) * t291 + Ifges(7,6) * t316;
t226 = mrSges(7,2) * t237 - mrSges(7,3) * t233 + Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t274 + t291 * t259 - t316 * t260;
t268 = Ifges(6,5) * t320 + Ifges(6,6) * t319 + Ifges(6,3) * t329;
t269 = Ifges(6,4) * t320 + Ifges(6,2) * t319 + Ifges(6,6) * t329;
t203 = mrSges(6,2) * t246 - mrSges(6,3) * t240 + Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t310 - pkin(11) * t221 - t372 * t225 + t376 * t226 + t319 * t268 - t329 * t269;
t270 = Ifges(6,1) * t320 + Ifges(6,4) * t319 + Ifges(6,5) * t329;
t381 = mrSges(7,1) * t233 - mrSges(7,2) * t234 + Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t274 + t292 * t260 - t291 * t261;
t204 = -mrSges(6,1) * t246 + mrSges(6,3) * t241 + Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t310 - pkin(5) * t221 - t320 * t268 + t329 * t270 - t381;
t296 = -Ifges(5,5) * t348 - Ifges(5,6) * t332 + Ifges(5,3) * t331;
t386 = mrSges(5,2) * t252 - mrSges(5,3) * t249 + Ifges(5,1) * t345 - Ifges(5,4) * t315 + Ifges(5,5) * t314 - pkin(10) * t211 + t377 * t203 - t373 * t204 - t332 * t296;
t187 = Ifges(4,3) * t345 + t332 * t299 - Ifges(4,6) * t314 + Ifges(4,5) * t315 - mrSges(4,2) * t255 + mrSges(4,1) * t254 + qJ(4) * (-t314 * mrSges(5,1) + t384) + t386 + pkin(3) * (-t345 * mrSges(5,2) + t348 * t322 + t391) + (-qJ(4) * t313 - t424) * t331;
t210 = -t314 * mrSges(5,2) - t331 * t322 + t397;
t383 = -mrSges(5,1) * t249 + mrSges(5,2) * t253 - pkin(4) * t219 - pkin(10) * t212 - t373 * t203 - t377 * t204;
t188 = -mrSges(4,1) * t263 + mrSges(4,3) * t255 - pkin(3) * t210 + t424 * t348 + (Ifges(4,6) - Ifges(5,5)) * t345 + t425 * t332 - t433 * t315 + (-Ifges(4,2) - Ifges(5,3)) * t314 + t383;
t385 = mrSges(6,1) * t240 - mrSges(6,2) * t241 + Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t310 + pkin(5) * t392 + pkin(11) * t408 + t376 * t225 + t372 * t226 + t320 * t269 - t319 * t270;
t380 = mrSges(5,1) * t252 - mrSges(5,3) * t253 + pkin(4) * t211 + t385;
t192 = mrSges(4,2) * t263 - mrSges(4,3) * t254 + t380 + (-t296 + t299) * t348 + (Ifges(4,5) - Ifges(5,4)) * t345 + t425 * t331 + (Ifges(4,1) + Ifges(5,2)) * t315 + t433 * t314 - qJ(4) * t210;
t388 = t432 * Ifges(3,6) + (Ifges(3,4) * t368 + Ifges(3,2) * t371) * t370;
t337 = t388 * qJD(1);
t389 = t432 * Ifges(3,5) + (Ifges(3,1) * t368 + Ifges(3,4) * t371) * t370;
t338 = t389 * qJD(1);
t178 = Ifges(3,3) * t407 + mrSges(3,1) * t317 - mrSges(3,2) * t318 + t431 * t187 + pkin(2) * t197 + (pkin(9) * t201 + t188 * t436 + t192 * t374) * t369 + (t403 * qJDD(1) + (t337 * t368 - t338 * t371) * qJD(1)) * t370;
t336 = (Ifges(3,3) * t432 + t370 * t403) * qJD(1);
t180 = t188 * t406 + t201 * t416 + t192 * t411 - mrSges(3,1) * t333 + mrSges(3,3) * t318 - pkin(2) * t196 - t369 * t187 + (-t336 * t428 + t338 * t432) * qJD(1) + t388 * qJDD(1);
t182 = mrSges(3,2) * t333 - mrSges(3,3) * t317 + t436 * t192 - t374 * t188 + (t336 * t426 - t337 * t432) * qJD(1) + (-t196 * t369 - t197 * t431) * pkin(9) + t389 * qJDD(1);
t390 = mrSges(2,1) * t363 - mrSges(2,2) * t364 + Ifges(2,3) * qJDD(1) + pkin(1) * t186 + t432 * t178 + t180 * t426 + t182 * t428 + t191 * t430;
t189 = m(2) * t364 - t379 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t191;
t185 = t432 * t195 + (t194 * t371 + t200 * t368) * t370;
t183 = m(2) * t363 + qJDD(1) * mrSges(2,1) - t379 * mrSges(2,2) + t186;
t176 = -mrSges(2,2) * g(3) - mrSges(2,3) * t363 + Ifges(2,5) * qJDD(1) - t379 * Ifges(2,6) - t368 * t180 + t371 * t182 + (-t185 * t370 - t186 * t432) * qJ(2);
t175 = qJ(2) * t191 * t432 + mrSges(2,1) * g(3) + mrSges(2,3) * t364 + t379 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t185 - t370 * t178 + t180 * t412 + t182 * t414;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t378 * t176 - t375 * t175 - pkin(8) * (t378 * t183 + t375 * t189), t176, t182, t192, -t331 * t298 + t386, t203, t226; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t375 * t176 + t378 * t175 + pkin(8) * (-t375 * t183 + t378 * t189), t175, t180, t188, Ifges(5,4) * t345 - Ifges(5,2) * t315 + Ifges(5,6) * t314 + t348 * t296 + t331 * t300 - t380, t204, t225; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t390, t390, t178, t187, Ifges(5,5) * t345 - Ifges(5,6) * t315 + Ifges(5,3) * t314 - t348 * t298 + t332 * t300 - t383, t385, t381;];
m_new  = t1;
