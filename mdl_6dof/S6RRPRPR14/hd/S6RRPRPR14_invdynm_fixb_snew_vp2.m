% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-05-06 17:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:54:45
% EndTime: 2019-05-06 16:55:13
% DurationCPUTime: 11.43s
% Computational Cost: add. (174898->407), mult. (392958->486), div. (0->0), fcn. (277961->10), ass. (0->159)
t353 = sin(qJ(1));
t357 = cos(qJ(1));
t337 = t353 * g(1) - t357 * g(2);
t358 = qJD(1) ^ 2;
t348 = sin(pkin(6));
t405 = pkin(8) * t348;
t319 = qJDD(1) * pkin(1) + t358 * t405 + t337;
t349 = cos(pkin(6));
t342 = t349 * qJD(1) + qJD(2);
t352 = sin(qJ(2));
t387 = qJD(1) * t348;
t380 = t352 * t387;
t324 = pkin(3) * t380 - t342 * pkin(9);
t356 = cos(qJ(2));
t386 = qJD(1) * t356;
t325 = (qJD(2) * t386 + qJDD(1) * t352) * t348;
t384 = qJDD(1) * t348;
t326 = -qJD(2) * t380 + t356 * t384;
t398 = t348 ^ 2 * t358;
t383 = t356 ^ 2 * t398;
t404 = t349 * g(3);
t406 = -pkin(2) - pkin(9);
t411 = -2 * qJD(3);
t409 = (pkin(2) * t342 + t411) * t380;
t226 = -pkin(3) * t383 - t404 - t325 * qJ(3) + t406 * t326 + (-t319 + (-qJ(3) * t342 * t356 - t324 * t352) * qJD(1)) * t348 + t409;
t341 = t349 * qJDD(1) + qJDD(2);
t321 = (-pkin(2) * t356 - qJ(3) * t352) * t387;
t340 = t342 ^ 2;
t338 = -t357 * g(1) - t353 * g(2);
t320 = -t358 * pkin(1) + pkin(8) * t384 + t338;
t396 = t348 * t356;
t388 = g(3) * t396 + t352 * t320;
t376 = -t340 * qJ(3) + t321 * t380 + qJDD(3) + t388;
t228 = t325 * pkin(3) + t406 * t341 + (-pkin(3) * t342 * t387 - pkin(9) * t352 * t398 - t319 * t349) * t356 + t376;
t351 = sin(qJ(4));
t355 = cos(qJ(4));
t220 = -t351 * t226 + t355 * t228;
t221 = t355 * t226 + t351 * t228;
t381 = t348 * t386;
t302 = t351 * t342 + t355 * t381;
t303 = t355 * t342 - t351 * t381;
t332 = qJD(4) + t380;
t255 = Ifges(5,4) * t303 - Ifges(5,2) * t302 + Ifges(5,6) * t332;
t268 = t303 * qJD(4) + t355 * t326 + t351 * t341;
t269 = -t302 * qJD(4) - t351 * t326 + t355 * t341;
t274 = -t302 * mrSges(6,2) - t303 * mrSges(6,3);
t278 = t302 * mrSges(6,1) - t332 * mrSges(6,3);
t314 = qJDD(4) + t325;
t272 = t302 * pkin(4) - t303 * qJ(5);
t329 = t332 ^ 2;
t216 = -t314 * pkin(4) - t329 * qJ(5) + t303 * t272 + qJDD(5) - t220;
t399 = t302 * t332;
t210 = (t302 * t303 - t314) * pkin(10) + (t269 + t399) * pkin(5) + t216;
t282 = t303 * pkin(5) - t332 * pkin(10);
t301 = t302 ^ 2;
t395 = t349 * t352;
t397 = t348 * t352;
t271 = -g(3) * t397 + t319 * t395 + t356 * t320;
t236 = t340 * pkin(2) - t341 * qJ(3) - t321 * t381 + t342 * t411 - t271;
t225 = t326 * pkin(3) - pkin(9) * t383 + t342 * t324 - t236;
t407 = -2 * qJD(5);
t360 = (-t269 + t399) * qJ(5) + t225 + (pkin(4) * t332 + t407) * t303;
t213 = -t303 * t282 - t301 * pkin(5) + (pkin(4) + pkin(10)) * t268 + t360;
t350 = sin(qJ(6));
t354 = cos(qJ(6));
t207 = t354 * t210 - t350 * t213;
t276 = t354 * t302 - t350 * t332;
t234 = t276 * qJD(6) + t350 * t268 + t354 * t314;
t277 = t350 * t302 + t354 * t332;
t244 = -t276 * mrSges(7,1) + t277 * mrSges(7,2);
t300 = qJD(6) + t303;
t249 = -t300 * mrSges(7,2) + t276 * mrSges(7,3);
t265 = qJDD(6) + t269;
t204 = m(7) * t207 + t265 * mrSges(7,1) - t234 * mrSges(7,3) - t277 * t244 + t300 * t249;
t208 = t350 * t210 + t354 * t213;
t233 = -t277 * qJD(6) + t354 * t268 - t350 * t314;
t250 = t300 * mrSges(7,1) - t277 * mrSges(7,3);
t205 = m(7) * t208 - t265 * mrSges(7,2) + t233 * mrSges(7,3) + t276 * t244 - t300 * t250;
t194 = t354 * t204 + t350 * t205;
t373 = -t329 * pkin(4) + t314 * qJ(5) - t302 * t272 + t221;
t212 = -t268 * pkin(5) - t301 * pkin(10) + ((2 * qJD(5)) + t282) * t332 + t373;
t238 = Ifges(7,5) * t277 + Ifges(7,6) * t276 + Ifges(7,3) * t300;
t240 = Ifges(7,1) * t277 + Ifges(7,4) * t276 + Ifges(7,5) * t300;
t197 = -mrSges(7,1) * t212 + mrSges(7,3) * t208 + Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t265 - t277 * t238 + t300 * t240;
t239 = Ifges(7,4) * t277 + Ifges(7,2) * t276 + Ifges(7,6) * t300;
t198 = mrSges(7,2) * t212 - mrSges(7,3) * t207 + Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t265 + t276 * t238 - t300 * t239;
t214 = t332 * t407 - t373;
t252 = Ifges(6,5) * t332 - Ifges(6,6) * t303 + Ifges(6,3) * t302;
t366 = -mrSges(6,2) * t216 + mrSges(6,3) * t214 - Ifges(6,1) * t314 + Ifges(6,4) * t269 - Ifges(6,5) * t268 + pkin(10) * t194 + t350 * t197 - t354 * t198 + t303 * t252;
t209 = -m(7) * t212 + t233 * mrSges(7,1) - t234 * mrSges(7,2) + t276 * t249 - t277 * t250;
t279 = t303 * mrSges(6,1) + t332 * mrSges(6,2);
t368 = -m(6) * t214 + t314 * mrSges(6,3) + t332 * t279 - t209;
t374 = -m(6) * t216 - t269 * mrSges(6,1) - t303 * t274 - t194;
t254 = Ifges(6,4) * t332 - Ifges(6,2) * t303 + Ifges(6,6) * t302;
t392 = Ifges(5,1) * t303 - Ifges(5,4) * t302 + Ifges(5,5) * t332 - t254;
t412 = -mrSges(5,2) * t221 + pkin(4) * (-t314 * mrSges(6,2) - t332 * t278 + t374) + qJ(5) * (-t268 * mrSges(6,1) - t302 * t274 + t368) + mrSges(5,1) * t220 + t303 * t255 - Ifges(5,6) * t268 + Ifges(5,5) * t269 + Ifges(5,3) * t314 - t366 + t392 * t302;
t273 = t302 * mrSges(5,1) + t303 * mrSges(5,2);
t391 = t332 * mrSges(5,2) + t302 * mrSges(5,3) + t278;
t402 = mrSges(5,1) - mrSges(6,2);
t189 = m(5) * t220 - t269 * mrSges(5,3) - t303 * t273 + t402 * t314 - t391 * t332 + t374;
t281 = t332 * mrSges(5,1) - t303 * mrSges(5,3);
t200 = m(5) * t221 - t314 * mrSges(5,2) - t332 * t281 + (-t273 - t274) * t302 + (-mrSges(5,3) - mrSges(6,1)) * t268 + t368;
t184 = t355 * t189 + t351 * t200;
t291 = -t348 * t319 - t404;
t237 = -t326 * pkin(2) + (-t342 * t381 - t325) * qJ(3) + t291 + t409;
t394 = t349 * t356;
t382 = t319 * t394;
t243 = -t341 * pkin(2) + t376 - t382;
t410 = mrSges(4,1) * t243 - mrSges(4,3) * t237 + pkin(3) * t184 + t412;
t195 = -t350 * t204 + t354 * t205;
t218 = t268 * pkin(4) + t360;
t377 = -m(6) * t218 + t269 * mrSges(6,3) + t303 * t279 - t195;
t408 = -m(5) * t225 - t269 * mrSges(5,2) - t402 * t268 - t303 * t281 + t391 * t302 + t377;
t403 = mrSges(3,1) - mrSges(4,2);
t401 = Ifges(3,4) + Ifges(4,6);
t400 = Ifges(5,4) + Ifges(6,6);
t256 = Ifges(6,1) * t332 - Ifges(6,4) * t303 + Ifges(6,5) * t302;
t393 = -Ifges(5,5) * t303 + Ifges(5,6) * t302 - Ifges(5,3) * t332 - t256;
t290 = Ifges(4,1) * t342 + (-Ifges(4,4) * t352 - Ifges(4,5) * t356) * t387;
t390 = Ifges(3,3) * t342 + (Ifges(3,5) * t352 + Ifges(3,6) * t356) * t387 + t290;
t288 = Ifges(4,5) * t342 + (-Ifges(4,6) * t352 - Ifges(4,3) * t356) * t387;
t389 = -Ifges(3,6) * t342 - (Ifges(3,4) * t352 + Ifges(3,2) * t356) * t387 + t288;
t270 = t382 - t388;
t316 = -t342 * mrSges(3,2) + mrSges(3,3) * t381;
t317 = -mrSges(4,1) * t381 - t342 * mrSges(4,3);
t322 = (mrSges(4,2) * t356 - mrSges(4,3) * t352) * t387;
t323 = (-mrSges(3,1) * t356 + mrSges(3,2) * t352) * t387;
t375 = -m(4) * t243 - t325 * mrSges(4,1) - t184;
t182 = m(3) * t270 - t325 * mrSges(3,3) + (t316 - t317) * t342 + t403 * t341 + (-t322 - t323) * t380 + t375;
t315 = t342 * mrSges(3,1) - mrSges(3,3) * t380;
t318 = mrSges(4,1) * t380 + t342 * mrSges(4,2);
t361 = -m(4) * t236 + t341 * mrSges(4,3) + t342 * t318 + t322 * t381 - t408;
t188 = -t342 * t315 + t361 - t341 * mrSges(3,2) + m(3) * t271 + t323 * t381 + (mrSges(3,3) + mrSges(4,1)) * t326;
t179 = -t352 * t182 + t356 * t188;
t185 = -t351 * t189 + t355 * t200;
t378 = m(4) * t237 - t325 * mrSges(4,3) + t317 * t381 + t185;
t181 = m(3) * t291 + t325 * mrSges(3,2) - t403 * t326 + (-t316 * t356 + (t315 - t318) * t352) * t387 + t378;
t173 = -t348 * t181 + t182 * t394 + t188 * t395;
t287 = Ifges(3,5) * t342 + (Ifges(3,1) * t352 + Ifges(3,4) * t356) * t387;
t193 = -t268 * mrSges(6,2) - t302 * t278 - t377;
t364 = -mrSges(6,1) * t214 + mrSges(6,2) * t218 - pkin(5) * t209 - pkin(10) * t195 - t354 * t197 - t350 * t198;
t174 = -mrSges(5,1) * t225 + mrSges(5,3) * t221 - pkin(4) * t193 + t392 * t332 + (Ifges(5,6) - Ifges(6,5)) * t314 + t393 * t303 + t400 * t269 + (-Ifges(5,2) - Ifges(6,3)) * t268 + t364;
t369 = mrSges(7,1) * t207 - mrSges(7,2) * t208 + Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t265 + t277 * t239 - t276 * t240;
t363 = mrSges(6,1) * t216 - mrSges(6,3) * t218 + pkin(5) * t194 + t369;
t176 = t363 + mrSges(5,2) * t225 - mrSges(5,3) * t220 - qJ(5) * t193 - t400 * t268 + (Ifges(5,5) - Ifges(6,4)) * t314 + t393 * t302 + (Ifges(5,1) + Ifges(6,2)) * t269 + (-t255 + t252) * t332;
t289 = Ifges(4,4) * t342 + (-Ifges(4,2) * t352 - Ifges(4,6) * t356) * t387;
t367 = mrSges(4,2) * t243 - mrSges(4,3) * t236 + Ifges(4,1) * t341 - Ifges(4,4) * t325 - Ifges(4,5) * t326 - pkin(9) * t184 - t351 * t174 + t355 * t176 + t289 * t381;
t165 = (-t287 * t356 + (-pkin(2) * t322 - t389) * t352) * t387 + t367 + Ifges(3,3) * t341 + Ifges(3,5) * t325 + Ifges(3,6) * t326 + mrSges(3,1) * t270 - mrSges(3,2) * t271 + pkin(2) * (-t341 * mrSges(4,2) - t342 * t317 + t375) + qJ(3) * (t326 * mrSges(4,1) + t361);
t183 = t326 * mrSges(4,2) - t318 * t380 + t378;
t365 = -mrSges(4,1) * t236 + mrSges(4,2) * t237 - pkin(3) * t408 - pkin(9) * t185 - t355 * t174 - t351 * t176;
t167 = -mrSges(3,1) * t291 + mrSges(3,3) * t271 - pkin(2) * t183 + (t287 - t289) * t342 + (Ifges(3,6) - Ifges(4,5)) * t341 + (Ifges(3,2) + Ifges(4,3)) * t326 + t401 * t325 - t390 * t380 + t365;
t169 = mrSges(3,2) * t291 - mrSges(3,3) * t270 - qJ(3) * t183 + t390 * t381 + (Ifges(3,1) + Ifges(4,2)) * t325 + t401 * t326 + (Ifges(3,5) - Ifges(4,4)) * t341 + t389 * t342 + t410;
t372 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t173 + t349 * t165 + t167 * t396 + t169 * t397 + t179 * t405;
t177 = m(2) * t338 - t358 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t179;
t172 = t349 * t181 + (t182 * t356 + t188 * t352) * t348;
t170 = m(2) * t337 + qJDD(1) * mrSges(2,1) - t358 * mrSges(2,2) + t173;
t163 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) - t352 * t167 + t356 * t169 + (-t172 * t348 - t173 * t349) * pkin(8);
t162 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t172 - t348 * t165 + (pkin(8) * t179 + t167 * t356 + t169 * t352) * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t163 - t353 * t162 - pkin(7) * (t357 * t170 + t353 * t177), t163, t169, -t288 * t380 + t367, t176, -t302 * t254 - t366, t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t353 * t163 + t357 * t162 + pkin(7) * (-t353 * t170 + t357 * t177), t162, t167, Ifges(4,4) * t341 - Ifges(4,2) * t325 - Ifges(4,6) * t326 - t342 * t288 - t290 * t381 - t410, t174, Ifges(6,4) * t314 - Ifges(6,2) * t269 + Ifges(6,6) * t268 - t332 * t252 + t302 * t256 - t363, t197; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372, t372, t165, Ifges(4,5) * t341 - Ifges(4,6) * t325 - Ifges(4,3) * t326 + t342 * t289 + t290 * t380 - t365, t412, Ifges(6,5) * t314 - Ifges(6,6) * t269 + Ifges(6,3) * t268 + t332 * t254 + t303 * t256 - t364, t369;];
m_new  = t1;
