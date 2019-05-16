% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 01:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 01:08:37
% EndTime: 2019-05-07 01:09:39
% DurationCPUTime: 29.51s
% Computational Cost: add. (515953->406), mult. (1143752->505), div. (0->0), fcn. (850118->12), ass. (0->164)
t408 = -2 * qJD(3);
t358 = cos(pkin(6));
t351 = t358 * qJD(1) + qJD(2);
t362 = sin(qJ(2));
t357 = sin(pkin(6));
t393 = qJD(1) * t357;
t386 = t362 * t393;
t407 = (pkin(2) * t351 + t408) * t386;
t363 = sin(qJ(1));
t368 = cos(qJ(1));
t346 = t363 * g(1) - t368 * g(2);
t369 = qJD(1) ^ 2;
t405 = pkin(8) * t357;
t326 = qJDD(1) * pkin(1) + t369 * t405 + t346;
t347 = -t368 * g(1) - t363 * g(2);
t390 = qJDD(1) * t357;
t327 = -t369 * pkin(1) + pkin(8) * t390 + t347;
t367 = cos(qJ(2));
t398 = t358 * t362;
t400 = t357 * t362;
t286 = -g(3) * t400 + t326 * t398 + t367 * t327;
t328 = (-pkin(2) * t367 - qJ(3) * t362) * t393;
t349 = t351 ^ 2;
t350 = t358 * qJDD(1) + qJDD(2);
t392 = qJD(1) * t367;
t387 = t357 * t392;
t257 = t349 * pkin(2) - t350 * qJ(3) - t328 * t387 + t351 * t408 - t286;
t406 = -pkin(2) - pkin(9);
t404 = t358 * g(3);
t403 = mrSges(3,1) - mrSges(4,2);
t402 = Ifges(3,4) + Ifges(4,6);
t401 = t357 ^ 2 * t369;
t399 = t357 * t367;
t397 = t358 * t367;
t331 = pkin(3) * t386 - t351 * pkin(9);
t332 = (qJD(2) * t392 + qJDD(1) * t362) * t357;
t333 = -qJD(2) * t386 + t367 * t390;
t389 = t367 ^ 2 * t401;
t248 = -pkin(3) * t389 - t404 - t332 * qJ(3) + t406 * t333 + (-t326 + (-qJ(3) * t351 * t367 - t331 * t362) * qJD(1)) * t357 + t407;
t394 = g(3) * t399 + t362 * t327;
t382 = -t349 * qJ(3) + t328 * t386 + qJDD(3) + t394;
t250 = t332 * pkin(3) + t406 * t350 + (-pkin(3) * t351 * t393 - pkin(9) * t362 * t401 - t326 * t358) * t367 + t382;
t361 = sin(qJ(4));
t366 = cos(qJ(4));
t237 = t366 * t248 + t361 * t250;
t312 = -t361 * t351 - t366 * t387;
t313 = t366 * t351 - t361 * t387;
t288 = -pkin(4) * t312 - pkin(10) * t313;
t321 = qJDD(4) + t332;
t341 = qJD(4) + t386;
t338 = t341 ^ 2;
t226 = -pkin(4) * t338 + pkin(10) * t321 + t288 * t312 + t237;
t247 = t333 * pkin(3) - pkin(9) * t389 + t351 * t331 - t257;
t283 = -t313 * qJD(4) - t333 * t366 - t361 * t350;
t284 = qJD(4) * t312 - t333 * t361 + t350 * t366;
t232 = (-t312 * t341 - t284) * pkin(10) + (t313 * t341 - t283) * pkin(4) + t247;
t360 = sin(qJ(5));
t365 = cos(qJ(5));
t221 = -t360 * t226 + t365 * t232;
t290 = -t313 * t360 + t341 * t365;
t253 = qJD(5) * t290 + t284 * t365 + t321 * t360;
t281 = qJDD(5) - t283;
t291 = t313 * t365 + t341 * t360;
t311 = qJD(5) - t312;
t219 = (t290 * t311 - t253) * pkin(11) + (t290 * t291 + t281) * pkin(5) + t221;
t222 = t365 * t226 + t360 * t232;
t252 = -qJD(5) * t291 - t284 * t360 + t321 * t365;
t272 = pkin(5) * t311 - pkin(11) * t291;
t289 = t290 ^ 2;
t220 = -pkin(5) * t289 + pkin(11) * t252 - t272 * t311 + t222;
t359 = sin(qJ(6));
t364 = cos(qJ(6));
t217 = t219 * t364 - t220 * t359;
t265 = t290 * t364 - t291 * t359;
t234 = qJD(6) * t265 + t252 * t359 + t253 * t364;
t266 = t290 * t359 + t291 * t364;
t243 = -mrSges(7,1) * t265 + mrSges(7,2) * t266;
t309 = qJD(6) + t311;
t254 = -mrSges(7,2) * t309 + mrSges(7,3) * t265;
t274 = qJDD(6) + t281;
t212 = m(7) * t217 + mrSges(7,1) * t274 - mrSges(7,3) * t234 - t243 * t266 + t254 * t309;
t218 = t219 * t359 + t220 * t364;
t233 = -qJD(6) * t266 + t252 * t364 - t253 * t359;
t255 = mrSges(7,1) * t309 - mrSges(7,3) * t266;
t213 = m(7) * t218 - mrSges(7,2) * t274 + mrSges(7,3) * t233 + t243 * t265 - t255 * t309;
t205 = t364 * t212 + t359 * t213;
t267 = -mrSges(6,1) * t290 + mrSges(6,2) * t291;
t270 = -mrSges(6,2) * t311 + mrSges(6,3) * t290;
t203 = m(6) * t221 + mrSges(6,1) * t281 - mrSges(6,3) * t253 - t267 * t291 + t270 * t311 + t205;
t271 = mrSges(6,1) * t311 - mrSges(6,3) * t291;
t384 = -t212 * t359 + t364 * t213;
t204 = m(6) * t222 - mrSges(6,2) * t281 + mrSges(6,3) * t252 + t267 * t290 - t271 * t311 + t384;
t198 = t365 * t203 + t360 * t204;
t300 = Ifges(4,1) * t351 + (-Ifges(4,4) * t362 - Ifges(4,5) * t367) * t393;
t396 = Ifges(3,3) * t351 + (Ifges(3,5) * t362 + Ifges(3,6) * t367) * t393 + t300;
t298 = Ifges(4,5) * t351 + (-Ifges(4,6) * t362 - Ifges(4,3) * t367) * t393;
t395 = -Ifges(3,6) * t351 - (Ifges(3,4) * t362 + Ifges(3,2) * t367) * t393 + t298;
t388 = t326 * t397;
t285 = t388 - t394;
t323 = -t351 * mrSges(3,2) + mrSges(3,3) * t387;
t324 = -mrSges(4,1) * t387 - t351 * mrSges(4,3);
t329 = (mrSges(4,2) * t367 - mrSges(4,3) * t362) * t393;
t330 = (-mrSges(3,1) * t367 + mrSges(3,2) * t362) * t393;
t287 = -mrSges(5,1) * t312 + mrSges(5,2) * t313;
t293 = mrSges(5,1) * t341 - mrSges(5,3) * t313;
t385 = -t360 * t203 + t365 * t204;
t195 = m(5) * t237 - mrSges(5,2) * t321 + mrSges(5,3) * t283 + t287 * t312 - t293 * t341 + t385;
t236 = -t361 * t248 + t250 * t366;
t292 = -mrSges(5,2) * t341 + mrSges(5,3) * t312;
t225 = -pkin(4) * t321 - pkin(10) * t338 + t313 * t288 - t236;
t223 = -pkin(5) * t252 - pkin(11) * t289 + t272 * t291 + t225;
t380 = m(7) * t223 - t233 * mrSges(7,1) + mrSges(7,2) * t234 - t265 * t254 + t255 * t266;
t372 = -m(6) * t225 + t252 * mrSges(6,1) - mrSges(6,2) * t253 + t290 * t270 - t271 * t291 - t380;
t208 = m(5) * t236 + mrSges(5,1) * t321 - mrSges(5,3) * t284 - t287 * t313 + t292 * t341 + t372;
t185 = t361 * t195 + t366 * t208;
t264 = -t350 * pkin(2) + t382 - t388;
t381 = -m(4) * t264 - t332 * mrSges(4,1) - t185;
t183 = m(3) * t285 - t332 * mrSges(3,3) + (t323 - t324) * t351 + t403 * t350 + (-t329 - t330) * t386 + t381;
t322 = t351 * mrSges(3,1) - mrSges(3,3) * t386;
t196 = -m(5) * t247 + t283 * mrSges(5,1) - t284 * mrSges(5,2) + t312 * t292 - t313 * t293 - t198;
t325 = mrSges(4,1) * t386 + t351 * mrSges(4,2);
t373 = -m(4) * t257 + t350 * mrSges(4,3) + t351 * t325 + t329 * t387 - t196;
t193 = t373 + (mrSges(3,3) + mrSges(4,1)) * t333 + t330 * t387 + m(3) * t286 - t350 * mrSges(3,2) - t351 * t322;
t179 = -t183 * t362 + t367 * t193;
t186 = t366 * t195 - t361 * t208;
t301 = -t357 * t326 - t404;
t258 = -t333 * pkin(2) + (-t351 * t387 - t332) * qJ(3) + t301 + t407;
t383 = m(4) * t258 - t332 * mrSges(4,3) + t324 * t387 + t186;
t182 = m(3) * t301 + t332 * mrSges(3,2) - t403 * t333 + (-t323 * t367 + (t322 - t325) * t362) * t393 + t383;
t174 = -t182 * t357 + t183 * t397 + t193 * t398;
t297 = Ifges(3,5) * t351 + (Ifges(3,1) * t362 + Ifges(3,4) * t367) * t393;
t239 = Ifges(7,5) * t266 + Ifges(7,6) * t265 + Ifges(7,3) * t309;
t241 = Ifges(7,1) * t266 + Ifges(7,4) * t265 + Ifges(7,5) * t309;
t206 = -mrSges(7,1) * t223 + mrSges(7,3) * t218 + Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t274 - t239 * t266 + t241 * t309;
t240 = Ifges(7,4) * t266 + Ifges(7,2) * t265 + Ifges(7,6) * t309;
t207 = mrSges(7,2) * t223 - mrSges(7,3) * t217 + Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t274 + t239 * t265 - t240 * t309;
t259 = Ifges(6,5) * t291 + Ifges(6,6) * t290 + Ifges(6,3) * t311;
t261 = Ifges(6,1) * t291 + Ifges(6,4) * t290 + Ifges(6,5) * t311;
t188 = -mrSges(6,1) * t225 + mrSges(6,3) * t222 + Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t281 - pkin(5) * t380 + pkin(11) * t384 + t364 * t206 + t359 * t207 - t291 * t259 + t311 * t261;
t260 = Ifges(6,4) * t291 + Ifges(6,2) * t290 + Ifges(6,6) * t311;
t190 = mrSges(6,2) * t225 - mrSges(6,3) * t221 + Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t281 - pkin(11) * t205 - t206 * t359 + t207 * t364 + t259 * t290 - t260 * t311;
t275 = Ifges(5,5) * t313 + Ifges(5,6) * t312 + Ifges(5,3) * t341;
t276 = Ifges(5,4) * t313 + Ifges(5,2) * t312 + Ifges(5,6) * t341;
t176 = mrSges(5,2) * t247 - mrSges(5,3) * t236 + Ifges(5,1) * t284 + Ifges(5,4) * t283 + Ifges(5,5) * t321 - pkin(10) * t198 - t188 * t360 + t190 * t365 + t275 * t312 - t276 * t341;
t277 = Ifges(5,1) * t313 + Ifges(5,4) * t312 + Ifges(5,5) * t341;
t377 = -mrSges(7,1) * t217 + mrSges(7,2) * t218 - Ifges(7,5) * t234 - Ifges(7,6) * t233 - Ifges(7,3) * t274 - t266 * t240 + t265 * t241;
t370 = mrSges(6,1) * t221 - mrSges(6,2) * t222 + Ifges(6,5) * t253 + Ifges(6,6) * t252 + Ifges(6,3) * t281 + pkin(5) * t205 + t291 * t260 - t290 * t261 - t377;
t180 = -mrSges(5,1) * t247 + mrSges(5,3) * t237 + Ifges(5,4) * t284 + Ifges(5,2) * t283 + Ifges(5,6) * t321 - pkin(4) * t198 - t313 * t275 + t341 * t277 - t370;
t299 = Ifges(4,4) * t351 + (-Ifges(4,2) * t362 - Ifges(4,6) * t367) * t393;
t376 = mrSges(4,2) * t264 - mrSges(4,3) * t257 + Ifges(4,1) * t350 - Ifges(4,4) * t332 - Ifges(4,5) * t333 - pkin(9) * t185 + t366 * t176 - t361 * t180 + t299 * t387;
t166 = pkin(2) * (-t350 * mrSges(4,2) - t351 * t324 + t381) + t376 + qJ(3) * (mrSges(4,1) * t333 + t373) + (-t297 * t367 + (-pkin(2) * t329 - t395) * t362) * t393 + mrSges(3,1) * t285 - mrSges(3,2) * t286 + Ifges(3,5) * t332 + Ifges(3,6) * t333 + Ifges(3,3) * t350;
t184 = t333 * mrSges(4,2) - t325 * t386 + t383;
t374 = -mrSges(4,1) * t257 + mrSges(4,2) * t258 - pkin(3) * t196 - pkin(9) * t186 - t361 * t176 - t366 * t180;
t168 = -mrSges(3,1) * t301 + mrSges(3,3) * t286 - pkin(2) * t184 + (t297 - t299) * t351 + (Ifges(3,6) - Ifges(4,5)) * t350 + (Ifges(3,2) + Ifges(4,3)) * t333 + t402 * t332 - t396 * t386 + t374;
t375 = -mrSges(5,1) * t236 + mrSges(5,2) * t237 - Ifges(5,5) * t284 - Ifges(5,6) * t283 - Ifges(5,3) * t321 - pkin(4) * t372 - pkin(10) * t385 - t365 * t188 - t360 * t190 - t313 * t276 + t312 * t277;
t371 = -mrSges(4,1) * t264 + mrSges(4,3) * t258 - pkin(3) * t185 + t375;
t170 = -qJ(3) * t184 + t396 * t387 - t371 + t395 * t351 + (Ifges(3,5) - Ifges(4,4)) * t350 + t402 * t333 + (Ifges(3,1) + Ifges(4,2)) * t332 - mrSges(3,3) * t285 + mrSges(3,2) * t301;
t379 = mrSges(2,1) * t346 - mrSges(2,2) * t347 + Ifges(2,3) * qJDD(1) + pkin(1) * t174 + t358 * t166 + t168 * t399 + t170 * t400 + t179 * t405;
t177 = m(2) * t347 - mrSges(2,1) * t369 - qJDD(1) * mrSges(2,2) + t179;
t173 = t358 * t182 + (t183 * t367 + t193 * t362) * t357;
t171 = m(2) * t346 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t369 + t174;
t164 = -mrSges(2,2) * g(3) - mrSges(2,3) * t346 + Ifges(2,5) * qJDD(1) - t369 * Ifges(2,6) - t362 * t168 + t367 * t170 + (-t173 * t357 - t174 * t358) * pkin(8);
t163 = mrSges(2,1) * g(3) + mrSges(2,3) * t347 + t369 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t173 - t357 * t166 + (pkin(8) * t179 + t168 * t367 + t170 * t362) * t358;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t368 * t164 - t363 * t163 - pkin(7) * (t171 * t368 + t177 * t363), t164, t170, -t298 * t386 + t376, t176, t190, t207; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t363 * t164 + t368 * t163 + pkin(7) * (-t363 * t171 + t177 * t368), t163, t168, Ifges(4,4) * t350 - Ifges(4,2) * t332 - Ifges(4,6) * t333 - t351 * t298 - t300 * t387 + t371, t180, t188, t206; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t379, t379, t166, Ifges(4,5) * t350 - Ifges(4,6) * t332 - Ifges(4,3) * t333 + t351 * t299 + t300 * t386 - t374, -t375, t370, -t377;];
m_new  = t1;
