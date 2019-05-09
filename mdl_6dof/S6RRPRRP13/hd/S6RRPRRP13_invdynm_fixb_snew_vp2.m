% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:06:18
% EndTime: 2019-05-06 19:07:01
% DurationCPUTime: 13.81s
% Computational Cost: add. (227270->401), mult. (506561->487), div. (0->0), fcn. (365666->10), ass. (0->154)
t405 = -2 * qJD(3);
t350 = cos(pkin(6));
t343 = qJD(1) * t350 + qJD(2);
t353 = sin(qJ(2));
t349 = sin(pkin(6));
t386 = qJD(1) * t349;
t377 = t353 * t386;
t404 = (pkin(2) * t343 + t405) * t377;
t354 = sin(qJ(1));
t358 = cos(qJ(1));
t338 = g(1) * t354 - g(2) * t358;
t359 = qJD(1) ^ 2;
t401 = pkin(8) * t349;
t321 = qJDD(1) * pkin(1) + t359 * t401 + t338;
t339 = -g(1) * t358 - g(2) * t354;
t383 = qJDD(1) * t349;
t322 = -pkin(1) * t359 + pkin(8) * t383 + t339;
t357 = cos(qJ(2));
t393 = t350 * t353;
t395 = t349 * t353;
t283 = -g(3) * t395 + t321 * t393 + t322 * t357;
t323 = (-pkin(2) * t357 - qJ(3) * t353) * t386;
t341 = t343 ^ 2;
t342 = qJDD(1) * t350 + qJDD(2);
t385 = qJD(1) * t357;
t378 = t349 * t385;
t246 = t341 * pkin(2) - qJ(3) * t342 - t323 * t378 + t343 * t405 - t283;
t352 = sin(qJ(4));
t356 = cos(qJ(4));
t307 = -t343 * t352 - t356 * t378;
t328 = -qJD(2) * t377 + t357 * t383;
t281 = qJD(4) * t307 - t328 * t352 + t342 * t356;
t308 = t343 * t356 - t352 * t378;
t333 = qJD(4) + t377;
t351 = sin(qJ(5));
t355 = cos(qJ(5));
t287 = -t308 * t351 + t333 * t355;
t327 = (qJD(2) * t385 + qJDD(1) * t353) * t349;
t316 = qJDD(4) + t327;
t243 = qJD(5) * t287 + t281 * t355 + t316 * t351;
t288 = t308 * t355 + t333 * t351;
t259 = -mrSges(7,1) * t287 + mrSges(7,2) * t288;
t326 = pkin(3) * t377 - pkin(9) * t343;
t396 = t349 ^ 2 * t359;
t382 = t357 ^ 2 * t396;
t400 = t350 * g(3);
t402 = -pkin(2) - pkin(9);
t227 = -pkin(3) * t382 - t400 - t327 * qJ(3) + t402 * t328 + (-t321 + (-qJ(3) * t343 * t357 - t326 * t353) * qJD(1)) * t349 + t404;
t394 = t349 * t357;
t387 = g(3) * t394 + t322 * t353;
t373 = -t341 * qJ(3) + t323 * t377 + qJDD(3) + t387;
t229 = t327 * pkin(3) + t402 * t342 + (-pkin(3) * t343 * t386 - pkin(9) * t353 * t396 - t321 * t350) * t357 + t373;
t222 = t227 * t356 + t229 * t352;
t285 = -pkin(4) * t307 - pkin(10) * t308;
t331 = t333 ^ 2;
t216 = -pkin(4) * t331 + pkin(10) * t316 + t285 * t307 + t222;
t226 = t328 * pkin(3) - pkin(9) * t382 + t326 * t343 - t246;
t280 = -qJD(4) * t308 - t328 * t356 - t342 * t352;
t219 = (-t307 * t333 - t281) * pkin(10) + (t308 * t333 - t280) * pkin(4) + t226;
t210 = -t351 * t216 + t219 * t355;
t278 = qJDD(5) - t280;
t306 = qJD(5) - t307;
t206 = -0.2e1 * qJD(6) * t288 + (t287 * t306 - t243) * qJ(6) + (t287 * t288 + t278) * pkin(5) + t210;
t263 = -mrSges(7,2) * t306 + mrSges(7,3) * t287;
t380 = m(7) * t206 + mrSges(7,1) * t278 + t263 * t306;
t203 = -t243 * mrSges(7,3) - t288 * t259 + t380;
t211 = t216 * t355 + t219 * t351;
t242 = -qJD(5) * t288 - t281 * t351 + t316 * t355;
t251 = Ifges(6,4) * t288 + Ifges(6,2) * t287 + Ifges(6,6) * t306;
t252 = Ifges(7,1) * t288 + Ifges(7,4) * t287 + Ifges(7,5) * t306;
t253 = Ifges(6,1) * t288 + Ifges(6,4) * t287 + Ifges(6,5) * t306;
t265 = pkin(5) * t306 - qJ(6) * t288;
t286 = t287 ^ 2;
t209 = -pkin(5) * t286 + qJ(6) * t242 + 0.2e1 * qJD(6) * t287 - t265 * t306 + t211;
t250 = Ifges(7,4) * t288 + Ifges(7,2) * t287 + Ifges(7,6) * t306;
t371 = -mrSges(7,1) * t206 + mrSges(7,2) * t209 - Ifges(7,5) * t243 - Ifges(7,6) * t242 - Ifges(7,3) * t278 - t250 * t288;
t403 = mrSges(6,1) * t210 - mrSges(6,2) * t211 + Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t278 + pkin(5) * t203 + t288 * t251 - (t253 + t252) * t287 - t371;
t399 = mrSges(3,1) - mrSges(4,2);
t398 = -mrSges(6,2) - mrSges(7,2);
t397 = Ifges(3,4) + Ifges(4,6);
t392 = t350 * t357;
t260 = -mrSges(6,1) * t287 + mrSges(6,2) * t288;
t264 = -mrSges(6,2) * t306 + mrSges(6,3) * t287;
t196 = m(6) * t210 + t278 * mrSges(6,1) + t306 * t264 + (-t259 - t260) * t288 + (-mrSges(6,3) - mrSges(7,3)) * t243 + t380;
t379 = m(7) * t209 + mrSges(7,3) * t242 + t259 * t287;
t266 = mrSges(7,1) * t306 - mrSges(7,3) * t288;
t390 = -mrSges(6,1) * t306 + mrSges(6,3) * t288 - t266;
t199 = m(6) * t211 + t242 * mrSges(6,3) + t287 * t260 + t278 * t398 + t306 * t390 + t379;
t193 = t196 * t355 + t199 * t351;
t297 = Ifges(4,1) * t343 + (-Ifges(4,4) * t353 - Ifges(4,5) * t357) * t386;
t389 = Ifges(3,3) * t343 + (Ifges(3,5) * t353 + Ifges(3,6) * t357) * t386 + t297;
t295 = Ifges(4,5) * t343 + (-Ifges(4,6) * t353 - Ifges(4,3) * t357) * t386;
t388 = -Ifges(3,6) * t343 - (Ifges(3,4) * t353 + Ifges(3,2) * t357) * t386 + t295;
t381 = t321 * t392;
t282 = t381 - t387;
t318 = -mrSges(3,2) * t343 + mrSges(3,3) * t378;
t319 = -mrSges(4,1) * t378 - mrSges(4,3) * t343;
t324 = (mrSges(4,2) * t357 - mrSges(4,3) * t353) * t386;
t325 = (-mrSges(3,1) * t357 + mrSges(3,2) * t353) * t386;
t284 = -mrSges(5,1) * t307 + mrSges(5,2) * t308;
t290 = mrSges(5,1) * t333 - mrSges(5,3) * t308;
t376 = -t196 * t351 + t199 * t355;
t188 = m(5) * t222 - mrSges(5,2) * t316 + mrSges(5,3) * t280 + t284 * t307 - t290 * t333 + t376;
t221 = -t227 * t352 + t229 * t356;
t289 = -mrSges(5,2) * t333 + mrSges(5,3) * t307;
t215 = -pkin(4) * t316 - pkin(10) * t331 + t285 * t308 - t221;
t213 = -pkin(5) * t242 - qJ(6) * t286 + t265 * t288 + qJDD(6) + t215;
t375 = -m(7) * t213 + mrSges(7,1) * t242 + t263 * t287;
t362 = -m(6) * t215 + mrSges(6,1) * t242 + t243 * t398 + t264 * t287 + t288 * t390 + t375;
t200 = m(5) * t221 + t316 * mrSges(5,1) - t281 * mrSges(5,3) - t308 * t284 + t333 * t289 + t362;
t180 = t352 * t188 + t356 * t200;
t257 = -t342 * pkin(2) + t373 - t381;
t369 = -m(4) * t257 - mrSges(4,1) * t327 - t180;
t178 = m(3) * t282 - t327 * mrSges(3,3) + (t318 - t319) * t343 + t399 * t342 + (-t324 - t325) * t377 + t369;
t317 = mrSges(3,1) * t343 - mrSges(3,3) * t377;
t191 = -m(5) * t226 + t280 * mrSges(5,1) - mrSges(5,2) * t281 + t307 * t289 - t290 * t308 - t193;
t320 = mrSges(4,1) * t377 + mrSges(4,2) * t343;
t363 = -m(4) * t246 + mrSges(4,3) * t342 + t320 * t343 + t324 * t378 - t191;
t186 = t363 + t325 * t378 - t342 * mrSges(3,2) - t343 * t317 + m(3) * t283 + (mrSges(3,3) + mrSges(4,1)) * t328;
t174 = -t178 * t353 + t186 * t357;
t181 = t188 * t356 - t352 * t200;
t298 = -t349 * t321 - t400;
t247 = -t328 * pkin(2) + (-t343 * t378 - t327) * qJ(3) + t298 + t404;
t374 = m(4) * t247 - mrSges(4,3) * t327 + t319 * t378 + t181;
t177 = m(3) * t298 + t327 * mrSges(3,2) - t399 * t328 + (-t318 * t357 + (t317 - t320) * t353) * t386 + t374;
t169 = -t177 * t349 + t178 * t392 + t186 * t393;
t372 = -mrSges(7,1) * t213 + mrSges(7,3) * t209 + Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t278 + t252 * t306;
t248 = Ifges(7,5) * t288 + Ifges(7,6) * t287 + Ifges(7,3) * t306;
t370 = mrSges(7,2) * t213 - mrSges(7,3) * t206 + Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t278 + t248 * t287;
t294 = Ifges(3,5) * t343 + (Ifges(3,1) * t353 + Ifges(3,4) * t357) * t386;
t249 = Ifges(6,5) * t288 + Ifges(6,6) * t287 + Ifges(6,3) * t306;
t183 = Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t278 + t306 * t253 - mrSges(6,1) * t215 + mrSges(6,3) * t211 - pkin(5) * (t243 * mrSges(7,2) - t375) + qJ(6) * (-t278 * mrSges(7,2) - t306 * t266 + t379) + (-pkin(5) * t266 - t248 - t249) * t288 + t372;
t190 = mrSges(6,2) * t215 - mrSges(6,3) * t210 + Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t278 - qJ(6) * t203 + t287 * t249 + (-t250 - t251) * t306 + t370;
t272 = Ifges(5,5) * t308 + Ifges(5,6) * t307 + Ifges(5,3) * t333;
t273 = Ifges(5,4) * t308 + Ifges(5,2) * t307 + Ifges(5,6) * t333;
t171 = mrSges(5,2) * t226 - mrSges(5,3) * t221 + Ifges(5,1) * t281 + Ifges(5,4) * t280 + Ifges(5,5) * t316 - pkin(10) * t193 - t183 * t351 + t190 * t355 + t272 * t307 - t273 * t333;
t274 = Ifges(5,1) * t308 + Ifges(5,4) * t307 + Ifges(5,5) * t333;
t175 = -mrSges(5,1) * t226 + mrSges(5,3) * t222 + Ifges(5,4) * t281 + Ifges(5,2) * t280 + Ifges(5,6) * t316 - pkin(4) * t193 - t308 * t272 + t333 * t274 - t403;
t296 = Ifges(4,4) * t343 + (-Ifges(4,2) * t353 - Ifges(4,6) * t357) * t386;
t366 = mrSges(4,2) * t257 - mrSges(4,3) * t246 + Ifges(4,1) * t342 - Ifges(4,4) * t327 - Ifges(4,5) * t328 - pkin(9) * t180 + t171 * t356 - t352 * t175 + t296 * t378;
t161 = qJ(3) * (mrSges(4,1) * t328 + t363) + (-t294 * t357 + (-pkin(2) * t324 - t388) * t353) * t386 + t366 + Ifges(3,3) * t342 + Ifges(3,5) * t327 + Ifges(3,6) * t328 + mrSges(3,1) * t282 - mrSges(3,2) * t283 + pkin(2) * (-t342 * mrSges(4,2) - t343 * t319 + t369);
t179 = t328 * mrSges(4,2) - t320 * t377 + t374;
t364 = -mrSges(4,1) * t246 + mrSges(4,2) * t247 - pkin(3) * t191 - pkin(9) * t181 - t352 * t171 - t356 * t175;
t163 = -mrSges(3,1) * t298 + mrSges(3,3) * t283 - pkin(2) * t179 + (t294 - t296) * t343 + (Ifges(3,6) - Ifges(4,5)) * t342 + (Ifges(3,2) + Ifges(4,3)) * t328 + t397 * t327 - t389 * t377 + t364;
t365 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t281 - Ifges(5,6) * t280 - Ifges(5,3) * t316 - pkin(4) * t362 - pkin(10) * t376 - t183 * t355 - t190 * t351 - t273 * t308 + t307 * t274;
t361 = -mrSges(4,1) * t257 + mrSges(4,3) * t247 - pkin(3) * t180 + t365;
t165 = t388 * t343 + (Ifges(3,5) - Ifges(4,4)) * t342 + t397 * t328 + (Ifges(3,1) + Ifges(4,2)) * t327 + t389 * t378 + mrSges(3,2) * t298 - t361 - mrSges(3,3) * t282 - qJ(3) * t179;
t368 = mrSges(2,1) * t338 - mrSges(2,2) * t339 + Ifges(2,3) * qJDD(1) + pkin(1) * t169 + t161 * t350 + t163 * t394 + t165 * t395 + t174 * t401;
t172 = m(2) * t339 - mrSges(2,1) * t359 - qJDD(1) * mrSges(2,2) + t174;
t168 = t350 * t177 + (t178 * t357 + t186 * t353) * t349;
t166 = m(2) * t338 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t359 + t169;
t159 = -mrSges(2,2) * g(3) - mrSges(2,3) * t338 + Ifges(2,5) * qJDD(1) - t359 * Ifges(2,6) - t353 * t163 + t357 * t165 + (-t168 * t349 - t169 * t350) * pkin(8);
t158 = mrSges(2,1) * g(3) + mrSges(2,3) * t339 + t359 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t168 - t349 * t161 + (pkin(8) * t174 + t163 * t357 + t165 * t353) * t350;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t358 * t159 - t354 * t158 - pkin(7) * (t166 * t358 + t172 * t354), t159, t165, -t295 * t377 + t366, t171, t190, -t250 * t306 + t370; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t354 * t159 + t358 * t158 + pkin(7) * (-t166 * t354 + t172 * t358), t158, t163, Ifges(4,4) * t342 - Ifges(4,2) * t327 - Ifges(4,6) * t328 - t343 * t295 - t297 * t378 + t361, t175, t183, -t288 * t248 + t372; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t368, t368, t161, Ifges(4,5) * t342 - Ifges(4,6) * t327 - Ifges(4,3) * t328 + t343 * t296 + t297 * t377 - t364, -t365, t403, -t287 * t252 - t371;];
m_new  = t1;
