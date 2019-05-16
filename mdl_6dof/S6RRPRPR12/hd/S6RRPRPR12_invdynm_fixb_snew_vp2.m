% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:14:08
% EndTime: 2019-05-06 16:14:59
% DurationCPUTime: 27.62s
% Computational Cost: add. (476025->406), mult. (1093152->507), div. (0->0), fcn. (811164->12), ass. (0->164)
t402 = -2 * qJD(3);
t352 = cos(pkin(6));
t343 = qJD(1) * t352 + qJD(2);
t355 = sin(qJ(2));
t350 = sin(pkin(6));
t386 = qJD(1) * t350;
t380 = t355 * t386;
t401 = (pkin(2) * t343 + t402) * t380;
t356 = sin(qJ(1));
t360 = cos(qJ(1));
t338 = t356 * g(1) - g(2) * t360;
t361 = qJD(1) ^ 2;
t398 = pkin(8) * t350;
t321 = qJDD(1) * pkin(1) + t361 * t398 + t338;
t339 = -g(1) * t360 - g(2) * t356;
t383 = qJDD(1) * t350;
t322 = -pkin(1) * t361 + pkin(8) * t383 + t339;
t359 = cos(qJ(2));
t391 = t352 * t355;
t393 = t350 * t355;
t281 = -g(3) * t393 + t321 * t391 + t359 * t322;
t323 = (-pkin(2) * t359 - qJ(3) * t355) * t386;
t341 = t343 ^ 2;
t342 = qJDD(1) * t352 + qJDD(2);
t385 = qJD(1) * t359;
t379 = t350 * t385;
t259 = t341 * pkin(2) - t342 * qJ(3) - t323 * t379 + t343 * t402 - t281;
t400 = 2 * qJD(5);
t399 = -pkin(2) - pkin(9);
t397 = t352 * g(3);
t396 = mrSges(3,1) - mrSges(4,2);
t395 = Ifges(3,4) + Ifges(4,6);
t394 = t350 ^ 2 * t361;
t392 = t350 * t359;
t390 = t352 * t359;
t326 = pkin(3) * t380 - pkin(9) * t343;
t327 = (qJD(2) * t385 + qJDD(1) * t355) * t350;
t328 = -qJD(2) * t380 + t359 * t383;
t382 = t359 ^ 2 * t394;
t241 = -pkin(3) * t382 - t397 - t327 * qJ(3) + t399 * t328 + (-t321 + (-qJ(3) * t343 * t359 - t326 * t355) * qJD(1)) * t350 + t401;
t387 = g(3) * t392 + t355 * t322;
t374 = -t341 * qJ(3) + t323 * t380 + qJDD(3) + t387;
t244 = t327 * pkin(3) + t399 * t342 + (-pkin(3) * t343 * t386 - pkin(9) * t355 * t394 - t321 * t352) * t359 + t374;
t354 = sin(qJ(4));
t358 = cos(qJ(4));
t226 = -t354 * t241 + t358 * t244;
t306 = -t343 * t354 - t358 * t379;
t279 = qJD(4) * t306 - t328 * t354 + t342 * t358;
t307 = t343 * t358 - t354 * t379;
t316 = qJDD(4) + t327;
t333 = qJD(4) + t380;
t222 = (t306 * t333 - t279) * qJ(5) + (t306 * t307 + t316) * pkin(4) + t226;
t227 = t358 * t241 + t354 * t244;
t278 = -qJD(4) * t307 - t328 * t358 - t342 * t354;
t288 = pkin(4) * t333 - qJ(5) * t307;
t305 = t306 ^ 2;
t224 = -pkin(4) * t305 + qJ(5) * t278 - t288 * t333 + t227;
t349 = sin(pkin(11));
t351 = cos(pkin(11));
t284 = t306 * t351 - t307 * t349;
t219 = t349 * t222 + t351 * t224 + t284 * t400;
t252 = t278 * t351 - t279 * t349;
t285 = t306 * t349 + t307 * t351;
t261 = -mrSges(6,1) * t284 + mrSges(6,2) * t285;
t269 = mrSges(6,1) * t333 - mrSges(6,3) * t285;
t262 = -pkin(5) * t284 - pkin(10) * t285;
t331 = t333 ^ 2;
t216 = -pkin(5) * t331 + pkin(10) * t316 + t262 * t284 + t219;
t240 = t328 * pkin(3) - pkin(9) * t382 + t343 * t326 - t259;
t229 = -t278 * pkin(4) - t305 * qJ(5) + t307 * t288 + qJDD(5) + t240;
t253 = t278 * t349 + t279 * t351;
t220 = t229 + (-t284 * t333 - t253) * pkin(10) + (t285 * t333 - t252) * pkin(5);
t353 = sin(qJ(6));
t357 = cos(qJ(6));
t213 = -t216 * t353 + t220 * t357;
t266 = -t285 * t353 + t333 * t357;
t232 = qJD(6) * t266 + t253 * t357 + t316 * t353;
t267 = t285 * t357 + t333 * t353;
t245 = -mrSges(7,1) * t266 + mrSges(7,2) * t267;
t283 = qJD(6) - t284;
t246 = -mrSges(7,2) * t283 + mrSges(7,3) * t266;
t251 = qJDD(6) - t252;
t209 = m(7) * t213 + mrSges(7,1) * t251 - mrSges(7,3) * t232 - t245 * t267 + t246 * t283;
t214 = t216 * t357 + t220 * t353;
t231 = -qJD(6) * t267 - t253 * t353 + t316 * t357;
t247 = mrSges(7,1) * t283 - mrSges(7,3) * t267;
t210 = m(7) * t214 - mrSges(7,2) * t251 + mrSges(7,3) * t231 + t245 * t266 - t247 * t283;
t377 = -t209 * t353 + t357 * t210;
t195 = m(6) * t219 - mrSges(6,2) * t316 + mrSges(6,3) * t252 + t261 * t284 - t269 * t333 + t377;
t376 = -t351 * t222 + t349 * t224;
t218 = -0.2e1 * qJD(5) * t285 - t376;
t268 = -mrSges(6,2) * t333 + mrSges(6,3) * t284;
t215 = -t316 * pkin(5) - t331 * pkin(10) + (t400 + t262) * t285 + t376;
t372 = -m(7) * t215 + t231 * mrSges(7,1) - mrSges(7,2) * t232 + t266 * t246 - t247 * t267;
t205 = m(6) * t218 + mrSges(6,1) * t316 - mrSges(6,3) * t253 - t261 * t285 + t268 * t333 + t372;
t189 = t349 * t195 + t351 * t205;
t198 = t357 * t209 + t353 * t210;
t296 = Ifges(4,1) * t343 + (-Ifges(4,4) * t355 - Ifges(4,5) * t359) * t386;
t389 = Ifges(3,3) * t343 + (Ifges(3,5) * t355 + Ifges(3,6) * t359) * t386 + t296;
t294 = Ifges(4,5) * t343 + (-Ifges(4,6) * t355 - Ifges(4,3) * t359) * t386;
t388 = -Ifges(3,6) * t343 - (Ifges(3,4) * t355 + Ifges(3,2) * t359) * t386 + t294;
t381 = t321 * t390;
t280 = t381 - t387;
t318 = -mrSges(3,2) * t343 + mrSges(3,3) * t379;
t319 = -mrSges(4,1) * t379 - mrSges(4,3) * t343;
t324 = (mrSges(4,2) * t359 - mrSges(4,3) * t355) * t386;
t325 = (-mrSges(3,1) * t359 + mrSges(3,2) * t355) * t386;
t286 = -mrSges(5,1) * t306 + mrSges(5,2) * t307;
t287 = -mrSges(5,2) * t333 + mrSges(5,3) * t306;
t186 = m(5) * t226 + mrSges(5,1) * t316 - mrSges(5,3) * t279 - t286 * t307 + t287 * t333 + t189;
t289 = mrSges(5,1) * t333 - mrSges(5,3) * t307;
t378 = t351 * t195 - t205 * t349;
t187 = m(5) * t227 - mrSges(5,2) * t316 + mrSges(5,3) * t278 + t286 * t306 - t289 * t333 + t378;
t181 = t358 * t186 + t354 * t187;
t264 = -t342 * pkin(2) + t374 - t381;
t373 = -m(4) * t264 - t327 * mrSges(4,1) - t181;
t179 = m(3) * t280 - t327 * mrSges(3,3) + (t318 - t319) * t343 + t396 * t342 + (-t324 - t325) * t380 + t373;
t317 = mrSges(3,1) * t343 - mrSges(3,3) * t380;
t369 = m(6) * t229 - t252 * mrSges(6,1) + t253 * mrSges(6,2) - t284 * t268 + t285 * t269 + t198;
t196 = -m(5) * t240 + t278 * mrSges(5,1) - t279 * mrSges(5,2) + t306 * t287 - t307 * t289 - t369;
t320 = mrSges(4,1) * t380 + mrSges(4,2) * t343;
t364 = -m(4) * t259 + t342 * mrSges(4,3) + t343 * t320 + t324 * t379 - t196;
t192 = -t343 * t317 - t342 * mrSges(3,2) + m(3) * t281 + t325 * t379 + t364 + (mrSges(3,3) + mrSges(4,1)) * t328;
t176 = -t179 * t355 + t359 * t192;
t182 = -t354 * t186 + t358 * t187;
t297 = -t350 * t321 - t397;
t260 = -t328 * pkin(2) + (-t343 * t379 - t327) * qJ(3) + t297 + t401;
t375 = m(4) * t260 - t327 * mrSges(4,3) + t319 * t379 + t182;
t178 = m(3) * t297 + t327 * mrSges(3,2) - t396 * t328 + (-t318 * t359 + (t317 - t320) * t355) * t386 + t375;
t170 = -t178 * t350 + t179 * t390 + t192 * t391;
t293 = Ifges(3,5) * t343 + (Ifges(3,1) * t355 + Ifges(3,4) * t359) * t386;
t233 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t283;
t235 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t283;
t202 = -mrSges(7,1) * t215 + mrSges(7,3) * t214 + Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t251 - t233 * t267 + t235 * t283;
t234 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t283;
t203 = mrSges(7,2) * t215 - mrSges(7,3) * t213 + Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t251 + t233 * t266 - t234 * t283;
t254 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t333;
t255 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t333;
t183 = mrSges(6,2) * t229 - mrSges(6,3) * t218 + Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t316 - pkin(10) * t198 - t202 * t353 + t203 * t357 + t254 * t284 - t255 * t333;
t256 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t333;
t365 = mrSges(7,1) * t213 - mrSges(7,2) * t214 + Ifges(7,5) * t232 + Ifges(7,6) * t231 + Ifges(7,3) * t251 + t234 * t267 - t235 * t266;
t184 = -mrSges(6,1) * t229 + mrSges(6,3) * t219 + Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t316 - pkin(5) * t198 - t254 * t285 + t256 * t333 - t365;
t270 = Ifges(5,5) * t307 + Ifges(5,6) * t306 + Ifges(5,3) * t333;
t272 = Ifges(5,1) * t307 + Ifges(5,4) * t306 + Ifges(5,5) * t333;
t171 = -mrSges(5,1) * t240 + mrSges(5,3) * t227 + Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t316 - pkin(4) * t369 + qJ(5) * t378 + t349 * t183 + t351 * t184 - t307 * t270 + t333 * t272;
t271 = Ifges(5,4) * t307 + Ifges(5,2) * t306 + Ifges(5,6) * t333;
t173 = mrSges(5,2) * t240 - mrSges(5,3) * t226 + Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t316 - qJ(5) * t189 + t183 * t351 - t184 * t349 + t270 * t306 - t271 * t333;
t295 = Ifges(4,4) * t343 + (-Ifges(4,2) * t355 - Ifges(4,6) * t359) * t386;
t368 = mrSges(4,2) * t264 - mrSges(4,3) * t259 + Ifges(4,1) * t342 - Ifges(4,4) * t327 - Ifges(4,5) * t328 - pkin(9) * t181 - t354 * t171 + t358 * t173 + t295 * t379;
t162 = qJ(3) * (t328 * mrSges(4,1) + t364) + Ifges(3,3) * t342 + Ifges(3,5) * t327 + Ifges(3,6) * t328 + mrSges(3,1) * t280 - mrSges(3,2) * t281 + (-t293 * t359 + (-pkin(2) * t324 - t388) * t355) * t386 + t368 + pkin(2) * (-t342 * mrSges(4,2) - t343 * t319 + t373);
t180 = t328 * mrSges(4,2) - t320 * t380 + t375;
t366 = -mrSges(4,1) * t259 + mrSges(4,2) * t260 - pkin(3) * t196 - pkin(9) * t182 - t358 * t171 - t354 * t173;
t164 = -mrSges(3,1) * t297 + mrSges(3,3) * t281 - pkin(2) * t180 + (t293 - t295) * t343 + (Ifges(3,6) - Ifges(4,5)) * t342 + (Ifges(3,2) + Ifges(4,3)) * t328 + t395 * t327 - t389 * t380 + t366;
t367 = -mrSges(6,1) * t218 + mrSges(6,2) * t219 - Ifges(6,5) * t253 - Ifges(6,6) * t252 - Ifges(6,3) * t316 - pkin(5) * t372 - pkin(10) * t377 - t357 * t202 - t353 * t203 - t285 * t255 + t284 * t256;
t363 = -mrSges(5,1) * t226 + mrSges(5,2) * t227 - Ifges(5,5) * t279 - Ifges(5,6) * t278 - Ifges(5,3) * t316 - pkin(4) * t189 - t307 * t271 + t306 * t272 + t367;
t362 = -mrSges(4,1) * t264 + mrSges(4,3) * t260 - pkin(3) * t181 + t363;
t166 = mrSges(3,2) * t297 - mrSges(3,3) * t280 + t389 * t379 - qJ(3) * t180 - t362 + (Ifges(3,1) + Ifges(4,2)) * t327 + t395 * t328 + (Ifges(3,5) - Ifges(4,4)) * t342 + t388 * t343;
t371 = mrSges(2,1) * t338 - mrSges(2,2) * t339 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t352 * t162 + t164 * t392 + t166 * t393 + t176 * t398;
t174 = m(2) * t339 - mrSges(2,1) * t361 - qJDD(1) * mrSges(2,2) + t176;
t169 = t352 * t178 + (t179 * t359 + t192 * t355) * t350;
t167 = m(2) * t338 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t361 + t170;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t338 + Ifges(2,5) * qJDD(1) - t361 * Ifges(2,6) - t355 * t164 + t359 * t166 + (-t169 * t350 - t170 * t352) * pkin(8);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t339 + t361 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t350 * t162 + (pkin(8) * t176 + t164 * t359 + t166 * t355) * t352;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t360 * t160 - t356 * t159 - pkin(7) * (t167 * t360 + t174 * t356), t160, t166, -t294 * t380 + t368, t173, t183, t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t356 * t160 + t360 * t159 + pkin(7) * (-t167 * t356 + t174 * t360), t159, t164, Ifges(4,4) * t342 - Ifges(4,2) * t327 - Ifges(4,6) * t328 - t343 * t294 - t296 * t379 + t362, t171, t184, t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t371, t371, t162, Ifges(4,5) * t342 - Ifges(4,6) * t327 - Ifges(4,3) * t328 + t343 * t295 + t296 * t380 - t366, -t363, -t367, t365;];
m_new  = t1;
