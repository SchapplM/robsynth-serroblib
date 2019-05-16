% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:15:27
% EndTime: 2019-05-06 02:17:11
% DurationCPUTime: 82.30s
% Computational Cost: add. (1294861->398), mult. (4017519->532), div. (0->0), fcn. (3417092->14), ass. (0->174)
t342 = sin(pkin(12));
t344 = sin(pkin(6));
t345 = cos(pkin(12));
t347 = cos(pkin(6));
t350 = sin(qJ(3));
t346 = cos(pkin(7));
t353 = cos(qJ(3));
t391 = t346 * t353;
t343 = sin(pkin(7));
t396 = t343 * t353;
t359 = t344 * (-t342 * t350 + t345 * t391) + t347 * t396;
t315 = t359 * qJD(1);
t392 = t346 * t350;
t397 = t343 * t350;
t361 = t347 * t397 + (t342 * t353 + t345 * t392) * t344;
t316 = t361 * qJD(1);
t302 = -qJD(3) * t316 + qJDD(1) * t359;
t394 = t344 * t346;
t327 = (t343 * t347 + t345 * t394) * qJD(1) * pkin(9);
t351 = sin(qJ(1));
t354 = cos(qJ(1));
t339 = -g(1) * t354 - g(2) * t351;
t355 = qJD(1) ^ 2;
t400 = qJ(2) * t344;
t331 = -pkin(1) * t355 + qJDD(1) * t400 + t339;
t404 = pkin(9) * t342;
t373 = -pkin(2) * t345 - t343 * t404;
t387 = qJD(1) * t344;
t401 = pkin(9) * qJDD(1);
t368 = qJD(1) * t373 * t387 + t346 * t401;
t338 = t351 * g(1) - g(2) * t354;
t330 = qJDD(1) * pkin(1) + t355 * t400 + t338;
t382 = qJD(2) * t387;
t393 = t345 * t347;
t395 = t344 * t345;
t374 = -g(3) * t395 + t330 * t393 - 0.2e1 * t342 * t382;
t283 = (pkin(2) * qJDD(1) + qJD(1) * t327) * t347 + (-t344 * t368 - t331) * t342 + t374;
t332 = (pkin(2) * t347 - t394 * t404) * qJD(1);
t398 = t342 * t347;
t383 = t330 * t398 + (t331 + 0.2e1 * t382) * t345;
t284 = (-qJD(1) * t332 + t343 * t401) * t347 + (-g(3) * t342 + t345 * t368) * t344 + t383;
t381 = -g(3) * t347 + qJDD(2);
t292 = (-t330 + t373 * qJDD(1) + (-t327 * t345 + t332 * t342) * qJD(1)) * t344 + t381;
t237 = -t350 * t284 + (t283 * t346 + t292 * t343) * t353;
t238 = t283 * t392 + t353 * t284 + t292 * t397;
t301 = -pkin(3) * t315 - pkin(10) * t316;
t369 = -t343 * t395 + t346 * t347;
t328 = qJD(1) * t369 + qJD(3);
t324 = t328 ^ 2;
t325 = qJDD(1) * t369 + qJDD(3);
t234 = -pkin(3) * t324 + pkin(10) * t325 + t301 * t315 + t238;
t257 = -t283 * t343 + t346 * t292;
t303 = qJD(3) * t315 + qJDD(1) * t361;
t236 = (-t315 * t328 - t303) * pkin(10) + (t316 * t328 - t302) * pkin(3) + t257;
t349 = sin(qJ(4));
t352 = cos(qJ(4));
t230 = t352 * t234 + t349 * t236;
t308 = -t316 * t349 + t328 * t352;
t309 = t316 * t352 + t328 * t349;
t286 = -pkin(4) * t308 - pkin(11) * t309;
t299 = qJDD(4) - t302;
t314 = qJD(4) - t315;
t313 = t314 ^ 2;
t226 = -pkin(4) * t313 + pkin(11) * t299 + t286 * t308 + t230;
t233 = -pkin(3) * t325 - pkin(10) * t324 + t316 * t301 - t237;
t278 = -qJD(4) * t309 - t303 * t349 + t325 * t352;
t279 = qJD(4) * t308 + t303 * t352 + t325 * t349;
t228 = (-t308 * t314 - t279) * pkin(11) + (t309 * t314 - t278) * pkin(4) + t233;
t348 = sin(qJ(5));
t405 = cos(qJ(5));
t222 = -t348 * t226 + t405 * t228;
t223 = t405 * t226 + t348 * t228;
t291 = t405 * t309 + t348 * t314;
t246 = qJD(5) * t291 + t279 * t348 - t405 * t299;
t290 = t309 * t348 - t405 * t314;
t247 = -t290 * qJD(5) + t405 * t279 + t348 * t299;
t305 = qJD(5) - t308;
t249 = Ifges(7,5) * t291 + Ifges(7,6) * t305 + Ifges(7,3) * t290;
t252 = Ifges(6,4) * t291 - Ifges(6,2) * t290 + Ifges(6,6) * t305;
t254 = Ifges(6,1) * t291 - Ifges(6,4) * t290 + Ifges(6,5) * t305;
t261 = mrSges(7,1) * t290 - mrSges(7,3) * t291;
t276 = qJDD(5) - t278;
t260 = pkin(5) * t290 - qJ(6) * t291;
t304 = t305 ^ 2;
t218 = -pkin(5) * t304 + qJ(6) * t276 + 0.2e1 * qJD(6) * t305 - t260 * t290 + t223;
t220 = -t276 * pkin(5) - t304 * qJ(6) + t291 * t260 + qJDD(6) - t222;
t253 = Ifges(7,1) * t291 + Ifges(7,4) * t305 + Ifges(7,5) * t290;
t366 = mrSges(7,1) * t220 - mrSges(7,3) * t218 - Ifges(7,4) * t247 - Ifges(7,2) * t276 - Ifges(7,6) * t246 - t290 * t253;
t264 = -mrSges(7,2) * t290 + mrSges(7,3) * t305;
t379 = -m(7) * t220 + t276 * mrSges(7,1) + t305 * t264;
t267 = -mrSges(7,1) * t305 + mrSges(7,2) * t291;
t384 = m(7) * t218 + t276 * mrSges(7,3) + t305 * t267;
t406 = -(-t252 + t249) * t291 + mrSges(6,1) * t222 - mrSges(6,2) * t223 + Ifges(6,5) * t247 - Ifges(6,6) * t246 + Ifges(6,3) * t276 + pkin(5) * (-mrSges(7,2) * t247 - t261 * t291 + t379) + qJ(6) * (-mrSges(7,2) * t246 - t261 * t290 + t384) + t290 * t254 - t366;
t403 = -mrSges(6,3) - mrSges(7,2);
t402 = Ifges(3,3) * t347;
t399 = t342 * t344;
t266 = mrSges(6,1) * t305 - mrSges(6,3) * t291;
t388 = -mrSges(6,1) * t290 - mrSges(6,2) * t291 - t261;
t210 = m(6) * t223 - mrSges(6,2) * t276 + t403 * t246 - t266 * t305 + t388 * t290 + t384;
t265 = -mrSges(6,2) * t305 - mrSges(6,3) * t290;
t211 = m(6) * t222 + mrSges(6,1) * t276 + t403 * t247 + t265 * t305 + t388 * t291 + t379;
t206 = t405 * t210 - t211 * t348;
t285 = -mrSges(5,1) * t308 + mrSges(5,2) * t309;
t294 = mrSges(5,1) * t314 - mrSges(5,3) * t309;
t202 = m(5) * t230 - mrSges(5,2) * t299 + mrSges(5,3) * t278 + t285 * t308 - t294 * t314 + t206;
t229 = -t349 * t234 + t236 * t352;
t225 = -pkin(4) * t299 - pkin(11) * t313 + t309 * t286 - t229;
t221 = -0.2e1 * qJD(6) * t291 + (t290 * t305 - t247) * qJ(6) + (t291 * t305 + t246) * pkin(5) + t225;
t215 = m(7) * t221 + mrSges(7,1) * t246 - t247 * mrSges(7,3) + t264 * t290 - t291 * t267;
t212 = -m(6) * t225 - t246 * mrSges(6,1) - mrSges(6,2) * t247 - t290 * t265 - t266 * t291 - t215;
t293 = -mrSges(5,2) * t314 + mrSges(5,3) * t308;
t208 = m(5) * t229 + mrSges(5,1) * t299 - mrSges(5,3) * t279 - t285 * t309 + t293 * t314 + t212;
t196 = t349 * t202 + t352 * t208;
t251 = Ifges(7,4) * t291 + Ifges(7,2) * t305 + Ifges(7,6) * t290;
t390 = -Ifges(6,5) * t291 + Ifges(6,6) * t290 - Ifges(6,3) * t305 - t251;
t300 = -mrSges(4,1) * t315 + mrSges(4,2) * t316;
t311 = mrSges(4,1) * t328 - mrSges(4,3) * t316;
t380 = t352 * t202 - t349 * t208;
t193 = m(4) * t238 - mrSges(4,2) * t325 + mrSges(4,3) * t302 + t300 * t315 - t311 * t328 + t380;
t310 = -mrSges(4,2) * t328 + mrSges(4,3) * t315;
t195 = m(4) * t257 - mrSges(4,1) * t302 + mrSges(4,2) * t303 - t310 * t315 + t311 * t316 + t196;
t205 = t348 * t210 + t405 * t211;
t358 = -m(5) * t233 + t278 * mrSges(5,1) - t279 * mrSges(5,2) + t308 * t293 - t309 * t294 - t205;
t199 = m(4) * t237 + t325 * mrSges(4,1) - t303 * mrSges(4,3) - t316 * t300 + t328 * t310 + t358;
t182 = t193 * t397 + t346 * t195 + t199 * t396;
t183 = t193 * t392 - t195 * t343 + t199 * t391;
t306 = -t331 * t342 + t374;
t378 = -mrSges(3,1) * t345 + mrSges(3,2) * t342;
t329 = t378 * t387;
t371 = -mrSges(3,2) * t347 + mrSges(3,3) * t395;
t334 = t371 * qJD(1);
t372 = mrSges(3,1) * t347 - mrSges(3,3) * t399;
t180 = m(3) * t306 + t372 * qJDD(1) + (-t329 * t399 + t334 * t347) * qJD(1) + t183;
t188 = t353 * t193 - t199 * t350;
t307 = -g(3) * t399 + t383;
t333 = t372 * qJD(1);
t187 = m(3) * t307 + t371 * qJDD(1) + (t329 * t395 - t333 * t347) * qJD(1) + t188;
t177 = -t180 * t342 + t345 * t187;
t317 = -t330 * t344 + t381;
t181 = m(3) * t317 + (t378 * qJDD(1) + (t333 * t342 - t334 * t345) * qJD(1)) * t344 + t182;
t172 = t180 * t393 - t181 * t344 + t187 * t398;
t377 = -mrSges(7,1) * t221 + mrSges(7,2) * t218;
t376 = Ifges(3,5) * t342 + Ifges(3,6) * t345;
t203 = -mrSges(6,1) * t225 + mrSges(6,3) * t223 - pkin(5) * t215 + (t253 + t254) * t305 + t390 * t291 + (Ifges(6,6) - Ifges(7,6)) * t276 + (Ifges(6,4) - Ifges(7,5)) * t247 + (-Ifges(6,2) - Ifges(7,3)) * t246 + t377;
t365 = mrSges(7,2) * t220 - mrSges(7,3) * t221 + Ifges(7,1) * t247 + Ifges(7,4) * t276 + Ifges(7,5) * t246 + t305 * t249;
t204 = mrSges(6,2) * t225 - mrSges(6,3) * t222 + Ifges(6,1) * t247 - Ifges(6,4) * t246 + Ifges(6,5) * t276 - qJ(6) * t215 - t252 * t305 + t390 * t290 + t365;
t272 = Ifges(5,5) * t309 + Ifges(5,6) * t308 + Ifges(5,3) * t314;
t273 = Ifges(5,4) * t309 + Ifges(5,2) * t308 + Ifges(5,6) * t314;
t184 = mrSges(5,2) * t233 - mrSges(5,3) * t229 + Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t299 - pkin(11) * t205 - t348 * t203 + t405 * t204 + t308 * t272 - t314 * t273;
t274 = Ifges(5,1) * t309 + Ifges(5,4) * t308 + Ifges(5,5) * t314;
t189 = -mrSges(5,1) * t233 + mrSges(5,3) * t230 + Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t299 - pkin(4) * t205 - t309 * t272 + t314 * t274 - t406;
t295 = Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t328;
t296 = Ifges(4,4) * t316 + Ifges(4,2) * t315 + Ifges(4,6) * t328;
t174 = mrSges(4,2) * t257 - mrSges(4,3) * t237 + Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t325 - pkin(10) * t196 + t184 * t352 - t189 * t349 + t295 * t315 - t296 * t328;
t297 = Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t328;
t356 = mrSges(5,1) * t229 - mrSges(5,2) * t230 + Ifges(5,5) * t279 + Ifges(5,6) * t278 + Ifges(5,3) * t299 + pkin(4) * t212 + pkin(11) * t206 + t405 * t203 + t348 * t204 + t309 * t273 - t308 * t274;
t178 = -mrSges(4,1) * t257 + mrSges(4,3) * t238 + Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t325 - pkin(3) * t196 - t316 * t295 + t328 * t297 - t356;
t367 = pkin(9) * t188 + t174 * t350 + t178 * t353;
t364 = Ifges(3,5) * t347 + (Ifges(3,1) * t342 + Ifges(3,4) * t345) * t344;
t363 = Ifges(3,6) * t347 + (Ifges(3,4) * t342 + Ifges(3,2) * t345) * t344;
t173 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * t325 + pkin(3) * t358 + pkin(10) * t380 + t349 * t184 + t352 * t189 + t316 * t296 - t315 * t297;
t321 = t363 * qJD(1);
t322 = t364 * qJD(1);
t164 = qJDD(1) * t402 + mrSges(3,1) * t306 - mrSges(3,2) * t307 + pkin(2) * t183 + t173 * t346 + t367 * t343 + (t376 * qJDD(1) + (t321 * t342 - t322 * t345) * qJD(1)) * t344;
t320 = (t376 * t344 + t402) * qJD(1);
t166 = -mrSges(3,1) * t317 + mrSges(3,3) * t307 - pkin(2) * t182 - t173 * t343 + (-t320 * t399 + t322 * t347) * qJD(1) + t367 * t346 + t363 * qJDD(1);
t168 = mrSges(3,2) * t317 - mrSges(3,3) * t306 + t174 * t353 - t178 * t350 + (t320 * t395 - t321 * t347) * qJD(1) + (-t182 * t343 - t183 * t346) * pkin(9) + t364 * qJDD(1);
t362 = mrSges(2,1) * t338 - mrSges(2,2) * t339 + Ifges(2,3) * qJDD(1) + pkin(1) * t172 + t347 * t164 + t166 * t395 + t168 * t399 + t177 * t400;
t175 = m(2) * t339 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t177;
t171 = t181 * t347 + (t180 * t345 + t187 * t342) * t344;
t169 = m(2) * t338 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t172;
t162 = -mrSges(2,2) * g(3) - mrSges(2,3) * t338 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t355 - t166 * t342 + t168 * t345 + (-t171 * t344 - t172 * t347) * qJ(2);
t161 = mrSges(2,1) * g(3) + mrSges(2,3) * t339 + Ifges(2,5) * t355 + Ifges(2,6) * qJDD(1) - pkin(1) * t171 - t344 * t164 + (qJ(2) * t177 + t166 * t345 + t168 * t342) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t162 - t351 * t161 - pkin(8) * (t169 * t354 + t175 * t351), t162, t168, t174, t184, t204, -t251 * t290 + t365; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t162 + t354 * t161 + pkin(8) * (-t351 * t169 + t175 * t354), t161, t166, t178, t189, t203, -t291 * t249 - t366; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t164, t173, t356, t406, Ifges(7,5) * t247 + Ifges(7,6) * t276 + Ifges(7,3) * t246 + t251 * t291 - t253 * t305 - t377;];
m_new  = t1;
