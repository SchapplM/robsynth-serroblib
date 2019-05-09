% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 16:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 15:42:38
% EndTime: 2019-05-08 15:49:45
% DurationCPUTime: 225.32s
% Computational Cost: add. (4039619->416), mult. (9936071->555), div. (0->0), fcn. (8461182->16), ass. (0->182)
t354 = cos(pkin(6));
t348 = qJD(1) * t354 + qJD(2);
t351 = sin(pkin(7));
t353 = cos(pkin(7));
t352 = sin(pkin(6));
t365 = cos(qJ(2));
t386 = qJD(1) * t365;
t382 = t352 * t386;
t332 = (t348 * t351 + t353 * t382) * pkin(10);
t359 = sin(qJ(2));
t388 = qJD(1) * t352;
t400 = pkin(10) * t351;
t336 = (-pkin(2) * t365 - t359 * t400) * t388;
t385 = qJD(1) * qJD(2);
t342 = (qJDD(1) * t359 + t365 * t385) * t352;
t347 = qJDD(1) * t354 + qJDD(2);
t360 = sin(qJ(1));
t366 = cos(qJ(1));
t345 = t360 * g(1) - g(2) * t366;
t367 = qJD(1) ^ 2;
t401 = pkin(9) * t352;
t339 = qJDD(1) * pkin(1) + t367 * t401 + t345;
t346 = -g(1) * t366 - g(2) * t360;
t340 = -pkin(1) * t367 + qJDD(1) * t401 + t346;
t390 = t354 * t365;
t379 = t339 * t390 - t359 * t340;
t387 = qJD(1) * t359;
t399 = pkin(10) * t353;
t288 = -t342 * t399 + t347 * pkin(2) + t348 * t332 + (-g(3) * t365 - t336 * t387) * t352 + t379;
t383 = t352 * t387;
t335 = pkin(2) * t348 - t383 * t399;
t343 = (qJDD(1) * t365 - t359 * t385) * t352;
t377 = t343 * t353 + t347 * t351;
t391 = t354 * t359;
t389 = t339 * t391 + t365 * t340;
t289 = -t348 * t335 + (-g(3) * t359 + t336 * t386) * t352 + t377 * pkin(10) + t389;
t398 = t354 * g(3);
t294 = -t342 * t400 - t343 * pkin(2) - t398 + (-t339 + (-t332 * t365 + t335 * t359) * qJD(1)) * t352;
t358 = sin(qJ(3));
t364 = cos(qJ(3));
t259 = -t358 * t289 + (t288 * t353 + t294 * t351) * t364;
t392 = t353 * t365;
t397 = t351 * t358;
t323 = t348 * t397 + (t358 * t392 + t359 * t364) * t388;
t306 = -t323 * qJD(3) - t358 * t342 + t364 * t377;
t396 = t351 * t364;
t322 = (-t358 * t359 + t364 * t392) * t388 + t348 * t396;
t395 = t352 * t359;
t394 = t352 * t365;
t393 = t353 * t358;
t260 = t288 * t393 + t364 * t289 + t294 * t397;
t309 = -pkin(3) * t322 - pkin(11) * t323;
t324 = -t343 * t351 + t347 * t353 + qJDD(3);
t333 = t348 * t353 - t351 * t382 + qJD(3);
t331 = t333 ^ 2;
t248 = -pkin(3) * t331 + pkin(11) * t324 + t309 * t322 + t260;
t267 = -t351 * t288 + t353 * t294;
t307 = t322 * qJD(3) + t364 * t342 + t358 * t377;
t250 = (-t322 * t333 - t307) * pkin(11) + (t323 * t333 - t306) * pkin(3) + t267;
t357 = sin(qJ(4));
t363 = cos(qJ(4));
t238 = t363 * t248 + t357 * t250;
t313 = -t357 * t323 + t333 * t363;
t314 = t323 * t363 + t333 * t357;
t291 = -pkin(4) * t313 - pkin(12) * t314;
t305 = qJDD(4) - t306;
t321 = qJD(4) - t322;
t320 = t321 ^ 2;
t233 = -pkin(4) * t320 + pkin(12) * t305 + t291 * t313 + t238;
t247 = -t324 * pkin(3) - t331 * pkin(11) + t323 * t309 - t259;
t275 = -t314 * qJD(4) - t357 * t307 + t324 * t363;
t276 = qJD(4) * t313 + t307 * t363 + t324 * t357;
t236 = (-t313 * t321 - t276) * pkin(12) + (t314 * t321 - t275) * pkin(4) + t247;
t356 = sin(qJ(5));
t362 = cos(qJ(5));
t228 = -t356 * t233 + t362 * t236;
t297 = -t314 * t356 + t321 * t362;
t258 = qJD(5) * t297 + t276 * t362 + t305 * t356;
t274 = qJDD(5) - t275;
t298 = t314 * t362 + t321 * t356;
t312 = qJD(5) - t313;
t226 = (t297 * t312 - t258) * pkin(13) + (t297 * t298 + t274) * pkin(5) + t228;
t229 = t362 * t233 + t356 * t236;
t257 = -qJD(5) * t298 - t276 * t356 + t305 * t362;
t280 = pkin(5) * t312 - pkin(13) * t298;
t296 = t297 ^ 2;
t227 = -pkin(5) * t296 + pkin(13) * t257 - t280 * t312 + t229;
t355 = sin(qJ(6));
t361 = cos(qJ(6));
t224 = t226 * t361 - t227 * t355;
t268 = t297 * t361 - t298 * t355;
t243 = qJD(6) * t268 + t257 * t355 + t258 * t361;
t269 = t297 * t355 + t298 * t361;
t255 = -mrSges(7,1) * t268 + mrSges(7,2) * t269;
t310 = qJD(6) + t312;
t261 = -mrSges(7,2) * t310 + mrSges(7,3) * t268;
t272 = qJDD(6) + t274;
t220 = m(7) * t224 + mrSges(7,1) * t272 - mrSges(7,3) * t243 - t255 * t269 + t261 * t310;
t225 = t226 * t355 + t227 * t361;
t242 = -qJD(6) * t269 + t257 * t361 - t258 * t355;
t262 = mrSges(7,1) * t310 - mrSges(7,3) * t269;
t221 = m(7) * t225 - mrSges(7,2) * t272 + mrSges(7,3) * t242 + t255 * t268 - t262 * t310;
t212 = t361 * t220 + t355 * t221;
t270 = -mrSges(6,1) * t297 + mrSges(6,2) * t298;
t278 = -mrSges(6,2) * t312 + mrSges(6,3) * t297;
t210 = m(6) * t228 + mrSges(6,1) * t274 - mrSges(6,3) * t258 - t270 * t298 + t278 * t312 + t212;
t279 = mrSges(6,1) * t312 - mrSges(6,3) * t298;
t380 = -t220 * t355 + t361 * t221;
t211 = m(6) * t229 - mrSges(6,2) * t274 + mrSges(6,3) * t257 + t270 * t297 - t279 * t312 + t380;
t208 = -t210 * t356 + t362 * t211;
t290 = -mrSges(5,1) * t313 + mrSges(5,2) * t314;
t300 = mrSges(5,1) * t321 - mrSges(5,3) * t314;
t206 = m(5) * t238 - mrSges(5,2) * t305 + mrSges(5,3) * t275 + t290 * t313 - t300 * t321 + t208;
t237 = -t357 * t248 + t250 * t363;
t232 = -pkin(4) * t305 - pkin(12) * t320 + t314 * t291 - t237;
t230 = -pkin(5) * t257 - pkin(13) * t296 + t280 * t298 + t232;
t373 = m(7) * t230 - t242 * mrSges(7,1) + mrSges(7,2) * t243 - t268 * t261 + t262 * t269;
t222 = -m(6) * t232 + t257 * mrSges(6,1) - mrSges(6,2) * t258 + t297 * t278 - t279 * t298 - t373;
t299 = -mrSges(5,2) * t321 + mrSges(5,3) * t313;
t216 = m(5) * t237 + mrSges(5,1) * t305 - mrSges(5,3) * t276 - t290 * t314 + t299 * t321 + t222;
t198 = t357 * t206 + t363 * t216;
t308 = -mrSges(4,1) * t322 + mrSges(4,2) * t323;
t316 = mrSges(4,1) * t333 - mrSges(4,3) * t323;
t381 = t363 * t206 - t216 * t357;
t195 = m(4) * t260 - mrSges(4,2) * t324 + mrSges(4,3) * t306 + t308 * t322 - t316 * t333 + t381;
t315 = -mrSges(4,2) * t333 + mrSges(4,3) * t322;
t197 = m(4) * t267 - mrSges(4,1) * t306 + mrSges(4,2) * t307 - t315 * t322 + t316 * t323 + t198;
t207 = t210 * t362 + t211 * t356;
t370 = -m(5) * t247 + t275 * mrSges(5,1) - mrSges(5,2) * t276 + t313 * t299 - t300 * t314 - t207;
t203 = m(4) * t259 + mrSges(4,1) * t324 - mrSges(4,3) * t307 - t308 * t323 + t315 * t333 + t370;
t184 = t195 * t397 + t353 * t197 + t203 * t396;
t185 = t353 * t364 * t203 + t195 * t393 - t197 * t351;
t317 = -g(3) * t394 + t379;
t338 = -mrSges(3,2) * t348 + mrSges(3,3) * t382;
t341 = (-mrSges(3,1) * t365 + mrSges(3,2) * t359) * t388;
t182 = m(3) * t317 + mrSges(3,1) * t347 - mrSges(3,3) * t342 + t338 * t348 - t341 * t383 + t185;
t190 = t364 * t195 - t203 * t358;
t318 = -g(3) * t395 + t389;
t337 = mrSges(3,1) * t348 - mrSges(3,3) * t383;
t189 = m(3) * t318 - mrSges(3,2) * t347 + mrSges(3,3) * t343 - t337 * t348 + t341 * t382 + t190;
t179 = -t182 * t359 + t365 * t189;
t328 = -t352 * t339 - t398;
t183 = m(3) * t328 - t343 * mrSges(3,1) + t342 * mrSges(3,2) + (t337 * t359 - t338 * t365) * t388 + t184;
t174 = t182 * t390 - t183 * t352 + t189 * t391;
t251 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t310;
t253 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t310;
t213 = -mrSges(7,1) * t230 + mrSges(7,3) * t225 + Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t272 - t251 * t269 + t253 * t310;
t252 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t310;
t214 = mrSges(7,2) * t230 - mrSges(7,3) * t224 + Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t272 + t251 * t268 - t252 * t310;
t263 = Ifges(6,5) * t298 + Ifges(6,6) * t297 + Ifges(6,3) * t312;
t265 = Ifges(6,1) * t298 + Ifges(6,4) * t297 + Ifges(6,5) * t312;
t199 = -mrSges(6,1) * t232 + mrSges(6,3) * t229 + Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t274 - pkin(5) * t373 + pkin(13) * t380 + t361 * t213 + t355 * t214 - t298 * t263 + t312 * t265;
t264 = Ifges(6,4) * t298 + Ifges(6,2) * t297 + Ifges(6,6) * t312;
t200 = mrSges(6,2) * t232 - mrSges(6,3) * t228 + Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t274 - pkin(13) * t212 - t213 * t355 + t214 * t361 + t263 * t297 - t264 * t312;
t281 = Ifges(5,5) * t314 + Ifges(5,6) * t313 + Ifges(5,3) * t321;
t282 = Ifges(5,4) * t314 + Ifges(5,2) * t313 + Ifges(5,6) * t321;
t186 = mrSges(5,2) * t247 - mrSges(5,3) * t237 + Ifges(5,1) * t276 + Ifges(5,4) * t275 + Ifges(5,5) * t305 - pkin(12) * t207 - t199 * t356 + t200 * t362 + t281 * t313 - t282 * t321;
t283 = Ifges(5,1) * t314 + Ifges(5,4) * t313 + Ifges(5,5) * t321;
t371 = -mrSges(7,1) * t224 + mrSges(7,2) * t225 - Ifges(7,5) * t243 - Ifges(7,6) * t242 - Ifges(7,3) * t272 - t269 * t252 + t268 * t253;
t368 = mrSges(6,1) * t228 - mrSges(6,2) * t229 + Ifges(6,5) * t258 + Ifges(6,6) * t257 + Ifges(6,3) * t274 + pkin(5) * t212 + t298 * t264 - t297 * t265 - t371;
t191 = -mrSges(5,1) * t247 + mrSges(5,3) * t238 + Ifges(5,4) * t276 + Ifges(5,2) * t275 + Ifges(5,6) * t305 - pkin(4) * t207 - t314 * t281 + t321 * t283 - t368;
t301 = Ifges(4,5) * t323 + Ifges(4,6) * t322 + Ifges(4,3) * t333;
t302 = Ifges(4,4) * t323 + Ifges(4,2) * t322 + Ifges(4,6) * t333;
t176 = mrSges(4,2) * t267 - mrSges(4,3) * t259 + Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t324 - pkin(11) * t198 + t186 * t363 - t191 * t357 + t301 * t322 - t302 * t333;
t303 = Ifges(4,1) * t323 + Ifges(4,4) * t322 + Ifges(4,5) * t333;
t369 = mrSges(5,1) * t237 - mrSges(5,2) * t238 + Ifges(5,5) * t276 + Ifges(5,6) * t275 + Ifges(5,3) * t305 + pkin(4) * t222 + pkin(12) * t208 + t362 * t199 + t356 * t200 + t314 * t282 - t313 * t283;
t180 = -mrSges(4,1) * t267 + mrSges(4,3) * t260 + Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t324 - pkin(3) * t198 - t323 * t301 + t333 * t303 - t369;
t374 = pkin(10) * t190 + t176 * t358 + t180 * t364;
t175 = mrSges(4,1) * t259 - mrSges(4,2) * t260 + Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t324 + pkin(3) * t370 + pkin(11) * t381 + t357 * t186 + t363 * t191 + t323 * t302 - t322 * t303;
t326 = Ifges(3,6) * t348 + (Ifges(3,4) * t359 + Ifges(3,2) * t365) * t388;
t327 = Ifges(3,5) * t348 + (Ifges(3,1) * t359 + Ifges(3,4) * t365) * t388;
t166 = mrSges(3,1) * t317 - mrSges(3,2) * t318 + Ifges(3,5) * t342 + Ifges(3,6) * t343 + Ifges(3,3) * t347 + pkin(2) * t185 + t353 * t175 + (t326 * t359 - t327 * t365) * t388 + t374 * t351;
t325 = Ifges(3,3) * t348 + (Ifges(3,5) * t359 + Ifges(3,6) * t365) * t388;
t168 = -mrSges(3,1) * t328 + mrSges(3,3) * t318 + Ifges(3,4) * t342 + Ifges(3,2) * t343 + Ifges(3,6) * t347 - pkin(2) * t184 - t351 * t175 - t325 * t383 + t348 * t327 + t353 * t374;
t170 = t325 * t382 + mrSges(3,2) * t328 - mrSges(3,3) * t317 + Ifges(3,1) * t342 + Ifges(3,4) * t343 + Ifges(3,5) * t347 + t364 * t176 - t358 * t180 - t348 * t326 + (-t184 * t351 - t185 * t353) * pkin(10);
t372 = mrSges(2,1) * t345 - mrSges(2,2) * t346 + Ifges(2,3) * qJDD(1) + pkin(1) * t174 + t354 * t166 + t168 * t394 + t170 * t395 + t179 * t401;
t177 = m(2) * t346 - mrSges(2,1) * t367 - qJDD(1) * mrSges(2,2) + t179;
t173 = t354 * t183 + (t182 * t365 + t189 * t359) * t352;
t171 = m(2) * t345 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t367 + t174;
t164 = -mrSges(2,2) * g(3) - mrSges(2,3) * t345 + Ifges(2,5) * qJDD(1) - t367 * Ifges(2,6) - t359 * t168 + t365 * t170 + (-t173 * t352 - t174 * t354) * pkin(9);
t163 = mrSges(2,1) * g(3) + mrSges(2,3) * t346 + t367 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t173 - t352 * t166 + (pkin(9) * t179 + t168 * t365 + t170 * t359) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t366 * t164 - t360 * t163 - pkin(8) * (t171 * t366 + t177 * t360), t164, t170, t176, t186, t200, t214; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t360 * t164 + t366 * t163 + pkin(8) * (-t171 * t360 + t177 * t366), t163, t168, t180, t191, t199, t213; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372, t372, t166, t175, t369, t368, -t371;];
m_new  = t1;
