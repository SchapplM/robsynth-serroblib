% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 02:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 02:08:11
% EndTime: 2019-05-08 02:13:57
% DurationCPUTime: 219.77s
% Computational Cost: add. (3909821->416), mult. (9694371->557), div. (0->0), fcn. (8221468->16), ass. (0->180)
t353 = cos(pkin(6));
t345 = t353 * qJD(1) + qJD(2);
t349 = sin(pkin(7));
t352 = cos(pkin(7));
t350 = sin(pkin(6));
t361 = cos(qJ(2));
t382 = qJD(1) * t361;
t378 = t350 * t382;
t329 = (t345 * t349 + t352 * t378) * pkin(10);
t357 = sin(qJ(2));
t384 = qJD(1) * t350;
t396 = pkin(10) * t349;
t333 = (-pkin(2) * t361 - t357 * t396) * t384;
t381 = qJD(1) * qJD(2);
t339 = (qJDD(1) * t357 + t361 * t381) * t350;
t344 = t353 * qJDD(1) + qJDD(2);
t358 = sin(qJ(1));
t362 = cos(qJ(1));
t342 = t358 * g(1) - t362 * g(2);
t363 = qJD(1) ^ 2;
t397 = pkin(9) * t350;
t336 = qJDD(1) * pkin(1) + t363 * t397 + t342;
t343 = -t362 * g(1) - t358 * g(2);
t337 = -t363 * pkin(1) + qJDD(1) * t397 + t343;
t386 = t353 * t361;
t375 = t336 * t386 - t357 * t337;
t383 = qJD(1) * t357;
t395 = pkin(10) * t352;
t285 = -t339 * t395 + t344 * pkin(2) + t345 * t329 + (-g(3) * t361 - t333 * t383) * t350 + t375;
t379 = t350 * t383;
t332 = t345 * pkin(2) - t379 * t395;
t340 = (qJDD(1) * t361 - t357 * t381) * t350;
t373 = t340 * t352 + t344 * t349;
t387 = t353 * t357;
t385 = t336 * t387 + t361 * t337;
t286 = -t345 * t332 + (-g(3) * t357 + t333 * t382) * t350 + t373 * pkin(10) + t385;
t394 = t353 * g(3);
t291 = -t339 * t396 - t340 * pkin(2) - t394 + (-t336 + (-t329 * t361 + t332 * t357) * qJD(1)) * t350;
t356 = sin(qJ(3));
t360 = cos(qJ(3));
t254 = -t356 * t286 + (t285 * t352 + t291 * t349) * t360;
t388 = t352 * t361;
t393 = t349 * t356;
t319 = t345 * t393 + (t356 * t388 + t357 * t360) * t384;
t304 = -t319 * qJD(3) - t356 * t339 + t360 * t373;
t392 = t349 * t360;
t318 = (-t356 * t357 + t360 * t388) * t384 + t345 * t392;
t398 = cos(qJ(4));
t391 = t350 * t357;
t390 = t350 * t361;
t389 = t352 * t356;
t255 = t285 * t389 + t360 * t286 + t291 * t393;
t307 = -t318 * pkin(3) - t319 * pkin(11);
t320 = -t349 * t340 + t352 * t344 + qJDD(3);
t330 = t352 * t345 - t349 * t378 + qJD(3);
t328 = t330 ^ 2;
t246 = -t328 * pkin(3) + t320 * pkin(11) + t318 * t307 + t255;
t265 = -t349 * t285 + t352 * t291;
t305 = t318 * qJD(3) + t360 * t339 + t356 * t373;
t248 = (-t318 * t330 - t305) * pkin(11) + (t319 * t330 - t304) * pkin(3) + t265;
t355 = sin(qJ(4));
t236 = t246 * t398 + t355 * t248;
t309 = t355 * t319 - t330 * t398;
t310 = t319 * t398 + t355 * t330;
t287 = t309 * pkin(4) - t310 * qJ(5);
t303 = qJDD(4) - t304;
t316 = qJD(4) - t318;
t315 = t316 ^ 2;
t231 = -t315 * pkin(4) + t303 * qJ(5) - t309 * t287 + t236;
t245 = -t320 * pkin(3) - t328 * pkin(11) + t319 * t307 - t254;
t272 = t310 * qJD(4) + t355 * t305 - t320 * t398;
t273 = -t309 * qJD(4) + t305 * t398 + t355 * t320;
t234 = (t309 * t316 - t273) * qJ(5) + (t310 * t316 + t272) * pkin(4) + t245;
t348 = sin(pkin(13));
t351 = cos(pkin(13));
t297 = t351 * t310 + t348 * t316;
t226 = -0.2e1 * qJD(5) * t297 - t348 * t231 + t351 * t234;
t258 = t351 * t273 + t348 * t303;
t296 = -t348 * t310 + t351 * t316;
t224 = (t296 * t309 - t258) * pkin(12) + (t296 * t297 + t272) * pkin(5) + t226;
t227 = 0.2e1 * qJD(5) * t296 + t351 * t231 + t348 * t234;
t257 = -t348 * t273 + t351 * t303;
t277 = t309 * pkin(5) - t297 * pkin(12);
t295 = t296 ^ 2;
t225 = -t295 * pkin(5) + t257 * pkin(12) - t309 * t277 + t227;
t354 = sin(qJ(6));
t359 = cos(qJ(6));
t222 = t359 * t224 - t354 * t225;
t266 = t359 * t296 - t354 * t297;
t241 = t266 * qJD(6) + t354 * t257 + t359 * t258;
t267 = t354 * t296 + t359 * t297;
t253 = -t266 * mrSges(7,1) + t267 * mrSges(7,2);
t308 = qJD(6) + t309;
t259 = -t308 * mrSges(7,2) + t266 * mrSges(7,3);
t271 = qJDD(6) + t272;
t218 = m(7) * t222 + t271 * mrSges(7,1) - t241 * mrSges(7,3) - t267 * t253 + t308 * t259;
t223 = t354 * t224 + t359 * t225;
t240 = -t267 * qJD(6) + t359 * t257 - t354 * t258;
t260 = t308 * mrSges(7,1) - t267 * mrSges(7,3);
t219 = m(7) * t223 - t271 * mrSges(7,2) + t240 * mrSges(7,3) + t266 * t253 - t308 * t260;
t210 = t359 * t218 + t354 * t219;
t268 = -t296 * mrSges(6,1) + t297 * mrSges(6,2);
t275 = -t309 * mrSges(6,2) + t296 * mrSges(6,3);
t208 = m(6) * t226 + t272 * mrSges(6,1) - t258 * mrSges(6,3) - t297 * t268 + t309 * t275 + t210;
t276 = t309 * mrSges(6,1) - t297 * mrSges(6,3);
t376 = -t354 * t218 + t359 * t219;
t209 = m(6) * t227 - t272 * mrSges(6,2) + t257 * mrSges(6,3) + t296 * t268 - t309 * t276 + t376;
t206 = -t348 * t208 + t351 * t209;
t288 = t309 * mrSges(5,1) + t310 * mrSges(5,2);
t299 = t316 * mrSges(5,1) - t310 * mrSges(5,3);
t204 = m(5) * t236 - t303 * mrSges(5,2) - t272 * mrSges(5,3) - t309 * t288 - t316 * t299 + t206;
t235 = -t355 * t246 + t248 * t398;
t230 = -t303 * pkin(4) - t315 * qJ(5) + t310 * t287 + qJDD(5) - t235;
t228 = -t257 * pkin(5) - t295 * pkin(12) + t297 * t277 + t230;
t369 = m(7) * t228 - t240 * mrSges(7,1) + t241 * mrSges(7,2) - t266 * t259 + t267 * t260;
t220 = -m(6) * t230 + t257 * mrSges(6,1) - t258 * mrSges(6,2) + t296 * t275 - t297 * t276 - t369;
t298 = -t316 * mrSges(5,2) - t309 * mrSges(5,3);
t214 = m(5) * t235 + t303 * mrSges(5,1) - t273 * mrSges(5,3) - t310 * t288 + t316 * t298 + t220;
t196 = t355 * t204 + t214 * t398;
t306 = -t318 * mrSges(4,1) + t319 * mrSges(4,2);
t312 = t330 * mrSges(4,1) - t319 * mrSges(4,3);
t377 = t204 * t398 - t355 * t214;
t193 = m(4) * t255 - t320 * mrSges(4,2) + t304 * mrSges(4,3) + t318 * t306 - t330 * t312 + t377;
t311 = -t330 * mrSges(4,2) + t318 * mrSges(4,3);
t195 = m(4) * t265 - t304 * mrSges(4,1) + t305 * mrSges(4,2) - t318 * t311 + t319 * t312 + t196;
t205 = t351 * t208 + t348 * t209;
t366 = -m(5) * t245 - t272 * mrSges(5,1) - t273 * mrSges(5,2) - t309 * t298 - t310 * t299 - t205;
t201 = m(4) * t254 + t320 * mrSges(4,1) - t305 * mrSges(4,3) - t319 * t306 + t330 * t311 + t366;
t182 = t193 * t393 + t352 * t195 + t201 * t392;
t183 = t352 * t360 * t201 + t193 * t389 - t349 * t195;
t313 = -g(3) * t390 + t375;
t335 = -t345 * mrSges(3,2) + mrSges(3,3) * t378;
t338 = (-mrSges(3,1) * t361 + mrSges(3,2) * t357) * t384;
t180 = m(3) * t313 + t344 * mrSges(3,1) - t339 * mrSges(3,3) + t345 * t335 - t338 * t379 + t183;
t188 = t360 * t193 - t356 * t201;
t314 = -g(3) * t391 + t385;
t334 = t345 * mrSges(3,1) - mrSges(3,3) * t379;
t187 = m(3) * t314 - t344 * mrSges(3,2) + t340 * mrSges(3,3) - t345 * t334 + t338 * t378 + t188;
t177 = -t357 * t180 + t361 * t187;
t324 = -t350 * t336 - t394;
t181 = m(3) * t324 - t340 * mrSges(3,1) + t339 * mrSges(3,2) + (t334 * t357 - t335 * t361) * t384 + t182;
t172 = t180 * t386 - t350 * t181 + t187 * t387;
t249 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t308;
t251 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t308;
t211 = -mrSges(7,1) * t228 + mrSges(7,3) * t223 + Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t271 - t267 * t249 + t308 * t251;
t250 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t308;
t212 = mrSges(7,2) * t228 - mrSges(7,3) * t222 + Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t271 + t266 * t249 - t308 * t250;
t261 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t309;
t263 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t309;
t197 = -mrSges(6,1) * t230 + mrSges(6,3) * t227 + Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t272 - pkin(5) * t369 + pkin(12) * t376 + t359 * t211 + t354 * t212 - t297 * t261 + t309 * t263;
t262 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t309;
t198 = mrSges(6,2) * t230 - mrSges(6,3) * t226 + Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t272 - pkin(12) * t210 - t354 * t211 + t359 * t212 + t296 * t261 - t309 * t262;
t278 = Ifges(5,5) * t310 - Ifges(5,6) * t309 + Ifges(5,3) * t316;
t279 = Ifges(5,4) * t310 - Ifges(5,2) * t309 + Ifges(5,6) * t316;
t184 = mrSges(5,2) * t245 - mrSges(5,3) * t235 + Ifges(5,1) * t273 - Ifges(5,4) * t272 + Ifges(5,5) * t303 - qJ(5) * t205 - t348 * t197 + t351 * t198 - t309 * t278 - t316 * t279;
t280 = Ifges(5,1) * t310 - Ifges(5,4) * t309 + Ifges(5,5) * t316;
t367 = -mrSges(7,1) * t222 + mrSges(7,2) * t223 - Ifges(7,5) * t241 - Ifges(7,6) * t240 - Ifges(7,3) * t271 - t267 * t250 + t266 * t251;
t365 = -mrSges(6,1) * t226 + mrSges(6,2) * t227 - Ifges(6,5) * t258 - Ifges(6,6) * t257 - pkin(5) * t210 - t297 * t262 + t296 * t263 + t367;
t189 = (-Ifges(5,2) - Ifges(6,3)) * t272 + t365 + t316 * t280 - t310 * t278 + Ifges(5,6) * t303 + Ifges(5,4) * t273 - mrSges(5,1) * t245 + mrSges(5,3) * t236 - pkin(4) * t205;
t300 = Ifges(4,5) * t319 + Ifges(4,6) * t318 + Ifges(4,3) * t330;
t301 = Ifges(4,4) * t319 + Ifges(4,2) * t318 + Ifges(4,6) * t330;
t174 = mrSges(4,2) * t265 - mrSges(4,3) * t254 + Ifges(4,1) * t305 + Ifges(4,4) * t304 + Ifges(4,5) * t320 - pkin(11) * t196 + t184 * t398 - t355 * t189 + t318 * t300 - t330 * t301;
t302 = Ifges(4,1) * t319 + Ifges(4,4) * t318 + Ifges(4,5) * t330;
t364 = mrSges(5,1) * t235 - mrSges(5,2) * t236 + Ifges(5,5) * t273 - Ifges(5,6) * t272 + Ifges(5,3) * t303 + pkin(4) * t220 + qJ(5) * t206 + t351 * t197 + t348 * t198 + t310 * t279 + t309 * t280;
t178 = -mrSges(4,1) * t265 + mrSges(4,3) * t255 + Ifges(4,4) * t305 + Ifges(4,2) * t304 + Ifges(4,6) * t320 - pkin(3) * t196 - t319 * t300 + t330 * t302 - t364;
t370 = pkin(10) * t188 + t174 * t356 + t178 * t360;
t173 = mrSges(4,1) * t254 - mrSges(4,2) * t255 + Ifges(4,5) * t305 + Ifges(4,6) * t304 + Ifges(4,3) * t320 + pkin(3) * t366 + pkin(11) * t377 + t355 * t184 + t189 * t398 + t319 * t301 - t318 * t302;
t322 = Ifges(3,6) * t345 + (Ifges(3,4) * t357 + Ifges(3,2) * t361) * t384;
t323 = Ifges(3,5) * t345 + (Ifges(3,1) * t357 + Ifges(3,4) * t361) * t384;
t164 = mrSges(3,1) * t313 - mrSges(3,2) * t314 + Ifges(3,5) * t339 + Ifges(3,6) * t340 + Ifges(3,3) * t344 + pkin(2) * t183 + t352 * t173 + (t322 * t357 - t323 * t361) * t384 + t370 * t349;
t321 = Ifges(3,3) * t345 + (Ifges(3,5) * t357 + Ifges(3,6) * t361) * t384;
t166 = -mrSges(3,1) * t324 + mrSges(3,3) * t314 + Ifges(3,4) * t339 + Ifges(3,2) * t340 + Ifges(3,6) * t344 - pkin(2) * t182 - t349 * t173 - t321 * t379 + t345 * t323 + t352 * t370;
t168 = t321 * t378 + mrSges(3,2) * t324 - mrSges(3,3) * t313 + Ifges(3,1) * t339 + Ifges(3,4) * t340 + Ifges(3,5) * t344 + t360 * t174 - t356 * t178 - t345 * t322 + (-t182 * t349 - t183 * t352) * pkin(10);
t368 = mrSges(2,1) * t342 - mrSges(2,2) * t343 + Ifges(2,3) * qJDD(1) + pkin(1) * t172 + t353 * t164 + t166 * t390 + t168 * t391 + t177 * t397;
t175 = m(2) * t343 - t363 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t177;
t171 = t353 * t181 + (t180 * t361 + t187 * t357) * t350;
t169 = m(2) * t342 + qJDD(1) * mrSges(2,1) - t363 * mrSges(2,2) + t172;
t162 = -mrSges(2,2) * g(3) - mrSges(2,3) * t342 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - t357 * t166 + t361 * t168 + (-t171 * t350 - t172 * t353) * pkin(9);
t161 = mrSges(2,1) * g(3) + mrSges(2,3) * t343 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t171 - t350 * t164 + (pkin(9) * t177 + t166 * t361 + t168 * t357) * t353;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t362 * t162 - t358 * t161 - pkin(8) * (t362 * t169 + t358 * t175), t162, t168, t174, t184, t198, t212; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t358 * t162 + t362 * t161 + pkin(8) * (-t358 * t169 + t362 * t175), t161, t166, t178, t189, t197, t211; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t368, t368, t164, t173, t364, Ifges(6,3) * t272 - t365, -t367;];
m_new  = t1;
