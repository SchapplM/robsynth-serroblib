% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 21:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:17:03
% EndTime: 2019-05-06 21:18:51
% DurationCPUTime: 74.90s
% Computational Cost: add. (1264986->396), mult. (3308279->517), div. (0->0), fcn. (2654355->14), ass. (0->165)
t387 = -2 * qJD(3);
t346 = sin(pkin(12));
t348 = cos(pkin(12));
t353 = sin(qJ(2));
t358 = cos(qJ(2));
t347 = sin(pkin(6));
t380 = qJD(1) * t347;
t324 = (t346 * t353 - t348 * t358) * t380;
t378 = qJD(1) * qJD(2);
t333 = (qJDD(1) * t353 + t358 * t378) * t347;
t349 = cos(pkin(6));
t340 = qJDD(1) * t349 + qJDD(2);
t341 = qJD(1) * t349 + qJD(2);
t354 = sin(qJ(1));
t359 = cos(qJ(1));
t337 = t354 * g(1) - g(2) * t359;
t360 = qJD(1) ^ 2;
t386 = pkin(8) * t347;
t330 = qJDD(1) * pkin(1) + t360 * t386 + t337;
t338 = -g(1) * t359 - g(2) * t354;
t331 = -pkin(1) * t360 + qJDD(1) * t386 + t338;
t381 = t349 * t358;
t370 = t330 * t381 - t353 * t331;
t385 = t347 ^ 2 * t360;
t264 = t340 * pkin(2) - t333 * qJ(3) + (pkin(2) * t353 * t385 + (qJ(3) * qJD(1) * t341 - g(3)) * t347) * t358 + t370;
t382 = t349 * t353;
t384 = t347 * t353;
t295 = -g(3) * t384 + t330 * t382 + t358 * t331;
t376 = t353 * t380;
t327 = pkin(2) * t341 - qJ(3) * t376;
t334 = (qJDD(1) * t358 - t353 * t378) * t347;
t377 = t358 ^ 2 * t385;
t268 = -pkin(2) * t377 + qJ(3) * t334 - t327 * t341 + t295;
t325 = (t346 * t358 + t348 * t353) * t380;
t244 = t348 * t264 - t346 * t268 + t325 * t387;
t383 = t347 * t358;
t245 = t346 * t264 + t348 * t268 + t324 * t387;
t296 = mrSges(4,1) * t324 + mrSges(4,2) * t325;
t303 = -t333 * t346 + t334 * t348;
t311 = mrSges(4,1) * t341 - mrSges(4,3) * t325;
t297 = pkin(3) * t324 - pkin(9) * t325;
t339 = t341 ^ 2;
t237 = -pkin(3) * t339 + pkin(9) * t340 - t297 * t324 + t245;
t315 = -t349 * g(3) - t347 * t330;
t280 = -t334 * pkin(2) - qJ(3) * t377 + t327 * t376 + qJDD(3) + t315;
t304 = t333 * t348 + t334 * t346;
t247 = (t324 * t341 - t304) * pkin(9) + (t325 * t341 - t303) * pkin(3) + t280;
t352 = sin(qJ(4));
t357 = cos(qJ(4));
t232 = t357 * t237 + t352 * t247;
t308 = -t352 * t325 + t341 * t357;
t309 = t325 * t357 + t341 * t352;
t283 = -pkin(4) * t308 - pkin(10) * t309;
t302 = qJDD(4) - t303;
t323 = qJD(4) + t324;
t322 = t323 ^ 2;
t222 = -pkin(4) * t322 + pkin(10) * t302 + t283 * t308 + t232;
t236 = -t340 * pkin(3) - t339 * pkin(9) + t325 * t297 - t244;
t277 = -t309 * qJD(4) - t352 * t304 + t340 * t357;
t278 = qJD(4) * t308 + t304 * t357 + t340 * t352;
t225 = (-t308 * t323 - t278) * pkin(10) + (t309 * t323 - t277) * pkin(4) + t236;
t351 = sin(qJ(5));
t356 = cos(qJ(5));
t217 = -t351 * t222 + t356 * t225;
t286 = -t309 * t351 + t323 * t356;
t250 = qJD(5) * t286 + t278 * t356 + t302 * t351;
t276 = qJDD(5) - t277;
t287 = t309 * t356 + t323 * t351;
t307 = qJD(5) - t308;
t215 = (t286 * t307 - t250) * pkin(11) + (t286 * t287 + t276) * pkin(5) + t217;
t218 = t356 * t222 + t351 * t225;
t249 = -qJD(5) * t287 - t278 * t351 + t302 * t356;
t267 = pkin(5) * t307 - pkin(11) * t287;
t285 = t286 ^ 2;
t216 = -pkin(5) * t285 + pkin(11) * t249 - t267 * t307 + t218;
t350 = sin(qJ(6));
t355 = cos(qJ(6));
t213 = t215 * t355 - t216 * t350;
t257 = t286 * t355 - t287 * t350;
t230 = qJD(6) * t257 + t249 * t350 + t250 * t355;
t258 = t286 * t350 + t287 * t355;
t242 = -mrSges(7,1) * t257 + mrSges(7,2) * t258;
t305 = qJD(6) + t307;
t251 = -mrSges(7,2) * t305 + mrSges(7,3) * t257;
t274 = qJDD(6) + t276;
t209 = m(7) * t213 + mrSges(7,1) * t274 - mrSges(7,3) * t230 - t242 * t258 + t251 * t305;
t214 = t215 * t350 + t216 * t355;
t229 = -qJD(6) * t258 + t249 * t355 - t250 * t350;
t252 = mrSges(7,1) * t305 - mrSges(7,3) * t258;
t210 = m(7) * t214 - mrSges(7,2) * t274 + mrSges(7,3) * t229 + t242 * t257 - t252 * t305;
t201 = t355 * t209 + t350 * t210;
t259 = -mrSges(6,1) * t286 + mrSges(6,2) * t287;
t265 = -mrSges(6,2) * t307 + mrSges(6,3) * t286;
t199 = m(6) * t217 + mrSges(6,1) * t276 - mrSges(6,3) * t250 - t259 * t287 + t265 * t307 + t201;
t266 = mrSges(6,1) * t307 - mrSges(6,3) * t287;
t372 = -t209 * t350 + t355 * t210;
t200 = m(6) * t218 - mrSges(6,2) * t276 + mrSges(6,3) * t249 + t259 * t286 - t266 * t307 + t372;
t197 = -t199 * t351 + t356 * t200;
t282 = -mrSges(5,1) * t308 + mrSges(5,2) * t309;
t289 = mrSges(5,1) * t323 - mrSges(5,3) * t309;
t194 = m(5) * t232 - mrSges(5,2) * t302 + mrSges(5,3) * t277 + t282 * t308 - t289 * t323 + t197;
t231 = -t352 * t237 + t247 * t357;
t221 = -pkin(4) * t302 - pkin(10) * t322 + t309 * t283 - t231;
t219 = -pkin(5) * t249 - pkin(11) * t285 + t267 * t287 + t221;
t368 = m(7) * t219 - t229 * mrSges(7,1) + mrSges(7,2) * t230 - t257 * t251 + t252 * t258;
t211 = -m(6) * t221 + t249 * mrSges(6,1) - mrSges(6,2) * t250 + t286 * t265 - t266 * t287 - t368;
t288 = -mrSges(5,2) * t323 + mrSges(5,3) * t308;
t205 = m(5) * t231 + mrSges(5,1) * t302 - mrSges(5,3) * t278 - t282 * t309 + t288 * t323 + t211;
t373 = t357 * t194 - t205 * t352;
t184 = m(4) * t245 - mrSges(4,2) * t340 + mrSges(4,3) * t303 - t296 * t324 - t311 * t341 + t373;
t310 = -mrSges(4,2) * t341 - mrSges(4,3) * t324;
t196 = t199 * t356 + t200 * t351;
t363 = -m(5) * t236 + t277 * mrSges(5,1) - mrSges(5,2) * t278 + t308 * t288 - t289 * t309 - t196;
t191 = m(4) * t244 + mrSges(4,1) * t340 - mrSges(4,3) * t304 - t296 * t325 + t310 * t341 + t363;
t179 = t346 * t184 + t348 * t191;
t187 = t352 * t194 + t357 * t205;
t375 = t358 * t380;
t294 = -g(3) * t383 + t370;
t329 = -mrSges(3,2) * t341 + mrSges(3,3) * t375;
t332 = (-mrSges(3,1) * t358 + mrSges(3,2) * t353) * t380;
t177 = m(3) * t294 + mrSges(3,1) * t340 - mrSges(3,3) * t333 + t329 * t341 - t332 * t376 + t179;
t328 = mrSges(3,1) * t341 - mrSges(3,3) * t376;
t374 = t348 * t184 - t191 * t346;
t178 = m(3) * t295 - mrSges(3,2) * t340 + mrSges(3,3) * t334 - t328 * t341 + t332 * t375 + t374;
t171 = -t177 * t353 + t358 * t178;
t366 = m(4) * t280 - t303 * mrSges(4,1) + t304 * mrSges(4,2) + t324 * t310 + t325 * t311 + t187;
t185 = m(3) * t315 - t334 * mrSges(3,1) + t333 * mrSges(3,2) + (t328 * t353 - t329 * t358) * t380 + t366;
t167 = t177 * t381 + t178 * t382 - t185 * t347;
t238 = Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t305;
t240 = Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t305;
t202 = -mrSges(7,1) * t219 + mrSges(7,3) * t214 + Ifges(7,4) * t230 + Ifges(7,2) * t229 + Ifges(7,6) * t274 - t238 * t258 + t240 * t305;
t239 = Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t305;
t203 = mrSges(7,2) * t219 - mrSges(7,3) * t213 + Ifges(7,1) * t230 + Ifges(7,4) * t229 + Ifges(7,5) * t274 + t238 * t257 - t239 * t305;
t253 = Ifges(6,5) * t287 + Ifges(6,6) * t286 + Ifges(6,3) * t307;
t255 = Ifges(6,1) * t287 + Ifges(6,4) * t286 + Ifges(6,5) * t307;
t188 = -mrSges(6,1) * t221 + mrSges(6,3) * t218 + Ifges(6,4) * t250 + Ifges(6,2) * t249 + Ifges(6,6) * t276 - pkin(5) * t368 + pkin(11) * t372 + t355 * t202 + t350 * t203 - t287 * t253 + t307 * t255;
t254 = Ifges(6,4) * t287 + Ifges(6,2) * t286 + Ifges(6,6) * t307;
t189 = mrSges(6,2) * t221 - mrSges(6,3) * t217 + Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t276 - pkin(11) * t201 - t202 * t350 + t203 * t355 + t253 * t286 - t254 * t307;
t269 = Ifges(5,5) * t309 + Ifges(5,6) * t308 + Ifges(5,3) * t323;
t270 = Ifges(5,4) * t309 + Ifges(5,2) * t308 + Ifges(5,6) * t323;
t173 = mrSges(5,2) * t236 - mrSges(5,3) * t231 + Ifges(5,1) * t278 + Ifges(5,4) * t277 + Ifges(5,5) * t302 - pkin(10) * t196 - t188 * t351 + t189 * t356 + t269 * t308 - t270 * t323;
t271 = Ifges(5,1) * t309 + Ifges(5,4) * t308 + Ifges(5,5) * t323;
t365 = -mrSges(7,1) * t213 + mrSges(7,2) * t214 - Ifges(7,5) * t230 - Ifges(7,6) * t229 - Ifges(7,3) * t274 - t258 * t239 + t257 * t240;
t361 = mrSges(6,1) * t217 - mrSges(6,2) * t218 + Ifges(6,5) * t250 + Ifges(6,6) * t249 + Ifges(6,3) * t276 + pkin(5) * t201 + t287 * t254 - t286 * t255 - t365;
t181 = -mrSges(5,1) * t236 + mrSges(5,3) * t232 + Ifges(5,4) * t278 + Ifges(5,2) * t277 + Ifges(5,6) * t302 - pkin(4) * t196 - t309 * t269 + t323 * t271 - t361;
t290 = Ifges(4,5) * t325 - Ifges(4,6) * t324 + Ifges(4,3) * t341;
t291 = Ifges(4,4) * t325 - Ifges(4,2) * t324 + Ifges(4,6) * t341;
t163 = mrSges(4,2) * t280 - mrSges(4,3) * t244 + Ifges(4,1) * t304 + Ifges(4,4) * t303 + Ifges(4,5) * t340 - pkin(9) * t187 + t173 * t357 - t181 * t352 - t290 * t324 - t291 * t341;
t292 = Ifges(4,1) * t325 - Ifges(4,4) * t324 + Ifges(4,5) * t341;
t362 = mrSges(5,1) * t231 - mrSges(5,2) * t232 + Ifges(5,5) * t278 + Ifges(5,6) * t277 + Ifges(5,3) * t302 + pkin(4) * t211 + pkin(10) * t197 + t356 * t188 + t351 * t189 + t309 * t270 - t308 * t271;
t168 = -mrSges(4,1) * t280 + mrSges(4,3) * t245 + Ifges(4,4) * t304 + Ifges(4,2) * t303 + Ifges(4,6) * t340 - pkin(3) * t187 - t325 * t290 + t341 * t292 - t362;
t312 = Ifges(3,3) * t341 + (Ifges(3,5) * t353 + Ifges(3,6) * t358) * t380;
t314 = Ifges(3,5) * t341 + (Ifges(3,1) * t353 + Ifges(3,4) * t358) * t380;
t158 = -mrSges(3,1) * t315 + mrSges(3,3) * t295 + Ifges(3,4) * t333 + Ifges(3,2) * t334 + Ifges(3,6) * t340 - pkin(2) * t366 + qJ(3) * t374 + t346 * t163 + t348 * t168 - t312 * t376 + t341 * t314;
t313 = Ifges(3,6) * t341 + (Ifges(3,4) * t353 + Ifges(3,2) * t358) * t380;
t160 = mrSges(3,2) * t315 - mrSges(3,3) * t294 + Ifges(3,1) * t333 + Ifges(3,4) * t334 + Ifges(3,5) * t340 - qJ(3) * t179 + t163 * t348 - t168 * t346 + t312 * t375 - t313 * t341;
t364 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t304 + Ifges(4,6) * t303 + Ifges(4,3) * t340 + pkin(3) * t363 + pkin(9) * t373 + t352 * t173 + t357 * t181 + t325 * t291 + t324 * t292;
t162 = pkin(2) * t179 + (t313 * t353 - t314 * t358) * t380 + Ifges(3,3) * t340 + Ifges(3,5) * t333 + Ifges(3,6) * t334 + mrSges(3,1) * t294 - mrSges(3,2) * t295 + t364;
t367 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t167 + t158 * t383 + t160 * t384 + t349 * t162 + t171 * t386;
t169 = m(2) * t338 - mrSges(2,1) * t360 - qJDD(1) * mrSges(2,2) + t171;
t166 = t349 * t185 + (t177 * t358 + t178 * t353) * t347;
t164 = m(2) * t337 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t360 + t167;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - t360 * Ifges(2,6) - t353 * t158 + t358 * t160 + (-t166 * t347 - t167 * t349) * pkin(8);
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + t360 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t347 * t162 + (pkin(8) * t171 + t158 * t358 + t160 * t353) * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t359 * t156 - t354 * t155 - pkin(7) * (t164 * t359 + t169 * t354), t156, t160, t163, t173, t189, t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t354 * t156 + t359 * t155 + pkin(7) * (-t164 * t354 + t169 * t359), t155, t158, t168, t181, t188, t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t367, t367, t162, t364, t362, t361, -t365;];
m_new  = t1;
