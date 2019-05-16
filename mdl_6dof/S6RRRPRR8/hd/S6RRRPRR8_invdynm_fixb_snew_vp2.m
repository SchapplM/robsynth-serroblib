% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 12:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:13:51
% EndTime: 2019-05-07 12:16:33
% DurationCPUTime: 83.78s
% Computational Cost: add. (1521440->398), mult. (3350810->518), div. (0->0), fcn. (2727061->14), ass. (0->165)
t345 = sin(pkin(6));
t351 = sin(qJ(2));
t356 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t332 = (-qJDD(1) * t356 + t351 * t376) * t345;
t379 = qJD(1) * t345;
t330 = (-pkin(2) * t356 - pkin(9) * t351) * t379;
t347 = cos(pkin(6));
t340 = t347 * qJD(1) + qJD(2);
t338 = t340 ^ 2;
t339 = t347 * qJDD(1) + qJDD(2);
t378 = qJD(1) * t356;
t352 = sin(qJ(1));
t357 = cos(qJ(1));
t336 = t352 * g(1) - t357 * g(2);
t358 = qJD(1) ^ 2;
t386 = pkin(8) * t345;
t327 = qJDD(1) * pkin(1) + t358 * t386 + t336;
t337 = -t357 * g(1) - t352 * g(2);
t328 = -t358 * pkin(1) + qJDD(1) * t386 + t337;
t382 = t347 * t351;
t380 = t327 * t382 + t356 * t328;
t283 = -t338 * pkin(2) + t339 * pkin(9) + (-g(3) * t351 + t330 * t378) * t345 + t380;
t331 = (qJDD(1) * t351 + t356 * t376) * t345;
t385 = t347 * g(3);
t284 = t332 * pkin(2) - t331 * pkin(9) - t385 + (-t327 + (pkin(2) * t351 - pkin(9) * t356) * t340 * qJD(1)) * t345;
t350 = sin(qJ(3));
t355 = cos(qJ(3));
t251 = -t350 * t283 + t355 * t284;
t375 = t351 * t379;
t319 = t355 * t340 - t350 * t375;
t298 = t319 * qJD(3) + t355 * t331 + t350 * t339;
t320 = t350 * t340 + t355 * t375;
t324 = qJDD(3) + t332;
t374 = t345 * t378;
t335 = qJD(3) - t374;
t236 = (t319 * t335 - t298) * qJ(4) + (t319 * t320 + t324) * pkin(3) + t251;
t252 = t355 * t283 + t350 * t284;
t297 = -t320 * qJD(3) - t350 * t331 + t355 * t339;
t309 = t335 * pkin(3) - t320 * qJ(4);
t318 = t319 ^ 2;
t243 = -t318 * pkin(3) + t297 * qJ(4) - t335 * t309 + t252;
t344 = sin(pkin(12));
t346 = cos(pkin(12));
t306 = t344 * t319 + t346 * t320;
t224 = -0.2e1 * qJD(4) * t306 + t346 * t236 - t344 * t243;
t384 = t345 * t351;
t383 = t345 * t356;
t381 = t347 * t356;
t305 = t346 * t319 - t344 * t320;
t225 = 0.2e1 * qJD(4) * t305 + t344 * t236 + t346 * t243;
t271 = t346 * t297 - t344 * t298;
t277 = -t305 * mrSges(5,1) + t306 * mrSges(5,2);
t289 = t335 * mrSges(5,1) - t306 * mrSges(5,3);
t278 = -t305 * pkin(4) - t306 * pkin(10);
t334 = t335 ^ 2;
t222 = -t334 * pkin(4) + t324 * pkin(10) + t305 * t278 + t225;
t300 = -g(3) * t383 + t327 * t381 - t351 * t328;
t282 = -t339 * pkin(2) - t338 * pkin(9) + t330 * t375 - t300;
t245 = -t297 * pkin(3) - t318 * qJ(4) + t320 * t309 + qJDD(4) + t282;
t272 = t344 * t297 + t346 * t298;
t228 = (-t305 * t335 - t272) * pkin(10) + (t306 * t335 - t271) * pkin(4) + t245;
t349 = sin(qJ(5));
t354 = cos(qJ(5));
t217 = -t349 * t222 + t354 * t228;
t286 = -t349 * t306 + t354 * t335;
t248 = t286 * qJD(5) + t354 * t272 + t349 * t324;
t270 = qJDD(5) - t271;
t287 = t354 * t306 + t349 * t335;
t304 = qJD(5) - t305;
t215 = (t286 * t304 - t248) * pkin(11) + (t286 * t287 + t270) * pkin(5) + t217;
t218 = t354 * t222 + t349 * t228;
t247 = -t287 * qJD(5) - t349 * t272 + t354 * t324;
t265 = t304 * pkin(5) - t287 * pkin(11);
t285 = t286 ^ 2;
t216 = -t285 * pkin(5) + t247 * pkin(11) - t304 * t265 + t218;
t348 = sin(qJ(6));
t353 = cos(qJ(6));
t213 = t353 * t215 - t348 * t216;
t259 = t353 * t286 - t348 * t287;
t233 = t259 * qJD(6) + t348 * t247 + t353 * t248;
t260 = t348 * t286 + t353 * t287;
t244 = -t259 * mrSges(7,1) + t260 * mrSges(7,2);
t299 = qJD(6) + t304;
t249 = -t299 * mrSges(7,2) + t259 * mrSges(7,3);
t266 = qJDD(6) + t270;
t208 = m(7) * t213 + t266 * mrSges(7,1) - t233 * mrSges(7,3) - t260 * t244 + t299 * t249;
t214 = t348 * t215 + t353 * t216;
t232 = -t260 * qJD(6) + t353 * t247 - t348 * t248;
t250 = t299 * mrSges(7,1) - t260 * mrSges(7,3);
t209 = m(7) * t214 - t266 * mrSges(7,2) + t232 * mrSges(7,3) + t259 * t244 - t299 * t250;
t200 = t353 * t208 + t348 * t209;
t261 = -t286 * mrSges(6,1) + t287 * mrSges(6,2);
t263 = -t304 * mrSges(6,2) + t286 * mrSges(6,3);
t198 = m(6) * t217 + t270 * mrSges(6,1) - t248 * mrSges(6,3) - t287 * t261 + t304 * t263 + t200;
t264 = t304 * mrSges(6,1) - t287 * mrSges(6,3);
t370 = -t348 * t208 + t353 * t209;
t199 = m(6) * t218 - t270 * mrSges(6,2) + t247 * mrSges(6,3) + t286 * t261 - t304 * t264 + t370;
t371 = -t349 * t198 + t354 * t199;
t191 = m(5) * t225 - t324 * mrSges(5,2) + t271 * mrSges(5,3) + t305 * t277 - t335 * t289 + t371;
t288 = -t335 * mrSges(5,2) + t305 * mrSges(5,3);
t221 = -t324 * pkin(4) - t334 * pkin(10) + t306 * t278 - t224;
t219 = -t247 * pkin(5) - t285 * pkin(11) + t287 * t265 + t221;
t367 = m(7) * t219 - t232 * mrSges(7,1) + t233 * mrSges(7,2) - t259 * t249 + t260 * t250;
t362 = -m(6) * t221 + t247 * mrSges(6,1) - t248 * mrSges(6,2) + t286 * t263 - t287 * t264 - t367;
t204 = m(5) * t224 + t324 * mrSges(5,1) - t272 * mrSges(5,3) - t306 * t277 + t335 * t288 + t362;
t182 = t344 * t191 + t346 * t204;
t307 = -t319 * mrSges(4,1) + t320 * mrSges(4,2);
t308 = -t335 * mrSges(4,2) + t319 * mrSges(4,3);
t180 = m(4) * t251 + t324 * mrSges(4,1) - t298 * mrSges(4,3) - t320 * t307 + t335 * t308 + t182;
t310 = t335 * mrSges(4,1) - t320 * mrSges(4,3);
t372 = t346 * t191 - t344 * t204;
t181 = m(4) * t252 - t324 * mrSges(4,2) + t297 * mrSges(4,3) + t319 * t307 - t335 * t310 + t372;
t175 = t355 * t180 + t350 * t181;
t193 = t354 * t198 + t349 * t199;
t301 = -g(3) * t384 + t380;
t325 = t340 * mrSges(3,1) - mrSges(3,3) * t375;
t329 = (-mrSges(3,1) * t356 + mrSges(3,2) * t351) * t379;
t373 = -t350 * t180 + t355 * t181;
t173 = m(3) * t301 - t339 * mrSges(3,2) - t332 * mrSges(3,3) - t340 * t325 + t329 * t374 + t373;
t326 = -t340 * mrSges(3,2) + mrSges(3,3) * t374;
t364 = m(5) * t245 - t271 * mrSges(5,1) + t272 * mrSges(5,2) - t305 * t288 + t306 * t289 + t193;
t361 = -m(4) * t282 + t297 * mrSges(4,1) - t298 * mrSges(4,2) + t319 * t308 - t320 * t310 - t364;
t188 = m(3) * t300 + t339 * mrSges(3,1) - t331 * mrSges(3,3) + t340 * t326 - t329 * t375 + t361;
t169 = t356 * t173 - t351 * t188;
t314 = -t345 * t327 - t385;
t174 = m(3) * t314 + t332 * mrSges(3,1) + t331 * mrSges(3,2) + (t325 * t351 - t326 * t356) * t379 + t175;
t166 = t173 * t382 - t345 * t174 + t188 * t381;
t237 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t299;
t239 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t299;
t201 = -mrSges(7,1) * t219 + mrSges(7,3) * t214 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t266 - t260 * t237 + t299 * t239;
t238 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t299;
t202 = mrSges(7,2) * t219 - mrSges(7,3) * t213 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t266 + t259 * t237 - t299 * t238;
t253 = Ifges(6,5) * t287 + Ifges(6,6) * t286 + Ifges(6,3) * t304;
t255 = Ifges(6,1) * t287 + Ifges(6,4) * t286 + Ifges(6,5) * t304;
t184 = -mrSges(6,1) * t221 + mrSges(6,3) * t218 + Ifges(6,4) * t248 + Ifges(6,2) * t247 + Ifges(6,6) * t270 - pkin(5) * t367 + pkin(11) * t370 + t353 * t201 + t348 * t202 - t287 * t253 + t304 * t255;
t254 = Ifges(6,4) * t287 + Ifges(6,2) * t286 + Ifges(6,6) * t304;
t186 = mrSges(6,2) * t221 - mrSges(6,3) * t217 + Ifges(6,1) * t248 + Ifges(6,4) * t247 + Ifges(6,5) * t270 - pkin(11) * t200 - t348 * t201 + t353 * t202 + t286 * t253 - t304 * t254;
t273 = Ifges(5,5) * t306 + Ifges(5,6) * t305 + Ifges(5,3) * t335;
t274 = Ifges(5,4) * t306 + Ifges(5,2) * t305 + Ifges(5,6) * t335;
t170 = mrSges(5,2) * t245 - mrSges(5,3) * t224 + Ifges(5,1) * t272 + Ifges(5,4) * t271 + Ifges(5,5) * t324 - pkin(10) * t193 - t349 * t184 + t354 * t186 + t305 * t273 - t335 * t274;
t275 = Ifges(5,1) * t306 + Ifges(5,4) * t305 + Ifges(5,5) * t335;
t365 = -mrSges(7,1) * t213 + mrSges(7,2) * t214 - Ifges(7,5) * t233 - Ifges(7,6) * t232 - Ifges(7,3) * t266 - t260 * t238 + t259 * t239;
t360 = mrSges(6,1) * t217 - mrSges(6,2) * t218 + Ifges(6,5) * t248 + Ifges(6,6) * t247 + Ifges(6,3) * t270 + pkin(5) * t200 + t287 * t254 - t286 * t255 - t365;
t176 = -mrSges(5,1) * t245 + mrSges(5,3) * t225 + Ifges(5,4) * t272 + Ifges(5,2) * t271 + Ifges(5,6) * t324 - pkin(4) * t193 - t306 * t273 + t335 * t275 - t360;
t291 = Ifges(4,5) * t320 + Ifges(4,6) * t319 + Ifges(4,3) * t335;
t293 = Ifges(4,1) * t320 + Ifges(4,4) * t319 + Ifges(4,5) * t335;
t159 = -mrSges(4,1) * t282 + mrSges(4,3) * t252 + Ifges(4,4) * t298 + Ifges(4,2) * t297 + Ifges(4,6) * t324 - pkin(3) * t364 + qJ(4) * t372 + t344 * t170 + t346 * t176 - t320 * t291 + t335 * t293;
t292 = Ifges(4,4) * t320 + Ifges(4,2) * t319 + Ifges(4,6) * t335;
t162 = mrSges(4,2) * t282 - mrSges(4,3) * t251 + Ifges(4,1) * t298 + Ifges(4,4) * t297 + Ifges(4,5) * t324 - qJ(4) * t182 + t346 * t170 - t344 * t176 + t319 * t291 - t335 * t292;
t312 = Ifges(3,6) * t340 + (Ifges(3,4) * t351 + Ifges(3,2) * t356) * t379;
t313 = Ifges(3,5) * t340 + (Ifges(3,1) * t351 + Ifges(3,4) * t356) * t379;
t156 = Ifges(3,5) * t331 - Ifges(3,6) * t332 + Ifges(3,3) * t339 + mrSges(3,1) * t300 - mrSges(3,2) * t301 + t350 * t162 + t355 * t159 + pkin(2) * t361 + pkin(9) * t373 + (t312 * t351 - t313 * t356) * t379;
t311 = Ifges(3,3) * t340 + (Ifges(3,5) * t351 + Ifges(3,6) * t356) * t379;
t158 = mrSges(3,2) * t314 - mrSges(3,3) * t300 + Ifges(3,1) * t331 - Ifges(3,4) * t332 + Ifges(3,5) * t339 - pkin(9) * t175 - t350 * t159 + t355 * t162 + t311 * t374 - t340 * t312;
t363 = -mrSges(5,1) * t224 + mrSges(5,2) * t225 - Ifges(5,5) * t272 - Ifges(5,6) * t271 - Ifges(5,3) * t324 - pkin(4) * t362 - pkin(10) * t371 - t354 * t184 - t349 * t186 - t306 * t274 + t305 * t275;
t359 = mrSges(4,1) * t251 - mrSges(4,2) * t252 + Ifges(4,5) * t298 + Ifges(4,6) * t297 + Ifges(4,3) * t324 + pkin(3) * t182 + t320 * t292 - t319 * t293 - t363;
t161 = -mrSges(3,1) * t314 + mrSges(3,3) * t301 + Ifges(3,4) * t331 - Ifges(3,2) * t332 + Ifges(3,6) * t339 - pkin(2) * t175 - t311 * t375 + t340 * t313 - t359;
t366 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t166 + t347 * t156 + t158 * t384 + t161 * t383 + t169 * t386;
t167 = m(2) * t337 - t358 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t169;
t165 = t347 * t174 + (t173 * t351 + t188 * t356) * t345;
t163 = m(2) * t336 + qJDD(1) * mrSges(2,1) - t358 * mrSges(2,2) + t166;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) + t356 * t158 - t351 * t161 + (-t165 * t345 - t166 * t347) * pkin(8);
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t345 * t156 + (pkin(8) * t169 + t158 * t351 + t161 * t356) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t154 - t352 * t153 - pkin(7) * (t357 * t163 + t352 * t167), t154, t158, t162, t170, t186, t202; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t154 + t357 * t153 + pkin(7) * (-t352 * t163 + t357 * t167), t153, t161, t159, t176, t184, t201; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t156, t359, -t363, t360, -t365;];
m_new  = t1;
