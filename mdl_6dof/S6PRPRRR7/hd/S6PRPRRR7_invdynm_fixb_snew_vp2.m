% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 02:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:49:47
% EndTime: 2019-05-05 01:52:21
% DurationCPUTime: 151.24s
% Computational Cost: add. (2540699->355), mult. (7174559->501), div. (0->0), fcn. (6122453->18), ass. (0->178)
t323 = sin(pkin(14));
t326 = sin(pkin(7));
t328 = cos(pkin(14));
t331 = cos(pkin(7));
t335 = sin(qJ(4));
t330 = cos(pkin(8));
t339 = cos(qJ(4));
t371 = t330 * t339;
t325 = sin(pkin(8));
t378 = t325 * t339;
t345 = t326 * (-t323 * t335 + t328 * t371) + t331 * t378;
t294 = t345 * qJD(2);
t372 = t330 * t335;
t379 = t325 * t335;
t347 = t331 * t379 + (t323 * t339 + t328 * t372) * t326;
t295 = t347 * qJD(2);
t278 = -t295 * qJD(4) + qJDD(2) * t345;
t324 = sin(pkin(13));
t329 = cos(pkin(13));
t318 = g(1) * t324 - g(2) * t329;
t319 = -g(1) * t329 - g(2) * t324;
t322 = -g(3) + qJDD(1);
t340 = cos(qJ(2));
t332 = cos(pkin(6));
t336 = sin(qJ(2));
t370 = t332 * t336;
t327 = sin(pkin(6));
t375 = t327 * t336;
t291 = t318 * t370 + t340 * t319 + t322 * t375;
t341 = qJD(2) ^ 2;
t382 = qJ(3) * t326;
t289 = -pkin(2) * t341 + qJDD(2) * t382 + t291;
t384 = pkin(10) * t323;
t358 = -pkin(3) * t328 - t325 * t384;
t368 = qJD(2) * t326;
t304 = t358 * t368;
t376 = t326 * t330;
t353 = pkin(10) * (t325 * t331 + t328 * t376);
t305 = qJD(2) * t353;
t369 = t332 * t340;
t374 = t327 * t340;
t290 = t318 * t369 - t319 * t336 + t322 * t374;
t288 = qJDD(2) * pkin(2) + t341 * t382 + t290;
t307 = -t318 * t327 + t322 * t332;
t363 = qJD(3) * t368;
t373 = t328 * t331;
t377 = t326 * t328;
t364 = t288 * t373 + t307 * t377 - 0.2e1 * t323 * t363;
t247 = (pkin(3) * qJDD(2) + qJD(2) * t305) * t331 + (-t289 + (-pkin(10) * qJDD(2) * t330 - qJD(2) * t304) * t326) * t323 + t364;
t380 = t323 * t331;
t381 = t323 * t326;
t263 = t288 * t380 + t307 * t381 + (t289 + 0.2e1 * t363) * t328;
t309 = (pkin(3) * t331 - t376 * t384) * qJD(2);
t248 = (t304 * t377 - t309 * t331) * qJD(2) + qJDD(2) * t353 + t263;
t366 = t331 * t307 + qJDD(3);
t261 = (-t288 + t358 * qJDD(2) + (-t305 * t328 + t309 * t323) * qJD(2)) * t326 + t366;
t233 = -t335 * t248 + (t247 * t330 + t261 * t325) * t339;
t234 = t247 * t372 + t339 * t248 + t261 * t379;
t276 = -mrSges(5,1) * t294 + mrSges(5,2) * t295;
t354 = -t325 * t377 + t330 * t331;
t306 = qJD(2) * t354 + qJD(4);
t284 = mrSges(5,1) * t306 - mrSges(5,3) * t295;
t303 = qJDD(2) * t354 + qJDD(4);
t277 = -pkin(4) * t294 - pkin(11) * t295;
t302 = t306 ^ 2;
t230 = -pkin(4) * t302 + pkin(11) * t303 + t277 * t294 + t234;
t235 = -t325 * t247 + t330 * t261;
t279 = t294 * qJD(4) + qJDD(2) * t347;
t232 = (-t294 * t306 - t279) * pkin(11) + (t295 * t306 - t278) * pkin(4) + t235;
t334 = sin(qJ(5));
t338 = cos(qJ(5));
t226 = t338 * t230 + t334 * t232;
t281 = -t295 * t334 + t306 * t338;
t282 = t295 * t338 + t306 * t334;
t265 = -pkin(5) * t281 - pkin(12) * t282;
t275 = qJDD(5) - t278;
t293 = qJD(5) - t294;
t292 = t293 ^ 2;
t224 = -pkin(5) * t292 + pkin(12) * t275 + t265 * t281 + t226;
t229 = -t303 * pkin(4) - t302 * pkin(11) + t295 * t277 - t233;
t257 = -qJD(5) * t282 - t279 * t334 + t303 * t338;
t258 = qJD(5) * t281 + t279 * t338 + t303 * t334;
t227 = (-t281 * t293 - t258) * pkin(12) + (t282 * t293 - t257) * pkin(5) + t229;
t333 = sin(qJ(6));
t337 = cos(qJ(6));
t220 = -t224 * t333 + t227 * t337;
t267 = -t282 * t333 + t293 * t337;
t238 = qJD(6) * t267 + t258 * t337 + t275 * t333;
t268 = t282 * t337 + t293 * t333;
t246 = -mrSges(7,1) * t267 + mrSges(7,2) * t268;
t280 = qJD(6) - t281;
t249 = -mrSges(7,2) * t280 + mrSges(7,3) * t267;
t255 = qJDD(6) - t257;
t218 = m(7) * t220 + mrSges(7,1) * t255 - mrSges(7,3) * t238 - t246 * t268 + t249 * t280;
t221 = t224 * t337 + t227 * t333;
t237 = -qJD(6) * t268 - t258 * t333 + t275 * t337;
t250 = mrSges(7,1) * t280 - mrSges(7,3) * t268;
t219 = m(7) * t221 - mrSges(7,2) * t255 + mrSges(7,3) * t237 + t246 * t267 - t250 * t280;
t212 = -t218 * t333 + t337 * t219;
t264 = -mrSges(6,1) * t281 + mrSges(6,2) * t282;
t270 = mrSges(6,1) * t293 - mrSges(6,3) * t282;
t210 = m(6) * t226 - mrSges(6,2) * t275 + mrSges(6,3) * t257 + t264 * t281 - t270 * t293 + t212;
t225 = -t230 * t334 + t232 * t338;
t223 = -pkin(5) * t275 - pkin(12) * t292 + t265 * t282 - t225;
t222 = -m(7) * t223 + t237 * mrSges(7,1) - mrSges(7,2) * t238 + t267 * t249 - t250 * t268;
t269 = -mrSges(6,2) * t293 + mrSges(6,3) * t281;
t216 = m(6) * t225 + mrSges(6,1) * t275 - mrSges(6,3) * t258 - t264 * t282 + t269 * t293 + t222;
t362 = t338 * t210 - t216 * t334;
t201 = m(5) * t234 - mrSges(5,2) * t303 + mrSges(5,3) * t278 + t276 * t294 - t284 * t306 + t362;
t204 = t334 * t210 + t338 * t216;
t283 = -mrSges(5,2) * t306 + mrSges(5,3) * t294;
t203 = m(5) * t235 - mrSges(5,1) * t278 + mrSges(5,2) * t279 - t283 * t294 + t284 * t295 + t204;
t211 = t218 * t337 + t219 * t333;
t344 = -m(6) * t229 + t257 * mrSges(6,1) - mrSges(6,2) * t258 + t281 * t269 - t270 * t282 - t211;
t207 = m(5) * t233 + mrSges(5,1) * t303 - mrSges(5,3) * t279 - t276 * t295 + t283 * t306 + t344;
t190 = t201 * t372 - t325 * t203 + t207 * t371;
t262 = -t323 * t289 + t364;
t361 = -mrSges(4,1) * t328 + mrSges(4,2) * t323;
t308 = t361 * t368;
t356 = -mrSges(4,2) * t331 + mrSges(4,3) * t377;
t313 = t356 * qJD(2);
t357 = mrSges(4,1) * t331 - mrSges(4,3) * t381;
t186 = m(4) * t262 + t357 * qJDD(2) + (-t308 * t381 + t313 * t331) * qJD(2) + t190;
t189 = t201 * t379 + t330 * t203 + t207 * t378;
t274 = -t326 * t288 + t366;
t312 = t357 * qJD(2);
t188 = m(4) * t274 + (t361 * qJDD(2) + (t312 * t323 - t313 * t328) * qJD(2)) * t326 + t189;
t195 = t339 * t201 - t335 * t207;
t194 = m(4) * t263 + t356 * qJDD(2) + (t308 * t377 - t312 * t331) * qJD(2) + t195;
t176 = t186 * t373 - t188 * t326 + t194 * t380;
t173 = m(3) * t290 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t341 + t176;
t180 = -t186 * t323 + t328 * t194;
t179 = m(3) * t291 - mrSges(3,1) * t341 - qJDD(2) * mrSges(3,2) + t180;
t171 = -t173 * t336 + t340 * t179;
t385 = pkin(9) * t171;
t383 = Ifges(4,3) * t331;
t175 = t186 * t377 + t331 * t188 + t194 * t381;
t174 = m(3) * t307 + t175;
t165 = t173 * t369 - t174 * t327 + t179 * t370;
t360 = Ifges(4,5) * t323 + Ifges(4,6) * t328;
t239 = Ifges(7,5) * t268 + Ifges(7,6) * t267 + Ifges(7,3) * t280;
t241 = Ifges(7,1) * t268 + Ifges(7,4) * t267 + Ifges(7,5) * t280;
t213 = -mrSges(7,1) * t223 + mrSges(7,3) * t221 + Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t255 - t239 * t268 + t241 * t280;
t240 = Ifges(7,4) * t268 + Ifges(7,2) * t267 + Ifges(7,6) * t280;
t214 = mrSges(7,2) * t223 - mrSges(7,3) * t220 + Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t255 + t239 * t267 - t240 * t280;
t251 = Ifges(6,5) * t282 + Ifges(6,6) * t281 + Ifges(6,3) * t293;
t252 = Ifges(6,4) * t282 + Ifges(6,2) * t281 + Ifges(6,6) * t293;
t196 = mrSges(6,2) * t229 - mrSges(6,3) * t225 + Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t275 - pkin(12) * t211 - t213 * t333 + t214 * t337 + t251 * t281 - t252 * t293;
t253 = Ifges(6,1) * t282 + Ifges(6,4) * t281 + Ifges(6,5) * t293;
t343 = mrSges(7,1) * t220 - mrSges(7,2) * t221 + Ifges(7,5) * t238 + Ifges(7,6) * t237 + Ifges(7,3) * t255 + t240 * t268 - t241 * t267;
t197 = -mrSges(6,1) * t229 + mrSges(6,3) * t226 + Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t275 - pkin(5) * t211 - t251 * t282 + t253 * t293 - t343;
t271 = Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t306;
t272 = Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t306;
t182 = mrSges(5,2) * t235 - mrSges(5,3) * t233 + Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t303 - pkin(11) * t204 + t196 * t338 - t197 * t334 + t271 * t294 - t272 * t306;
t273 = Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t306;
t342 = mrSges(6,1) * t225 - mrSges(6,2) * t226 + Ifges(6,5) * t258 + Ifges(6,6) * t257 + Ifges(6,3) * t275 + pkin(5) * t222 + pkin(12) * t212 + t337 * t213 + t333 * t214 + t282 * t252 - t281 * t253;
t183 = -mrSges(5,1) * t235 + mrSges(5,3) * t234 + Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t303 - pkin(4) * t204 - t295 * t271 + t306 * t273 - t342;
t352 = pkin(10) * t195 + t182 * t335 + t183 * t339;
t181 = mrSges(5,1) * t233 - mrSges(5,2) * t234 + Ifges(5,5) * t279 + Ifges(5,6) * t278 + Ifges(5,3) * t303 + pkin(4) * t344 + pkin(11) * t362 + t334 * t196 + t338 * t197 + t295 * t272 - t294 * t273;
t298 = (t326 * t360 + t383) * qJD(2);
t349 = Ifges(4,5) * t331 + (Ifges(4,1) * t323 + Ifges(4,4) * t328) * t326;
t300 = t349 * qJD(2);
t348 = Ifges(4,6) * t331 + (Ifges(4,4) * t323 + Ifges(4,2) * t328) * t326;
t167 = -mrSges(4,1) * t274 + mrSges(4,3) * t263 - pkin(3) * t189 - t325 * t181 + (-t298 * t381 + t300 * t331) * qJD(2) + t352 * t330 + t348 * qJDD(2);
t299 = t348 * qJD(2);
t168 = mrSges(4,2) * t274 - mrSges(4,3) * t262 + t339 * t182 - t335 * t183 + (t298 * t377 - t299 * t331) * qJD(2) + (-t189 * t325 - t190 * t330) * pkin(10) + t349 * qJDD(2);
t351 = qJ(3) * t180 + t167 * t328 + t168 * t323;
t166 = qJDD(2) * t383 + mrSges(4,1) * t262 - mrSges(4,2) * t263 + pkin(3) * t190 + t330 * t181 + t352 * t325 + (t360 * qJDD(2) + (t299 * t323 - t300 * t328) * qJD(2)) * t326;
t157 = mrSges(3,1) * t290 - mrSges(3,2) * t291 + Ifges(3,3) * qJDD(2) + pkin(2) * t176 + t331 * t166 + t326 * t351;
t159 = -mrSges(3,1) * t307 + mrSges(3,3) * t291 + t341 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t175 - t326 * t166 + t331 * t351;
t161 = mrSges(3,2) * t307 - mrSges(3,3) * t290 + Ifges(3,5) * qJDD(2) - t341 * Ifges(3,6) - t323 * t167 + t328 * t168 + (-t175 * t326 - t176 * t331) * qJ(3);
t350 = mrSges(2,1) * t318 - mrSges(2,2) * t319 + pkin(1) * t165 + t332 * t157 + t159 * t374 + t161 * t375 + t327 * t385;
t169 = m(2) * t319 + t171;
t164 = t332 * t174 + (t173 * t340 + t179 * t336) * t327;
t162 = m(2) * t318 + t165;
t155 = mrSges(2,2) * t322 - mrSges(2,3) * t318 - t336 * t159 + t340 * t161 + (-t164 * t327 - t165 * t332) * pkin(9);
t154 = -mrSges(2,1) * t322 + mrSges(2,3) * t319 - pkin(1) * t164 - t327 * t157 + (t159 * t340 + t161 * t336 + t385) * t332;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t329 * t155 - t324 * t154 - qJ(1) * (t162 * t329 + t169 * t324), t155, t161, t168, t182, t196, t214; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t324 * t155 + t329 * t154 + qJ(1) * (-t162 * t324 + t169 * t329), t154, t159, t167, t183, t197, t213; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t157, t166, t181, t342, t343;];
m_new  = t1;
