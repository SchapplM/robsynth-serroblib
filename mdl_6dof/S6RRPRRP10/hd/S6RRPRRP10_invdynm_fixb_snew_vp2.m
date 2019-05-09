% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:39:33
% EndTime: 2019-05-06 18:40:38
% DurationCPUTime: 36.91s
% Computational Cost: add. (641422->395), mult. (1453046->499), div. (0->0), fcn. (1164468->12), ass. (0->154)
t342 = sin(qJ(2));
t345 = cos(qJ(2));
t337 = sin(pkin(6));
t367 = qJD(1) * t337;
t319 = (-pkin(2) * t345 - qJ(3) * t342) * t367;
t339 = cos(pkin(6));
t332 = qJD(1) * t339 + qJD(2);
t330 = t332 ^ 2;
t331 = qJDD(1) * t339 + qJDD(2);
t366 = qJD(1) * t345;
t343 = sin(qJ(1));
t346 = cos(qJ(1));
t327 = g(1) * t343 - g(2) * t346;
t347 = qJD(1) ^ 2;
t378 = pkin(8) * t337;
t317 = qJDD(1) * pkin(1) + t347 * t378 + t327;
t328 = -g(1) * t346 - g(2) * t343;
t365 = qJDD(1) * t337;
t318 = -pkin(1) * t347 + pkin(8) * t365 + t328;
t373 = t339 * t342;
t368 = t317 * t373 + t318 * t345;
t275 = -t330 * pkin(2) + t331 * qJ(3) + (-g(3) * t342 + t319 * t366) * t337 + t368;
t321 = (qJD(2) * t366 + qJDD(1) * t342) * t337;
t363 = t342 * t367;
t322 = -qJD(2) * t363 + t345 * t365;
t377 = t339 * g(3);
t276 = -t322 * pkin(2) - t377 - t321 * qJ(3) + (-t317 + (pkin(2) * t342 - qJ(3) * t345) * t332 * qJD(1)) * t337;
t336 = sin(pkin(11));
t338 = cos(pkin(11));
t310 = t332 * t336 + t338 * t363;
t225 = -0.2e1 * qJD(3) * t310 - t336 * t275 + t276 * t338;
t297 = t321 * t338 + t331 * t336;
t309 = t332 * t338 - t336 * t363;
t362 = t337 * t366;
t219 = (-t309 * t362 - t297) * pkin(9) + (t309 * t310 - t322) * pkin(3) + t225;
t226 = 0.2e1 * qJD(3) * t309 + t275 * t338 + t276 * t336;
t296 = -t321 * t336 + t331 * t338;
t298 = -pkin(3) * t362 - pkin(9) * t310;
t307 = t309 ^ 2;
t222 = -pkin(3) * t307 + pkin(9) * t296 + t298 * t362 + t226;
t341 = sin(qJ(4));
t344 = cos(qJ(4));
t215 = t219 * t341 + t222 * t344;
t290 = t309 * t344 - t310 * t341;
t291 = t309 * t341 + t310 * t344;
t270 = -pkin(4) * t290 - pkin(10) * t291;
t314 = qJDD(4) - t322;
t326 = qJD(4) - t362;
t325 = t326 ^ 2;
t212 = -pkin(4) * t325 + pkin(10) * t314 + t270 * t290 + t215;
t372 = t339 * t345;
t374 = t337 * t345;
t286 = -g(3) * t374 + t317 * t372 - t318 * t342;
t274 = -t331 * pkin(2) - t330 * qJ(3) + t319 * t363 + qJDD(3) - t286;
t235 = -t296 * pkin(3) - t307 * pkin(9) + t298 * t310 + t274;
t259 = -qJD(4) * t291 + t296 * t344 - t297 * t341;
t260 = qJD(4) * t290 + t296 * t341 + t297 * t344;
t217 = (-t290 * t326 - t260) * pkin(10) + (t291 * t326 - t259) * pkin(4) + t235;
t340 = sin(qJ(5));
t379 = cos(qJ(5));
t208 = -t212 * t340 + t217 * t379;
t209 = t212 * t379 + t217 * t340;
t278 = t291 * t379 + t326 * t340;
t233 = qJD(5) * t278 + t260 * t340 - t314 * t379;
t277 = t291 * t340 - t326 * t379;
t234 = -qJD(5) * t277 + t260 * t379 + t314 * t340;
t289 = qJD(5) - t290;
t236 = Ifges(7,5) * t278 + Ifges(7,6) * t289 + Ifges(7,3) * t277;
t239 = Ifges(6,4) * t278 - Ifges(6,2) * t277 + Ifges(6,6) * t289;
t241 = Ifges(6,1) * t278 - Ifges(6,4) * t277 + Ifges(6,5) * t289;
t249 = mrSges(7,1) * t277 - mrSges(7,3) * t278;
t258 = qJDD(5) - t259;
t248 = pkin(5) * t277 - qJ(6) * t278;
t288 = t289 ^ 2;
t204 = -pkin(5) * t288 + qJ(6) * t258 + 0.2e1 * qJD(6) * t289 - t248 * t277 + t209;
t206 = -pkin(5) * t258 - qJ(6) * t288 + t248 * t278 + qJDD(6) - t208;
t240 = Ifges(7,1) * t278 + Ifges(7,4) * t289 + Ifges(7,5) * t277;
t356 = mrSges(7,1) * t206 - mrSges(7,3) * t204 - Ifges(7,4) * t234 - Ifges(7,2) * t258 - Ifges(7,6) * t233 - t240 * t277;
t261 = -mrSges(7,2) * t277 + mrSges(7,3) * t289;
t358 = -m(7) * t206 + mrSges(7,1) * t258 + t261 * t289;
t264 = -mrSges(7,1) * t289 + mrSges(7,2) * t278;
t364 = m(7) * t204 + mrSges(7,3) * t258 + t264 * t289;
t380 = -(-t239 + t236) * t278 + mrSges(6,1) * t208 - mrSges(6,2) * t209 + Ifges(6,5) * t234 - Ifges(6,6) * t233 + Ifges(6,3) * t258 + pkin(5) * (-mrSges(7,2) * t234 - t249 * t278 + t358) + qJ(6) * (-mrSges(7,2) * t233 - t249 * t277 + t364) + t277 * t241 - t356;
t376 = -mrSges(6,3) - mrSges(7,2);
t375 = t337 * t342;
t269 = -mrSges(5,1) * t290 + mrSges(5,2) * t291;
t280 = mrSges(5,1) * t326 - mrSges(5,3) * t291;
t263 = mrSges(6,1) * t289 - mrSges(6,3) * t278;
t369 = -mrSges(6,1) * t277 - mrSges(6,2) * t278 - t249;
t194 = m(6) * t209 - mrSges(6,2) * t258 + t233 * t376 - t263 * t289 + t277 * t369 + t364;
t262 = -mrSges(6,2) * t289 - mrSges(6,3) * t277;
t196 = m(6) * t208 + mrSges(6,1) * t258 + t234 * t376 + t262 * t289 + t278 * t369 + t358;
t359 = t194 * t379 - t196 * t340;
t182 = m(5) * t215 - mrSges(5,2) * t314 + mrSges(5,3) * t259 + t269 * t290 - t280 * t326 + t359;
t214 = t344 * t219 - t222 * t341;
t279 = -mrSges(5,2) * t326 + mrSges(5,3) * t290;
t211 = -t314 * pkin(4) - t325 * pkin(10) + t270 * t291 - t214;
t207 = -0.2e1 * qJD(6) * t278 + (t277 * t289 - t234) * qJ(6) + (t278 * t289 + t233) * pkin(5) + t211;
t201 = m(7) * t207 + mrSges(7,1) * t233 - mrSges(7,3) * t234 + t261 * t277 - t264 * t278;
t351 = -m(6) * t211 - mrSges(6,1) * t233 - mrSges(6,2) * t234 - t262 * t277 - t263 * t278 - t201;
t191 = m(5) * t214 + mrSges(5,1) * t314 - mrSges(5,3) * t260 - t269 * t291 + t279 * t326 + t351;
t177 = t182 * t341 + t191 * t344;
t292 = -mrSges(4,1) * t309 + mrSges(4,2) * t310;
t294 = mrSges(4,2) * t362 + mrSges(4,3) * t309;
t175 = m(4) * t225 - mrSges(4,1) * t322 - mrSges(4,3) * t297 - t292 * t310 - t294 * t362 + t177;
t295 = -mrSges(4,1) * t362 - mrSges(4,3) * t310;
t360 = t182 * t344 - t191 * t341;
t176 = m(4) * t226 + mrSges(4,2) * t322 + mrSges(4,3) * t296 + t292 * t309 + t295 * t362 + t360;
t169 = t175 * t338 + t176 * t336;
t188 = t194 * t340 + t196 * t379;
t238 = Ifges(7,4) * t278 + Ifges(7,2) * t289 + Ifges(7,6) * t277;
t371 = -Ifges(6,5) * t278 + Ifges(6,6) * t277 - Ifges(6,3) * t289 - t238;
t287 = -g(3) * t375 + t368;
t315 = mrSges(3,1) * t332 - mrSges(3,3) * t363;
t320 = (-mrSges(3,1) * t345 + mrSges(3,2) * t342) * t367;
t361 = -t175 * t336 + t176 * t338;
t167 = m(3) * t287 - mrSges(3,2) * t331 + mrSges(3,3) * t322 - t315 * t332 + t320 * t362 + t361;
t316 = -mrSges(3,2) * t332 + mrSges(3,3) * t362;
t353 = m(5) * t235 - mrSges(5,1) * t259 + mrSges(5,2) * t260 - t279 * t290 + t280 * t291 + t188;
t350 = -m(4) * t274 + mrSges(4,1) * t296 - mrSges(4,2) * t297 + t294 * t309 - t295 * t310 - t353;
t179 = m(3) * t286 + mrSges(3,1) * t331 - mrSges(3,3) * t321 + t316 * t332 - t320 * t363 + t350;
t164 = t167 * t345 - t179 * t342;
t302 = -t337 * t317 - t377;
t168 = m(3) * t302 - t322 * mrSges(3,1) + t321 * mrSges(3,2) + (t315 * t342 - t316 * t345) * t367 + t169;
t160 = t167 * t373 - t168 * t337 + t179 * t372;
t357 = -mrSges(7,1) * t207 + mrSges(7,2) * t204;
t355 = mrSges(7,2) * t206 - mrSges(7,3) * t207 + Ifges(7,1) * t234 + Ifges(7,4) * t258 + Ifges(7,5) * t233 + t236 * t289;
t184 = -mrSges(6,1) * t211 + mrSges(6,3) * t209 - pkin(5) * t201 + (t240 + t241) * t289 + t371 * t278 + (Ifges(6,6) - Ifges(7,6)) * t258 + (Ifges(6,4) - Ifges(7,5)) * t234 + (-Ifges(6,2) - Ifges(7,3)) * t233 + t357;
t186 = mrSges(6,2) * t211 - mrSges(6,3) * t208 + Ifges(6,1) * t234 - Ifges(6,4) * t233 + Ifges(6,5) * t258 - qJ(6) * t201 - t239 * t289 + t277 * t371 + t355;
t265 = Ifges(5,5) * t291 + Ifges(5,6) * t290 + Ifges(5,3) * t326;
t266 = Ifges(5,4) * t291 + Ifges(5,2) * t290 + Ifges(5,6) * t326;
t170 = mrSges(5,2) * t235 - mrSges(5,3) * t214 + Ifges(5,1) * t260 + Ifges(5,4) * t259 + Ifges(5,5) * t314 - pkin(10) * t188 - t184 * t340 + t186 * t379 + t265 * t290 - t266 * t326;
t267 = Ifges(5,1) * t291 + Ifges(5,4) * t290 + Ifges(5,5) * t326;
t171 = -mrSges(5,1) * t235 + mrSges(5,3) * t215 + Ifges(5,4) * t260 + Ifges(5,2) * t259 + Ifges(5,6) * t314 - pkin(4) * t188 - t291 * t265 + t326 * t267 - t380;
t281 = Ifges(4,5) * t310 + Ifges(4,6) * t309 - Ifges(4,3) * t362;
t283 = Ifges(4,1) * t310 + Ifges(4,4) * t309 - Ifges(4,5) * t362;
t156 = -mrSges(4,1) * t274 + mrSges(4,3) * t226 + Ifges(4,4) * t297 + Ifges(4,2) * t296 - Ifges(4,6) * t322 - pkin(3) * t353 + pkin(9) * t360 + t341 * t170 + t344 * t171 - t310 * t281 - t283 * t362;
t282 = Ifges(4,4) * t310 + Ifges(4,2) * t309 - Ifges(4,6) * t362;
t161 = mrSges(4,2) * t274 - mrSges(4,3) * t225 + Ifges(4,1) * t297 + Ifges(4,4) * t296 - Ifges(4,5) * t322 - pkin(9) * t177 + t170 * t344 - t171 * t341 + t281 * t309 + t282 * t362;
t300 = Ifges(3,6) * t332 + (Ifges(3,4) * t342 + Ifges(3,2) * t345) * t367;
t301 = Ifges(3,5) * t332 + (Ifges(3,1) * t342 + Ifges(3,4) * t345) * t367;
t151 = Ifges(3,5) * t321 + Ifges(3,6) * t322 + Ifges(3,3) * t331 + mrSges(3,1) * t286 - mrSges(3,2) * t287 + t336 * t161 + t338 * t156 + pkin(2) * t350 + qJ(3) * t361 + (t300 * t342 - t301 * t345) * t367;
t299 = Ifges(3,3) * t332 + (Ifges(3,5) * t342 + Ifges(3,6) * t345) * t367;
t153 = mrSges(3,2) * t302 - mrSges(3,3) * t286 + Ifges(3,1) * t321 + Ifges(3,4) * t322 + Ifges(3,5) * t331 - qJ(3) * t169 - t156 * t336 + t161 * t338 + t299 * t362 - t300 * t332;
t352 = -mrSges(5,1) * t214 + mrSges(5,2) * t215 - Ifges(5,5) * t260 - Ifges(5,6) * t259 - Ifges(5,3) * t314 - pkin(4) * t351 - pkin(10) * t359 - t184 * t379 - t186 * t340 - t266 * t291 + t290 * t267;
t348 = -mrSges(4,1) * t225 + mrSges(4,2) * t226 - Ifges(4,5) * t297 - Ifges(4,6) * t296 - pkin(3) * t177 - t310 * t282 + t309 * t283 + t352;
t155 = Ifges(3,4) * t321 - t299 * t363 + Ifges(3,6) * t331 + t332 * t301 + t348 - mrSges(3,1) * t302 - pkin(2) * t169 + mrSges(3,3) * t287 + (Ifges(3,2) + Ifges(4,3)) * t322;
t354 = mrSges(2,1) * t327 - mrSges(2,2) * t328 + Ifges(2,3) * qJDD(1) + pkin(1) * t160 + t151 * t339 + t153 * t375 + t155 * t374 + t164 * t378;
t162 = m(2) * t328 - mrSges(2,1) * t347 - qJDD(1) * mrSges(2,2) + t164;
t159 = t339 * t168 + (t167 * t342 + t179 * t345) * t337;
t157 = m(2) * t327 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t347 + t160;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - t347 * Ifges(2,6) + t345 * t153 - t342 * t155 + (-t159 * t337 - t160 * t339) * pkin(8);
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t159 - t337 * t151 + (pkin(8) * t164 + t153 * t342 + t155 * t345) * t339;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t149 - t343 * t148 - pkin(7) * (t157 * t346 + t162 * t343), t149, t153, t161, t170, t186, -t238 * t277 + t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t343 * t149 + t346 * t148 + pkin(7) * (-t157 * t343 + t162 * t346), t148, t155, t156, t171, t184, -t278 * t236 - t356; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t354, t354, t151, -Ifges(4,3) * t322 - t348, -t352, t380, Ifges(7,5) * t234 + Ifges(7,6) * t258 + Ifges(7,3) * t233 + t238 * t278 - t240 * t289 - t357;];
m_new  = t1;
