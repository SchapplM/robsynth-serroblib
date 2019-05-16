% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:30:15
% EndTime: 2019-05-05 21:30:45
% DurationCPUTime: 17.05s
% Computational Cost: add. (277296->361), mult. (664425->438), div. (0->0), fcn. (499055->10), ass. (0->144)
t328 = qJD(1) ^ 2;
t319 = sin(pkin(9));
t320 = cos(pkin(9));
t322 = sin(qJ(3));
t325 = cos(qJ(3));
t342 = t319 * t322 - t320 * t325;
t298 = t342 * qJD(1);
t323 = sin(qJ(1));
t326 = cos(qJ(1));
t307 = -t326 * g(1) - t323 * g(2);
t300 = -t328 * pkin(1) + qJDD(1) * qJ(2) + t307;
t358 = qJD(1) * qJD(2);
t355 = -t320 * g(3) - 0.2e1 * t319 * t358;
t368 = pkin(2) * t320;
t266 = (-pkin(7) * qJDD(1) + t328 * t368 - t300) * t319 + t355;
t287 = -t319 * g(3) + (t300 + 0.2e1 * t358) * t320;
t357 = qJDD(1) * t320;
t315 = t320 ^ 2;
t364 = t315 * t328;
t267 = -pkin(2) * t364 + pkin(7) * t357 + t287;
t240 = t322 * t266 + t325 * t267;
t343 = t319 * t325 + t320 * t322;
t299 = t343 * qJD(1);
t277 = t298 * mrSges(4,1) + t299 * mrSges(4,2);
t359 = t299 * qJD(3);
t284 = -qJDD(1) * t342 - t359;
t293 = qJD(3) * mrSges(4,1) - t299 * mrSges(4,3);
t282 = t298 * pkin(3) - t299 * pkin(8);
t327 = qJD(3) ^ 2;
t209 = -t327 * pkin(3) + qJDD(3) * pkin(8) - t298 * t282 + t240;
t314 = t319 ^ 2;
t306 = t323 * g(1) - t326 * g(2);
t349 = qJDD(2) - t306;
t283 = (-pkin(1) - t368) * qJDD(1) + (-qJ(2) + (-t314 - t315) * pkin(7)) * t328 + t349;
t360 = t298 * qJD(3);
t285 = qJDD(1) * t343 - t360;
t223 = (-t285 + t360) * pkin(8) + (-t284 + t359) * pkin(3) + t283;
t321 = sin(qJ(4));
t324 = cos(qJ(4));
t202 = -t321 * t209 + t324 * t223;
t290 = t324 * qJD(3) - t321 * t299;
t253 = t290 * qJD(4) + t321 * qJDD(3) + t324 * t285;
t281 = qJDD(4) - t284;
t291 = t321 * qJD(3) + t324 * t299;
t296 = qJD(4) + t298;
t198 = (t290 * t296 - t253) * qJ(5) + (t290 * t291 + t281) * pkin(4) + t202;
t203 = t324 * t209 + t321 * t223;
t252 = -t291 * qJD(4) + t324 * qJDD(3) - t321 * t285;
t262 = t296 * pkin(4) - t291 * qJ(5);
t289 = t290 ^ 2;
t200 = -t289 * pkin(4) + t252 * qJ(5) - t296 * t262 + t203;
t318 = sin(pkin(10));
t365 = cos(pkin(10));
t255 = -t290 * t365 + t318 * t291;
t369 = -2 * qJD(5);
t194 = t318 * t198 + t365 * t200 + t255 * t369;
t219 = -t252 * t365 + t318 * t253;
t256 = t318 * t290 + t291 * t365;
t243 = t296 * mrSges(6,1) - t256 * mrSges(6,3);
t233 = t255 * pkin(5) - t256 * qJ(6);
t295 = t296 ^ 2;
t189 = -t295 * pkin(5) + t281 * qJ(6) + 0.2e1 * qJD(6) * t296 - t255 * t233 + t194;
t244 = -t296 * mrSges(7,1) + t256 * mrSges(7,2);
t356 = m(7) * t189 + t281 * mrSges(7,3) + t296 * t244;
t234 = t255 * mrSges(7,1) - t256 * mrSges(7,3);
t362 = -t255 * mrSges(6,1) - t256 * mrSges(6,2) - t234;
t367 = -mrSges(6,3) - mrSges(7,2);
t175 = m(6) * t194 - t281 * mrSges(6,2) + t219 * t367 - t296 * t243 + t255 * t362 + t356;
t340 = t198 * t365 - t318 * t200;
t193 = t256 * t369 + t340;
t220 = t318 * t252 + t253 * t365;
t242 = -t296 * mrSges(6,2) - t255 * mrSges(6,3);
t191 = -t281 * pkin(5) - t295 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t233) * t256 - t340;
t241 = -t255 * mrSges(7,2) + t296 * mrSges(7,3);
t350 = -m(7) * t191 + t281 * mrSges(7,1) + t296 * t241;
t177 = m(6) * t193 + t281 * mrSges(6,1) + t220 * t367 + t296 * t242 + t256 * t362 + t350;
t170 = t318 * t175 + t365 * t177;
t258 = -t290 * mrSges(5,1) + t291 * mrSges(5,2);
t261 = -t296 * mrSges(5,2) + t290 * mrSges(5,3);
t168 = m(5) * t202 + t281 * mrSges(5,1) - t253 * mrSges(5,3) - t291 * t258 + t296 * t261 + t170;
t263 = t296 * mrSges(5,1) - t291 * mrSges(5,3);
t351 = t365 * t175 - t318 * t177;
t169 = m(5) * t203 - t281 * mrSges(5,2) + t252 * mrSges(5,3) + t290 * t258 - t296 * t263 + t351;
t352 = -t321 * t168 + t324 * t169;
t161 = m(4) * t240 - qJDD(3) * mrSges(4,2) + t284 * mrSges(4,3) - qJD(3) * t293 - t298 * t277 + t352;
t239 = t325 * t266 - t322 * t267;
t292 = -qJD(3) * mrSges(4,2) - t298 * mrSges(4,3);
t208 = -qJDD(3) * pkin(3) - t327 * pkin(8) + t299 * t282 - t239;
t201 = -t252 * pkin(4) - t289 * qJ(5) + t291 * t262 + qJDD(5) + t208;
t196 = -0.2e1 * qJD(6) * t256 + (t255 * t296 - t220) * qJ(6) + (t256 * t296 + t219) * pkin(5) + t201;
t186 = m(7) * t196 + t219 * mrSges(7,1) - t220 * mrSges(7,3) + t255 * t241 - t256 * t244;
t334 = m(6) * t201 + t219 * mrSges(6,1) + t220 * mrSges(6,2) + t255 * t242 + t256 * t243 + t186;
t330 = -m(5) * t208 + t252 * mrSges(5,1) - t253 * mrSges(5,2) + t290 * t261 - t291 * t263 - t334;
t179 = m(4) * t239 + qJDD(3) * mrSges(4,1) - t285 * mrSges(4,3) + qJD(3) * t292 - t299 * t277 + t330;
t156 = t322 * t161 + t325 * t179;
t286 = -t319 * t300 + t355;
t228 = Ifges(7,1) * t256 + Ifges(7,4) * t296 + Ifges(7,5) * t255;
t229 = Ifges(6,1) * t256 - Ifges(6,4) * t255 + Ifges(6,5) * t296;
t348 = -mrSges(7,1) * t196 + mrSges(7,2) * t189;
t226 = Ifges(7,4) * t256 + Ifges(7,2) * t296 + Ifges(7,6) * t255;
t363 = -Ifges(6,5) * t256 + Ifges(6,6) * t255 - Ifges(6,3) * t296 - t226;
t171 = -mrSges(6,1) * t201 + mrSges(6,3) * t194 - pkin(5) * t186 + (t228 + t229) * t296 + (Ifges(6,6) - Ifges(7,6)) * t281 + t363 * t256 + (Ifges(6,4) - Ifges(7,5)) * t220 + (-Ifges(6,2) - Ifges(7,3)) * t219 + t348;
t227 = Ifges(6,4) * t256 - Ifges(6,2) * t255 + Ifges(6,6) * t296;
t224 = Ifges(7,5) * t256 + Ifges(7,6) * t296 + Ifges(7,3) * t255;
t339 = mrSges(7,2) * t191 - mrSges(7,3) * t196 + Ifges(7,1) * t220 + Ifges(7,4) * t281 + Ifges(7,5) * t219 + t296 * t224;
t172 = mrSges(6,2) * t201 - mrSges(6,3) * t193 + Ifges(6,1) * t220 - Ifges(6,4) * t219 + Ifges(6,5) * t281 - qJ(6) * t186 - t296 * t227 + t255 * t363 + t339;
t246 = Ifges(5,5) * t291 + Ifges(5,6) * t290 + Ifges(5,3) * t296;
t248 = Ifges(5,1) * t291 + Ifges(5,4) * t290 + Ifges(5,5) * t296;
t150 = -mrSges(5,1) * t208 + mrSges(5,3) * t203 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t281 - pkin(4) * t334 + qJ(5) * t351 + t171 * t365 + t318 * t172 - t291 * t246 + t296 * t248;
t247 = Ifges(5,4) * t291 + Ifges(5,2) * t290 + Ifges(5,6) * t296;
t152 = mrSges(5,2) * t208 - mrSges(5,3) * t202 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t281 - qJ(5) * t170 - t318 * t171 + t172 * t365 + t290 * t246 - t296 * t247;
t269 = Ifges(4,4) * t299 - Ifges(4,2) * t298 + Ifges(4,6) * qJD(3);
t270 = Ifges(4,1) * t299 - Ifges(4,4) * t298 + Ifges(4,5) * qJD(3);
t335 = -mrSges(4,1) * t239 + mrSges(4,2) * t240 - Ifges(4,5) * t285 - Ifges(4,6) * t284 - Ifges(4,3) * qJDD(3) - pkin(3) * t330 - pkin(8) * t352 - t324 * t150 - t321 * t152 - t299 * t269 - t298 * t270;
t346 = Ifges(3,4) * t319 + Ifges(3,2) * t320;
t347 = Ifges(3,1) * t319 + Ifges(3,4) * t320;
t370 = -mrSges(3,1) * t286 + mrSges(3,2) * t287 - pkin(2) * t156 - (t319 * t346 - t320 * t347) * t328 + t335;
t366 = mrSges(3,2) * t319;
t163 = t324 * t168 + t321 * t169;
t345 = Ifges(3,5) * t319 + Ifges(3,6) * t320;
t361 = t328 * t345;
t341 = mrSges(3,3) * qJDD(1) + t328 * (-mrSges(3,1) * t320 + t366);
t154 = m(3) * t286 - t319 * t341 + t156;
t353 = t325 * t161 - t322 * t179;
t155 = m(3) * t287 + t320 * t341 + t353;
t354 = -t319 * t154 + t320 * t155;
t268 = Ifges(4,5) * t299 - Ifges(4,6) * t298 + Ifges(4,3) * qJD(3);
t144 = mrSges(4,2) * t283 - mrSges(4,3) * t239 + Ifges(4,1) * t285 + Ifges(4,4) * t284 + Ifges(4,5) * qJDD(3) - pkin(8) * t163 - qJD(3) * t269 - t321 * t150 + t324 * t152 - t298 * t268;
t337 = mrSges(7,1) * t191 - mrSges(7,3) * t189 - Ifges(7,4) * t220 - Ifges(7,2) * t281 - Ifges(7,6) * t219 + t256 * t224 - t255 * t228;
t331 = mrSges(6,2) * t194 - t255 * t229 - qJ(6) * (-t219 * mrSges(7,2) - t255 * t234 + t356) - pkin(5) * (-t220 * mrSges(7,2) - t256 * t234 + t350) - mrSges(6,1) * t193 - t256 * t227 + Ifges(6,6) * t219 - Ifges(6,5) * t220 - Ifges(6,3) * t281 + t337;
t329 = mrSges(5,1) * t202 - mrSges(5,2) * t203 + Ifges(5,5) * t253 + Ifges(5,6) * t252 + Ifges(5,3) * t281 + pkin(4) * t170 + t291 * t247 - t290 * t248 - t331;
t148 = -mrSges(4,1) * t283 + mrSges(4,3) * t240 + Ifges(4,4) * t285 + Ifges(4,2) * t284 + Ifges(4,6) * qJDD(3) - pkin(3) * t163 + qJD(3) * t270 - t299 * t268 - t329;
t297 = -qJDD(1) * pkin(1) - t328 * qJ(2) + t349;
t336 = m(4) * t283 - t284 * mrSges(4,1) + t285 * mrSges(4,2) + t298 * t292 + t299 * t293 + t163;
t140 = -mrSges(3,1) * t297 + mrSges(3,3) * t287 - pkin(2) * t336 + pkin(7) * t353 + qJDD(1) * t346 + t322 * t144 + t325 * t148 - t319 * t361;
t143 = mrSges(3,2) * t297 - mrSges(3,3) * t286 - pkin(7) * t156 + qJDD(1) * t347 + t325 * t144 - t322 * t148 + t320 * t361;
t333 = -m(3) * t297 + mrSges(3,1) * t357 - t336 + (t314 * t328 + t364) * mrSges(3,3);
t338 = -mrSges(2,2) * t307 + qJ(2) * t354 + t320 * t140 + t319 * t143 + pkin(1) * (-qJDD(1) * t366 + t333) + mrSges(2,1) * t306 + Ifges(2,3) * qJDD(1);
t157 = (mrSges(2,1) - t366) * qJDD(1) + t333 - t328 * mrSges(2,2) + m(2) * t306;
t147 = t320 * t154 + t319 * t155;
t145 = m(2) * t307 - t328 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t354;
t141 = -pkin(1) * t147 + mrSges(2,1) * g(3) + (Ifges(2,6) - t345) * qJDD(1) + t328 * Ifges(2,5) + mrSges(2,3) * t307 + t370;
t138 = -mrSges(2,2) * g(3) - mrSges(2,3) * t306 + Ifges(2,5) * qJDD(1) - t328 * Ifges(2,6) - qJ(2) * t147 - t319 * t140 + t320 * t143;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t326 * t138 - t323 * t141 - pkin(6) * (t323 * t145 + t326 * t157), t138, t143, t144, t152, t172, -t255 * t226 + t339; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t323 * t138 + t326 * t141 + pkin(6) * (t326 * t145 - t323 * t157), t141, t140, t148, t150, t171, -t337; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t338, t338, qJDD(1) * t345 - t370, -t335, t329, -t331, Ifges(7,5) * t220 + Ifges(7,6) * t281 + Ifges(7,3) * t219 + t256 * t226 - t296 * t228 - t348;];
m_new  = t1;
