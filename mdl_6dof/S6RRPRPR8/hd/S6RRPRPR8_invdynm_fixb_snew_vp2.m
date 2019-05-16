% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:54:46
% EndTime: 2019-05-06 14:55:22
% DurationCPUTime: 16.17s
% Computational Cost: add. (260970->383), mult. (564089->464), div. (0->0), fcn. (398538->10), ass. (0->144)
t334 = sin(pkin(10));
t335 = cos(pkin(10));
t338 = sin(qJ(2));
t366 = qJD(1) * t338;
t313 = qJD(2) * t335 - t334 * t366;
t314 = qJD(2) * t334 + t335 * t366;
t337 = sin(qJ(4));
t371 = cos(qJ(4));
t280 = -t371 * t313 + t337 * t314;
t341 = cos(qJ(2));
t364 = qJD(1) * qJD(2);
t363 = t341 * t364;
t319 = qJDD(1) * t338 + t363;
t289 = qJDD(2) * t335 - t319 * t334;
t290 = qJDD(2) * t334 + t319 * t335;
t246 = -t280 * qJD(4) + t337 * t289 + t371 * t290;
t365 = qJD(1) * t341;
t291 = -pkin(3) * t365 - pkin(8) * t314;
t312 = t313 ^ 2;
t339 = sin(qJ(1));
t342 = cos(qJ(1));
t326 = -g(1) * t342 - g(2) * t339;
t344 = qJD(1) ^ 2;
t303 = -pkin(1) * t344 + qJDD(1) * pkin(7) + t326;
t284 = -t341 * g(3) - t338 * t303;
t317 = (-pkin(2) * t341 - qJ(3) * t338) * qJD(1);
t343 = qJD(2) ^ 2;
t354 = qJDD(2) * pkin(2) + t343 * qJ(3) - t317 * t366 - qJDD(3) + t284;
t351 = t289 * pkin(3) + t312 * pkin(8) - t314 * t291 + t354;
t329 = qJD(4) - t365;
t369 = t280 * t329;
t374 = (-t246 + t369) * qJ(5) - t351;
t281 = t337 * t313 + t371 * t314;
t245 = t281 * qJD(4) - t371 * t289 + t337 * t290;
t273 = -pkin(5) * t329 - pkin(9) * t281;
t279 = t280 ^ 2;
t372 = 2 * qJD(5);
t196 = -t279 * pkin(9) + (-pkin(4) - pkin(5)) * t245 + (-pkin(4) * t329 + t273 + t372) * t281 - t374;
t336 = sin(qJ(6));
t340 = cos(qJ(6));
t257 = t280 * t336 + t281 * t340;
t211 = -qJD(6) * t257 + t245 * t340 - t246 * t336;
t256 = t280 * t340 - t281 * t336;
t212 = qJD(6) * t256 + t245 * t336 + t246 * t340;
t327 = qJD(6) - t329;
t234 = -mrSges(7,2) * t327 + mrSges(7,3) * t256;
t235 = mrSges(7,1) * t327 - mrSges(7,3) * t257;
t192 = -m(7) * t196 + t211 * mrSges(7,1) - t212 * mrSges(7,2) + t256 * t234 - t257 * t235;
t203 = -0.2e1 * qJD(5) * t281 + (t281 * t329 + t245) * pkin(4) + t374;
t271 = -mrSges(6,1) * t329 + mrSges(6,2) * t281;
t272 = -mrSges(6,2) * t280 + mrSges(6,3) * t329;
t184 = m(6) * t203 + t245 * mrSges(6,1) - t246 * mrSges(6,3) - t281 * t271 + t280 * t272 + t192;
t325 = t339 * g(1) - t342 * g(2);
t302 = -qJDD(1) * pkin(1) - t344 * pkin(7) - t325;
t330 = t338 * t364;
t320 = qJDD(1) * t341 - t330;
t263 = (-t319 - t363) * qJ(3) + (-t320 + t330) * pkin(2) + t302;
t285 = -g(3) * t338 + t341 * t303;
t268 = -pkin(2) * t343 + qJDD(2) * qJ(3) + t317 * t365 + t285;
t226 = -0.2e1 * qJD(3) * t314 + t335 * t263 - t334 * t268;
t216 = (-t313 * t365 - t290) * pkin(8) + (t313 * t314 - t320) * pkin(3) + t226;
t227 = 0.2e1 * qJD(3) * t313 + t334 * t263 + t335 * t268;
t219 = -pkin(3) * t312 + pkin(8) * t289 + t291 * t365 + t227;
t206 = t337 * t216 + t371 * t219;
t251 = Ifges(6,1) * t281 + Ifges(6,4) * t329 + Ifges(6,5) * t280;
t252 = Ifges(5,1) * t281 - Ifges(5,4) * t280 + Ifges(5,5) * t329;
t316 = qJDD(4) - t320;
t205 = t371 * t216 - t337 * t219;
t258 = pkin(4) * t280 - qJ(5) * t281;
t328 = t329 ^ 2;
t201 = -t316 * pkin(4) - t328 * qJ(5) + t281 * t258 + qJDD(5) - t205;
t193 = (-t246 - t369) * pkin(9) + (t280 * t281 - t316) * pkin(5) + t201;
t199 = -pkin(4) * t328 + t316 * qJ(5) - t280 * t258 + t329 * t372 + t206;
t194 = -pkin(5) * t279 + pkin(9) * t245 + t273 * t329 + t199;
t190 = t193 * t340 - t194 * t336;
t225 = -mrSges(7,1) * t256 + mrSges(7,2) * t257;
t310 = qJDD(6) - t316;
t186 = m(7) * t190 + mrSges(7,1) * t310 - mrSges(7,3) * t212 - t225 * t257 + t234 * t327;
t191 = t193 * t336 + t194 * t340;
t187 = m(7) * t191 - mrSges(7,2) * t310 + mrSges(7,3) * t211 + t225 * t256 - t235 * t327;
t178 = -t336 * t186 + t340 * t187;
t220 = Ifges(7,5) * t257 + Ifges(7,6) * t256 + Ifges(7,3) * t327;
t222 = Ifges(7,1) * t257 + Ifges(7,4) * t256 + Ifges(7,5) * t327;
t180 = -mrSges(7,1) * t196 + mrSges(7,3) * t191 + Ifges(7,4) * t212 + Ifges(7,2) * t211 + Ifges(7,6) * t310 - t220 * t257 + t222 * t327;
t221 = Ifges(7,4) * t257 + Ifges(7,2) * t256 + Ifges(7,6) * t327;
t181 = mrSges(7,2) * t196 - mrSges(7,3) * t190 + Ifges(7,1) * t212 + Ifges(7,4) * t211 + Ifges(7,5) * t310 + t220 * t256 - t221 * t327;
t352 = -mrSges(6,1) * t203 + mrSges(6,2) * t199 - pkin(5) * t192 - pkin(9) * t178 - t340 * t180 - t336 * t181;
t249 = Ifges(6,4) * t281 + Ifges(6,2) * t329 + Ifges(6,6) * t280;
t368 = -Ifges(5,5) * t281 + Ifges(5,6) * t280 - Ifges(5,3) * t329 - t249;
t162 = mrSges(5,1) * t351 + mrSges(5,3) * t206 - pkin(4) * t184 + (t252 + t251) * t329 + (Ifges(5,6) - Ifges(6,6)) * t316 + t368 * t281 + (Ifges(5,4) - Ifges(6,5)) * t246 + (-Ifges(5,2) - Ifges(6,3)) * t245 + t352;
t250 = Ifges(5,4) * t281 - Ifges(5,2) * t280 + Ifges(5,6) * t329;
t177 = t340 * t186 + t336 * t187;
t247 = Ifges(6,5) * t281 + Ifges(6,6) * t329 + Ifges(6,3) * t280;
t353 = mrSges(6,2) * t201 - mrSges(6,3) * t203 + Ifges(6,1) * t246 + Ifges(6,4) * t316 + Ifges(6,5) * t245 - pkin(9) * t177 - t336 * t180 + t340 * t181 + t329 * t247;
t163 = -mrSges(5,2) * t351 - mrSges(5,3) * t205 + Ifges(5,1) * t246 - Ifges(5,4) * t245 + Ifges(5,5) * t316 - qJ(5) * t184 - t329 * t250 + t368 * t280 + t353;
t274 = Ifges(4,5) * t314 + Ifges(4,6) * t313 - Ifges(4,3) * t365;
t276 = Ifges(4,1) * t314 + Ifges(4,4) * t313 - Ifges(4,5) * t365;
t269 = -mrSges(5,2) * t329 - mrSges(5,3) * t280;
t270 = mrSges(5,1) * t329 - mrSges(5,3) * t281;
t349 = -m(5) * t351 + t245 * mrSges(5,1) + t246 * mrSges(5,2) + t280 * t269 + t281 * t270 + t184;
t358 = m(6) * t199 + t316 * mrSges(6,3) + t329 * t271 + t178;
t259 = mrSges(6,1) * t280 - mrSges(6,3) * t281;
t367 = -mrSges(5,1) * t280 - mrSges(5,2) * t281 - t259;
t370 = -mrSges(5,3) - mrSges(6,2);
t170 = m(5) * t206 - t316 * mrSges(5,2) + t370 * t245 - t329 * t270 + t367 * t280 + t358;
t355 = -m(6) * t201 + t316 * mrSges(6,1) + t329 * t272 - t177;
t172 = m(5) * t205 + t316 * mrSges(5,1) + t370 * t246 + t329 * t269 + t367 * t281 + t355;
t361 = t371 * t170 - t172 * t337;
t151 = mrSges(4,1) * t354 + mrSges(4,3) * t227 + Ifges(4,4) * t290 + Ifges(4,2) * t289 - Ifges(4,6) * t320 - pkin(3) * t349 + pkin(8) * t361 + t371 * t162 + t337 * t163 - t314 * t274 - t276 * t365;
t167 = t337 * t170 + t371 * t172;
t275 = Ifges(4,4) * t314 + Ifges(4,2) * t313 - Ifges(4,6) * t365;
t152 = -mrSges(4,2) * t354 - mrSges(4,3) * t226 + Ifges(4,1) * t290 + Ifges(4,4) * t289 - Ifges(4,5) * t320 - pkin(8) * t167 - t337 * t162 + t371 * t163 + t313 * t274 + t275 * t365;
t282 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t287 = mrSges(4,2) * t365 + mrSges(4,3) * t313;
t165 = m(4) * t226 - mrSges(4,1) * t320 - mrSges(4,3) * t290 - t282 * t314 - t287 * t365 + t167;
t288 = -mrSges(4,1) * t365 - mrSges(4,3) * t314;
t166 = m(4) * t227 + mrSges(4,2) * t320 + mrSges(4,3) * t289 + t282 * t313 + t288 * t365 + t361;
t161 = -t165 * t334 + t335 * t166;
t183 = m(4) * t354 + t289 * mrSges(4,1) - t290 * mrSges(4,2) + t313 * t287 - t314 * t288 - t349;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t341) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t341) * qJD(1);
t373 = mrSges(3,1) * t284 - mrSges(3,2) * t285 + Ifges(3,5) * t319 + Ifges(3,6) * t320 + Ifges(3,3) * qJDD(2) + pkin(2) * t183 + qJ(3) * t161 + t335 * t151 + t334 * t152 + (t300 * t338 - t301 * t341) * qJD(1);
t318 = (-mrSges(3,1) * t341 + mrSges(3,2) * t338) * qJD(1);
t322 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t366;
t159 = m(3) * t285 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t320 - qJD(2) * t322 + t318 * t365 + t161;
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t365;
t182 = m(3) * t284 + qJDD(2) * mrSges(3,1) - t319 * mrSges(3,3) + qJD(2) * t323 - t318 * t366 + t183;
t362 = t341 * t159 - t182 * t338;
t160 = t165 * t335 + t166 * t334;
t357 = mrSges(7,1) * t190 - mrSges(7,2) * t191 + Ifges(7,5) * t212 + Ifges(7,6) * t211 + Ifges(7,3) * t310 + t257 * t221 - t256 * t222;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t341) * qJD(1);
t148 = mrSges(3,2) * t302 - mrSges(3,3) * t284 + Ifges(3,1) * t319 + Ifges(3,4) * t320 + Ifges(3,5) * qJDD(2) - qJ(3) * t160 - qJD(2) * t300 - t151 * t334 + t152 * t335 + t299 * t365;
t348 = mrSges(6,1) * t201 - mrSges(6,3) * t199 - Ifges(6,4) * t246 - Ifges(6,2) * t316 - Ifges(6,6) * t245 + pkin(5) * t177 + t281 * t247 - t280 * t251 + t357;
t346 = mrSges(5,2) * t206 - t280 * t252 - qJ(5) * (-t245 * mrSges(6,2) - t280 * t259 + t358) - pkin(4) * (-t246 * mrSges(6,2) - t281 * t259 + t355) - mrSges(5,1) * t205 - t281 * t250 + Ifges(5,6) * t245 - Ifges(5,5) * t246 - Ifges(5,3) * t316 + t348;
t345 = mrSges(4,1) * t226 - mrSges(4,2) * t227 + Ifges(4,5) * t290 + Ifges(4,6) * t289 + pkin(3) * t167 + t314 * t275 - t313 * t276 - t346;
t150 = Ifges(3,6) * qJDD(2) - t299 * t366 - t345 - pkin(2) * t160 + (Ifges(3,2) + Ifges(4,3)) * t320 + mrSges(3,3) * t285 + qJD(2) * t301 - mrSges(3,1) * t302 + Ifges(3,4) * t319;
t350 = -m(3) * t302 + t320 * mrSges(3,1) - mrSges(3,2) * t319 - t322 * t366 + t323 * t365 - t160;
t356 = mrSges(2,1) * t325 - mrSges(2,2) * t326 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t362 + t338 * t148 + t341 * t150;
t156 = m(2) * t325 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t344 + t350;
t155 = t159 * t338 + t182 * t341;
t153 = m(2) * t326 - mrSges(2,1) * t344 - qJDD(1) * mrSges(2,2) + t362;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t326 + t344 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t155 - t373;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t325 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t344 - pkin(7) * t155 + t148 * t341 - t150 * t338;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t145 - t339 * t146 - pkin(6) * (t153 * t339 + t156 * t342), t145, t148, t152, t163, -t249 * t280 + t353, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t145 + t342 * t146 + pkin(6) * (t153 * t342 - t156 * t339), t146, t150, t151, t162, -t348, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t373, -Ifges(4,3) * t320 + t345, -t346, Ifges(6,5) * t246 + Ifges(6,6) * t316 + Ifges(6,3) * t245 + t281 * t249 - t329 * t251 - t352, t357;];
m_new  = t1;
