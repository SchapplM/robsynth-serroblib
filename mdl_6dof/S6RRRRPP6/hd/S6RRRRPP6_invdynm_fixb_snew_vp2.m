% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:31:56
% EndTime: 2019-05-07 18:32:18
% DurationCPUTime: 9.06s
% Computational Cost: add. (145326->384), mult. (288531->442), div. (0->0), fcn. (199066->8), ass. (0->140)
t328 = sin(qJ(3));
t331 = cos(qJ(3));
t329 = sin(qJ(2));
t360 = qJD(1) * t329;
t309 = qJD(2) * t328 + t331 * t360;
t332 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t356 = t332 * t358;
t312 = qJDD(1) * t329 + t356;
t275 = -qJD(3) * t309 + qJDD(2) * t331 - t312 * t328;
t308 = qJD(2) * t331 - t328 * t360;
t276 = qJD(3) * t308 + qJDD(2) * t328 + t312 * t331;
t327 = sin(qJ(4));
t368 = cos(qJ(4));
t278 = -t368 * t308 + t309 * t327;
t226 = -t278 * qJD(4) + t327 * t275 + t368 * t276;
t279 = t327 * t308 + t368 * t309;
t246 = -mrSges(7,2) * t279 + mrSges(7,3) * t278;
t330 = sin(qJ(1));
t333 = cos(qJ(1));
t317 = g(1) * t330 - t333 * g(2);
t335 = qJD(1) ^ 2;
t301 = -qJDD(1) * pkin(1) - pkin(7) * t335 - t317;
t323 = t329 * t358;
t313 = qJDD(1) * t332 - t323;
t252 = (-t312 - t356) * pkin(8) + (-t313 + t323) * pkin(2) + t301;
t318 = -g(1) * t333 - g(2) * t330;
t302 = -pkin(1) * t335 + qJDD(1) * pkin(7) + t318;
t284 = -g(3) * t329 + t332 * t302;
t311 = (-pkin(2) * t332 - pkin(8) * t329) * qJD(1);
t334 = qJD(2) ^ 2;
t359 = qJD(1) * t332;
t258 = -pkin(2) * t334 + qJDD(2) * pkin(8) + t311 * t359 + t284;
t227 = t331 * t252 - t258 * t328;
t307 = qJDD(3) - t313;
t322 = qJD(3) - t359;
t200 = (t308 * t322 - t276) * pkin(9) + (t308 * t309 + t307) * pkin(3) + t227;
t228 = t328 * t252 + t331 * t258;
t285 = pkin(3) * t322 - pkin(9) * t309;
t306 = t308 ^ 2;
t203 = -pkin(3) * t306 + pkin(9) * t275 - t285 * t322 + t228;
t197 = t368 * t200 - t327 * t203;
t247 = pkin(4) * t278 - qJ(5) * t279;
t303 = qJDD(4) + t307;
t320 = -qJD(4) - t322;
t319 = t320 ^ 2;
t193 = -t303 * pkin(4) - t319 * qJ(5) + t279 * t247 + qJDD(5) - t197;
t364 = t278 * t320;
t369 = 2 * qJD(6);
t184 = t320 * t369 + (t278 * t279 - t303) * qJ(6) + (t226 - t364) * pkin(5) + t193;
t262 = -mrSges(7,1) * t278 - mrSges(7,2) * t320;
t352 = -m(7) * t184 + t303 * mrSges(7,3) - t320 * t262;
t181 = mrSges(7,1) * t226 + t246 * t279 - t352;
t198 = t327 * t200 + t368 * t203;
t225 = qJD(4) * t279 - t368 * t275 + t276 * t327;
t237 = Ifges(5,4) * t279 - Ifges(5,2) * t278 - Ifges(5,6) * t320;
t261 = mrSges(6,1) * t278 + mrSges(6,3) * t320;
t346 = -pkin(4) * t319 + qJ(5) * t303 - t247 * t278 + t198;
t191 = 0.2e1 * qJD(5) * t320 - t346;
t233 = -Ifges(6,5) * t320 - Ifges(6,6) * t279 + Ifges(6,3) * t278;
t259 = pkin(5) * t279 + qJ(6) * t320;
t277 = t278 ^ 2;
t370 = -0.2e1 * qJD(5);
t187 = -pkin(5) * t225 - qJ(6) * t277 + qJDD(6) + (t370 - t259) * t320 + t346;
t232 = -Ifges(7,5) * t320 + Ifges(7,6) * t278 + Ifges(7,3) * t279;
t235 = -Ifges(7,4) * t320 + Ifges(7,2) * t278 + Ifges(7,6) * t279;
t345 = -mrSges(7,2) * t187 + mrSges(7,3) * t184 - Ifges(7,1) * t303 - Ifges(7,4) * t225 - Ifges(7,5) * t226 - t278 * t232 + t279 * t235;
t340 = -mrSges(6,2) * t193 + mrSges(6,3) * t191 - Ifges(6,1) * t303 + Ifges(6,4) * t226 - Ifges(6,5) * t225 + qJ(6) * t181 + t279 * t233 + t345;
t263 = mrSges(6,1) * t279 - mrSges(6,2) * t320;
t260 = mrSges(7,1) * t279 + mrSges(7,3) * t320;
t357 = m(7) * t187 + t303 * mrSges(7,2) - t320 * t260;
t350 = -m(6) * t191 + t303 * mrSges(6,3) - t320 * t263 + t357;
t249 = -mrSges(6,2) * t278 - mrSges(6,3) * t279;
t361 = -t246 - t249;
t236 = -Ifges(6,4) * t320 - Ifges(6,2) * t279 + Ifges(6,6) * t278;
t362 = Ifges(5,1) * t279 - Ifges(5,4) * t278 - Ifges(5,5) * t320 - t236;
t372 = -m(6) * t193 - t226 * mrSges(6,1) - t279 * t249;
t374 = -mrSges(5,2) * t198 + pkin(4) * (-mrSges(6,2) * t303 + t261 * t320 - t181 + t372) + qJ(5) * (t361 * t278 + (-mrSges(6,1) - mrSges(7,1)) * t225 + t350) + mrSges(5,1) * t197 - Ifges(5,6) * t225 + Ifges(5,5) * t226 + t279 * t237 + Ifges(5,3) * t303 - t340 + t362 * t278;
t248 = mrSges(5,1) * t278 + mrSges(5,2) * t279;
t264 = mrSges(5,2) * t320 - mrSges(5,3) * t278;
t366 = -mrSges(7,1) - mrSges(5,3);
t171 = m(5) * t197 + (t261 - t264) * t320 + (mrSges(5,1) - mrSges(6,2)) * t303 + (-t246 - t248) * t279 + t366 * t226 + t352 + t372;
t265 = -mrSges(5,1) * t320 - mrSges(5,3) * t279;
t174 = m(5) * t198 - mrSges(5,2) * t303 + t265 * t320 + (-t248 + t361) * t278 + (-mrSges(6,1) + t366) * t225 + t350;
t167 = t368 * t171 + t327 * t174;
t269 = Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * t322;
t270 = Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * t322;
t373 = mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t276 + Ifges(4,6) * t275 + Ifges(4,3) * t307 + pkin(3) * t167 + t309 * t269 - t308 * t270 + t374;
t283 = -t332 * g(3) - t329 * t302;
t257 = -qJDD(2) * pkin(2) - pkin(8) * t334 + t311 * t360 - t283;
t204 = -pkin(3) * t275 - pkin(9) * t306 + t309 * t285 + t257;
t339 = (-t226 - t364) * qJ(5) + t204 + (-t320 * pkin(4) + t370) * t279;
t190 = t339 + (pkin(4) + qJ(6)) * t225 - pkin(5) * t277 + t278 * t369 - t259 * t279;
t180 = m(7) * t190 - t226 * mrSges(7,2) + t225 * mrSges(7,3) - t279 * t260 + t278 * t262;
t195 = pkin(4) * t225 + t339;
t175 = m(6) * t195 - t225 * mrSges(6,2) - t226 * mrSges(6,3) - t278 * t261 - t279 * t263 + t180;
t234 = Ifges(5,5) * t279 - Ifges(5,6) * t278 - Ifges(5,3) * t320;
t239 = -Ifges(6,1) * t320 - Ifges(6,4) * t279 + Ifges(6,5) * t278;
t238 = -Ifges(7,1) * t320 + Ifges(7,4) * t278 + Ifges(7,5) * t279;
t349 = mrSges(7,1) * t187 - mrSges(7,3) * t190 - Ifges(7,4) * t303 - Ifges(7,2) * t225 - Ifges(7,6) * t226 - t279 * t238;
t342 = mrSges(6,1) * t191 - mrSges(6,2) * t195 + pkin(5) * (mrSges(7,1) * t225 + t246 * t278 - t357) + qJ(6) * t180 - t349;
t365 = Ifges(5,4) + Ifges(6,6);
t162 = -mrSges(5,1) * t204 + mrSges(5,3) * t198 - pkin(4) * t175 + (-t232 - t362) * t320 + (Ifges(5,6) - Ifges(6,5)) * t303 + (-t234 - t239) * t279 + (-Ifges(5,2) - Ifges(6,3)) * t225 + t365 * t226 - t342;
t348 = -mrSges(7,1) * t184 + mrSges(7,2) * t190 - Ifges(7,5) * t303 - Ifges(7,6) * t225 - Ifges(7,3) * t226 + t320 * t235;
t344 = -mrSges(6,1) * t193 + mrSges(6,3) * t195 - pkin(5) * t181 + t348;
t363 = t238 + t239;
t163 = mrSges(5,2) * t204 - mrSges(5,3) * t197 - qJ(5) * t175 - t344 - t365 * t225 + (Ifges(5,1) + Ifges(6,2)) * t226 + (-t234 - t363) * t278 + (Ifges(5,5) - Ifges(6,4)) * t303 + (t237 - t233) * t320;
t268 = Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * t322;
t341 = m(5) * t204 + t225 * mrSges(5,1) + t226 * mrSges(5,2) + t278 * t264 + t279 * t265 + t175;
t353 = -t171 * t327 + t368 * t174;
t151 = -mrSges(4,1) * t257 + mrSges(4,3) * t228 + Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * t307 - pkin(3) * t341 + pkin(9) * t353 + t368 * t162 + t327 * t163 - t309 * t268 + t322 * t270;
t152 = mrSges(4,2) * t257 - mrSges(4,3) * t227 + Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * t307 - pkin(9) * t167 - t327 * t162 + t368 * t163 + t308 * t268 - t322 * t269;
t280 = -mrSges(4,1) * t308 + mrSges(4,2) * t309;
t281 = -mrSges(4,2) * t322 + mrSges(4,3) * t308;
t165 = m(4) * t227 + mrSges(4,1) * t307 - mrSges(4,3) * t276 - t280 * t309 + t281 * t322 + t167;
t282 = mrSges(4,1) * t322 - mrSges(4,3) * t309;
t166 = m(4) * t228 - mrSges(4,2) * t307 + mrSges(4,3) * t275 + t280 * t308 - t282 * t322 + t353;
t161 = -t165 * t328 + t331 * t166;
t169 = -m(4) * t257 + t275 * mrSges(4,1) - t276 * mrSges(4,2) + t308 * t281 - t309 * t282 - t341;
t299 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t329 + Ifges(3,2) * t332) * qJD(1);
t300 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t329 + Ifges(3,4) * t332) * qJD(1);
t371 = mrSges(3,1) * t283 - mrSges(3,2) * t284 + Ifges(3,5) * t312 + Ifges(3,6) * t313 + Ifges(3,3) * qJDD(2) + pkin(2) * t169 + pkin(8) * t161 + t331 * t151 + t328 * t152 + (t299 * t329 - t300 * t332) * qJD(1);
t310 = (-mrSges(3,1) * t332 + mrSges(3,2) * t329) * qJD(1);
t315 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t360;
t159 = m(3) * t284 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t313 - qJD(2) * t315 + t310 * t359 + t161;
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t168 = m(3) * t283 + qJDD(2) * mrSges(3,1) - t312 * mrSges(3,3) + qJD(2) * t316 - t310 * t360 + t169;
t354 = t332 * t159 - t168 * t329;
t160 = t165 * t331 + t166 * t328;
t298 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t329 + Ifges(3,6) * t332) * qJD(1);
t148 = mrSges(3,2) * t301 - mrSges(3,3) * t283 + Ifges(3,1) * t312 + Ifges(3,4) * t313 + Ifges(3,5) * qJDD(2) - pkin(8) * t160 - qJD(2) * t299 - t151 * t328 + t152 * t331 + t298 * t359;
t150 = -mrSges(3,1) * t301 + mrSges(3,3) * t284 + Ifges(3,4) * t312 + Ifges(3,2) * t313 + Ifges(3,6) * qJDD(2) - pkin(2) * t160 + qJD(2) * t300 - t298 * t360 - t373;
t343 = -m(3) * t301 + t313 * mrSges(3,1) - mrSges(3,2) * t312 - t315 * t360 + t316 * t359 - t160;
t347 = mrSges(2,1) * t317 - mrSges(2,2) * t318 + Ifges(2,3) * qJDD(1) + pkin(1) * t343 + pkin(7) * t354 + t329 * t148 + t332 * t150;
t156 = m(2) * t317 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t335 + t343;
t155 = t159 * t329 + t168 * t332;
t153 = m(2) * t318 - mrSges(2,1) * t335 - qJDD(1) * mrSges(2,2) + t354;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t318 + t335 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t155 - t371;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t317 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t335 - pkin(7) * t155 + t148 * t332 - t150 * t329;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t333 * t145 - t330 * t146 - pkin(6) * (t153 * t330 + t156 * t333), t145, t148, t152, t163, -t278 * t236 - t340, -t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t330 * t145 + t333 * t146 + pkin(6) * (t153 * t333 - t156 * t330), t146, t150, t151, t162, Ifges(6,4) * t303 - Ifges(6,2) * t226 + Ifges(6,6) * t225 + t320 * t233 + t363 * t278 + t344, t320 * t232 - t349; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t347, t347, t371, t373, t374, Ifges(6,5) * t303 + t279 * t239 + Ifges(6,3) * t225 - Ifges(6,6) * t226 + (t232 - t236) * t320 + t342, -t278 * t238 - t348;];
m_new  = t1;
