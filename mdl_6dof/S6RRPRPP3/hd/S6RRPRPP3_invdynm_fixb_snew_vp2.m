% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:29:06
% EndTime: 2019-05-06 12:29:24
% DurationCPUTime: 8.37s
% Computational Cost: add. (125791->384), mult. (270001->444), div. (0->0), fcn. (183865->8), ass. (0->138)
t327 = sin(pkin(9));
t328 = cos(pkin(9));
t330 = sin(qJ(2));
t360 = qJD(1) * t330;
t307 = qJD(2) * t328 - t327 * t360;
t308 = qJD(2) * t327 + t328 * t360;
t329 = sin(qJ(4));
t368 = cos(qJ(4));
t274 = -t368 * t307 + t308 * t329;
t332 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t356 = t332 * t358;
t313 = qJDD(1) * t330 + t356;
t283 = qJDD(2) * t328 - t313 * t327;
t284 = qJDD(2) * t327 + t313 * t328;
t231 = -t274 * qJD(4) + t329 * t283 + t368 * t284;
t275 = t329 * t307 + t368 * t308;
t246 = -mrSges(7,2) * t275 + mrSges(7,3) * t274;
t331 = sin(qJ(1));
t333 = cos(qJ(1));
t318 = g(1) * t331 - t333 * g(2);
t335 = qJD(1) ^ 2;
t295 = -qJDD(1) * pkin(1) - pkin(7) * t335 - t318;
t323 = t330 * t358;
t314 = qJDD(1) * t332 - t323;
t252 = (-t313 - t356) * qJ(3) + (-t314 + t323) * pkin(2) + t295;
t319 = -g(1) * t333 - g(2) * t331;
t296 = -pkin(1) * t335 + qJDD(1) * pkin(7) + t319;
t279 = -g(3) * t330 + t332 * t296;
t311 = (-pkin(2) * t332 - qJ(3) * t330) * qJD(1);
t334 = qJD(2) ^ 2;
t359 = qJD(1) * t332;
t258 = -pkin(2) * t334 + qJDD(2) * qJ(3) + t311 * t359 + t279;
t204 = -0.2e1 * qJD(3) * t308 + t328 * t252 - t258 * t327;
t200 = (-t307 * t359 - t284) * pkin(8) + (t307 * t308 - t314) * pkin(3) + t204;
t205 = 0.2e1 * qJD(3) * t307 + t327 * t252 + t328 * t258;
t285 = -pkin(3) * t359 - pkin(8) * t308;
t306 = t307 ^ 2;
t203 = -pkin(3) * t306 + pkin(8) * t283 + t285 * t359 + t205;
t197 = t368 * t200 - t329 * t203;
t247 = pkin(4) * t274 - qJ(5) * t275;
t310 = qJDD(4) - t314;
t321 = -qJD(4) + t359;
t320 = t321 ^ 2;
t193 = -t310 * pkin(4) - t320 * qJ(5) + t275 * t247 + qJDD(5) - t197;
t364 = t274 * t321;
t369 = 2 * qJD(6);
t184 = t321 * t369 + (t274 * t275 - t310) * qJ(6) + (t231 - t364) * pkin(5) + t193;
t264 = -mrSges(7,1) * t274 - mrSges(7,2) * t321;
t352 = -m(7) * t184 + t310 * mrSges(7,3) - t321 * t264;
t180 = mrSges(7,1) * t231 + t246 * t275 - t352;
t198 = t329 * t200 + t368 * t203;
t230 = qJD(4) * t275 - t368 * t283 + t284 * t329;
t237 = Ifges(5,4) * t275 - Ifges(5,2) * t274 - Ifges(5,6) * t321;
t263 = mrSges(6,1) * t274 + mrSges(6,3) * t321;
t346 = -pkin(4) * t320 + qJ(5) * t310 - t247 * t274 + t198;
t191 = 0.2e1 * qJD(5) * t321 - t346;
t233 = -Ifges(6,5) * t321 - Ifges(6,6) * t275 + Ifges(6,3) * t274;
t261 = pkin(5) * t275 + qJ(6) * t321;
t273 = t274 ^ 2;
t370 = -0.2e1 * qJD(5);
t187 = -pkin(5) * t230 - qJ(6) * t273 + qJDD(6) + (t370 - t261) * t321 + t346;
t232 = -Ifges(7,5) * t321 + Ifges(7,6) * t274 + Ifges(7,3) * t275;
t235 = -Ifges(7,4) * t321 + Ifges(7,2) * t274 + Ifges(7,6) * t275;
t345 = -mrSges(7,2) * t187 + mrSges(7,3) * t184 - Ifges(7,1) * t310 - Ifges(7,4) * t230 - Ifges(7,5) * t231 - t274 * t232 + t275 * t235;
t340 = -mrSges(6,2) * t193 + mrSges(6,3) * t191 - Ifges(6,1) * t310 + Ifges(6,4) * t231 - Ifges(6,5) * t230 + qJ(6) * t180 + t275 * t233 + t345;
t265 = mrSges(6,1) * t275 - mrSges(6,2) * t321;
t262 = mrSges(7,1) * t275 + mrSges(7,3) * t321;
t357 = m(7) * t187 + t310 * mrSges(7,2) - t321 * t262;
t350 = -m(6) * t191 + t310 * mrSges(6,3) - t321 * t265 + t357;
t249 = -mrSges(6,2) * t274 - mrSges(6,3) * t275;
t361 = -t246 - t249;
t236 = -Ifges(6,4) * t321 - Ifges(6,2) * t275 + Ifges(6,6) * t274;
t362 = Ifges(5,1) * t275 - Ifges(5,4) * t274 - Ifges(5,5) * t321 - t236;
t372 = -m(6) * t193 - t231 * mrSges(6,1) - t275 * t249;
t374 = t362 * t274 - mrSges(5,2) * t198 + pkin(4) * (-mrSges(6,2) * t310 + t263 * t321 - t180 + t372) + qJ(5) * (t361 * t274 + (-mrSges(6,1) - mrSges(7,1)) * t230 + t350) + mrSges(5,1) * t197 + t275 * t237 - Ifges(5,6) * t230 + Ifges(5,5) * t231 + Ifges(5,3) * t310 - t340;
t248 = mrSges(5,1) * t274 + mrSges(5,2) * t275;
t259 = mrSges(5,2) * t321 - mrSges(5,3) * t274;
t366 = -mrSges(7,1) - mrSges(5,3);
t169 = m(5) * t197 + (-t259 + t263) * t321 + (mrSges(5,1) - mrSges(6,2)) * t310 + (-t246 - t248) * t275 + t366 * t231 + t352 + t372;
t260 = -mrSges(5,1) * t321 - mrSges(5,3) * t275;
t172 = m(5) * t198 - mrSges(5,2) * t310 + t260 * t321 + (-t248 + t361) * t274 + (-mrSges(6,1) + t366) * t230 + t350;
t167 = t368 * t169 + t329 * t172;
t269 = Ifges(4,4) * t308 + Ifges(4,2) * t307 - Ifges(4,6) * t359;
t270 = Ifges(4,1) * t308 + Ifges(4,4) * t307 - Ifges(4,5) * t359;
t373 = -mrSges(4,1) * t204 + mrSges(4,2) * t205 - Ifges(4,5) * t284 - Ifges(4,6) * t283 - pkin(3) * t167 - t308 * t269 + t307 * t270 - t374;
t278 = -t332 * g(3) - t330 * t296;
t257 = -qJDD(2) * pkin(2) - qJ(3) * t334 + t311 * t360 + qJDD(3) - t278;
t206 = -pkin(3) * t283 - pkin(8) * t306 + t308 * t285 + t257;
t338 = (-t231 - t364) * qJ(5) + t206 + (-t321 * pkin(4) + t370) * t275;
t190 = t338 + (pkin(4) + qJ(6)) * t230 - t261 * t275 - pkin(5) * t273 + t274 * t369;
t181 = m(7) * t190 - t231 * mrSges(7,2) + t230 * mrSges(7,3) - t275 * t262 + t274 * t264;
t195 = pkin(4) * t230 + t338;
t179 = m(6) * t195 - t230 * mrSges(6,2) - t231 * mrSges(6,3) - t274 * t263 - t275 * t265 + t181;
t234 = Ifges(5,5) * t275 - Ifges(5,6) * t274 - Ifges(5,3) * t321;
t239 = -Ifges(6,1) * t321 - Ifges(6,4) * t275 + Ifges(6,5) * t274;
t238 = -Ifges(7,1) * t321 + Ifges(7,4) * t274 + Ifges(7,5) * t275;
t349 = mrSges(7,1) * t187 - mrSges(7,3) * t190 - Ifges(7,4) * t310 - Ifges(7,2) * t230 - Ifges(7,6) * t231 - t275 * t238;
t342 = mrSges(6,1) * t191 - mrSges(6,2) * t195 + pkin(5) * (mrSges(7,1) * t230 + t246 * t274 - t357) + qJ(6) * t181 - t349;
t365 = Ifges(5,4) + Ifges(6,6);
t162 = (-t232 - t362) * t321 + (Ifges(5,6) - Ifges(6,5)) * t310 + (-t234 - t239) * t275 + t365 * t231 + (-Ifges(5,2) - Ifges(6,3)) * t230 - mrSges(5,1) * t206 + mrSges(5,3) * t198 - pkin(4) * t179 - t342;
t348 = -mrSges(7,1) * t184 + mrSges(7,2) * t190 - Ifges(7,5) * t310 - Ifges(7,6) * t230 - Ifges(7,3) * t231 + t321 * t235;
t344 = -mrSges(6,1) * t193 + mrSges(6,3) * t195 - pkin(5) * t180 + t348;
t363 = t238 + t239;
t163 = (t237 - t233) * t321 + (Ifges(5,5) - Ifges(6,4)) * t310 + (-t234 - t363) * t274 + (Ifges(5,1) + Ifges(6,2)) * t231 - t365 * t230 - t344 + mrSges(5,2) * t206 - mrSges(5,3) * t197 - qJ(5) * t179;
t268 = Ifges(4,5) * t308 + Ifges(4,6) * t307 - Ifges(4,3) * t359;
t341 = m(5) * t206 + t230 * mrSges(5,1) + t231 * mrSges(5,2) + t274 * t259 + t275 * t260 + t179;
t353 = -t169 * t329 + t368 * t172;
t151 = -mrSges(4,1) * t257 + mrSges(4,3) * t205 + Ifges(4,4) * t284 + Ifges(4,2) * t283 - Ifges(4,6) * t314 - pkin(3) * t341 + pkin(8) * t353 + t368 * t162 + t329 * t163 - t308 * t268 - t270 * t359;
t152 = mrSges(4,2) * t257 - mrSges(4,3) * t204 + Ifges(4,1) * t284 + Ifges(4,4) * t283 - Ifges(4,5) * t314 - pkin(8) * t167 - t329 * t162 + t368 * t163 + t307 * t268 + t269 * t359;
t276 = -mrSges(4,1) * t307 + mrSges(4,2) * t308;
t281 = mrSges(4,2) * t359 + mrSges(4,3) * t307;
t165 = m(4) * t204 - mrSges(4,1) * t314 - mrSges(4,3) * t284 - t276 * t308 - t281 * t359 + t167;
t282 = -mrSges(4,1) * t359 - mrSges(4,3) * t308;
t166 = m(4) * t205 + mrSges(4,2) * t314 + mrSges(4,3) * t283 + t276 * t307 + t282 * t359 + t353;
t161 = -t165 * t327 + t328 * t166;
t174 = -m(4) * t257 + t283 * mrSges(4,1) - t284 * mrSges(4,2) + t307 * t281 - t308 * t282 - t341;
t293 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t330 + Ifges(3,2) * t332) * qJD(1);
t294 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t330 + Ifges(3,4) * t332) * qJD(1);
t371 = mrSges(3,1) * t278 - mrSges(3,2) * t279 + Ifges(3,5) * t313 + Ifges(3,6) * t314 + Ifges(3,3) * qJDD(2) + pkin(2) * t174 + qJ(3) * t161 + t328 * t151 + t327 * t152 + (t293 * t330 - t294 * t332) * qJD(1);
t312 = (-mrSges(3,1) * t332 + mrSges(3,2) * t330) * qJD(1);
t316 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t360;
t159 = m(3) * t279 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t314 - qJD(2) * t316 + t312 * t359 + t161;
t317 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t173 = m(3) * t278 + qJDD(2) * mrSges(3,1) - t313 * mrSges(3,3) + qJD(2) * t317 - t312 * t360 + t174;
t354 = t332 * t159 - t173 * t330;
t160 = t165 * t328 + t166 * t327;
t292 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t330 + Ifges(3,6) * t332) * qJD(1);
t148 = mrSges(3,2) * t295 - mrSges(3,3) * t278 + Ifges(3,1) * t313 + Ifges(3,4) * t314 + Ifges(3,5) * qJDD(2) - qJ(3) * t160 - qJD(2) * t293 - t151 * t327 + t152 * t328 + t292 * t359;
t150 = Ifges(3,6) * qJDD(2) + Ifges(3,4) * t313 + qJD(2) * t294 - mrSges(3,1) * t295 + mrSges(3,3) * t279 - pkin(2) * t160 + (Ifges(3,2) + Ifges(4,3)) * t314 - t292 * t360 + t373;
t343 = -m(3) * t295 + t314 * mrSges(3,1) - mrSges(3,2) * t313 - t316 * t360 + t317 * t359 - t160;
t347 = mrSges(2,1) * t318 - mrSges(2,2) * t319 + Ifges(2,3) * qJDD(1) + pkin(1) * t343 + pkin(7) * t354 + t330 * t148 + t332 * t150;
t156 = m(2) * t318 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t335 + t343;
t155 = t159 * t330 + t173 * t332;
t153 = m(2) * t319 - mrSges(2,1) * t335 - qJDD(1) * mrSges(2,2) + t354;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t319 + t335 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t155 - t371;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t318 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t335 - pkin(7) * t155 + t148 * t332 - t150 * t330;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t333 * t145 - t331 * t146 - pkin(6) * (t153 * t331 + t156 * t333), t145, t148, t152, t163, -t274 * t236 - t340, -t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t331 * t145 + t333 * t146 + pkin(6) * (t153 * t333 - t156 * t331), t146, t150, t151, t162, Ifges(6,4) * t310 - Ifges(6,2) * t231 + Ifges(6,6) * t230 + t321 * t233 + t363 * t274 + t344, t321 * t232 - t349; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t347, t347, t371, -Ifges(4,3) * t314 - t373, t374, Ifges(6,5) * t310 + t275 * t239 - Ifges(6,6) * t231 + Ifges(6,3) * t230 + t342 + (t232 - t236) * t321, -t274 * t238 - t348;];
m_new  = t1;
