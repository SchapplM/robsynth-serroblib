% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:46:40
% EndTime: 2019-05-06 03:47:43
% DurationCPUTime: 42.73s
% Computational Cost: add. (738816->366), mult. (1756575->457), div. (0->0), fcn. (1366244->12), ass. (0->157)
t341 = qJD(1) ^ 2;
t328 = sin(pkin(11));
t371 = qJD(1) * t328;
t329 = cos(pkin(11));
t370 = qJD(1) * t329;
t334 = sin(qJ(1));
t339 = cos(qJ(1));
t315 = -g(1) * t339 - g(2) * t334;
t308 = -pkin(1) * t341 + qJDD(1) * qJ(2) + t315;
t368 = qJD(1) * qJD(2);
t365 = -t329 * g(3) - 0.2e1 * t328 * t368;
t374 = pkin(2) * t329;
t274 = (-pkin(7) * qJDD(1) + t341 * t374 - t308) * t328 + t365;
t294 = -g(3) * t328 + (t308 + 0.2e1 * t368) * t329;
t366 = qJDD(1) * t329;
t325 = t329 ^ 2;
t372 = t325 * t341;
t275 = -pkin(2) * t372 + pkin(7) * t366 + t294;
t333 = sin(qJ(3));
t338 = cos(qJ(3));
t250 = t333 * t274 + t338 * t275;
t306 = -t333 * t371 + t338 * t370;
t353 = t328 * t338 + t329 * t333;
t307 = t353 * qJD(1);
t284 = -mrSges(4,1) * t306 + mrSges(4,2) * t307;
t303 = t307 * qJD(3);
t367 = qJDD(1) * t328;
t291 = -t333 * t367 + t338 * t366 - t303;
t299 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t307;
t289 = -pkin(3) * t306 - pkin(8) * t307;
t340 = qJD(3) ^ 2;
t235 = -pkin(3) * t340 + qJDD(3) * pkin(8) + t289 * t306 + t250;
t324 = t328 ^ 2;
t314 = t334 * g(1) - t339 * g(2);
t359 = qJDD(2) - t314;
t290 = (-pkin(1) - t374) * qJDD(1) + (-qJ(2) + (-t324 - t325) * pkin(7)) * t341 + t359;
t369 = t306 * qJD(3);
t292 = t353 * qJDD(1) + t369;
t240 = (-t292 - t369) * pkin(8) + (-t291 + t303) * pkin(3) + t290;
t332 = sin(qJ(4));
t337 = cos(qJ(4));
t219 = -t332 * t235 + t337 * t240;
t296 = qJD(3) * t337 - t307 * t332;
t261 = qJD(4) * t296 + qJDD(3) * t332 + t292 * t337;
t288 = qJDD(4) - t291;
t297 = qJD(3) * t332 + t307 * t337;
t304 = qJD(4) - t306;
t214 = (t296 * t304 - t261) * pkin(9) + (t296 * t297 + t288) * pkin(4) + t219;
t220 = t337 * t235 + t332 * t240;
t260 = -qJD(4) * t297 + qJDD(3) * t337 - t292 * t332;
t273 = pkin(4) * t304 - pkin(9) * t297;
t295 = t296 ^ 2;
t216 = -pkin(4) * t295 + pkin(9) * t260 - t273 * t304 + t220;
t331 = sin(qJ(5));
t336 = cos(qJ(5));
t202 = t336 * t214 - t331 * t216;
t263 = t296 * t336 - t297 * t331;
t231 = qJD(5) * t263 + t260 * t331 + t261 * t336;
t264 = t296 * t331 + t297 * t336;
t283 = qJDD(5) + t288;
t302 = qJD(5) + t304;
t199 = (t263 * t302 - t231) * pkin(10) + (t263 * t264 + t283) * pkin(5) + t202;
t203 = t331 * t214 + t336 * t216;
t230 = -qJD(5) * t264 + t260 * t336 - t261 * t331;
t253 = pkin(5) * t302 - pkin(10) * t264;
t262 = t263 ^ 2;
t200 = -pkin(5) * t262 + pkin(10) * t230 - t253 * t302 + t203;
t330 = sin(qJ(6));
t335 = cos(qJ(6));
t197 = t199 * t335 - t200 * t330;
t245 = t263 * t335 - t264 * t330;
t211 = qJD(6) * t245 + t230 * t330 + t231 * t335;
t246 = t263 * t330 + t264 * t335;
t226 = -mrSges(7,1) * t245 + mrSges(7,2) * t246;
t300 = qJD(6) + t302;
t238 = -mrSges(7,2) * t300 + mrSges(7,3) * t245;
t281 = qJDD(6) + t283;
t192 = m(7) * t197 + mrSges(7,1) * t281 - mrSges(7,3) * t211 - t226 * t246 + t238 * t300;
t198 = t199 * t330 + t200 * t335;
t210 = -qJD(6) * t246 + t230 * t335 - t231 * t330;
t239 = mrSges(7,1) * t300 - mrSges(7,3) * t246;
t193 = m(7) * t198 - mrSges(7,2) * t281 + mrSges(7,3) * t210 + t226 * t245 - t239 * t300;
t184 = t335 * t192 + t330 * t193;
t247 = -mrSges(6,1) * t263 + mrSges(6,2) * t264;
t251 = -mrSges(6,2) * t302 + mrSges(6,3) * t263;
t181 = m(6) * t202 + mrSges(6,1) * t283 - mrSges(6,3) * t231 - t247 * t264 + t251 * t302 + t184;
t252 = mrSges(6,1) * t302 - mrSges(6,3) * t264;
t360 = -t192 * t330 + t335 * t193;
t182 = m(6) * t203 - mrSges(6,2) * t283 + mrSges(6,3) * t230 + t247 * t263 - t252 * t302 + t360;
t177 = t336 * t181 + t331 * t182;
t266 = -mrSges(5,1) * t296 + mrSges(5,2) * t297;
t269 = -mrSges(5,2) * t304 + mrSges(5,3) * t296;
t175 = m(5) * t219 + mrSges(5,1) * t288 - mrSges(5,3) * t261 - t266 * t297 + t269 * t304 + t177;
t270 = mrSges(5,1) * t304 - mrSges(5,3) * t297;
t361 = -t181 * t331 + t336 * t182;
t176 = m(5) * t220 - mrSges(5,2) * t288 + mrSges(5,3) * t260 + t266 * t296 - t270 * t304 + t361;
t362 = -t175 * t332 + t337 * t176;
t166 = m(4) * t250 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t291 - qJD(3) * t299 + t284 * t306 + t362;
t249 = t274 * t338 - t333 * t275;
t298 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t306;
t234 = -qJDD(3) * pkin(3) - pkin(8) * t340 + t307 * t289 - t249;
t218 = -pkin(4) * t260 - pkin(9) * t295 + t297 * t273 + t234;
t205 = -pkin(5) * t230 - pkin(10) * t262 + t253 * t264 + t218;
t355 = m(7) * t205 - t210 * mrSges(7,1) + t211 * mrSges(7,2) - t245 * t238 + t246 * t239;
t347 = m(6) * t218 - t230 * mrSges(6,1) + mrSges(6,2) * t231 - t263 * t251 + t252 * t264 + t355;
t343 = -m(5) * t234 + t260 * mrSges(5,1) - mrSges(5,2) * t261 + t296 * t269 - t270 * t297 - t347;
t188 = m(4) * t249 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t292 + qJD(3) * t298 - t284 * t307 + t343;
t161 = t333 * t166 + t338 * t188;
t293 = -t328 * t308 + t365;
t221 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t300;
t223 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t300;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t211 + Ifges(7,2) * t210 + Ifges(7,6) * t281 - t221 * t246 + t223 * t300;
t222 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t300;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t211 + Ifges(7,4) * t210 + Ifges(7,5) * t281 + t221 * t245 - t222 * t300;
t241 = Ifges(6,5) * t264 + Ifges(6,6) * t263 + Ifges(6,3) * t302;
t243 = Ifges(6,1) * t264 + Ifges(6,4) * t263 + Ifges(6,5) * t302;
t170 = -mrSges(6,1) * t218 + mrSges(6,3) * t203 + Ifges(6,4) * t231 + Ifges(6,2) * t230 + Ifges(6,6) * t283 - pkin(5) * t355 + pkin(10) * t360 + t335 * t185 + t330 * t186 - t264 * t241 + t302 * t243;
t242 = Ifges(6,4) * t264 + Ifges(6,2) * t263 + Ifges(6,6) * t302;
t171 = mrSges(6,2) * t218 - mrSges(6,3) * t202 + Ifges(6,1) * t231 + Ifges(6,4) * t230 + Ifges(6,5) * t283 - pkin(10) * t184 - t185 * t330 + t186 * t335 + t241 * t263 - t242 * t302;
t254 = Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t304;
t256 = Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t304;
t155 = -mrSges(5,1) * t234 + mrSges(5,3) * t220 + Ifges(5,4) * t261 + Ifges(5,2) * t260 + Ifges(5,6) * t288 - pkin(4) * t347 + pkin(9) * t361 + t336 * t170 + t331 * t171 - t297 * t254 + t304 * t256;
t255 = Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t304;
t157 = mrSges(5,2) * t234 - mrSges(5,3) * t219 + Ifges(5,1) * t261 + Ifges(5,4) * t260 + Ifges(5,5) * t288 - pkin(9) * t177 - t170 * t331 + t171 * t336 + t254 * t296 - t255 * t304;
t278 = Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * qJD(3);
t279 = Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * qJD(3);
t348 = -mrSges(4,1) * t249 + mrSges(4,2) * t250 - Ifges(4,5) * t292 - Ifges(4,6) * t291 - Ifges(4,3) * qJDD(3) - pkin(3) * t343 - pkin(8) * t362 - t337 * t155 - t332 * t157 - t307 * t278 + t306 * t279;
t357 = Ifges(3,4) * t328 + Ifges(3,2) * t329;
t358 = Ifges(3,1) * t328 + Ifges(3,4) * t329;
t375 = -mrSges(3,1) * t293 + mrSges(3,2) * t294 - pkin(2) * t161 - (t357 * t371 - t358 * t370) * qJD(1) + t348;
t373 = mrSges(3,2) * t328;
t168 = t337 * t175 + t332 * t176;
t352 = mrSges(3,3) * qJDD(1) + t341 * (-mrSges(3,1) * t329 + t373);
t159 = m(3) * t293 - t352 * t328 + t161;
t363 = t338 * t166 - t333 * t188;
t160 = m(3) * t294 + t352 * t329 + t363;
t364 = -t159 * t328 + t329 * t160;
t356 = Ifges(3,5) * t328 + Ifges(3,6) * t329;
t277 = Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * qJD(3);
t149 = mrSges(4,2) * t290 - mrSges(4,3) * t249 + Ifges(4,1) * t292 + Ifges(4,4) * t291 + Ifges(4,5) * qJDD(3) - pkin(8) * t168 - qJD(3) * t278 - t155 * t332 + t157 * t337 + t277 * t306;
t350 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t211 - Ifges(7,6) * t210 - Ifges(7,3) * t281 - t246 * t222 + t245 * t223;
t345 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t231 - Ifges(6,6) * t230 - Ifges(6,3) * t283 - pkin(5) * t184 - t264 * t242 + t263 * t243 + t350;
t342 = mrSges(5,1) * t219 - mrSges(5,2) * t220 + Ifges(5,5) * t261 + Ifges(5,6) * t260 + Ifges(5,3) * t288 + pkin(4) * t177 + t297 * t255 - t296 * t256 - t345;
t153 = -mrSges(4,1) * t290 + mrSges(4,3) * t250 + Ifges(4,4) * t292 + Ifges(4,2) * t291 + Ifges(4,6) * qJDD(3) - pkin(3) * t168 + qJD(3) * t279 - t307 * t277 - t342;
t305 = -qJDD(1) * pkin(1) - t341 * qJ(2) + t359;
t310 = t356 * qJD(1);
t349 = m(4) * t290 - t291 * mrSges(4,1) + t292 * mrSges(4,2) - t306 * t298 + t307 * t299 + t168;
t145 = -mrSges(3,1) * t305 + mrSges(3,3) * t294 - pkin(2) * t349 + pkin(7) * t363 + t357 * qJDD(1) + t333 * t149 + t338 * t153 - t310 * t371;
t148 = mrSges(3,2) * t305 - mrSges(3,3) * t293 - pkin(7) * t161 + t358 * qJDD(1) + t338 * t149 - t333 * t153 + t310 * t370;
t346 = -m(3) * t305 + mrSges(3,1) * t366 - t349 + (t324 * t341 + t372) * mrSges(3,3);
t351 = -mrSges(2,2) * t315 + qJ(2) * t364 + t329 * t145 + t328 * t148 + pkin(1) * (-mrSges(3,2) * t367 + t346) + mrSges(2,1) * t314 + Ifges(2,3) * qJDD(1);
t162 = -t341 * mrSges(2,2) + m(2) * t314 + t346 + (mrSges(2,1) - t373) * qJDD(1);
t152 = t159 * t329 + t160 * t328;
t150 = m(2) * t315 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t364;
t146 = mrSges(2,1) * g(3) + (Ifges(2,6) - t356) * qJDD(1) + t341 * Ifges(2,5) + mrSges(2,3) * t315 - pkin(1) * t152 + t375;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t314 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t341 - qJ(2) * t152 - t145 * t328 + t148 * t329;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t143 - t334 * t146 - pkin(6) * (t150 * t334 + t162 * t339), t143, t148, t149, t157, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t143 + t339 * t146 + pkin(6) * (t150 * t339 - t162 * t334), t146, t145, t153, t155, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t351, t351, t356 * qJDD(1) - t375, -t348, t342, -t345, -t350;];
m_new  = t1;
