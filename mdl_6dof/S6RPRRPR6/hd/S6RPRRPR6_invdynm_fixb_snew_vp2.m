% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:48:20
% EndTime: 2019-05-05 22:49:15
% DurationCPUTime: 40.50s
% Computational Cost: add. (707131->365), mult. (1705990->459), div. (0->0), fcn. (1316306->12), ass. (0->155)
t341 = qJD(1) ^ 2;
t329 = sin(pkin(10));
t371 = qJD(1) * t329;
t331 = cos(pkin(10));
t370 = qJD(1) * t331;
t335 = sin(qJ(1));
t339 = cos(qJ(1));
t315 = -g(1) * t339 - g(2) * t335;
t308 = -pkin(1) * t341 + qJDD(1) * qJ(2) + t315;
t368 = qJD(1) * qJD(2);
t365 = -t331 * g(3) - 0.2e1 * t329 * t368;
t374 = pkin(2) * t331;
t276 = (-pkin(7) * qJDD(1) + t341 * t374 - t308) * t329 + t365;
t295 = -g(3) * t329 + (t308 + 0.2e1 * t368) * t331;
t366 = qJDD(1) * t331;
t325 = t331 ^ 2;
t372 = t325 * t341;
t277 = -pkin(2) * t372 + pkin(7) * t366 + t295;
t334 = sin(qJ(3));
t338 = cos(qJ(3));
t250 = t334 * t276 + t338 * t277;
t306 = -t334 * t371 + t338 * t370;
t353 = t329 * t338 + t331 * t334;
t307 = t353 * qJD(1);
t285 = -mrSges(4,1) * t306 + mrSges(4,2) * t307;
t303 = t307 * qJD(3);
t367 = qJDD(1) * t329;
t292 = -t334 * t367 + t338 * t366 - t303;
t300 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t307;
t290 = -pkin(3) * t306 - pkin(8) * t307;
t340 = qJD(3) ^ 2;
t230 = -pkin(3) * t340 + qJDD(3) * pkin(8) + t290 * t306 + t250;
t324 = t329 ^ 2;
t314 = t335 * g(1) - t339 * g(2);
t359 = qJDD(2) - t314;
t291 = (-pkin(1) - t374) * qJDD(1) + (-qJ(2) + (-t324 - t325) * pkin(7)) * t341 + t359;
t369 = t306 * qJD(3);
t293 = t353 * qJDD(1) + t369;
t240 = (-t293 - t369) * pkin(8) + (-t292 + t303) * pkin(3) + t291;
t333 = sin(qJ(4));
t337 = cos(qJ(4));
t219 = -t333 * t230 + t337 * t240;
t297 = qJD(3) * t337 - t307 * t333;
t261 = qJD(4) * t297 + qJDD(3) * t333 + t293 * t337;
t289 = qJDD(4) - t292;
t298 = qJD(3) * t333 + t307 * t337;
t304 = qJD(4) - t306;
t208 = (t297 * t304 - t261) * qJ(5) + (t297 * t298 + t289) * pkin(4) + t219;
t220 = t337 * t230 + t333 * t240;
t260 = -qJD(4) * t298 + qJDD(3) * t337 - t293 * t333;
t272 = pkin(4) * t304 - qJ(5) * t298;
t296 = t297 ^ 2;
t210 = -pkin(4) * t296 + qJ(5) * t260 - t272 * t304 + t220;
t328 = sin(pkin(11));
t330 = cos(pkin(11));
t266 = t297 * t328 + t298 * t330;
t202 = -0.2e1 * qJD(5) * t266 + t330 * t208 - t328 * t210;
t235 = t260 * t328 + t261 * t330;
t265 = t297 * t330 - t298 * t328;
t199 = (t265 * t304 - t235) * pkin(9) + (t265 * t266 + t289) * pkin(5) + t202;
t203 = 0.2e1 * qJD(5) * t265 + t328 * t208 + t330 * t210;
t234 = t260 * t330 - t261 * t328;
t253 = pkin(5) * t304 - pkin(9) * t266;
t264 = t265 ^ 2;
t200 = -pkin(5) * t264 + pkin(9) * t234 - t253 * t304 + t203;
t332 = sin(qJ(6));
t336 = cos(qJ(6));
t197 = t199 * t336 - t200 * t332;
t245 = t265 * t336 - t266 * t332;
t216 = qJD(6) * t245 + t234 * t332 + t235 * t336;
t246 = t265 * t332 + t266 * t336;
t226 = -mrSges(7,1) * t245 + mrSges(7,2) * t246;
t302 = qJD(6) + t304;
t238 = -mrSges(7,2) * t302 + mrSges(7,3) * t245;
t284 = qJDD(6) + t289;
t190 = m(7) * t197 + mrSges(7,1) * t284 - mrSges(7,3) * t216 - t226 * t246 + t238 * t302;
t198 = t199 * t332 + t200 * t336;
t215 = -qJD(6) * t246 + t234 * t336 - t235 * t332;
t239 = mrSges(7,1) * t302 - mrSges(7,3) * t246;
t191 = m(7) * t198 - mrSges(7,2) * t284 + mrSges(7,3) * t215 + t226 * t245 - t239 * t302;
t184 = t336 * t190 + t332 * t191;
t247 = -mrSges(6,1) * t265 + mrSges(6,2) * t266;
t251 = -mrSges(6,2) * t304 + mrSges(6,3) * t265;
t181 = m(6) * t202 + mrSges(6,1) * t289 - mrSges(6,3) * t235 - t247 * t266 + t251 * t304 + t184;
t252 = mrSges(6,1) * t304 - mrSges(6,3) * t266;
t360 = -t190 * t332 + t336 * t191;
t182 = m(6) * t203 - mrSges(6,2) * t289 + mrSges(6,3) * t234 + t247 * t265 - t252 * t304 + t360;
t177 = t330 * t181 + t328 * t182;
t268 = -mrSges(5,1) * t297 + mrSges(5,2) * t298;
t271 = -mrSges(5,2) * t304 + mrSges(5,3) * t297;
t175 = m(5) * t219 + mrSges(5,1) * t289 - mrSges(5,3) * t261 - t268 * t298 + t271 * t304 + t177;
t273 = mrSges(5,1) * t304 - mrSges(5,3) * t298;
t361 = -t181 * t328 + t330 * t182;
t176 = m(5) * t220 - mrSges(5,2) * t289 + mrSges(5,3) * t260 + t268 * t297 - t273 * t304 + t361;
t362 = -t175 * t333 + t337 * t176;
t166 = m(4) * t250 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t292 - qJD(3) * t300 + t285 * t306 + t362;
t249 = t276 * t338 - t334 * t277;
t299 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t306;
t229 = -qJDD(3) * pkin(3) - pkin(8) * t340 + t307 * t290 - t249;
t218 = -pkin(4) * t260 - qJ(5) * t296 + t298 * t272 + qJDD(5) + t229;
t205 = -pkin(5) * t234 - pkin(9) * t264 + t253 * t266 + t218;
t355 = m(7) * t205 - t215 * mrSges(7,1) + t216 * mrSges(7,2) - t245 * t238 + t246 * t239;
t347 = m(6) * t218 - t234 * mrSges(6,1) + mrSges(6,2) * t235 - t265 * t251 + t252 * t266 + t355;
t343 = -m(5) * t229 + t260 * mrSges(5,1) - mrSges(5,2) * t261 + t297 * t271 - t273 * t298 - t347;
t193 = m(4) * t249 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t293 + qJD(3) * t299 - t285 * t307 + t343;
t161 = t334 * t166 + t338 * t193;
t294 = -t329 * t308 + t365;
t221 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t302;
t223 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t302;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t284 - t221 * t246 + t223 * t302;
t222 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t302;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t284 + t221 * t245 - t222 * t302;
t241 = Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t304;
t243 = Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t304;
t170 = -mrSges(6,1) * t218 + mrSges(6,3) * t203 + Ifges(6,4) * t235 + Ifges(6,2) * t234 + Ifges(6,6) * t289 - pkin(5) * t355 + pkin(9) * t360 + t336 * t185 + t332 * t186 - t266 * t241 + t304 * t243;
t242 = Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t304;
t171 = mrSges(6,2) * t218 - mrSges(6,3) * t202 + Ifges(6,1) * t235 + Ifges(6,4) * t234 + Ifges(6,5) * t289 - pkin(9) * t184 - t185 * t332 + t186 * t336 + t241 * t265 - t242 * t304;
t254 = Ifges(5,5) * t298 + Ifges(5,6) * t297 + Ifges(5,3) * t304;
t256 = Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t304;
t155 = -mrSges(5,1) * t229 + mrSges(5,3) * t220 + Ifges(5,4) * t261 + Ifges(5,2) * t260 + Ifges(5,6) * t289 - pkin(4) * t347 + qJ(5) * t361 + t330 * t170 + t328 * t171 - t298 * t254 + t304 * t256;
t255 = Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t304;
t157 = mrSges(5,2) * t229 - mrSges(5,3) * t219 + Ifges(5,1) * t261 + Ifges(5,4) * t260 + Ifges(5,5) * t289 - qJ(5) * t177 - t170 * t328 + t171 * t330 + t254 * t297 - t255 * t304;
t279 = Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * qJD(3);
t280 = Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * qJD(3);
t348 = -mrSges(4,1) * t249 + mrSges(4,2) * t250 - Ifges(4,5) * t293 - Ifges(4,6) * t292 - Ifges(4,3) * qJDD(3) - pkin(3) * t343 - pkin(8) * t362 - t337 * t155 - t333 * t157 - t307 * t279 + t306 * t280;
t357 = Ifges(3,4) * t329 + Ifges(3,2) * t331;
t358 = Ifges(3,1) * t329 + Ifges(3,4) * t331;
t375 = -mrSges(3,1) * t294 + mrSges(3,2) * t295 - pkin(2) * t161 - (t357 * t371 - t358 * t370) * qJD(1) + t348;
t373 = mrSges(3,2) * t329;
t168 = t337 * t175 + t333 * t176;
t352 = mrSges(3,3) * qJDD(1) + t341 * (-mrSges(3,1) * t331 + t373);
t159 = m(3) * t294 - t352 * t329 + t161;
t363 = t338 * t166 - t334 * t193;
t160 = m(3) * t295 + t352 * t331 + t363;
t364 = -t159 * t329 + t331 * t160;
t356 = Ifges(3,5) * t329 + Ifges(3,6) * t331;
t278 = Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * qJD(3);
t149 = mrSges(4,2) * t291 - mrSges(4,3) * t249 + Ifges(4,1) * t293 + Ifges(4,4) * t292 + Ifges(4,5) * qJDD(3) - pkin(8) * t168 - qJD(3) * t279 - t155 * t333 + t157 * t337 + t278 * t306;
t350 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t216 - Ifges(7,6) * t215 - Ifges(7,3) * t284 - t246 * t222 + t245 * t223;
t345 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t235 - Ifges(6,6) * t234 - Ifges(6,3) * t289 - pkin(5) * t184 - t266 * t242 + t265 * t243 + t350;
t342 = mrSges(5,1) * t219 - mrSges(5,2) * t220 + Ifges(5,5) * t261 + Ifges(5,6) * t260 + Ifges(5,3) * t289 + pkin(4) * t177 + t298 * t255 - t297 * t256 - t345;
t153 = -mrSges(4,1) * t291 + mrSges(4,3) * t250 + Ifges(4,4) * t293 + Ifges(4,2) * t292 + Ifges(4,6) * qJDD(3) - pkin(3) * t168 + qJD(3) * t280 - t307 * t278 - t342;
t305 = -qJDD(1) * pkin(1) - t341 * qJ(2) + t359;
t310 = t356 * qJD(1);
t349 = m(4) * t291 - t292 * mrSges(4,1) + t293 * mrSges(4,2) - t306 * t299 + t307 * t300 + t168;
t145 = -mrSges(3,1) * t305 + mrSges(3,3) * t295 - pkin(2) * t349 + pkin(7) * t363 + t357 * qJDD(1) + t334 * t149 + t338 * t153 - t310 * t371;
t148 = mrSges(3,2) * t305 - mrSges(3,3) * t294 - pkin(7) * t161 + t358 * qJDD(1) + t338 * t149 - t334 * t153 + t310 * t370;
t346 = -m(3) * t305 + mrSges(3,1) * t366 - t349 + (t324 * t341 + t372) * mrSges(3,3);
t351 = -mrSges(2,2) * t315 + qJ(2) * t364 + t331 * t145 + t329 * t148 + pkin(1) * (-mrSges(3,2) * t367 + t346) + mrSges(2,1) * t314 + Ifges(2,3) * qJDD(1);
t162 = -t341 * mrSges(2,2) + m(2) * t314 + t346 + (mrSges(2,1) - t373) * qJDD(1);
t152 = t159 * t331 + t160 * t329;
t150 = m(2) * t315 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t364;
t146 = mrSges(2,1) * g(3) + (Ifges(2,6) - t356) * qJDD(1) - pkin(1) * t152 + t341 * Ifges(2,5) + mrSges(2,3) * t315 + t375;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t314 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t341 - qJ(2) * t152 - t145 * t329 + t148 * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t143 - t335 * t146 - pkin(6) * (t150 * t335 + t162 * t339), t143, t148, t149, t157, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t143 + t339 * t146 + pkin(6) * (t150 * t339 - t162 * t335), t146, t145, t153, t155, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t351, t351, t356 * qJDD(1) - t375, -t348, t342, -t345, -t350;];
m_new  = t1;
