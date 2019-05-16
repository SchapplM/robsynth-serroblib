% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR4
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
% Datum: 2019-05-05 22:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:26:56
% EndTime: 2019-05-05 22:27:47
% DurationCPUTime: 40.49s
% Computational Cost: add. (686550->365), mult. (1672903->457), div. (0->0), fcn. (1329207->12), ass. (0->152)
t338 = qJD(1) ^ 2;
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t312 = -t337 * g(1) - t334 * g(2);
t305 = -t338 * pkin(1) + qJDD(1) * qJ(2) + t312;
t328 = sin(pkin(10));
t330 = cos(pkin(10));
t365 = qJD(1) * qJD(2);
t363 = -t330 * g(3) - 0.2e1 * t328 * t365;
t370 = pkin(2) * t330;
t276 = (-pkin(7) * qJDD(1) + t338 * t370 - t305) * t328 + t363;
t295 = -t328 * g(3) + (t305 + 0.2e1 * t365) * t330;
t364 = qJDD(1) * t330;
t323 = t330 ^ 2;
t368 = t323 * t338;
t277 = -pkin(2) * t368 + pkin(7) * t364 + t295;
t333 = sin(qJ(3));
t336 = cos(qJ(3));
t253 = t336 * t276 - t333 * t277;
t352 = t328 * t336 + t330 * t333;
t351 = -t328 * t333 + t330 * t336;
t303 = t351 * qJD(1);
t366 = t303 * qJD(3);
t293 = t352 * qJDD(1) + t366;
t304 = t352 * qJD(1);
t224 = (-t293 + t366) * pkin(8) + (t303 * t304 + qJDD(3)) * pkin(3) + t253;
t254 = t333 * t276 + t336 * t277;
t292 = -t304 * qJD(3) + t351 * qJDD(1);
t298 = qJD(3) * pkin(3) - t304 * pkin(8);
t302 = t303 ^ 2;
t227 = -t302 * pkin(3) + t292 * pkin(8) - qJD(3) * t298 + t254;
t332 = sin(qJ(4));
t371 = cos(qJ(4));
t217 = t332 * t224 + t371 * t227;
t283 = t332 * t303 + t371 * t304;
t247 = t283 * qJD(4) - t371 * t292 + t332 * t293;
t282 = -t371 * t303 + t332 * t304;
t264 = t282 * mrSges(5,1) + t283 * mrSges(5,2);
t324 = qJD(3) + qJD(4);
t274 = t324 * mrSges(5,1) - t283 * mrSges(5,3);
t321 = qJDD(3) + qJDD(4);
t263 = t282 * pkin(4) - t283 * qJ(5);
t320 = t324 ^ 2;
t205 = -t320 * pkin(4) + t321 * qJ(5) - t282 * t263 + t217;
t322 = t328 ^ 2;
t311 = t334 * g(1) - t337 * g(2);
t357 = qJDD(2) - t311;
t291 = (-pkin(1) - t370) * qJDD(1) + (-qJ(2) + (-t322 - t323) * pkin(7)) * t338 + t357;
t239 = -t292 * pkin(3) - t302 * pkin(8) + t304 * t298 + t291;
t248 = -t282 * qJD(4) + t332 * t292 + t371 * t293;
t208 = (t282 * t324 - t248) * qJ(5) + (t283 * t324 + t247) * pkin(4) + t239;
t327 = sin(pkin(11));
t329 = cos(pkin(11));
t270 = t329 * t283 + t327 * t324;
t200 = -0.2e1 * qJD(5) * t270 - t327 * t205 + t329 * t208;
t236 = t329 * t248 + t327 * t321;
t269 = -t327 * t283 + t329 * t324;
t198 = (t269 * t282 - t236) * pkin(9) + (t269 * t270 + t247) * pkin(5) + t200;
t201 = 0.2e1 * qJD(5) * t269 + t329 * t205 + t327 * t208;
t235 = -t327 * t248 + t329 * t321;
t257 = t282 * pkin(5) - t270 * pkin(9);
t268 = t269 ^ 2;
t199 = -t268 * pkin(5) + t235 * pkin(9) - t282 * t257 + t201;
t331 = sin(qJ(6));
t335 = cos(qJ(6));
t196 = t335 * t198 - t331 * t199;
t249 = t335 * t269 - t331 * t270;
t213 = t249 * qJD(6) + t331 * t235 + t335 * t236;
t250 = t331 * t269 + t335 * t270;
t222 = -t249 * mrSges(7,1) + t250 * mrSges(7,2);
t278 = qJD(6) + t282;
t228 = -t278 * mrSges(7,2) + t249 * mrSges(7,3);
t246 = qJDD(6) + t247;
t191 = m(7) * t196 + t246 * mrSges(7,1) - t213 * mrSges(7,3) - t250 * t222 + t278 * t228;
t197 = t331 * t198 + t335 * t199;
t212 = -t250 * qJD(6) + t335 * t235 - t331 * t236;
t229 = t278 * mrSges(7,1) - t250 * mrSges(7,3);
t192 = m(7) * t197 - t246 * mrSges(7,2) + t212 * mrSges(7,3) + t249 * t222 - t278 * t229;
t183 = t335 * t191 + t331 * t192;
t251 = -t269 * mrSges(6,1) + t270 * mrSges(6,2);
t255 = -t282 * mrSges(6,2) + t269 * mrSges(6,3);
t181 = m(6) * t200 + t247 * mrSges(6,1) - t236 * mrSges(6,3) - t270 * t251 + t282 * t255 + t183;
t256 = t282 * mrSges(6,1) - t270 * mrSges(6,3);
t358 = -t331 * t191 + t335 * t192;
t182 = m(6) * t201 - t247 * mrSges(6,2) + t235 * mrSges(6,3) + t269 * t251 - t282 * t256 + t358;
t359 = -t327 * t181 + t329 * t182;
t174 = m(5) * t217 - t321 * mrSges(5,2) - t247 * mrSges(5,3) - t282 * t264 - t324 * t274 + t359;
t216 = t371 * t224 - t332 * t227;
t273 = -t324 * mrSges(5,2) - t282 * mrSges(5,3);
t204 = -t321 * pkin(4) - t320 * qJ(5) + t283 * t263 + qJDD(5) - t216;
t202 = -t235 * pkin(5) - t268 * pkin(9) + t270 * t257 + t204;
t347 = m(7) * t202 - t212 * mrSges(7,1) + t213 * mrSges(7,2) - t249 * t228 + t250 * t229;
t343 = -m(6) * t204 + t235 * mrSges(6,1) - t236 * mrSges(6,2) + t269 * t255 - t270 * t256 - t347;
t187 = m(5) * t216 + t321 * mrSges(5,1) - t248 * mrSges(5,3) - t283 * t264 + t324 * t273 + t343;
t165 = t332 * t174 + t371 * t187;
t287 = -t303 * mrSges(4,1) + t304 * mrSges(4,2);
t296 = -qJD(3) * mrSges(4,2) + t303 * mrSges(4,3);
t162 = m(4) * t253 + qJDD(3) * mrSges(4,1) - t293 * mrSges(4,3) + qJD(3) * t296 - t304 * t287 + t165;
t297 = qJD(3) * mrSges(4,1) - t304 * mrSges(4,3);
t360 = t371 * t174 - t332 * t187;
t163 = m(4) * t254 - qJDD(3) * mrSges(4,2) + t292 * mrSges(4,3) - qJD(3) * t297 + t303 * t287 + t360;
t157 = t336 * t162 + t333 * t163;
t294 = -t328 * t305 + t363;
t280 = Ifges(4,4) * t304 + Ifges(4,2) * t303 + Ifges(4,6) * qJD(3);
t281 = Ifges(4,1) * t304 + Ifges(4,4) * t303 + Ifges(4,5) * qJD(3);
t218 = Ifges(7,5) * t250 + Ifges(7,6) * t249 + Ifges(7,3) * t278;
t220 = Ifges(7,1) * t250 + Ifges(7,4) * t249 + Ifges(7,5) * t278;
t184 = -mrSges(7,1) * t202 + mrSges(7,3) * t197 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t246 - t250 * t218 + t278 * t220;
t219 = Ifges(7,4) * t250 + Ifges(7,2) * t249 + Ifges(7,6) * t278;
t185 = mrSges(7,2) * t202 - mrSges(7,3) * t196 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t246 + t249 * t218 - t278 * t219;
t230 = Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t282;
t232 = Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t282;
t167 = -mrSges(6,1) * t204 + mrSges(6,3) * t201 + Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * t247 - pkin(5) * t347 + pkin(9) * t358 + t335 * t184 + t331 * t185 - t270 * t230 + t282 * t232;
t231 = Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t282;
t169 = mrSges(6,2) * t204 - mrSges(6,3) * t200 + Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * t247 - pkin(9) * t183 - t331 * t184 + t335 * t185 + t269 * t230 - t282 * t231;
t259 = Ifges(5,4) * t283 - Ifges(5,2) * t282 + Ifges(5,6) * t324;
t260 = Ifges(5,1) * t283 - Ifges(5,4) * t282 + Ifges(5,5) * t324;
t345 = -mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t248 + Ifges(5,6) * t247 - Ifges(5,3) * t321 - pkin(4) * t343 - qJ(5) * t359 - t329 * t167 - t327 * t169 - t283 * t259 - t282 * t260;
t341 = -mrSges(4,1) * t253 + mrSges(4,2) * t254 - Ifges(4,5) * t293 - Ifges(4,6) * t292 - Ifges(4,3) * qJDD(3) - pkin(3) * t165 - t304 * t280 + t303 * t281 + t345;
t355 = Ifges(3,4) * t328 + Ifges(3,2) * t330;
t356 = Ifges(3,1) * t328 + Ifges(3,4) * t330;
t372 = -mrSges(3,1) * t294 + mrSges(3,2) * t295 - pkin(2) * t157 - (t328 * t355 - t330 * t356) * t338 + t341;
t369 = mrSges(3,2) * t328;
t176 = t329 * t181 + t327 * t182;
t354 = Ifges(3,5) * t328 + Ifges(3,6) * t330;
t367 = t338 * t354;
t350 = mrSges(3,3) * qJDD(1) + t338 * (-mrSges(3,1) * t330 + t369);
t155 = m(3) * t294 - t350 * t328 + t157;
t361 = -t333 * t162 + t336 * t163;
t156 = m(3) * t295 + t350 * t330 + t361;
t362 = -t328 * t155 + t330 * t156;
t349 = m(5) * t239 + t247 * mrSges(5,1) + t248 * mrSges(5,2) + t282 * t273 + t283 * t274 + t176;
t258 = Ifges(5,5) * t283 - Ifges(5,6) * t282 + Ifges(5,3) * t324;
t153 = mrSges(5,2) * t239 - mrSges(5,3) * t216 + Ifges(5,1) * t248 - Ifges(5,4) * t247 + Ifges(5,5) * t321 - qJ(5) * t176 - t327 * t167 + t329 * t169 - t282 * t258 - t324 * t259;
t346 = -mrSges(7,1) * t196 + mrSges(7,2) * t197 - Ifges(7,5) * t213 - Ifges(7,6) * t212 - Ifges(7,3) * t246 - t250 * t219 + t249 * t220;
t340 = -mrSges(6,1) * t200 + mrSges(6,2) * t201 - Ifges(6,5) * t236 - Ifges(6,6) * t235 - pkin(5) * t183 - t270 * t231 + t269 * t232 + t346;
t158 = t340 + (-Ifges(5,2) - Ifges(6,3)) * t247 + Ifges(5,6) * t321 + t324 * t260 - t283 * t258 + Ifges(5,4) * t248 - mrSges(5,1) * t239 + mrSges(5,3) * t217 - pkin(4) * t176;
t279 = Ifges(4,5) * t304 + Ifges(4,6) * t303 + Ifges(4,3) * qJD(3);
t148 = -mrSges(4,1) * t291 + mrSges(4,3) * t254 + Ifges(4,4) * t293 + Ifges(4,2) * t292 + Ifges(4,6) * qJDD(3) - pkin(3) * t349 + pkin(8) * t360 + qJD(3) * t281 + t332 * t153 + t371 * t158 - t304 * t279;
t149 = mrSges(4,2) * t291 - mrSges(4,3) * t253 + Ifges(4,1) * t293 + Ifges(4,4) * t292 + Ifges(4,5) * qJDD(3) - pkin(8) * t165 - qJD(3) * t280 + t371 * t153 - t332 * t158 + t303 * t279;
t301 = -qJDD(1) * pkin(1) - t338 * qJ(2) + t357;
t344 = m(4) * t291 - t292 * mrSges(4,1) + t293 * mrSges(4,2) - t303 * t296 + t304 * t297 + t349;
t144 = -mrSges(3,1) * t301 + mrSges(3,3) * t295 - pkin(2) * t344 + pkin(7) * t361 + t355 * qJDD(1) + t336 * t148 + t333 * t149 - t328 * t367;
t146 = mrSges(3,2) * t301 - mrSges(3,3) * t294 - pkin(7) * t157 + t356 * qJDD(1) - t333 * t148 + t336 * t149 + t330 * t367;
t342 = -m(3) * t301 + mrSges(3,1) * t364 - t344 + (t322 * t338 + t368) * mrSges(3,3);
t348 = -mrSges(2,2) * t312 + qJ(2) * t362 + t330 * t144 + t328 * t146 + pkin(1) * (-qJDD(1) * t369 + t342) + mrSges(2,1) * t311 + Ifges(2,3) * qJDD(1);
t170 = (mrSges(2,1) - t369) * qJDD(1) + t342 - t338 * mrSges(2,2) + m(2) * t311;
t152 = t330 * t155 + t328 * t156;
t150 = m(2) * t312 - t338 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t362;
t147 = mrSges(2,1) * g(3) - pkin(1) * t152 + (Ifges(2,6) - t354) * qJDD(1) + t338 * Ifges(2,5) + mrSges(2,3) * t312 + t372;
t142 = -mrSges(2,2) * g(3) - mrSges(2,3) * t311 + Ifges(2,5) * qJDD(1) - t338 * Ifges(2,6) - qJ(2) * t152 - t328 * t144 + t330 * t146;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t142 - t334 * t147 - pkin(6) * (t334 * t150 + t337 * t170), t142, t146, t149, t153, t169, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t142 + t337 * t147 + pkin(6) * (t337 * t150 - t334 * t170), t147, t144, t148, t158, t167, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t348, t348, t354 * qJDD(1) - t372, -t341, -t345, Ifges(6,3) * t247 - t340, -t346;];
m_new  = t1;
