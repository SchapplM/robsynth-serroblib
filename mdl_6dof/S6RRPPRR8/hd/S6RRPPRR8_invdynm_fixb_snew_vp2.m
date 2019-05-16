% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 11:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:17:49
% EndTime: 2019-05-06 11:18:19
% DurationCPUTime: 14.82s
% Computational Cost: add. (246634->381), mult. (543744->461), div. (0->0), fcn. (368938->10), ass. (0->144)
t394 = -2 * qJD(3);
t393 = -2 * qJD(4);
t355 = sin(qJ(1));
t359 = cos(qJ(1));
t337 = -g(1) * t359 - g(2) * t355;
t361 = qJD(1) ^ 2;
t315 = -pkin(1) * t361 + qJDD(1) * pkin(7) + t337;
t354 = sin(qJ(2));
t358 = cos(qJ(2));
t286 = -t358 * g(3) - t354 * t315;
t329 = (-pkin(2) * t358 - qJ(3) * t354) * qJD(1);
t360 = qJD(2) ^ 2;
t384 = qJD(1) * t354;
t264 = -qJDD(2) * pkin(2) - t360 * qJ(3) + t329 * t384 + qJDD(3) - t286;
t381 = qJD(1) * qJD(2);
t379 = t358 * t381;
t331 = qJDD(1) * t354 + t379;
t351 = sin(pkin(10));
t389 = cos(pkin(10));
t300 = -qJDD(2) * t389 + t331 * t351;
t301 = t351 * qJDD(2) + t331 * t389;
t322 = t351 * qJD(2) + t384 * t389;
t344 = t358 * qJD(1);
t321 = -qJD(2) * t389 + t351 * t384;
t380 = t321 * t344;
t238 = t264 + (-t301 - t380) * qJ(4) + (-t322 * t344 + t300) * pkin(3) + t322 * t393;
t336 = t355 * g(1) - t359 * g(2);
t314 = -qJDD(1) * pkin(1) - t361 * pkin(7) - t336;
t341 = t354 * t381;
t332 = qJDD(1) * t358 - t341;
t260 = (-t331 - t379) * qJ(3) + (-t332 + t341) * pkin(2) + t314;
t287 = -g(3) * t354 + t358 * t315;
t265 = -pkin(2) * t360 + qJDD(2) * qJ(3) + t329 * t344 + t287;
t239 = t260 * t389 - t351 * t265 + t322 * t394;
t240 = t351 * t260 + t389 * t265 + t321 * t394;
t276 = Ifges(4,1) * t322 - Ifges(4,4) * t321 - Ifges(4,5) * t344;
t284 = mrSges(5,1) * t321 - mrSges(5,3) * t322;
t296 = -mrSges(5,2) * t321 - mrSges(5,3) * t344;
t299 = mrSges(5,1) * t344 + mrSges(5,2) * t322;
t283 = pkin(3) * t321 - qJ(4) * t322;
t388 = t358 ^ 2 * t361;
t231 = t332 * pkin(3) - qJ(4) * t388 + t322 * t283 + qJDD(4) - t239;
t214 = (-t301 + t380) * pkin(8) + (t321 * t322 + t332) * pkin(4) + t231;
t229 = -pkin(3) * t388 - t332 * qJ(4) - t321 * t283 + t344 * t393 + t240;
t302 = pkin(4) * t344 - pkin(8) * t322;
t319 = t321 ^ 2;
t216 = -pkin(4) * t319 + pkin(8) * t300 - t302 * t344 + t229;
t353 = sin(qJ(5));
t357 = cos(qJ(5));
t208 = t357 * t214 - t353 * t216;
t281 = t321 * t357 - t322 * t353;
t250 = qJD(5) * t281 + t300 * t353 + t301 * t357;
t282 = t321 * t353 + t322 * t357;
t328 = qJDD(5) + t332;
t339 = t344 + qJD(5);
t204 = (t281 * t339 - t250) * pkin(9) + (t281 * t282 + t328) * pkin(5) + t208;
t209 = t353 * t214 + t357 * t216;
t249 = -qJD(5) * t282 + t300 * t357 - t301 * t353;
t270 = pkin(5) * t339 - pkin(9) * t282;
t280 = t281 ^ 2;
t205 = -pkin(5) * t280 + pkin(9) * t249 - t270 * t339 + t209;
t352 = sin(qJ(6));
t356 = cos(qJ(6));
t202 = t204 * t356 - t205 * t352;
t256 = t281 * t356 - t282 * t352;
t222 = qJD(6) * t256 + t249 * t352 + t250 * t356;
t257 = t281 * t352 + t282 * t356;
t237 = -mrSges(7,1) * t256 + mrSges(7,2) * t257;
t338 = qJD(6) + t339;
t243 = -mrSges(7,2) * t338 + mrSges(7,3) * t256;
t317 = qJDD(6) + t328;
t197 = m(7) * t202 + mrSges(7,1) * t317 - mrSges(7,3) * t222 - t237 * t257 + t243 * t338;
t203 = t204 * t352 + t205 * t356;
t221 = -qJD(6) * t257 + t249 * t356 - t250 * t352;
t244 = mrSges(7,1) * t338 - mrSges(7,3) * t257;
t198 = m(7) * t203 - mrSges(7,2) * t317 + mrSges(7,3) * t221 + t237 * t256 - t244 * t338;
t188 = t356 * t197 + t352 * t198;
t258 = -mrSges(6,1) * t281 + mrSges(6,2) * t282;
t266 = -mrSges(6,2) * t339 + mrSges(6,3) * t281;
t185 = m(6) * t208 + mrSges(6,1) * t328 - mrSges(6,3) * t250 - t258 * t282 + t266 * t339 + t188;
t267 = mrSges(6,1) * t339 - mrSges(6,3) * t282;
t377 = -t197 * t352 + t356 * t198;
t186 = m(6) * t209 - mrSges(6,2) * t328 + mrSges(6,3) * t249 + t258 * t281 - t267 * t339 + t377;
t182 = t357 * t185 + t353 * t186;
t275 = Ifges(5,1) * t322 - Ifges(5,4) * t344 + Ifges(5,5) * t321;
t252 = Ifges(6,4) * t282 + Ifges(6,2) * t281 + Ifges(6,6) * t339;
t253 = Ifges(6,1) * t282 + Ifges(6,4) * t281 + Ifges(6,5) * t339;
t233 = Ifges(7,4) * t257 + Ifges(7,2) * t256 + Ifges(7,6) * t338;
t234 = Ifges(7,1) * t257 + Ifges(7,4) * t256 + Ifges(7,5) * t338;
t372 = mrSges(7,1) * t202 - mrSges(7,2) * t203 + Ifges(7,5) * t222 + Ifges(7,6) * t221 + Ifges(7,3) * t317 + t257 * t233 - t256 * t234;
t366 = mrSges(6,1) * t208 - mrSges(6,2) * t209 + Ifges(6,5) * t250 + Ifges(6,6) * t249 + Ifges(6,3) * t328 + pkin(5) * t188 + t282 * t252 - t281 * t253 + t372;
t363 = mrSges(5,1) * t231 - mrSges(5,3) * t229 - Ifges(5,4) * t301 + Ifges(5,2) * t332 - Ifges(5,6) * t300 + pkin(4) * t182 - t321 * t275 + t366;
t371 = -m(5) * t231 - t332 * mrSges(5,1) - t182;
t183 = -t353 * t185 + t357 * t186;
t373 = m(5) * t229 - t332 * mrSges(5,3) + t183;
t271 = Ifges(5,5) * t322 - Ifges(5,6) * t344 + Ifges(5,3) * t321;
t387 = -Ifges(4,4) * t322 + Ifges(4,2) * t321 + Ifges(4,6) * t344 + t271;
t392 = t387 * t322 - mrSges(4,1) * t239 + mrSges(4,2) * t240 - Ifges(4,5) * t301 + Ifges(4,6) * t300 - pkin(3) * (-t301 * mrSges(5,2) - t322 * t284 - t296 * t344 + t371) - qJ(4) * (-t300 * mrSges(5,2) - t321 * t284 - t299 * t344 + t373) - t321 * t276 + t363;
t226 = -pkin(4) * t300 - pkin(8) * t319 + t322 * t302 - t238;
t211 = -pkin(5) * t249 - pkin(9) * t280 + t270 * t282 + t226;
t375 = m(7) * t211 - t221 * mrSges(7,1) + t222 * mrSges(7,2) - t256 * t243 + t257 * t244;
t199 = -m(6) * t226 + t249 * mrSges(6,1) - t250 * mrSges(6,2) + t281 * t266 - t282 * t267 - t375;
t193 = m(5) * t238 + t300 * mrSges(5,1) - t301 * mrSges(5,3) + t321 * t296 - t322 * t299 + t199;
t232 = Ifges(7,5) * t257 + Ifges(7,6) * t256 + Ifges(7,3) * t338;
t189 = -mrSges(7,1) * t211 + mrSges(7,3) * t203 + Ifges(7,4) * t222 + Ifges(7,2) * t221 + Ifges(7,6) * t317 - t232 * t257 + t234 * t338;
t190 = mrSges(7,2) * t211 - mrSges(7,3) * t202 + Ifges(7,1) * t222 + Ifges(7,4) * t221 + Ifges(7,5) * t317 + t232 * t256 - t233 * t338;
t251 = Ifges(6,5) * t282 + Ifges(6,6) * t281 + Ifges(6,3) * t339;
t174 = -mrSges(6,1) * t226 + mrSges(6,3) * t209 + Ifges(6,4) * t250 + Ifges(6,2) * t249 + Ifges(6,6) * t328 - pkin(5) * t375 + pkin(9) * t377 + t356 * t189 + t352 * t190 - t282 * t251 + t339 * t253;
t176 = mrSges(6,2) * t226 - mrSges(6,3) * t208 + Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t328 - pkin(9) * t188 - t189 * t352 + t190 * t356 + t251 * t281 - t252 * t339;
t367 = -mrSges(5,1) * t238 + mrSges(5,2) * t229 - pkin(4) * t199 - pkin(8) * t183 - t357 * t174 - t353 * t176;
t273 = Ifges(5,4) * t322 - Ifges(5,2) * t344 + Ifges(5,6) * t321;
t386 = -Ifges(4,5) * t322 + Ifges(4,6) * t321 + Ifges(4,3) * t344 - t273;
t163 = -mrSges(4,1) * t264 + mrSges(4,3) * t240 - pkin(3) * t193 + (-Ifges(4,6) + Ifges(5,6)) * t332 + t386 * t322 + (Ifges(4,4) - Ifges(5,5)) * t301 + (-Ifges(4,2) - Ifges(5,3)) * t300 + (-t275 - t276) * t344 + t367;
t369 = mrSges(5,2) * t231 - mrSges(5,3) * t238 + Ifges(5,1) * t301 - Ifges(5,4) * t332 + Ifges(5,5) * t300 - pkin(8) * t182 - t353 * t174 + t357 * t176;
t164 = mrSges(4,2) * t264 - mrSges(4,3) * t239 + Ifges(4,1) * t301 - Ifges(4,4) * t300 - Ifges(4,5) * t332 - qJ(4) * t193 + t321 * t386 - t344 * t387 + t369;
t298 = -mrSges(4,1) * t344 - mrSges(4,3) * t322;
t385 = -mrSges(4,1) * t321 - mrSges(4,2) * t322 - t284;
t390 = -mrSges(4,3) - mrSges(5,2);
t178 = m(4) * t240 + t332 * mrSges(4,2) + t385 * t321 + t390 * t300 + (t298 - t299) * t344 + t373;
t297 = mrSges(4,2) * t344 - mrSges(4,3) * t321;
t179 = m(4) * t239 - t332 * mrSges(4,1) + t385 * t322 + t390 * t301 + (-t296 - t297) * t344 + t371;
t173 = t389 * t178 - t179 * t351;
t192 = -m(4) * t264 - t300 * mrSges(4,1) - t301 * mrSges(4,2) - t321 * t297 - t322 * t298 - t193;
t312 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t354 + Ifges(3,2) * t358) * qJD(1);
t313 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t354 + Ifges(3,4) * t358) * qJD(1);
t391 = mrSges(3,1) * t286 - mrSges(3,2) * t287 + Ifges(3,5) * t331 + Ifges(3,6) * t332 + Ifges(3,3) * qJDD(2) + pkin(2) * t192 + qJ(3) * t173 + (t312 * t354 - t313 * t358) * qJD(1) + t163 * t389 + t351 * t164;
t330 = (-mrSges(3,1) * t358 + mrSges(3,2) * t354) * qJD(1);
t334 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t384;
t171 = m(3) * t287 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t332 - qJD(2) * t334 + t330 * t344 + t173;
t335 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t344;
t191 = m(3) * t286 + qJDD(2) * mrSges(3,1) - t331 * mrSges(3,3) + qJD(2) * t335 - t330 * t384 + t192;
t378 = t358 * t171 - t191 * t354;
t172 = t351 * t178 + t179 * t389;
t311 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t354 + Ifges(3,6) * t358) * qJD(1);
t160 = mrSges(3,2) * t314 - mrSges(3,3) * t286 + Ifges(3,1) * t331 + Ifges(3,4) * t332 + Ifges(3,5) * qJDD(2) - qJ(3) * t172 - qJD(2) * t312 - t351 * t163 + t164 * t389 + t311 * t344;
t162 = Ifges(3,6) * qJDD(2) - pkin(2) * t172 + mrSges(3,3) * t287 + (Ifges(3,2) + Ifges(4,3)) * t332 - t311 * t384 + qJD(2) * t313 - mrSges(3,1) * t314 + Ifges(3,4) * t331 + t392;
t365 = -m(3) * t314 + t332 * mrSges(3,1) - t331 * mrSges(3,2) - t334 * t384 + t335 * t344 - t172;
t370 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t365 + pkin(7) * t378 + t354 * t160 + t358 * t162;
t168 = m(2) * t336 + qJDD(1) * mrSges(2,1) - t361 * mrSges(2,2) + t365;
t167 = t171 * t354 + t191 * t358;
t165 = m(2) * t337 - mrSges(2,1) * t361 - qJDD(1) * mrSges(2,2) + t378;
t158 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + t361 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t391;
t157 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t361 - pkin(7) * t167 + t160 * t358 - t162 * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t359 * t157 - t355 * t158 - pkin(6) * (t165 * t355 + t168 * t359), t157, t160, t164, -t271 * t344 - t321 * t273 + t369, t176, t190; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t355 * t157 + t359 * t158 + pkin(6) * (t165 * t359 - t168 * t355), t158, t162, t163, -t322 * t271 - t363, t174, t189; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t370, t370, t391, -Ifges(4,3) * t332 - t392, Ifges(5,5) * t301 - Ifges(5,6) * t332 + Ifges(5,3) * t300 + t322 * t273 + t275 * t344 - t367, t366, t372;];
m_new  = t1;
