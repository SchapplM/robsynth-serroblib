% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:23:42
% EndTime: 2019-05-06 01:24:16
% DurationCPUTime: 19.18s
% Computational Cost: add. (317418->360), mult. (761474->437), div. (0->0), fcn. (594503->10), ass. (0->144)
t333 = qJD(1) ^ 2;
t323 = sin(pkin(10));
t324 = cos(pkin(10));
t327 = sin(qJ(3));
t331 = cos(qJ(3));
t347 = -t323 * t327 + t324 * t331;
t299 = t347 * qJD(1);
t348 = t323 * t331 + t324 * t327;
t300 = t348 * qJD(1);
t326 = sin(qJ(4));
t330 = cos(qJ(4));
t280 = t299 * t330 - t300 * t326;
t289 = -t300 * qJD(3) + qJDD(1) * t347;
t364 = t299 * qJD(3);
t290 = qJDD(1) * t348 + t364;
t247 = qJD(4) * t280 + t289 * t326 + t290 * t330;
t281 = t299 * t326 + t300 * t330;
t320 = qJD(3) + qJD(4);
t325 = sin(qJ(5));
t329 = cos(qJ(5));
t267 = -t281 * t325 + t320 * t329;
t317 = qJDD(3) + qJDD(4);
t219 = qJD(5) * t267 + t247 * t329 + t317 * t325;
t268 = t281 * t329 + t320 * t325;
t248 = -mrSges(7,1) * t267 + mrSges(7,2) * t268;
t328 = sin(qJ(1));
t332 = cos(qJ(1));
t308 = -g(1) * t332 - g(2) * t328;
t301 = -pkin(1) * t333 + qJDD(1) * qJ(2) + t308;
t363 = qJD(1) * qJD(2);
t359 = -t324 * g(3) - 0.2e1 * t323 * t363;
t371 = pkin(2) * t324;
t274 = (-pkin(7) * qJDD(1) + t333 * t371 - t301) * t323 + t359;
t292 = -g(3) * t323 + (t301 + 0.2e1 * t363) * t324;
t362 = qJDD(1) * t324;
t319 = t324 ^ 2;
t368 = t319 * t333;
t275 = -pkin(2) * t368 + pkin(7) * t362 + t292;
t251 = t331 * t274 - t327 * t275;
t208 = (-t290 + t364) * pkin(8) + (t299 * t300 + qJDD(3)) * pkin(3) + t251;
t252 = t327 * t274 + t331 * t275;
t295 = qJD(3) * pkin(3) - pkin(8) * t300;
t298 = t299 ^ 2;
t222 = -pkin(3) * t298 + pkin(8) * t289 - qJD(3) * t295 + t252;
t203 = t326 * t208 + t330 * t222;
t264 = -pkin(4) * t280 - pkin(9) * t281;
t316 = t320 ^ 2;
t197 = -pkin(4) * t316 + pkin(9) * t317 + t264 * t280 + t203;
t318 = t323 ^ 2;
t307 = t328 * g(1) - t332 * g(2);
t353 = qJDD(2) - t307;
t288 = (-pkin(1) - t371) * qJDD(1) + (-qJ(2) + (-t318 - t319) * pkin(7)) * t333 + t353;
t235 = -t289 * pkin(3) - t298 * pkin(8) + t300 * t295 + t288;
t246 = -qJD(4) * t281 + t289 * t330 - t290 * t326;
t200 = (-t280 * t320 - t247) * pkin(9) + (t281 * t320 - t246) * pkin(4) + t235;
t191 = -t325 * t197 + t329 * t200;
t245 = qJDD(5) - t246;
t276 = qJD(5) - t280;
t187 = -0.2e1 * qJD(6) * t268 + (t267 * t276 - t219) * qJ(6) + (t267 * t268 + t245) * pkin(5) + t191;
t253 = -mrSges(7,2) * t276 + mrSges(7,3) * t267;
t361 = m(7) * t187 + t245 * mrSges(7,1) + t276 * t253;
t184 = -t219 * mrSges(7,3) - t268 * t248 + t361;
t192 = t329 * t197 + t325 * t200;
t218 = -qJD(5) * t268 - t247 * t325 + t317 * t329;
t227 = Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * t276;
t228 = Ifges(7,1) * t268 + Ifges(7,4) * t267 + Ifges(7,5) * t276;
t229 = Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t276;
t255 = pkin(5) * t276 - qJ(6) * t268;
t266 = t267 ^ 2;
t190 = -pkin(5) * t266 + qJ(6) * t218 + 0.2e1 * qJD(6) * t267 - t255 * t276 + t192;
t226 = Ifges(7,4) * t268 + Ifges(7,2) * t267 + Ifges(7,6) * t276;
t343 = -mrSges(7,1) * t187 + mrSges(7,2) * t190 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t245 - t268 * t226;
t373 = mrSges(6,1) * t191 - mrSges(6,2) * t192 + Ifges(6,5) * t219 + Ifges(6,6) * t218 + Ifges(6,3) * t245 + pkin(5) * t184 + t268 * t227 - (t229 + t228) * t267 - t343;
t263 = -mrSges(5,1) * t280 + mrSges(5,2) * t281;
t272 = mrSges(5,1) * t320 - mrSges(5,3) * t281;
t249 = -mrSges(6,1) * t267 + mrSges(6,2) * t268;
t254 = -mrSges(6,2) * t276 + mrSges(6,3) * t267;
t176 = m(6) * t191 + t245 * mrSges(6,1) + t276 * t254 + (-t248 - t249) * t268 + (-mrSges(6,3) - mrSges(7,3)) * t219 + t361;
t360 = m(7) * t190 + t218 * mrSges(7,3) + t267 * t248;
t256 = mrSges(7,1) * t276 - mrSges(7,3) * t268;
t366 = -mrSges(6,1) * t276 + mrSges(6,3) * t268 - t256;
t370 = -mrSges(6,2) - mrSges(7,2);
t179 = m(6) * t192 + t218 * mrSges(6,3) + t370 * t245 + t267 * t249 + t366 * t276 + t360;
t355 = -t176 * t325 + t329 * t179;
t169 = m(5) * t203 - mrSges(5,2) * t317 + mrSges(5,3) * t246 + t263 * t280 - t272 * t320 + t355;
t202 = t208 * t330 - t326 * t222;
t271 = -mrSges(5,2) * t320 + mrSges(5,3) * t280;
t196 = -pkin(4) * t317 - pkin(9) * t316 + t281 * t264 - t202;
t194 = -pkin(5) * t218 - qJ(6) * t266 + t255 * t268 + qJDD(6) + t196;
t354 = -m(7) * t194 + t218 * mrSges(7,1) + t267 * t253;
t338 = -m(6) * t196 + t218 * mrSges(6,1) + t370 * t219 + t267 * t254 + t366 * t268 + t354;
t181 = m(5) * t202 + t317 * mrSges(5,1) - t247 * mrSges(5,3) - t281 * t263 + t320 * t271 + t338;
t162 = t326 * t169 + t330 * t181;
t284 = -mrSges(4,1) * t299 + mrSges(4,2) * t300;
t293 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t299;
t159 = m(4) * t251 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t290 + qJD(3) * t293 - t284 * t300 + t162;
t294 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t300;
t356 = t330 * t169 - t181 * t326;
t160 = m(4) * t252 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t289 - qJD(3) * t294 + t284 * t299 + t356;
t154 = t331 * t159 + t327 * t160;
t291 = -t323 * t301 + t359;
t278 = Ifges(4,4) * t300 + Ifges(4,2) * t299 + Ifges(4,6) * qJD(3);
t279 = Ifges(4,1) * t300 + Ifges(4,4) * t299 + Ifges(4,5) * qJD(3);
t224 = Ifges(7,5) * t268 + Ifges(7,6) * t267 + Ifges(7,3) * t276;
t225 = Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t276;
t344 = -mrSges(7,1) * t194 + mrSges(7,3) * t190 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t245 + t276 * t228;
t164 = Ifges(6,4) * t219 + Ifges(6,2) * t218 + Ifges(6,6) * t245 + t276 * t229 - mrSges(6,1) * t196 + mrSges(6,3) * t192 - pkin(5) * (t219 * mrSges(7,2) - t354) + qJ(6) * (-t245 * mrSges(7,2) - t276 * t256 + t360) + (-pkin(5) * t256 - t224 - t225) * t268 + t344;
t342 = mrSges(7,2) * t194 - mrSges(7,3) * t187 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t245 + t267 * t224;
t171 = mrSges(6,2) * t196 - mrSges(6,3) * t191 + Ifges(6,1) * t219 + Ifges(6,4) * t218 + Ifges(6,5) * t245 - qJ(6) * t184 + t267 * t225 + (-t226 - t227) * t276 + t342;
t259 = Ifges(5,4) * t281 + Ifges(5,2) * t280 + Ifges(5,6) * t320;
t260 = Ifges(5,1) * t281 + Ifges(5,4) * t280 + Ifges(5,5) * t320;
t340 = -mrSges(5,1) * t202 + mrSges(5,2) * t203 - Ifges(5,5) * t247 - Ifges(5,6) * t246 - Ifges(5,3) * t317 - pkin(4) * t338 - pkin(9) * t355 - t329 * t164 - t325 * t171 - t281 * t259 + t280 * t260;
t335 = -mrSges(4,1) * t251 + mrSges(4,2) * t252 - Ifges(4,5) * t290 - Ifges(4,6) * t289 - Ifges(4,3) * qJDD(3) - pkin(3) * t162 - t300 * t278 + t299 * t279 + t340;
t351 = Ifges(3,4) * t323 + Ifges(3,2) * t324;
t352 = Ifges(3,1) * t323 + Ifges(3,4) * t324;
t372 = -mrSges(3,1) * t291 + mrSges(3,2) * t292 - pkin(2) * t154 - (t323 * t351 - t324 * t352) * t333 + t335;
t369 = mrSges(3,2) * t323;
t173 = t329 * t176 + t325 * t179;
t350 = Ifges(3,5) * t323 + Ifges(3,6) * t324;
t365 = t333 * t350;
t346 = mrSges(3,3) * qJDD(1) + t333 * (-mrSges(3,1) * t324 + t369);
t152 = m(3) * t291 - t323 * t346 + t154;
t357 = -t327 * t159 + t331 * t160;
t153 = m(3) * t292 + t324 * t346 + t357;
t358 = -t152 * t323 + t324 * t153;
t345 = m(5) * t235 - t246 * mrSges(5,1) + t247 * mrSges(5,2) - t280 * t271 + t281 * t272 + t173;
t258 = Ifges(5,5) * t281 + Ifges(5,6) * t280 + Ifges(5,3) * t320;
t150 = mrSges(5,2) * t235 - mrSges(5,3) * t202 + Ifges(5,1) * t247 + Ifges(5,4) * t246 + Ifges(5,5) * t317 - pkin(9) * t173 - t164 * t325 + t171 * t329 + t258 * t280 - t259 * t320;
t155 = -mrSges(5,1) * t235 + mrSges(5,3) * t203 + Ifges(5,4) * t247 + Ifges(5,2) * t246 + Ifges(5,6) * t317 - pkin(4) * t173 - t281 * t258 + t320 * t260 - t373;
t277 = Ifges(4,5) * t300 + Ifges(4,6) * t299 + Ifges(4,3) * qJD(3);
t145 = -mrSges(4,1) * t288 + mrSges(4,3) * t252 + Ifges(4,4) * t290 + Ifges(4,2) * t289 + Ifges(4,6) * qJDD(3) - pkin(3) * t345 + pkin(8) * t356 + qJD(3) * t279 + t326 * t150 + t330 * t155 - t300 * t277;
t146 = mrSges(4,2) * t288 - mrSges(4,3) * t251 + Ifges(4,1) * t290 + Ifges(4,4) * t289 + Ifges(4,5) * qJDD(3) - pkin(8) * t162 - qJD(3) * t278 + t150 * t330 - t155 * t326 + t277 * t299;
t297 = -qJDD(1) * pkin(1) - t333 * qJ(2) + t353;
t339 = m(4) * t288 - t289 * mrSges(4,1) + t290 * mrSges(4,2) - t299 * t293 + t300 * t294 + t345;
t141 = -mrSges(3,1) * t297 + mrSges(3,3) * t292 - pkin(2) * t339 + pkin(7) * t357 + qJDD(1) * t351 + t331 * t145 + t327 * t146 - t323 * t365;
t143 = mrSges(3,2) * t297 - mrSges(3,3) * t291 - pkin(7) * t154 + qJDD(1) * t352 - t327 * t145 + t331 * t146 + t324 * t365;
t337 = -m(3) * t297 + mrSges(3,1) * t362 - t339 + (t318 * t333 + t368) * mrSges(3,3);
t341 = -mrSges(2,2) * t308 + qJ(2) * t358 + t324 * t141 + t323 * t143 + pkin(1) * (-qJDD(1) * t369 + t337) + mrSges(2,1) * t307 + Ifges(2,3) * qJDD(1);
t165 = t337 - t333 * mrSges(2,2) + m(2) * t307 + (mrSges(2,1) - t369) * qJDD(1);
t149 = t152 * t324 + t153 * t323;
t147 = m(2) * t308 - mrSges(2,1) * t333 - qJDD(1) * mrSges(2,2) + t358;
t144 = -pkin(1) * t149 + mrSges(2,1) * g(3) + (Ifges(2,6) - t350) * qJDD(1) + t333 * Ifges(2,5) + mrSges(2,3) * t308 + t372;
t139 = -mrSges(2,2) * g(3) - mrSges(2,3) * t307 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t333 - qJ(2) * t149 - t141 * t323 + t143 * t324;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t139 - t328 * t144 - pkin(6) * (t147 * t328 + t165 * t332), t139, t143, t146, t150, t171, -t226 * t276 + t342; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t328 * t139 + t332 * t144 + pkin(6) * (t147 * t332 - t165 * t328), t144, t141, t145, t155, t164, -t268 * t224 + t344; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t341, t341, qJDD(1) * t350 - t372, -t335, -t340, t373, -t267 * t228 - t343;];
m_new  = t1;
