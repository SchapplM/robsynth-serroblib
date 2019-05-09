% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:37:13
% EndTime: 2019-05-07 07:37:50
% DurationCPUTime: 19.77s
% Computational Cost: add. (354440->383), mult. (732191->469), div. (0->0), fcn. (525471->10), ass. (0->142)
t332 = sin(qJ(2));
t334 = cos(qJ(2));
t357 = qJD(1) * qJD(2);
t312 = qJDD(1) * t332 + t334 * t357;
t333 = sin(qJ(1));
t335 = cos(qJ(1));
t319 = -g(1) * t335 - g(2) * t333;
t336 = qJD(1) ^ 2;
t306 = -pkin(1) * t336 + qJDD(1) * pkin(7) + t319;
t362 = t306 * t332;
t364 = pkin(2) * t336;
t263 = qJDD(2) * pkin(2) - pkin(8) * t312 - t362 + (pkin(8) * t357 + t332 * t364 - g(3)) * t334;
t294 = -g(3) * t332 + t334 * t306;
t313 = qJDD(1) * t334 - t332 * t357;
t359 = qJD(1) * t332;
t317 = qJD(2) * pkin(2) - pkin(8) * t359;
t327 = t334 ^ 2;
t264 = pkin(8) * t313 - qJD(2) * t317 - t327 * t364 + t294;
t331 = sin(qJ(3));
t366 = cos(qJ(3));
t240 = t331 * t263 + t366 * t264;
t304 = (t331 * t334 + t366 * t332) * qJD(1);
t275 = qJD(3) * t304 + t312 * t331 - t366 * t313;
t358 = qJD(1) * t334;
t303 = t331 * t359 - t366 * t358;
t286 = mrSges(4,1) * t303 + mrSges(4,2) * t304;
t325 = qJD(2) + qJD(3);
t296 = mrSges(4,1) * t325 - mrSges(4,3) * t304;
t324 = qJDD(2) + qJDD(3);
t276 = -t303 * qJD(3) + t366 * t312 + t331 * t313;
t318 = g(1) * t333 - t335 * g(2);
t348 = -qJDD(1) * pkin(1) - t318;
t277 = -pkin(2) * t313 + t317 * t359 + (-pkin(8) * t327 - pkin(7)) * t336 + t348;
t221 = (t303 * t325 - t276) * qJ(4) + (t304 * t325 + t275) * pkin(3) + t277;
t285 = pkin(3) * t303 - qJ(4) * t304;
t323 = t325 ^ 2;
t224 = -pkin(3) * t323 + qJ(4) * t324 - t285 * t303 + t240;
t328 = sin(pkin(10));
t329 = cos(pkin(10));
t292 = t304 * t329 + t325 * t328;
t203 = -0.2e1 * qJD(4) * t292 + t329 * t221 - t224 * t328;
t257 = t276 * t329 + t324 * t328;
t291 = -t304 * t328 + t325 * t329;
t200 = (t291 * t303 - t257) * pkin(9) + (t291 * t292 + t275) * pkin(4) + t203;
t204 = 0.2e1 * qJD(4) * t291 + t328 * t221 + t329 * t224;
t256 = -t276 * t328 + t324 * t329;
t280 = pkin(4) * t303 - pkin(9) * t292;
t290 = t291 ^ 2;
t202 = -pkin(4) * t290 + pkin(9) * t256 - t280 * t303 + t204;
t330 = sin(qJ(5));
t365 = cos(qJ(5));
t196 = t330 * t200 + t365 * t202;
t255 = t330 * t291 + t365 * t292;
t219 = qJD(5) * t255 - t365 * t256 + t257 * t330;
t299 = qJD(5) + t303;
t245 = mrSges(6,1) * t299 - mrSges(6,3) * t255;
t254 = -t365 * t291 + t292 * t330;
t274 = qJDD(5) + t275;
t235 = pkin(5) * t254 - qJ(6) * t255;
t298 = t299 ^ 2;
t191 = -pkin(5) * t298 + qJ(6) * t274 + 0.2e1 * qJD(6) * t299 - t235 * t254 + t196;
t246 = -mrSges(7,1) * t299 + mrSges(7,2) * t255;
t356 = m(7) * t191 + t274 * mrSges(7,3) + t299 * t246;
t236 = mrSges(7,1) * t254 - mrSges(7,3) * t255;
t360 = -mrSges(6,1) * t254 - mrSges(6,2) * t255 - t236;
t363 = -mrSges(6,3) - mrSges(7,2);
t177 = m(6) * t196 - mrSges(6,2) * t274 + t363 * t219 - t245 * t299 + t360 * t254 + t356;
t195 = t365 * t200 - t330 * t202;
t220 = -t254 * qJD(5) + t330 * t256 + t365 * t257;
t244 = -mrSges(6,2) * t299 - mrSges(6,3) * t254;
t193 = -t274 * pkin(5) - t298 * qJ(6) + t255 * t235 + qJDD(6) - t195;
t243 = -mrSges(7,2) * t254 + mrSges(7,3) * t299;
t351 = -m(7) * t193 + t274 * mrSges(7,1) + t299 * t243;
t179 = m(6) * t195 + mrSges(6,1) * t274 + t363 * t220 + t244 * t299 + t360 * t255 + t351;
t172 = t330 * t177 + t365 * t179;
t259 = -mrSges(5,1) * t291 + mrSges(5,2) * t292;
t278 = -mrSges(5,2) * t303 + mrSges(5,3) * t291;
t170 = m(5) * t203 + mrSges(5,1) * t275 - mrSges(5,3) * t257 - t259 * t292 + t278 * t303 + t172;
t279 = mrSges(5,1) * t303 - mrSges(5,3) * t292;
t352 = t365 * t177 - t179 * t330;
t171 = m(5) * t204 - mrSges(5,2) * t275 + mrSges(5,3) * t256 + t259 * t291 - t279 * t303 + t352;
t353 = -t170 * t328 + t329 * t171;
t163 = m(4) * t240 - mrSges(4,2) * t324 - mrSges(4,3) * t275 - t286 * t303 - t296 * t325 + t353;
t239 = t366 * t263 - t331 * t264;
t295 = -mrSges(4,2) * t325 - mrSges(4,3) * t303;
t223 = -t324 * pkin(3) - t323 * qJ(4) + t304 * t285 + qJDD(4) - t239;
t205 = -t256 * pkin(4) - t290 * pkin(9) + t292 * t280 + t223;
t198 = -0.2e1 * qJD(6) * t255 + (t254 * t299 - t220) * qJ(6) + (t255 * t299 + t219) * pkin(5) + t205;
t188 = m(7) * t198 + t219 * mrSges(7,1) - t220 * mrSges(7,3) + t254 * t243 - t255 * t246;
t342 = m(6) * t205 + t219 * mrSges(6,1) + mrSges(6,2) * t220 + t254 * t244 + t245 * t255 + t188;
t339 = -m(5) * t223 + t256 * mrSges(5,1) - mrSges(5,2) * t257 + t291 * t278 - t279 * t292 - t342;
t181 = m(4) * t239 + mrSges(4,1) * t324 - mrSges(4,3) * t276 - t286 * t304 + t295 * t325 + t339;
t158 = t331 * t163 + t366 * t181;
t293 = -g(3) * t334 - t362;
t301 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t332 + Ifges(3,2) * t334) * qJD(1);
t302 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t332 + Ifges(3,4) * t334) * qJD(1);
t230 = Ifges(7,1) * t255 + Ifges(7,4) * t299 + Ifges(7,5) * t254;
t231 = Ifges(6,1) * t255 - Ifges(6,4) * t254 + Ifges(6,5) * t299;
t350 = -mrSges(7,1) * t198 + mrSges(7,2) * t191;
t228 = Ifges(7,4) * t255 + Ifges(7,2) * t299 + Ifges(7,6) * t254;
t361 = -Ifges(6,5) * t255 + Ifges(6,6) * t254 - Ifges(6,3) * t299 - t228;
t173 = -mrSges(6,1) * t205 + mrSges(6,3) * t196 - pkin(5) * t188 + (t230 + t231) * t299 + (Ifges(6,6) - Ifges(7,6)) * t274 + t361 * t255 + (Ifges(6,4) - Ifges(7,5)) * t220 + (-Ifges(6,2) - Ifges(7,3)) * t219 + t350;
t229 = Ifges(6,4) * t255 - Ifges(6,2) * t254 + Ifges(6,6) * t299;
t226 = Ifges(7,5) * t255 + Ifges(7,6) * t299 + Ifges(7,3) * t254;
t347 = mrSges(7,2) * t193 - mrSges(7,3) * t198 + Ifges(7,1) * t220 + Ifges(7,4) * t274 + Ifges(7,5) * t219 + t299 * t226;
t174 = mrSges(6,2) * t205 - mrSges(6,3) * t195 + Ifges(6,1) * t220 - Ifges(6,4) * t219 + Ifges(6,5) * t274 - qJ(6) * t188 - t229 * t299 + t361 * t254 + t347;
t247 = Ifges(5,5) * t292 + Ifges(5,6) * t291 + Ifges(5,3) * t303;
t249 = Ifges(5,1) * t292 + Ifges(5,4) * t291 + Ifges(5,5) * t303;
t152 = -mrSges(5,1) * t223 + mrSges(5,3) * t204 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * t275 - pkin(4) * t342 + pkin(9) * t352 + t365 * t173 + t330 * t174 - t292 * t247 + t303 * t249;
t248 = Ifges(5,4) * t292 + Ifges(5,2) * t291 + Ifges(5,6) * t303;
t154 = mrSges(5,2) * t223 - mrSges(5,3) * t203 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * t275 - pkin(9) * t172 - t330 * t173 + t365 * t174 + t291 * t247 - t303 * t248;
t282 = Ifges(4,4) * t304 - Ifges(4,2) * t303 + Ifges(4,6) * t325;
t283 = Ifges(4,1) * t304 - Ifges(4,4) * t303 + Ifges(4,5) * t325;
t343 = -mrSges(4,1) * t239 + mrSges(4,2) * t240 - Ifges(4,5) * t276 + Ifges(4,6) * t275 - Ifges(4,3) * t324 - pkin(3) * t339 - qJ(4) * t353 - t329 * t152 - t328 * t154 - t304 * t282 - t303 * t283;
t367 = mrSges(3,1) * t293 - mrSges(3,2) * t294 + Ifges(3,5) * t312 + Ifges(3,6) * t313 + Ifges(3,3) * qJDD(2) + pkin(2) * t158 + (t332 * t301 - t334 * t302) * qJD(1) - t343;
t165 = t329 * t170 + t328 * t171;
t311 = (-mrSges(3,1) * t334 + mrSges(3,2) * t332) * qJD(1);
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t358;
t156 = m(3) * t293 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t312 + qJD(2) * t316 - t311 * t359 + t158;
t315 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t359;
t354 = t366 * t163 - t181 * t331;
t157 = m(3) * t294 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t313 - qJD(2) * t315 + t311 * t358 + t354;
t355 = -t156 * t332 + t334 * t157;
t281 = Ifges(4,5) * t304 - Ifges(4,6) * t303 + Ifges(4,3) * t325;
t146 = mrSges(4,2) * t277 - mrSges(4,3) * t239 + Ifges(4,1) * t276 - Ifges(4,4) * t275 + Ifges(4,5) * t324 - qJ(4) * t165 - t152 * t328 + t154 * t329 - t281 * t303 - t282 * t325;
t345 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t220 - Ifges(7,2) * t274 - Ifges(7,6) * t219 + t255 * t226 - t254 * t230;
t340 = mrSges(6,2) * t196 - t254 * t231 - qJ(6) * (-mrSges(7,2) * t219 - t236 * t254 + t356) - pkin(5) * (-mrSges(7,2) * t220 - t236 * t255 + t351) - mrSges(6,1) * t195 - t255 * t229 + Ifges(6,6) * t219 - Ifges(6,5) * t220 - Ifges(6,3) * t274 + t345;
t337 = mrSges(5,1) * t203 - mrSges(5,2) * t204 + Ifges(5,5) * t257 + Ifges(5,6) * t256 + pkin(4) * t172 + t292 * t248 - t291 * t249 - t340;
t150 = -t337 + (-Ifges(4,2) - Ifges(5,3)) * t275 + Ifges(4,6) * t324 + t325 * t283 - t304 * t281 + Ifges(4,4) * t276 - mrSges(4,1) * t277 + mrSges(4,3) * t240 - pkin(3) * t165;
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t332 + Ifges(3,6) * t334) * qJD(1);
t305 = -pkin(7) * t336 + t348;
t344 = m(4) * t277 + t275 * mrSges(4,1) + mrSges(4,2) * t276 + t303 * t295 + t296 * t304 + t165;
t142 = -mrSges(3,1) * t305 + mrSges(3,3) * t294 + Ifges(3,4) * t312 + Ifges(3,2) * t313 + Ifges(3,6) * qJDD(2) - pkin(2) * t344 + pkin(8) * t354 + qJD(2) * t302 + t331 * t146 + t366 * t150 - t300 * t359;
t145 = mrSges(3,2) * t305 - mrSges(3,3) * t293 + Ifges(3,1) * t312 + Ifges(3,4) * t313 + Ifges(3,5) * qJDD(2) - pkin(8) * t158 - qJD(2) * t301 + t366 * t146 - t331 * t150 + t300 * t358;
t341 = -m(3) * t305 + t313 * mrSges(3,1) - mrSges(3,2) * t312 - t315 * t359 + t316 * t358 - t344;
t346 = mrSges(2,1) * t318 - mrSges(2,2) * t319 + Ifges(2,3) * qJDD(1) + pkin(1) * t341 + pkin(7) * t355 + t334 * t142 + t332 * t145;
t159 = m(2) * t318 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t336 + t341;
t149 = t156 * t334 + t157 * t332;
t147 = m(2) * t319 - mrSges(2,1) * t336 - qJDD(1) * mrSges(2,2) + t355;
t143 = mrSges(2,1) * g(3) + mrSges(2,3) * t319 + t336 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t149 - t367;
t140 = -mrSges(2,2) * g(3) - mrSges(2,3) * t318 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t336 - pkin(7) * t149 - t142 * t332 + t145 * t334;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t335 * t140 - t333 * t143 - pkin(6) * (t147 * t333 + t159 * t335), t140, t145, t146, t154, t174, -t228 * t254 + t347; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t333 * t140 + t335 * t143 + pkin(6) * (t147 * t335 - t333 * t159), t143, t142, t150, t152, t173, -t345; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t346, t346, t367, -t343, Ifges(5,3) * t275 + t337, -t340, Ifges(7,5) * t220 + Ifges(7,6) * t274 + Ifges(7,3) * t219 + t228 * t255 - t230 * t299 - t350;];
m_new  = t1;
