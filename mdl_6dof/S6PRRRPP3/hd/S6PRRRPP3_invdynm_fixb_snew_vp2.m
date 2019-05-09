% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:49:47
% EndTime: 2019-05-05 06:50:04
% DurationCPUTime: 7.77s
% Computational Cost: add. (102481->343), mult. (193274->396), div. (0->0), fcn. (125323->10), ass. (0->134)
t295 = sin(pkin(10));
t297 = cos(pkin(10));
t284 = g(1) * t295 - g(2) * t297;
t285 = -g(1) * t297 - g(2) * t295;
t294 = -g(3) + qJDD(1);
t303 = cos(qJ(2));
t298 = cos(pkin(6));
t301 = sin(qJ(2));
t334 = t298 * t301;
t296 = sin(pkin(6));
t335 = t296 * t301;
t212 = t284 * t334 + t303 * t285 + t294 * t335;
t305 = qJD(2) ^ 2;
t208 = -pkin(2) * t305 + qJDD(2) * pkin(8) + t212;
t257 = -t284 * t296 + t294 * t298;
t300 = sin(qJ(3));
t302 = cos(qJ(3));
t200 = t302 * t208 + t300 * t257;
t280 = (-pkin(3) * t302 - pkin(9) * t300) * qJD(2);
t304 = qJD(3) ^ 2;
t328 = qJD(2) * t302;
t196 = -pkin(3) * t304 + qJDD(3) * pkin(9) + t280 * t328 + t200;
t211 = -t301 * t285 + (t284 * t298 + t294 * t296) * t303;
t207 = -qJDD(2) * pkin(2) - t305 * pkin(8) - t211;
t327 = qJD(2) * qJD(3);
t324 = t302 * t327;
t281 = qJDD(2) * t300 + t324;
t325 = t300 * t327;
t282 = qJDD(2) * t302 - t325;
t198 = (-t281 - t324) * pkin(9) + (-t282 + t325) * pkin(3) + t207;
t299 = sin(qJ(4));
t344 = cos(qJ(4));
t191 = -t299 * t196 + t344 * t198;
t329 = qJD(2) * t300;
t277 = -t344 * qJD(3) + t299 * t329;
t240 = -t277 * qJD(4) + t299 * qJDD(3) + t344 * t281;
t278 = t299 * qJD(3) + t344 * t329;
t244 = -mrSges(7,2) * t278 + mrSges(7,3) * t277;
t246 = mrSges(5,1) * t277 + mrSges(5,2) * t278;
t290 = -qJD(4) + t328;
t249 = mrSges(5,2) * t290 - mrSges(5,3) * t277;
t253 = mrSges(6,1) * t277 + mrSges(6,3) * t290;
t274 = qJDD(4) - t282;
t245 = pkin(4) * t277 - qJ(5) * t278;
t289 = t290 ^ 2;
t189 = -t274 * pkin(4) - t289 * qJ(5) + t278 * t245 + qJDD(5) - t191;
t337 = t277 * t290;
t345 = 2 * qJD(6);
t180 = t290 * t345 + (t277 * t278 - t274) * qJ(6) + (t240 - t337) * pkin(5) + t189;
t254 = -mrSges(7,1) * t277 - mrSges(7,2) * t290;
t322 = -m(7) * t180 + t274 * mrSges(7,3) - t290 * t254;
t341 = -mrSges(7,1) - mrSges(5,3);
t247 = -mrSges(6,2) * t277 - mrSges(6,3) * t278;
t349 = -m(6) * t189 - t240 * mrSges(6,1) - t278 * t247;
t169 = m(5) * t191 + (-t249 + t253) * t290 + (-t244 - t246) * t278 + (mrSges(5,1) - mrSges(6,2)) * t274 + t341 * t240 + t322 + t349;
t192 = t344 * t196 + t299 * t198;
t239 = qJD(4) * t278 - t344 * qJDD(3) + t281 * t299;
t250 = -mrSges(5,1) * t290 - mrSges(5,3) * t278;
t313 = -t289 * pkin(4) + t274 * qJ(5) - t277 * t245 + t192;
t187 = 0.2e1 * qJD(5) * t290 - t313;
t255 = mrSges(6,1) * t278 - mrSges(6,2) * t290;
t251 = pkin(5) * t278 + qJ(6) * t290;
t273 = t277 ^ 2;
t346 = -0.2e1 * qJD(5);
t183 = -t239 * pkin(5) - t273 * qJ(6) + qJDD(6) + (t346 - t251) * t290 + t313;
t252 = mrSges(7,1) * t278 + mrSges(7,3) * t290;
t326 = m(7) * t183 + t274 * mrSges(7,2) - t290 * t252;
t319 = -m(6) * t187 + t274 * mrSges(6,3) - t290 * t255 + t326;
t330 = -t244 - t247;
t171 = m(5) * t192 - t274 * mrSges(5,2) + t290 * t250 + (-t246 + t330) * t277 + (-mrSges(6,1) + t341) * t239 + t319;
t166 = -t169 * t299 + t344 * t171;
t279 = (-mrSges(4,1) * t302 + mrSges(4,2) * t300) * qJD(2);
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t329;
t164 = m(4) * t200 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t282 - qJD(3) * t286 + t279 * t328 + t166;
t199 = -t300 * t208 + t302 * t257;
t195 = -qJDD(3) * pkin(3) - t304 * pkin(9) + t280 * t329 - t199;
t310 = (-t240 - t337) * qJ(5) + t195 + (-pkin(4) * t290 + t346) * t278;
t186 = -t273 * pkin(5) + t277 * t345 - t278 * t251 + (pkin(4) + qJ(6)) * t239 + t310;
t177 = m(7) * t186 - t240 * mrSges(7,2) + t239 * mrSges(7,3) - t278 * t252 + t277 * t254;
t190 = t239 * pkin(4) + t310;
t314 = -m(6) * t190 + t239 * mrSges(6,2) + t277 * t253 - t177;
t172 = -m(5) * t195 - t239 * mrSges(5,1) - t277 * t249 + (-t250 + t255) * t278 + (-mrSges(5,2) + mrSges(6,3)) * t240 + t314;
t287 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t328;
t168 = m(4) * t199 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,3) + qJD(3) * t287 - t279 * t329 + t172;
t158 = t300 * t164 + t302 * t168;
t175 = -t240 * mrSges(6,3) - t278 * t255 - t314;
t213 = -Ifges(7,5) * t290 + Ifges(7,6) * t277 + Ifges(7,3) * t278;
t215 = Ifges(5,5) * t278 - Ifges(5,6) * t277 - Ifges(5,3) * t290;
t220 = -Ifges(6,1) * t290 - Ifges(6,4) * t278 + Ifges(6,5) * t277;
t219 = -Ifges(7,1) * t290 + Ifges(7,4) * t277 + Ifges(7,5) * t278;
t318 = mrSges(7,1) * t183 - mrSges(7,3) * t186 - Ifges(7,4) * t274 - Ifges(7,2) * t239 - Ifges(7,6) * t240 - t278 * t219;
t309 = mrSges(6,1) * t187 - mrSges(6,2) * t190 + pkin(5) * (t239 * mrSges(7,1) + t277 * t244 - t326) + qJ(6) * t177 - t318;
t217 = -Ifges(6,4) * t290 - Ifges(6,2) * t278 + Ifges(6,6) * t277;
t331 = Ifges(5,1) * t278 - Ifges(5,4) * t277 - Ifges(5,5) * t290 - t217;
t340 = Ifges(5,4) + Ifges(6,6);
t153 = mrSges(5,3) * t192 - mrSges(5,1) * t195 - t309 - pkin(4) * t175 + (-t213 - t331) * t290 + (-Ifges(5,2) - Ifges(6,3)) * t239 + t340 * t240 + (Ifges(5,6) - Ifges(6,5)) * t274 + (-t215 - t220) * t278;
t214 = -Ifges(6,5) * t290 - Ifges(6,6) * t278 + Ifges(6,3) * t277;
t218 = Ifges(5,4) * t278 - Ifges(5,2) * t277 - Ifges(5,6) * t290;
t176 = t240 * mrSges(7,1) + t278 * t244 - t322;
t216 = -Ifges(7,4) * t290 + Ifges(7,2) * t277 + Ifges(7,6) * t278;
t317 = -mrSges(7,1) * t180 + mrSges(7,2) * t186 - Ifges(7,5) * t274 - Ifges(7,6) * t239 - Ifges(7,3) * t240 + t290 * t216;
t312 = -mrSges(6,1) * t189 + mrSges(6,3) * t190 - pkin(5) * t176 + t317;
t332 = t219 + t220;
t159 = -t312 - mrSges(5,3) * t191 + mrSges(5,2) * t195 - qJ(5) * t175 - t340 * t239 + (Ifges(5,1) + Ifges(6,2)) * t240 + (Ifges(5,5) - Ifges(6,4)) * t274 + (-t215 - t332) * t277 + (t218 - t214) * t290;
t262 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t300 + Ifges(4,2) * t302) * qJD(2);
t263 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t300 + Ifges(4,4) * t302) * qJD(2);
t348 = mrSges(4,1) * t199 - mrSges(4,2) * t200 + Ifges(4,5) * t281 + Ifges(4,6) * t282 + Ifges(4,3) * qJDD(3) + pkin(3) * t172 + pkin(9) * t166 + (t262 * t300 - t263 * t302) * qJD(2) + t344 * t153 + t299 * t159;
t143 = -mrSges(3,1) * t257 + mrSges(3,3) * t212 + t305 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t158 - t348;
t323 = t302 * t164 - t168 * t300;
t156 = m(3) * t212 - mrSges(3,1) * t305 - qJDD(2) * mrSges(3,2) + t323;
t165 = t344 * t169 + t299 * t171;
t311 = -m(4) * t207 + t282 * mrSges(4,1) - t281 * mrSges(4,2) - t286 * t329 + t287 * t328 - t165;
t161 = m(3) * t211 + qJDD(2) * mrSges(3,1) - t305 * mrSges(3,2) + t311;
t152 = t303 * t156 - t161 * t301;
t350 = pkin(7) * t152 + t143 * t303;
t316 = -mrSges(7,2) * t183 + mrSges(7,3) * t180 - Ifges(7,1) * t274 - Ifges(7,4) * t239 - Ifges(7,5) * t240 - t277 * t213;
t308 = -mrSges(6,2) * t189 + mrSges(6,3) * t187 - Ifges(6,1) * t274 + Ifges(6,4) * t240 - Ifges(6,5) * t239 + qJ(6) * t176 + t278 * t214 + t316;
t347 = t331 * t277 + (t218 - t216) * t278 + mrSges(5,1) * t191 - mrSges(5,2) * t192 + Ifges(5,5) * t240 - Ifges(5,6) * t239 + Ifges(5,3) * t274 + pkin(4) * (-t274 * mrSges(6,2) + t290 * t253 - t176 + t349) + qJ(5) * (t330 * t277 + (-mrSges(6,1) - mrSges(7,1)) * t239 + t319) - t308;
t338 = t161 * t303;
t336 = t278 * t216;
t157 = m(3) * t257 + t158;
t147 = t156 * t334 - t157 * t296 + t298 * t338;
t261 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t300 + Ifges(4,6) * t302) * qJD(2);
t148 = mrSges(4,2) * t207 - mrSges(4,3) * t199 + Ifges(4,1) * t281 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - pkin(9) * t165 - qJD(3) * t262 - t299 * t153 + t344 * t159 + t261 * t328;
t149 = -mrSges(4,1) * t207 + mrSges(4,3) * t200 + Ifges(4,4) * t281 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t165 + qJD(3) * t263 - t261 * t329 - t347;
t139 = mrSges(3,1) * t211 - mrSges(3,2) * t212 + Ifges(3,3) * qJDD(2) + pkin(2) * t311 + pkin(8) * t323 + t300 * t148 + t302 * t149;
t141 = mrSges(3,2) * t257 - mrSges(3,3) * t211 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t305 - pkin(8) * t158 + t148 * t302 - t149 * t300;
t315 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + pkin(1) * t147 + t298 * t139 + t141 * t335 + t350 * t296;
t150 = m(2) * t285 + t152;
t146 = t298 * t157 + (t156 * t301 + t338) * t296;
t144 = m(2) * t284 + t147;
t137 = mrSges(2,2) * t294 - mrSges(2,3) * t284 + t303 * t141 - t301 * t143 + (-t146 * t296 - t147 * t298) * pkin(7);
t136 = -mrSges(2,1) * t294 + mrSges(2,3) * t285 - pkin(1) * t146 - t296 * t139 + (t141 * t301 + t350) * t298;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t297 * t137 - t295 * t136 - qJ(1) * (t144 * t297 + t150 * t295), t137, t141, t148, t159, -t277 * t217 - t308 - t336, -t316 - t336; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t295 * t137 + t297 * t136 + qJ(1) * (-t144 * t295 + t150 * t297), t136, t143, t149, t153, Ifges(6,4) * t274 - Ifges(6,2) * t240 + Ifges(6,6) * t239 + t290 * t214 + t332 * t277 + t312, t290 * t213 - t318; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, t139, t348, t347, t278 * t220 + Ifges(6,5) * t274 + Ifges(6,3) * t239 - Ifges(6,6) * t240 + t309 + (t213 - t217) * t290, -t277 * t219 - t317;];
m_new  = t1;
