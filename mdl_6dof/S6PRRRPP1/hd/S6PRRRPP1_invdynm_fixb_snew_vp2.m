% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-05-05 06:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:33:33
% EndTime: 2019-05-05 06:34:00
% DurationCPUTime: 14.92s
% Computational Cost: add. (274693->336), mult. (543998->418), div. (0->0), fcn. (382513->12), ass. (0->135)
t293 = sin(pkin(10));
t295 = cos(pkin(10));
t284 = t293 * g(1) - t295 * g(2);
t285 = -t295 * g(1) - t293 * g(2);
t291 = -g(3) + qJDD(1);
t302 = cos(qJ(2));
t296 = cos(pkin(6));
t299 = sin(qJ(2));
t328 = t296 * t299;
t294 = sin(pkin(6));
t329 = t294 * t299;
t241 = t284 * t328 + t302 * t285 + t291 * t329;
t304 = qJD(2) ^ 2;
t235 = -t304 * pkin(2) + qJDD(2) * pkin(8) + t241;
t259 = -t294 * t284 + t296 * t291;
t298 = sin(qJ(3));
t301 = cos(qJ(3));
t226 = t301 * t235 + t298 * t259;
t280 = (-pkin(3) * t301 - pkin(9) * t298) * qJD(2);
t303 = qJD(3) ^ 2;
t324 = t301 * qJD(2);
t200 = -t303 * pkin(3) + qJDD(3) * pkin(9) + t280 * t324 + t226;
t240 = -t299 * t285 + (t284 * t296 + t291 * t294) * t302;
t234 = -qJDD(2) * pkin(2) - t304 * pkin(8) - t240;
t323 = qJD(2) * qJD(3);
t320 = t301 * t323;
t281 = t298 * qJDD(2) + t320;
t321 = t298 * t323;
t282 = t301 * qJDD(2) - t321;
t204 = (-t281 - t320) * pkin(9) + (-t282 + t321) * pkin(3) + t234;
t297 = sin(qJ(4));
t300 = cos(qJ(4));
t193 = -t297 * t200 + t300 * t204;
t325 = qJD(2) * t298;
t277 = t300 * qJD(3) - t297 * t325;
t251 = t277 * qJD(4) + t297 * qJDD(3) + t300 * t281;
t274 = qJDD(4) - t282;
t278 = t297 * qJD(3) + t300 * t325;
t290 = qJD(4) - t324;
t190 = (t277 * t290 - t251) * qJ(5) + (t277 * t278 + t274) * pkin(4) + t193;
t194 = t300 * t200 + t297 * t204;
t250 = -t278 * qJD(4) + t300 * qJDD(3) - t297 * t281;
t257 = t290 * pkin(4) - t278 * qJ(5);
t273 = t277 ^ 2;
t192 = -t273 * pkin(4) + t250 * qJ(5) - t290 * t257 + t194;
t292 = sin(pkin(11));
t332 = cos(pkin(11));
t252 = -t332 * t277 + t292 * t278;
t335 = -2 * qJD(5);
t186 = t292 * t190 + t332 * t192 + t252 * t335;
t220 = -t332 * t250 + t292 * t251;
t253 = t292 * t277 + t332 * t278;
t237 = t290 * mrSges(6,1) - t253 * mrSges(6,3);
t227 = t252 * pkin(5) - t253 * qJ(6);
t289 = t290 ^ 2;
t181 = -t289 * pkin(5) + t274 * qJ(6) + 0.2e1 * qJD(6) * t290 - t252 * t227 + t186;
t238 = -t290 * mrSges(7,1) + t253 * mrSges(7,2);
t322 = m(7) * t181 + t274 * mrSges(7,3) + t290 * t238;
t228 = t252 * mrSges(7,1) - t253 * mrSges(7,3);
t326 = -t252 * mrSges(6,1) - t253 * mrSges(6,2) - t228;
t333 = -mrSges(6,3) - mrSges(7,2);
t168 = m(6) * t186 - t274 * mrSges(6,2) + t333 * t220 - t290 * t237 + t326 * t252 + t322;
t313 = t332 * t190 - t292 * t192;
t185 = t253 * t335 + t313;
t221 = t292 * t250 + t332 * t251;
t236 = -t290 * mrSges(6,2) - t252 * mrSges(6,3);
t183 = -t274 * pkin(5) - t289 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t227) * t253 - t313;
t239 = -t252 * mrSges(7,2) + t290 * mrSges(7,3);
t317 = -m(7) * t183 + t274 * mrSges(7,1) + t290 * t239;
t170 = m(6) * t185 + t274 * mrSges(6,1) + t333 * t221 + t290 * t236 + t326 * t253 + t317;
t163 = t292 * t168 + t332 * t170;
t254 = -t277 * mrSges(5,1) + t278 * mrSges(5,2);
t256 = -t290 * mrSges(5,2) + t277 * mrSges(5,3);
t161 = m(5) * t193 + t274 * mrSges(5,1) - t251 * mrSges(5,3) - t278 * t254 + t290 * t256 + t163;
t258 = t290 * mrSges(5,1) - t278 * mrSges(5,3);
t318 = t332 * t168 - t292 * t170;
t162 = m(5) * t194 - t274 * mrSges(5,2) + t250 * mrSges(5,3) + t277 * t254 - t290 * t258 + t318;
t159 = -t297 * t161 + t300 * t162;
t279 = (-mrSges(4,1) * t301 + mrSges(4,2) * t298) * qJD(2);
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t325;
t157 = m(4) * t226 - qJDD(3) * mrSges(4,2) + t282 * mrSges(4,3) - qJD(3) * t286 + t279 * t324 + t159;
t225 = -t298 * t235 + t301 * t259;
t199 = -qJDD(3) * pkin(3) - t303 * pkin(9) + t280 * t325 - t225;
t195 = -t250 * pkin(4) - t273 * qJ(5) + t278 * t257 + qJDD(5) + t199;
t188 = -0.2e1 * qJD(6) * t253 + (t252 * t290 - t221) * qJ(6) + (t253 * t290 + t220) * pkin(5) + t195;
t178 = m(7) * t188 + t220 * mrSges(7,1) - t221 * mrSges(7,3) - t253 * t238 + t252 * t239;
t309 = m(6) * t195 + t220 * mrSges(6,1) + t221 * mrSges(6,2) + t252 * t236 + t253 * t237 + t178;
t173 = -m(5) * t199 + t250 * mrSges(5,1) - t251 * mrSges(5,2) + t277 * t256 - t278 * t258 - t309;
t287 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t324;
t172 = m(4) * t225 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,3) + qJD(3) * t287 - t279 * t325 + t173;
t152 = t298 * t157 + t301 * t172;
t218 = Ifges(7,1) * t253 + Ifges(7,4) * t290 + Ifges(7,5) * t252;
t219 = Ifges(6,1) * t253 - Ifges(6,4) * t252 + Ifges(6,5) * t290;
t316 = -mrSges(7,1) * t188 + mrSges(7,2) * t181;
t216 = Ifges(7,4) * t253 + Ifges(7,2) * t290 + Ifges(7,6) * t252;
t327 = -Ifges(6,5) * t253 + Ifges(6,6) * t252 - Ifges(6,3) * t290 - t216;
t164 = -mrSges(6,1) * t195 + mrSges(6,3) * t186 - pkin(5) * t178 + (t218 + t219) * t290 + (Ifges(6,6) - Ifges(7,6)) * t274 + t327 * t253 + (Ifges(6,4) - Ifges(7,5)) * t221 + (-Ifges(6,2) - Ifges(7,3)) * t220 + t316;
t217 = Ifges(6,4) * t253 - Ifges(6,2) * t252 + Ifges(6,6) * t290;
t214 = Ifges(7,5) * t253 + Ifges(7,6) * t290 + Ifges(7,3) * t252;
t312 = mrSges(7,2) * t183 - mrSges(7,3) * t188 + Ifges(7,1) * t221 + Ifges(7,4) * t274 + Ifges(7,5) * t220 + t290 * t214;
t165 = mrSges(6,2) * t195 - mrSges(6,3) * t185 + Ifges(6,1) * t221 - Ifges(6,4) * t220 + Ifges(6,5) * t274 - qJ(6) * t178 - t290 * t217 + t327 * t252 + t312;
t243 = Ifges(5,5) * t278 + Ifges(5,6) * t277 + Ifges(5,3) * t290;
t245 = Ifges(5,1) * t278 + Ifges(5,4) * t277 + Ifges(5,5) * t290;
t146 = -mrSges(5,1) * t199 + mrSges(5,3) * t194 + Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * t274 - pkin(4) * t309 + qJ(5) * t318 + t332 * t164 + t292 * t165 - t278 * t243 + t290 * t245;
t244 = Ifges(5,4) * t278 + Ifges(5,2) * t277 + Ifges(5,6) * t290;
t147 = mrSges(5,2) * t199 - mrSges(5,3) * t193 + Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * t274 - qJ(5) * t163 - t292 * t164 + t332 * t165 + t277 * t243 - t290 * t244;
t264 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t298 + Ifges(4,2) * t301) * qJD(2);
t265 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t298 + Ifges(4,4) * t301) * qJD(2);
t336 = mrSges(4,1) * t225 - mrSges(4,2) * t226 + Ifges(4,5) * t281 + Ifges(4,6) * t282 + Ifges(4,3) * qJDD(3) + pkin(3) * t173 + pkin(9) * t159 + t300 * t146 + t297 * t147 + (t298 * t264 - t301 * t265) * qJD(2);
t136 = -mrSges(3,1) * t259 + mrSges(3,3) * t241 + t304 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t152 - t336;
t319 = t301 * t157 - t298 * t172;
t150 = m(3) * t241 - t304 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t319;
t158 = t300 * t161 + t297 * t162;
t308 = -m(4) * t234 + t282 * mrSges(4,1) - t281 * mrSges(4,2) - t286 * t325 + t287 * t324 - t158;
t154 = m(3) * t240 + qJDD(2) * mrSges(3,1) - t304 * mrSges(3,2) + t308;
t144 = t302 * t150 - t299 * t154;
t337 = pkin(7) * t144 + t136 * t302;
t330 = t154 * t302;
t151 = m(3) * t259 + t152;
t141 = t150 * t328 - t294 * t151 + t296 * t330;
t263 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t298 + Ifges(4,6) * t301) * qJD(2);
t137 = mrSges(4,2) * t234 - mrSges(4,3) * t225 + Ifges(4,1) * t281 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - pkin(9) * t158 - qJD(3) * t264 - t297 * t146 + t300 * t147 + t263 * t324;
t310 = mrSges(7,1) * t183 - mrSges(7,3) * t181 - Ifges(7,4) * t221 - Ifges(7,2) * t274 - Ifges(7,6) * t220 + t253 * t214 - t252 * t218;
t307 = mrSges(6,2) * t186 - t252 * t219 - qJ(6) * (-t220 * mrSges(7,2) - t252 * t228 + t322) - pkin(5) * (-t221 * mrSges(7,2) - t253 * t228 + t317) - mrSges(6,1) * t185 - t253 * t217 + Ifges(6,6) * t220 - Ifges(6,5) * t221 - Ifges(6,3) * t274 + t310;
t305 = mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t251 + Ifges(5,6) * t250 + Ifges(5,3) * t274 + pkin(4) * t163 + t278 * t244 - t277 * t245 - t307;
t145 = -mrSges(4,1) * t234 + mrSges(4,3) * t226 + Ifges(4,4) * t281 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t158 + qJD(3) * t265 - t263 * t325 - t305;
t132 = mrSges(3,1) * t240 - mrSges(3,2) * t241 + Ifges(3,3) * qJDD(2) + pkin(2) * t308 + pkin(8) * t319 + t298 * t137 + t301 * t145;
t134 = mrSges(3,2) * t259 - mrSges(3,3) * t240 + Ifges(3,5) * qJDD(2) - t304 * Ifges(3,6) - pkin(8) * t152 + t301 * t137 - t298 * t145;
t311 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + pkin(1) * t141 + t296 * t132 + t134 * t329 + t337 * t294;
t142 = m(2) * t285 + t144;
t140 = t296 * t151 + (t150 * t299 + t330) * t294;
t138 = m(2) * t284 + t141;
t130 = mrSges(2,2) * t291 - mrSges(2,3) * t284 + t302 * t134 - t299 * t136 + (-t140 * t294 - t141 * t296) * pkin(7);
t129 = -mrSges(2,1) * t291 + mrSges(2,3) * t285 - pkin(1) * t140 - t294 * t132 + (t134 * t299 + t337) * t296;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t295 * t130 - t293 * t129 - qJ(1) * (t295 * t138 + t293 * t142), t130, t134, t137, t147, t165, -t252 * t216 + t312; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t293 * t130 + t295 * t129 + qJ(1) * (-t293 * t138 + t295 * t142), t129, t136, t145, t146, t164, -t310; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t311, t311, t132, t336, t305, -t307, Ifges(7,5) * t221 + Ifges(7,6) * t274 + Ifges(7,3) * t220 + t253 * t216 - t290 * t218 - t316;];
m_new  = t1;
