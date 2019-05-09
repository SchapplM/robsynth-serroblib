% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:05:23
% EndTime: 2019-05-05 16:05:41
% DurationCPUTime: 12.63s
% Computational Cost: add. (220412->320), mult. (507395->390), div. (0->0), fcn. (375871->10), ass. (0->138)
t303 = qJD(1) ^ 2;
t293 = sin(pkin(10));
t283 = t293 ^ 2;
t294 = cos(pkin(10));
t337 = t294 ^ 2 + t283;
t329 = t337 * mrSges(4,3);
t345 = t303 * t329;
t298 = sin(qJ(1));
t302 = cos(qJ(1));
t269 = t298 * g(1) - t302 * g(2);
t319 = -t303 * qJ(2) + qJDD(2) - t269;
t339 = -pkin(1) - qJ(3);
t344 = -(2 * qJD(1) * qJD(3)) + t339 * qJDD(1) + t319;
t270 = -t302 * g(1) - t298 * g(2);
t343 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t270;
t342 = pkin(3) * t303;
t341 = mrSges(2,1) - mrSges(3,2);
t340 = -Ifges(2,6) + Ifges(3,5);
t338 = Ifges(4,6) * t293;
t244 = t293 * g(3) + t344 * t294;
t223 = (-pkin(7) * qJDD(1) - t293 * t342) * t294 + t244;
t245 = -g(3) * t294 + t344 * t293;
t333 = qJDD(1) * t293;
t226 = -pkin(7) * t333 - t283 * t342 + t245;
t297 = sin(qJ(4));
t301 = cos(qJ(4));
t208 = t301 * t223 - t297 * t226;
t322 = -t293 * t297 + t294 * t301;
t323 = -t293 * t301 - t294 * t297;
t262 = t323 * qJD(1);
t335 = t262 * qJD(4);
t247 = t322 * qJDD(1) + t335;
t263 = t322 * qJD(1);
t188 = (-t247 + t335) * pkin(8) + (t262 * t263 + qJDD(4)) * pkin(4) + t208;
t209 = t297 * t223 + t301 * t226;
t246 = -t263 * qJD(4) + t323 * qJDD(1);
t254 = qJD(4) * pkin(4) - pkin(8) * t263;
t261 = t262 ^ 2;
t190 = -pkin(4) * t261 + pkin(8) * t246 - qJD(4) * t254 + t209;
t296 = sin(qJ(5));
t300 = cos(qJ(5));
t186 = t296 * t188 + t300 * t190;
t236 = t262 * t296 + t263 * t300;
t205 = -qJD(5) * t236 + t246 * t300 - t247 * t296;
t235 = t262 * t300 - t263 * t296;
t217 = -mrSges(6,1) * t235 + mrSges(6,2) * t236;
t285 = qJD(4) + qJD(5);
t228 = mrSges(6,1) * t285 - mrSges(6,3) * t236;
t282 = qJDD(4) + qJDD(5);
t218 = -pkin(5) * t235 - pkin(9) * t236;
t281 = t285 ^ 2;
t182 = -pkin(5) * t281 + pkin(9) * t282 + t218 * t235 + t186;
t318 = qJDD(3) + t343;
t230 = pkin(3) * t333 + (-t337 * pkin(7) + t339) * t303 + t318;
t199 = -t246 * pkin(4) - t261 * pkin(8) + t263 * t254 + t230;
t206 = qJD(5) * t235 + t246 * t296 + t247 * t300;
t183 = (-t235 * t285 - t206) * pkin(9) + (t236 * t285 - t205) * pkin(5) + t199;
t295 = sin(qJ(6));
t299 = cos(qJ(6));
t179 = -t182 * t295 + t183 * t299;
t224 = -t236 * t295 + t285 * t299;
t193 = qJD(6) * t224 + t206 * t299 + t282 * t295;
t204 = qJDD(6) - t205;
t225 = t236 * t299 + t285 * t295;
t210 = -mrSges(7,1) * t224 + mrSges(7,2) * t225;
t231 = qJD(6) - t235;
t211 = -mrSges(7,2) * t231 + mrSges(7,3) * t224;
t175 = m(7) * t179 + mrSges(7,1) * t204 - t193 * mrSges(7,3) - t210 * t225 + t211 * t231;
t180 = t182 * t299 + t183 * t295;
t192 = -qJD(6) * t225 - t206 * t295 + t282 * t299;
t212 = mrSges(7,1) * t231 - mrSges(7,3) * t225;
t176 = m(7) * t180 - mrSges(7,2) * t204 + t192 * mrSges(7,3) + t210 * t224 - t212 * t231;
t326 = -t175 * t295 + t299 * t176;
t162 = m(6) * t186 - mrSges(6,2) * t282 + mrSges(6,3) * t205 + t217 * t235 - t228 * t285 + t326;
t185 = t188 * t300 - t190 * t296;
t227 = -mrSges(6,2) * t285 + mrSges(6,3) * t235;
t181 = -pkin(5) * t282 - pkin(9) * t281 + t218 * t236 - t185;
t316 = -m(7) * t181 + t192 * mrSges(7,1) - t193 * mrSges(7,2) + t224 * t211 - t212 * t225;
t171 = m(6) * t185 + mrSges(6,1) * t282 - mrSges(6,3) * t206 - t217 * t236 + t227 * t285 + t316;
t156 = t296 * t162 + t300 * t171;
t240 = -mrSges(5,1) * t262 + mrSges(5,2) * t263;
t252 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t262;
t153 = m(5) * t208 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t247 + qJD(4) * t252 - t240 * t263 + t156;
t253 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t263;
t327 = t300 * t162 - t171 * t296;
t154 = m(5) * t209 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t246 - qJD(4) * t253 + t240 * t262 + t327;
t147 = t301 * t153 + t297 * t154;
t164 = t299 * t175 + t295 * t176;
t336 = t303 * (Ifges(4,5) * t294 - t338);
t332 = qJDD(1) * t294;
t330 = Ifges(3,4) + t338;
t321 = -mrSges(4,3) * qJDD(1) - t303 * (mrSges(4,1) * t293 + mrSges(4,2) * t294);
t144 = m(4) * t244 + t321 * t294 + t147;
t328 = -t297 * t153 + t301 * t154;
t145 = m(4) * t245 + t321 * t293 + t328;
t141 = -t144 * t293 + t294 * t145;
t325 = Ifges(4,1) * t294 - Ifges(4,4) * t293;
t324 = Ifges(4,4) * t294 - Ifges(4,2) * t293;
t140 = t294 * t144 + t293 * t145;
t260 = -qJDD(1) * pkin(1) + t319;
t317 = -m(3) * t260 + t303 * mrSges(3,3) - t140;
t315 = -m(6) * t199 + t205 * mrSges(6,1) - t206 * mrSges(6,2) + t235 * t227 - t236 * t228 - t164;
t194 = Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t231;
t196 = Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t231;
t168 = -mrSges(7,1) * t181 + mrSges(7,3) * t180 + Ifges(7,4) * t193 + Ifges(7,2) * t192 + Ifges(7,6) * t204 - t194 * t225 + t196 * t231;
t195 = Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t231;
t169 = mrSges(7,2) * t181 - mrSges(7,3) * t179 + Ifges(7,1) * t193 + Ifges(7,4) * t192 + Ifges(7,5) * t204 + t194 * t224 - t195 * t231;
t213 = Ifges(6,5) * t236 + Ifges(6,6) * t235 + Ifges(6,3) * t285;
t214 = Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * t285;
t148 = mrSges(6,2) * t199 - mrSges(6,3) * t185 + Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t282 - pkin(9) * t164 - t168 * t295 + t169 * t299 + t213 * t235 - t214 * t285;
t215 = Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * t285;
t310 = mrSges(7,1) * t179 - mrSges(7,2) * t180 + Ifges(7,5) * t193 + Ifges(7,6) * t192 + Ifges(7,3) * t204 + t195 * t225 - t196 * t224;
t149 = -mrSges(6,1) * t199 + mrSges(6,3) * t186 + Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t282 - pkin(5) * t164 - t213 * t236 + t215 * t285 - t310;
t232 = Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * qJD(4);
t234 = Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * qJD(4);
t136 = -mrSges(5,1) * t230 + mrSges(5,3) * t209 + Ifges(5,4) * t247 + Ifges(5,2) * t246 + Ifges(5,6) * qJDD(4) + pkin(4) * t315 + pkin(8) * t327 + qJD(4) * t234 + t296 * t148 + t300 * t149 - t263 * t232;
t233 = Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * qJD(4);
t142 = mrSges(5,2) * t230 - mrSges(5,3) * t208 + Ifges(5,1) * t247 + Ifges(5,4) * t246 + Ifges(5,5) * qJDD(4) - pkin(8) * t156 - qJD(4) * t233 + t148 * t300 - t149 * t296 + t232 * t262;
t251 = t339 * t303 + t318;
t309 = -m(5) * t230 + t246 * mrSges(5,1) - t247 * mrSges(5,2) + t262 * t252 - t263 * t253 + t315;
t133 = -mrSges(4,1) * t251 + mrSges(4,3) * t245 + pkin(3) * t309 + pkin(7) * t328 + t324 * qJDD(1) + t301 * t136 + t297 * t142 - t294 * t336;
t135 = mrSges(4,2) * t251 - mrSges(4,3) * t244 - pkin(7) * t147 + t325 * qJDD(1) - t297 * t136 + t301 * t142 - t293 * t336;
t256 = t303 * pkin(1) - t343;
t314 = mrSges(3,2) * t260 - mrSges(3,3) * t256 + Ifges(3,1) * qJDD(1) - qJ(3) * t140 - t133 * t293 + t294 * t135;
t308 = -m(4) * t251 - mrSges(4,1) * t333 - mrSges(4,2) * t332 + t309;
t313 = -mrSges(3,1) * t256 - pkin(2) * (t308 + t345) - qJ(3) * t141 - t294 * t133 - t293 * t135;
t312 = -mrSges(6,1) * t185 + mrSges(6,2) * t186 - Ifges(6,5) * t206 - Ifges(6,6) * t205 - Ifges(6,3) * t282 - pkin(5) * t316 - pkin(9) * t326 - t299 * t168 - t295 * t169 - t236 * t214 + t235 * t215;
t306 = -m(3) * t256 + t303 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t308;
t311 = -mrSges(2,2) * t270 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t317) + qJ(2) * (t306 - t345) + mrSges(2,1) * t269 + Ifges(2,3) * qJDD(1) + t314;
t307 = -mrSges(5,1) * t208 + mrSges(5,2) * t209 - Ifges(5,5) * t247 - Ifges(5,6) * t246 - Ifges(5,3) * qJDD(4) - pkin(4) * t156 - t263 * t233 + t262 * t234 + t312;
t305 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t332 + pkin(3) * t147 - t307 + (t293 * t325 + t294 * t324) * t303;
t304 = -mrSges(3,1) * t260 - pkin(2) * t140 - t305;
t157 = -qJDD(1) * mrSges(2,2) + t306 + m(2) * t270 + (-mrSges(2,1) - t329) * t303;
t139 = -m(3) * g(3) + t141;
t137 = m(2) * t269 - t303 * mrSges(2,2) + t341 * qJDD(1) + t317;
t132 = -t304 + t340 * t303 + (Ifges(2,5) - t330) * qJDD(1) - mrSges(2,3) * t269 - qJ(2) * t139 + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t131 = mrSges(2,3) * t270 - pkin(1) * t139 + (-Ifges(3,4) + Ifges(2,5)) * t303 - t340 * qJDD(1) + t341 * g(3) + t313;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t302 * t132 - t298 * t131 - pkin(6) * (t137 * t302 + t157 * t298), t132, t314, t135, t142, t148, t169; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t132 + t302 * t131 + pkin(6) * (-t137 * t298 + t157 * t302), t131, -mrSges(3,3) * g(3) - t303 * Ifges(3,5) + t330 * qJDD(1) + t304, t133, t136, t149, t168; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t311, t311, mrSges(3,2) * g(3) + t303 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t313, -Ifges(4,6) * t333 + t305, -t307, -t312, t310;];
m_new  = t1;
