% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:09:57
% EndTime: 2019-05-05 15:10:18
% DurationCPUTime: 18.77s
% Computational Cost: add. (342012->321), mult. (777457->401), div. (0->0), fcn. (574315->12), ass. (0->139)
t305 = qJD(1) ^ 2;
t300 = sin(qJ(1));
t304 = cos(qJ(1));
t273 = t300 * g(1) - g(2) * t304;
t270 = qJDD(1) * pkin(1) + t273;
t274 = -g(1) * t304 - g(2) * t300;
t271 = -pkin(1) * t305 + t274;
t294 = sin(pkin(10));
t296 = cos(pkin(10));
t256 = t294 * t270 + t296 * t271;
t245 = -pkin(2) * t305 + qJDD(1) * qJ(3) + t256;
t293 = sin(pkin(11));
t292 = -g(3) + qJDD(2);
t295 = cos(pkin(11));
t330 = qJD(1) * qJD(3);
t333 = t295 * t292 - 0.2e1 * t293 * t330;
t336 = pkin(3) * t295;
t227 = (-pkin(7) * qJDD(1) + t305 * t336 - t245) * t293 + t333;
t231 = t293 * t292 + (t245 + 0.2e1 * t330) * t295;
t329 = qJDD(1) * t295;
t287 = t295 ^ 2;
t334 = t287 * t305;
t228 = -pkin(3) * t334 + pkin(7) * t329 + t231;
t299 = sin(qJ(4));
t303 = cos(qJ(4));
t203 = t303 * t227 - t299 * t228;
t318 = t293 * t303 + t295 * t299;
t317 = -t293 * t299 + t295 * t303;
t261 = t317 * qJD(1);
t331 = t261 * qJD(4);
t253 = t318 * qJDD(1) + t331;
t262 = t318 * qJD(1);
t191 = (-t253 + t331) * pkin(8) + (t261 * t262 + qJDD(4)) * pkin(4) + t203;
t204 = t299 * t227 + t303 * t228;
t252 = -t262 * qJD(4) + t317 * qJDD(1);
t259 = qJD(4) * pkin(4) - pkin(8) * t262;
t260 = t261 ^ 2;
t193 = -pkin(4) * t260 + pkin(8) * t252 - qJD(4) * t259 + t204;
t298 = sin(qJ(5));
t302 = cos(qJ(5));
t189 = t298 * t191 + t302 * t193;
t244 = t261 * t298 + t262 * t302;
t212 = -qJD(5) * t244 + t252 * t302 - t253 * t298;
t243 = t261 * t302 - t262 * t298;
t222 = -mrSges(6,1) * t243 + mrSges(6,2) * t244;
t288 = qJD(4) + qJD(5);
t236 = mrSges(6,1) * t288 - mrSges(6,3) * t244;
t285 = qJDD(4) + qJDD(5);
t223 = -pkin(5) * t243 - pkin(9) * t244;
t284 = t288 ^ 2;
t185 = -pkin(5) * t284 + pkin(9) * t285 + t223 * t243 + t189;
t286 = t293 ^ 2;
t255 = t296 * t270 - t294 * t271;
t320 = qJDD(3) - t255;
t229 = (-pkin(2) - t336) * qJDD(1) + (-qJ(3) + (-t286 - t287) * pkin(7)) * t305 + t320;
t198 = -t252 * pkin(4) - t260 * pkin(8) + t262 * t259 + t229;
t213 = qJD(5) * t243 + t252 * t298 + t253 * t302;
t186 = (-t243 * t288 - t213) * pkin(9) + (t244 * t288 - t212) * pkin(5) + t198;
t297 = sin(qJ(6));
t301 = cos(qJ(6));
t182 = -t185 * t297 + t186 * t301;
t233 = -t244 * t297 + t288 * t301;
t196 = qJD(6) * t233 + t213 * t301 + t285 * t297;
t211 = qJDD(6) - t212;
t234 = t244 * t301 + t288 * t297;
t214 = -mrSges(7,1) * t233 + mrSges(7,2) * t234;
t238 = qJD(6) - t243;
t215 = -mrSges(7,2) * t238 + mrSges(7,3) * t233;
t178 = m(7) * t182 + mrSges(7,1) * t211 - mrSges(7,3) * t196 - t214 * t234 + t215 * t238;
t183 = t185 * t301 + t186 * t297;
t195 = -qJD(6) * t234 - t213 * t297 + t285 * t301;
t216 = mrSges(7,1) * t238 - mrSges(7,3) * t234;
t179 = m(7) * t183 - mrSges(7,2) * t211 + mrSges(7,3) * t195 + t214 * t233 - t216 * t238;
t324 = -t178 * t297 + t301 * t179;
t165 = m(6) * t189 - mrSges(6,2) * t285 + mrSges(6,3) * t212 + t222 * t243 - t236 * t288 + t324;
t188 = t191 * t302 - t193 * t298;
t235 = -mrSges(6,2) * t288 + mrSges(6,3) * t243;
t184 = -pkin(5) * t285 - pkin(9) * t284 + t223 * t244 - t188;
t313 = -m(7) * t184 + t195 * mrSges(7,1) - mrSges(7,2) * t196 + t233 * t215 - t216 * t234;
t174 = m(6) * t188 + mrSges(6,1) * t285 - mrSges(6,3) * t213 - t222 * t244 + t235 * t288 + t313;
t159 = t298 * t165 + t302 * t174;
t248 = -mrSges(5,1) * t261 + mrSges(5,2) * t262;
t257 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t261;
t156 = m(5) * t203 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t253 + qJD(4) * t257 - t248 * t262 + t159;
t258 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t262;
t325 = t302 * t165 - t174 * t298;
t157 = m(5) * t204 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t252 - qJD(4) * t258 + t248 * t261 + t325;
t150 = t303 * t156 + t299 * t157;
t230 = -t245 * t293 + t333;
t241 = Ifges(5,4) * t262 + Ifges(5,2) * t261 + Ifges(5,6) * qJD(4);
t242 = Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * qJD(4);
t200 = Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t238;
t202 = Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t238;
t171 = -mrSges(7,1) * t184 + mrSges(7,3) * t183 + Ifges(7,4) * t196 + Ifges(7,2) * t195 + Ifges(7,6) * t211 - t200 * t234 + t202 * t238;
t201 = Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t238;
t172 = mrSges(7,2) * t184 - mrSges(7,3) * t182 + Ifges(7,1) * t196 + Ifges(7,4) * t195 + Ifges(7,5) * t211 + t200 * t233 - t201 * t238;
t218 = Ifges(6,4) * t244 + Ifges(6,2) * t243 + Ifges(6,6) * t288;
t219 = Ifges(6,1) * t244 + Ifges(6,4) * t243 + Ifges(6,5) * t288;
t311 = -mrSges(6,1) * t188 + mrSges(6,2) * t189 - Ifges(6,5) * t213 - Ifges(6,6) * t212 - Ifges(6,3) * t285 - pkin(5) * t313 - pkin(9) * t324 - t301 * t171 - t297 * t172 - t244 * t218 + t243 * t219;
t307 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t253 - Ifges(5,6) * t252 - Ifges(5,3) * qJDD(4) - pkin(4) * t159 - t262 * t241 + t261 * t242 + t311;
t322 = Ifges(4,4) * t293 + Ifges(4,2) * t295;
t323 = Ifges(4,1) * t293 + Ifges(4,4) * t295;
t337 = -mrSges(4,1) * t230 + mrSges(4,2) * t231 - pkin(3) * t150 - (t293 * t322 - t295 * t323) * t305 + t307;
t335 = mrSges(4,2) * t293;
t316 = mrSges(4,3) * qJDD(1) + t305 * (-mrSges(4,1) * t295 + t335);
t148 = m(4) * t230 - t316 * t293 + t150;
t326 = -t299 * t156 + t303 * t157;
t149 = m(4) * t231 + t316 * t295 + t326;
t327 = -t148 * t293 + t295 * t149;
t141 = m(3) * t256 - mrSges(3,1) * t305 - qJDD(1) * mrSges(3,2) + t327;
t239 = -qJDD(1) * pkin(2) - t305 * qJ(3) + t320;
t167 = t301 * t178 + t297 * t179;
t315 = m(6) * t198 - t212 * mrSges(6,1) + t213 * mrSges(6,2) - t243 * t235 + t244 * t236 + t167;
t310 = m(5) * t229 - t252 * mrSges(5,1) + t253 * mrSges(5,2) - t261 * t257 + t262 * t258 + t315;
t308 = -m(4) * t239 + mrSges(4,1) * t329 - t310 + (t286 * t305 + t334) * mrSges(4,3);
t161 = (mrSges(3,1) - t335) * qJDD(1) + t308 - t305 * mrSges(3,2) + m(3) * t255;
t137 = t294 * t141 + t296 * t161;
t143 = t295 * t148 + t293 * t149;
t321 = Ifges(4,5) * t293 + Ifges(4,6) * t295;
t332 = t305 * t321;
t328 = t296 * t141 - t161 * t294;
t217 = Ifges(6,5) * t244 + Ifges(6,6) * t243 + Ifges(6,3) * t288;
t151 = mrSges(6,2) * t198 - mrSges(6,3) * t188 + Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t285 - pkin(9) * t167 - t171 * t297 + t172 * t301 + t217 * t243 - t218 * t288;
t309 = mrSges(7,1) * t182 - mrSges(7,2) * t183 + Ifges(7,5) * t196 + Ifges(7,6) * t195 + Ifges(7,3) * t211 + t201 * t234 - t202 * t233;
t152 = -mrSges(6,1) * t198 + mrSges(6,3) * t189 + Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t285 - pkin(5) * t167 - t217 * t244 + t219 * t288 - t309;
t240 = Ifges(5,5) * t262 + Ifges(5,6) * t261 + Ifges(5,3) * qJD(4);
t138 = -mrSges(5,1) * t229 + mrSges(5,3) * t204 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * qJDD(4) - pkin(4) * t315 + pkin(8) * t325 + qJD(4) * t242 + t298 * t151 + t302 * t152 - t262 * t240;
t144 = mrSges(5,2) * t229 - mrSges(5,3) * t203 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * qJDD(4) - pkin(8) * t159 - qJD(4) * t241 + t151 * t302 - t152 * t298 + t240 * t261;
t130 = -mrSges(4,1) * t239 + mrSges(4,3) * t231 - pkin(3) * t310 + pkin(7) * t326 + t322 * qJDD(1) + t303 * t138 + t299 * t144 - t293 * t332;
t132 = mrSges(4,2) * t239 - mrSges(4,3) * t230 - pkin(7) * t150 + t323 * qJDD(1) - t299 * t138 + t303 * t144 + t295 * t332;
t314 = -mrSges(3,2) * t256 + qJ(3) * t327 + t295 * t130 + t293 * t132 + pkin(2) * (-qJDD(1) * t335 + t308) + mrSges(3,1) * t255 + Ifges(3,3) * qJDD(1);
t312 = mrSges(2,1) * t273 - mrSges(2,2) * t274 + Ifges(2,3) * qJDD(1) + pkin(1) * t137 + t314;
t135 = m(2) * t274 - mrSges(2,1) * t305 - qJDD(1) * mrSges(2,2) + t328;
t134 = m(2) * t273 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t305 + t137;
t133 = (Ifges(3,6) - t321) * qJDD(1) + t305 * Ifges(3,5) - mrSges(3,1) * t292 + mrSges(3,3) * t256 - pkin(2) * t143 + t337;
t128 = mrSges(3,2) * t292 - mrSges(3,3) * t255 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t305 - qJ(3) * t143 - t130 * t293 + t132 * t295;
t127 = -mrSges(2,2) * g(3) - mrSges(2,3) * t273 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t305 - qJ(2) * t137 + t128 * t296 - t133 * t294;
t126 = Ifges(2,6) * qJDD(1) + t305 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t274 + t294 * t128 + t296 * t133 - pkin(1) * (m(3) * t292 + t143) + qJ(2) * t328;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t127 - t300 * t126 - pkin(6) * (t134 * t304 + t135 * t300), t127, t128, t132, t144, t151, t172; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t300 * t127 + t304 * t126 + pkin(6) * (-t134 * t300 + t135 * t304), t126, t133, t130, t138, t152, t171; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t312, t312, t314, t321 * qJDD(1) - t337, -t307, -t311, t309;];
m_new  = t1;
