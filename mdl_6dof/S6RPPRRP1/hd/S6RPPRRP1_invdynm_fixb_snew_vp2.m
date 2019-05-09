% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:43:40
% EndTime: 2019-05-05 14:43:52
% DurationCPUTime: 8.20s
% Computational Cost: add. (140156->316), mult. (310507->382), div. (0->0), fcn. (213127->10), ass. (0->130)
t301 = qJD(1) ^ 2;
t290 = sin(pkin(10));
t292 = cos(pkin(10));
t295 = sin(qJ(4));
t298 = cos(qJ(4));
t314 = t290 * t295 - t292 * t298;
t262 = t314 * qJD(1);
t315 = t290 * t298 + t292 * t295;
t331 = t262 * qJD(4);
t252 = t315 * qJDD(1) - t331;
t263 = t315 * qJD(1);
t294 = sin(qJ(5));
t297 = cos(qJ(5));
t257 = t297 * qJD(4) - t294 * t263;
t219 = t257 * qJD(5) + t294 * qJDD(4) + t297 * t252;
t258 = t294 * qJD(4) + t297 * t263;
t226 = -t257 * mrSges(7,1) + t258 * mrSges(7,2);
t296 = sin(qJ(1));
t299 = cos(qJ(1));
t274 = t296 * g(1) - t299 * g(2);
t271 = qJDD(1) * pkin(1) + t274;
t275 = -t299 * g(1) - t296 * g(2);
t272 = -t301 * pkin(1) + t275;
t291 = sin(pkin(9));
t293 = cos(pkin(9));
t255 = t291 * t271 + t293 * t272;
t239 = -t301 * pkin(2) + qJDD(1) * qJ(3) + t255;
t289 = -g(3) + qJDD(2);
t329 = qJD(1) * qJD(3);
t333 = t292 * t289 - 0.2e1 * t290 * t329;
t339 = pkin(3) * t292;
t217 = (-pkin(7) * qJDD(1) + t301 * t339 - t239) * t290 + t333;
t224 = t290 * t289 + (t239 + 0.2e1 * t329) * t292;
t328 = qJDD(1) * t292;
t285 = t292 ^ 2;
t336 = t285 * t301;
t220 = -pkin(3) * t336 + pkin(7) * t328 + t224;
t192 = t295 * t217 + t298 * t220;
t250 = t262 * pkin(4) - t263 * pkin(8);
t300 = qJD(4) ^ 2;
t186 = -t300 * pkin(4) + qJDD(4) * pkin(8) - t262 * t250 + t192;
t284 = t290 ^ 2;
t254 = t293 * t271 - t291 * t272;
t317 = qJDD(3) - t254;
t222 = (-pkin(2) - t339) * qJDD(1) + (-qJ(3) + (-t284 - t285) * pkin(7)) * t301 + t317;
t330 = t263 * qJD(4);
t251 = -t314 * qJDD(1) - t330;
t189 = (-t252 + t331) * pkin(8) + (-t251 + t330) * pkin(4) + t222;
t180 = -t294 * t186 + t297 * t189;
t249 = qJDD(5) - t251;
t261 = qJD(5) + t262;
t176 = -0.2e1 * qJD(6) * t258 + (t257 * t261 - t219) * qJ(6) + (t257 * t258 + t249) * pkin(5) + t180;
t229 = -t261 * mrSges(7,2) + t257 * mrSges(7,3);
t327 = m(7) * t176 + t249 * mrSges(7,1) + t261 * t229;
t173 = -t219 * mrSges(7,3) - t258 * t226 + t327;
t181 = t297 * t186 + t294 * t189;
t201 = Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t261;
t202 = Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t261;
t203 = Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t261;
t218 = -t258 * qJD(5) + t297 * qJDD(4) - t294 * t252;
t231 = t261 * pkin(5) - t258 * qJ(6);
t256 = t257 ^ 2;
t179 = -t256 * pkin(5) + t218 * qJ(6) + 0.2e1 * qJD(6) * t257 - t261 * t231 + t181;
t200 = Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t261;
t311 = -mrSges(7,1) * t176 + mrSges(7,2) * t179 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t249 - t258 * t200;
t341 = mrSges(6,1) * t180 - mrSges(6,2) * t181 + Ifges(6,5) * t219 + Ifges(6,6) * t218 + Ifges(6,3) * t249 + pkin(5) * t173 + t258 * t201 - (t203 + t202) * t257 - t311;
t245 = t262 * mrSges(5,1) + t263 * mrSges(5,2);
t260 = qJD(4) * mrSges(5,1) - t263 * mrSges(5,3);
t227 = -t257 * mrSges(6,1) + t258 * mrSges(6,2);
t230 = -t261 * mrSges(6,2) + t257 * mrSges(6,3);
t165 = m(6) * t180 + t249 * mrSges(6,1) + t261 * t230 + (-t226 - t227) * t258 + (-mrSges(6,3) - mrSges(7,3)) * t219 + t327;
t326 = m(7) * t179 + t218 * mrSges(7,3) + t257 * t226;
t232 = t261 * mrSges(7,1) - t258 * mrSges(7,3);
t334 = -t261 * mrSges(6,1) + t258 * mrSges(6,3) - t232;
t338 = -mrSges(6,2) - mrSges(7,2);
t168 = m(6) * t181 + t218 * mrSges(6,3) + t257 * t227 + t338 * t249 + t334 * t261 + t326;
t322 = -t294 * t165 + t297 * t168;
t158 = m(5) * t192 - qJDD(4) * mrSges(5,2) + t251 * mrSges(5,3) - qJD(4) * t260 - t262 * t245 + t322;
t191 = t298 * t217 - t295 * t220;
t259 = -qJD(4) * mrSges(5,2) - t262 * mrSges(5,3);
t185 = -qJDD(4) * pkin(4) - t300 * pkin(8) + t263 * t250 - t191;
t183 = -t218 * pkin(5) - t256 * qJ(6) + t258 * t231 + qJDD(6) + t185;
t321 = -m(7) * t183 + t218 * mrSges(7,1) + t257 * t229;
t305 = -m(6) * t185 + t218 * mrSges(6,1) + t338 * t219 + t257 * t230 + t334 * t258 + t321;
t170 = m(5) * t191 + qJDD(4) * mrSges(5,1) - t252 * mrSges(5,3) + qJD(4) * t259 - t263 * t245 + t305;
t150 = t295 * t158 + t298 * t170;
t223 = -t290 * t239 + t333;
t198 = Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t261;
t199 = Ifges(6,5) * t258 + Ifges(6,6) * t257 + Ifges(6,3) * t261;
t312 = -mrSges(7,1) * t183 + mrSges(7,3) * t179 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t249 + t261 * t202;
t152 = Ifges(6,4) * t219 + Ifges(6,2) * t218 + Ifges(6,6) * t249 + t261 * t203 - mrSges(6,1) * t185 + mrSges(6,3) * t181 - pkin(5) * (t219 * mrSges(7,2) - t321) + qJ(6) * (-t249 * mrSges(7,2) - t261 * t232 + t326) + (-pkin(5) * t232 - t198 - t199) * t258 + t312;
t310 = mrSges(7,2) * t183 - mrSges(7,3) * t176 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t249 + t257 * t198;
t160 = mrSges(6,2) * t185 - mrSges(6,3) * t180 + Ifges(6,1) * t219 + Ifges(6,4) * t218 + Ifges(6,5) * t249 - qJ(6) * t173 + t257 * t199 + (-t200 - t201) * t261 + t310;
t237 = Ifges(5,4) * t263 - Ifges(5,2) * t262 + Ifges(5,6) * qJD(4);
t238 = Ifges(5,1) * t263 - Ifges(5,4) * t262 + Ifges(5,5) * qJD(4);
t306 = -mrSges(5,1) * t191 + mrSges(5,2) * t192 - Ifges(5,5) * t252 - Ifges(5,6) * t251 - Ifges(5,3) * qJDD(4) - pkin(4) * t305 - pkin(8) * t322 - t297 * t152 - t294 * t160 - t263 * t237 - t262 * t238;
t319 = Ifges(4,4) * t290 + Ifges(4,2) * t292;
t320 = Ifges(4,1) * t290 + Ifges(4,4) * t292;
t340 = -mrSges(4,1) * t223 + mrSges(4,2) * t224 - pkin(3) * t150 - (t290 * t319 - t292 * t320) * t301 + t306;
t337 = mrSges(4,2) * t290;
t313 = mrSges(4,3) * qJDD(1) + t301 * (-mrSges(4,1) * t292 + t337);
t148 = m(4) * t223 - t313 * t290 + t150;
t323 = t298 * t158 - t295 * t170;
t149 = m(4) * t224 + t313 * t292 + t323;
t324 = -t290 * t148 + t292 * t149;
t141 = m(3) * t255 - t301 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t324;
t235 = -qJDD(1) * pkin(2) - t301 * qJ(3) + t317;
t162 = t297 * t165 + t294 * t168;
t308 = m(5) * t222 - t251 * mrSges(5,1) + t252 * mrSges(5,2) + t262 * t259 + t263 * t260 + t162;
t304 = -m(4) * t235 + mrSges(4,1) * t328 - t308 + (t284 * t301 + t336) * mrSges(4,3);
t154 = -t301 * mrSges(3,2) + m(3) * t254 + t304 + (mrSges(3,1) - t337) * qJDD(1);
t137 = t291 * t141 + t293 * t154;
t143 = t292 * t148 + t290 * t149;
t318 = Ifges(4,5) * t290 + Ifges(4,6) * t292;
t332 = t301 * t318;
t325 = t293 * t141 - t291 * t154;
t236 = Ifges(5,5) * t263 - Ifges(5,6) * t262 + Ifges(5,3) * qJD(4);
t138 = mrSges(5,2) * t222 - mrSges(5,3) * t191 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - pkin(8) * t162 - qJD(4) * t237 - t294 * t152 + t297 * t160 - t262 * t236;
t144 = -mrSges(5,1) * t222 + mrSges(5,3) * t192 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * qJDD(4) - pkin(4) * t162 + qJD(4) * t238 - t263 * t236 - t341;
t131 = -mrSges(4,1) * t235 + mrSges(4,3) * t224 - pkin(3) * t308 + pkin(7) * t323 + t319 * qJDD(1) + t295 * t138 + t298 * t144 - t290 * t332;
t133 = mrSges(4,2) * t235 - mrSges(4,3) * t223 - pkin(7) * t150 + t320 * qJDD(1) + t298 * t138 - t295 * t144 + t292 * t332;
t309 = -mrSges(3,2) * t255 + qJ(3) * t324 + t292 * t131 + t290 * t133 + pkin(2) * (-qJDD(1) * t337 + t304) + mrSges(3,1) * t254 + Ifges(3,3) * qJDD(1);
t307 = mrSges(2,1) * t274 - mrSges(2,2) * t275 + Ifges(2,3) * qJDD(1) + pkin(1) * t137 + t309;
t135 = m(2) * t275 - t301 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t325;
t134 = m(2) * t274 + qJDD(1) * mrSges(2,1) - t301 * mrSges(2,2) + t137;
t129 = t301 * Ifges(3,5) - mrSges(3,1) * t289 + mrSges(3,3) * t255 - pkin(2) * t143 + (Ifges(3,6) - t318) * qJDD(1) + t340;
t128 = mrSges(3,2) * t289 - mrSges(3,3) * t254 + Ifges(3,5) * qJDD(1) - t301 * Ifges(3,6) - qJ(3) * t143 - t290 * t131 + t292 * t133;
t127 = -mrSges(2,2) * g(3) - mrSges(2,3) * t274 + Ifges(2,5) * qJDD(1) - t301 * Ifges(2,6) - qJ(2) * t137 + t293 * t128 - t291 * t129;
t126 = Ifges(2,6) * qJDD(1) + t301 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t275 + t291 * t128 + t293 * t129 - pkin(1) * (m(3) * t289 + t143) + qJ(2) * t325;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t299 * t127 - t296 * t126 - pkin(6) * (t299 * t134 + t296 * t135), t127, t128, t133, t138, t160, -t261 * t200 + t310; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t296 * t127 + t299 * t126 + pkin(6) * (-t296 * t134 + t299 * t135), t126, t129, t131, t144, t152, -t258 * t198 + t312; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t307, t307, t309, t318 * qJDD(1) - t340, -t306, t341, -t257 * t202 - t311;];
m_new  = t1;
