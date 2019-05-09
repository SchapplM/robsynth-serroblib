% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:00:42
% EndTime: 2019-05-05 23:01:06
% DurationCPUTime: 14.06s
% Computational Cost: add. (254200->342), mult. (548541->422), div. (0->0), fcn. (382367->10), ass. (0->138)
t301 = sin(qJ(1));
t305 = cos(qJ(1));
t279 = -t305 * g(1) - t301 * g(2);
t321 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t279;
t336 = 2 * qJD(5);
t335 = -pkin(1) - pkin(7);
t334 = mrSges(2,1) - mrSges(3,2);
t333 = Ifges(2,5) - Ifges(3,4);
t332 = (-Ifges(2,6) + Ifges(3,5));
t278 = t301 * g(1) - t305 * g(2);
t306 = qJD(1) ^ 2;
t320 = -t306 * qJ(2) + qJDD(2) - t278;
t253 = qJDD(1) * t335 + t320;
t300 = sin(qJ(3));
t304 = cos(qJ(3));
t244 = t300 * g(3) + t304 * t253;
t329 = qJD(1) * qJD(3);
t327 = t300 * t329;
t273 = t304 * qJDD(1) - t327;
t219 = (-t273 - t327) * pkin(8) + (-t300 * t304 * t306 + qJDD(3)) * pkin(3) + t244;
t245 = -t304 * g(3) + t300 * t253;
t272 = -t300 * qJDD(1) - t304 * t329;
t330 = qJD(1) * t304;
t277 = qJD(3) * pkin(3) - pkin(8) * t330;
t293 = t300 ^ 2;
t220 = -t293 * t306 * pkin(3) + t272 * pkin(8) - qJD(3) * t277 + t245;
t299 = sin(qJ(4));
t303 = cos(qJ(4));
t199 = t303 * t219 - t299 * t220;
t263 = (-t299 * t304 - t300 * t303) * qJD(1);
t230 = t263 * qJD(4) + t299 * t272 + t303 * t273;
t264 = (-t299 * t300 + t303 * t304) * qJD(1);
t286 = qJDD(3) + qJDD(4);
t287 = qJD(3) + qJD(4);
t185 = (t263 * t287 - t230) * qJ(5) + (t263 * t264 + t286) * pkin(4) + t199;
t200 = t299 * t219 + t303 * t220;
t229 = -t264 * qJD(4) + t303 * t272 - t299 * t273;
t250 = t287 * pkin(4) - t264 * qJ(5);
t259 = t263 ^ 2;
t187 = -t259 * pkin(4) + t229 * qJ(5) - t287 * t250 + t200;
t296 = sin(pkin(10));
t297 = cos(pkin(10));
t240 = t297 * t263 - t296 * t264;
t182 = t296 * t185 + t297 * t187 + t240 * t336;
t205 = t297 * t229 - t296 * t230;
t241 = t296 * t263 + t297 * t264;
t214 = -t240 * mrSges(6,1) + t241 * mrSges(6,2);
t232 = t287 * mrSges(6,1) - t241 * mrSges(6,3);
t215 = -t240 * pkin(5) - t241 * pkin(9);
t285 = t287 ^ 2;
t179 = -t285 * pkin(5) + t286 * pkin(9) + t240 * t215 + t182;
t223 = -t272 * pkin(3) + t277 * t330 + (-pkin(8) * t293 + t335) * t306 + t321;
t192 = -t229 * pkin(4) - t259 * qJ(5) + t264 * t250 + qJDD(5) + t223;
t206 = t296 * t229 + t297 * t230;
t183 = t192 + (-t240 * t287 - t206) * pkin(9) + (t241 * t287 - t205) * pkin(5);
t298 = sin(qJ(6));
t302 = cos(qJ(6));
t176 = -t298 * t179 + t302 * t183;
t227 = -t298 * t241 + t302 * t287;
t190 = t227 * qJD(6) + t302 * t206 + t298 * t286;
t204 = qJDD(6) - t205;
t228 = t302 * t241 + t298 * t287;
t207 = -t227 * mrSges(7,1) + t228 * mrSges(7,2);
t234 = qJD(6) - t240;
t208 = -t234 * mrSges(7,2) + t227 * mrSges(7,3);
t172 = m(7) * t176 + t204 * mrSges(7,1) - t190 * mrSges(7,3) - t228 * t207 + t234 * t208;
t177 = t302 * t179 + t298 * t183;
t189 = -t228 * qJD(6) - t298 * t206 + t302 * t286;
t209 = t234 * mrSges(7,1) - t228 * mrSges(7,3);
t173 = m(7) * t177 - t204 * mrSges(7,2) + t189 * mrSges(7,3) + t227 * t207 - t234 * t209;
t324 = -t298 * t172 + t302 * t173;
t159 = m(6) * t182 - t286 * mrSges(6,2) + t205 * mrSges(6,3) + t240 * t214 - t287 * t232 + t324;
t323 = -t297 * t185 + t296 * t187;
t181 = -0.2e1 * qJD(5) * t241 - t323;
t231 = -t287 * mrSges(6,2) + t240 * mrSges(6,3);
t178 = -t286 * pkin(5) - t285 * pkin(9) + (t336 + t215) * t241 + t323;
t318 = -m(7) * t178 + t189 * mrSges(7,1) - t190 * mrSges(7,2) + t227 * t208 - t228 * t209;
t168 = m(6) * t181 + t286 * mrSges(6,1) - t206 * mrSges(6,3) - t241 * t214 + t287 * t231 + t318;
t153 = t296 * t159 + t297 * t168;
t242 = -t263 * mrSges(5,1) + t264 * mrSges(5,2);
t249 = -t287 * mrSges(5,2) + t263 * mrSges(5,3);
t150 = m(5) * t199 + t286 * mrSges(5,1) - t230 * mrSges(5,3) - t264 * t242 + t287 * t249 + t153;
t251 = t287 * mrSges(5,1) - t264 * mrSges(5,3);
t325 = t297 * t159 - t296 * t168;
t151 = m(5) * t200 - t286 * mrSges(5,2) + t229 * mrSges(5,3) + t263 * t242 - t287 * t251 + t325;
t144 = t303 * t150 + t299 * t151;
t161 = t302 * t172 + t298 * t173;
t331 = qJD(1) * t300;
t271 = (mrSges(4,1) * t300 + mrSges(4,2) * t304) * qJD(1);
t275 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t331;
t141 = m(4) * t244 + qJDD(3) * mrSges(4,1) - t273 * mrSges(4,3) + qJD(3) * t275 - t271 * t330 + t144;
t276 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t326 = -t299 * t150 + t303 * t151;
t142 = m(4) * t245 - qJDD(3) * mrSges(4,2) + t272 * mrSges(4,3) - qJD(3) * t276 - t271 * t331 + t326;
t138 = -t300 * t141 + t304 * t142;
t137 = t304 * t141 + t300 * t142;
t258 = -qJDD(1) * pkin(1) + t320;
t319 = -m(3) * t258 + (t306 * mrSges(3,3)) - t137;
t317 = m(6) * t192 - t205 * mrSges(6,1) + t206 * mrSges(6,2) - t240 * t231 + t241 * t232 + t161;
t193 = Ifges(7,5) * t228 + Ifges(7,6) * t227 + Ifges(7,3) * t234;
t195 = Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t234;
t165 = -mrSges(7,1) * t178 + mrSges(7,3) * t177 + Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t204 - t228 * t193 + t234 * t195;
t194 = Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t234;
t166 = mrSges(7,2) * t178 - mrSges(7,3) * t176 + Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t204 + t227 * t193 - t234 * t194;
t210 = Ifges(6,5) * t241 + Ifges(6,6) * t240 + Ifges(6,3) * t287;
t211 = Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t287;
t145 = mrSges(6,2) * t192 - mrSges(6,3) * t181 + Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t286 - pkin(9) * t161 - t298 * t165 + t302 * t166 + t240 * t210 - t287 * t211;
t212 = Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t287;
t312 = mrSges(7,1) * t176 - mrSges(7,2) * t177 + Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t204 + t228 * t194 - t227 * t195;
t146 = -mrSges(6,1) * t192 + mrSges(6,3) * t182 + Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t286 - pkin(5) * t161 - t241 * t210 + t287 * t212 - t312;
t235 = Ifges(5,5) * t264 + Ifges(5,6) * t263 + Ifges(5,3) * t287;
t237 = Ifges(5,1) * t264 + Ifges(5,4) * t263 + Ifges(5,5) * t287;
t133 = -mrSges(5,1) * t223 + mrSges(5,3) * t200 + Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t286 - pkin(4) * t317 + qJ(5) * t325 + t296 * t145 + t297 * t146 - t264 * t235 + t287 * t237;
t236 = Ifges(5,4) * t264 + Ifges(5,2) * t263 + Ifges(5,6) * t287;
t139 = mrSges(5,2) * t223 - mrSges(5,3) * t199 + Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t286 - qJ(5) * t153 + t297 * t145 - t296 * t146 + t263 * t235 - t287 * t236;
t252 = t306 * t335 + t321;
t260 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t304 - Ifges(4,6) * t300) * qJD(1);
t262 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t304 - Ifges(4,4) * t300) * qJD(1);
t311 = m(5) * t223 - t229 * mrSges(5,1) + t230 * mrSges(5,2) - t263 * t249 + t264 * t251 + t317;
t130 = -mrSges(4,1) * t252 + mrSges(4,3) * t245 + Ifges(4,4) * t273 + Ifges(4,2) * t272 + Ifges(4,6) * qJDD(3) - pkin(3) * t311 + pkin(8) * t326 + qJD(3) * t262 + t303 * t133 + t299 * t139 - t260 * t330;
t261 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t304 - Ifges(4,2) * t300) * qJD(1);
t132 = mrSges(4,2) * t252 - mrSges(4,3) * t244 + Ifges(4,1) * t273 + Ifges(4,4) * t272 + Ifges(4,5) * qJDD(3) - pkin(8) * t144 - qJD(3) * t261 - t299 * t133 + t303 * t139 - t260 * t331;
t256 = t306 * pkin(1) - t321;
t316 = mrSges(3,2) * t258 - mrSges(3,3) * t256 + Ifges(3,1) * qJDD(1) - pkin(7) * t137 - t300 * t130 + t304 * t132;
t156 = -m(4) * t252 + t272 * mrSges(4,1) - t273 * mrSges(4,2) - t275 * t331 - t276 * t330 - t311;
t315 = -mrSges(3,1) * t256 - pkin(2) * t156 - pkin(7) * t138 - t304 * t130 - t300 * t132;
t314 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t206 - Ifges(6,6) * t205 - Ifges(6,3) * t286 - pkin(5) * t318 - pkin(9) * t324 - t302 * t165 - t298 * t166 - t241 * t211 + t240 * t212;
t309 = -m(3) * t256 + t306 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t156;
t313 = -mrSges(2,2) * t279 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t319) + qJ(2) * t309 + mrSges(2,1) * t278 + Ifges(2,3) * qJDD(1) + t316;
t310 = -mrSges(5,1) * t199 + mrSges(5,2) * t200 - Ifges(5,5) * t230 - Ifges(5,6) * t229 - Ifges(5,3) * t286 - pkin(4) * t153 - t264 * t236 + t263 * t237 + t314;
t308 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t273 + Ifges(4,6) * t272 + Ifges(4,3) * qJDD(3) + pkin(3) * t144 + t261 * t330 + t262 * t331 - t310;
t307 = -mrSges(3,1) * t258 - pkin(2) * t137 - t308;
t154 = m(2) * t279 - t306 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t309;
t136 = -m(3) * g(3) + t138;
t134 = m(2) * t278 - t306 * mrSges(2,2) + qJDD(1) * t334 + t319;
t129 = -t307 - qJ(2) * t136 + (t332 * t306) + t333 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t278;
t128 = mrSges(2,3) * t279 - pkin(1) * t136 + g(3) * t334 - qJDD(1) * t332 + t306 * t333 + t315;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t305 * t129 - t301 * t128 - pkin(6) * (t305 * t134 + t301 * t154), t129, t316, t132, t139, t145, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t129 + t305 * t128 + pkin(6) * (-t301 * t134 + t305 * t154), t128, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t306 * Ifges(3,5)) + t307, t130, t133, t146, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t313, t313, mrSges(3,2) * g(3) + t306 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t315, t308, -t310, -t314, t312;];
m_new  = t1;
