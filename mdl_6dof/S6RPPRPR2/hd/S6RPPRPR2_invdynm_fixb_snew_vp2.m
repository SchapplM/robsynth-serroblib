% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-05-05 14:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:01:52
% EndTime: 2019-05-05 14:02:03
% DurationCPUTime: 6.94s
% Computational Cost: add. (108096->323), mult. (243414->386), div. (0->0), fcn. (161559->10), ass. (0->135)
t299 = sin(qJ(1));
t301 = cos(qJ(1));
t273 = t299 * g(1) - t301 * g(2);
t270 = qJDD(1) * pkin(1) + t273;
t274 = -t301 * g(1) - t299 * g(2);
t303 = qJD(1) ^ 2;
t271 = -t303 * pkin(1) + t274;
t294 = sin(pkin(9));
t296 = cos(pkin(9));
t248 = t294 * t270 + t296 * t271;
t227 = -t303 * pkin(2) + qJDD(1) * qJ(3) + t248;
t293 = sin(pkin(10));
t292 = -g(3) + qJDD(2);
t295 = cos(pkin(10));
t332 = qJD(1) * qJD(3);
t336 = t295 * t292 - 0.2e1 * t293 * t332;
t342 = pkin(3) * t295;
t205 = (-pkin(7) * qJDD(1) + t303 * t342 - t227) * t293 + t336;
t211 = t293 * t292 + (t227 + 0.2e1 * t332) * t295;
t330 = qJDD(1) * t295;
t286 = t295 ^ 2;
t339 = t286 * t303;
t208 = -pkin(3) * t339 + pkin(7) * t330 + t211;
t298 = sin(qJ(4));
t343 = cos(qJ(4));
t192 = t205 * t343 - t298 * t208;
t193 = t298 * t205 + t343 * t208;
t329 = t295 * t343;
t335 = qJD(1) * t293;
t261 = -qJD(1) * t329 + t298 * t335;
t318 = t293 * t343 + t295 * t298;
t262 = t318 * qJD(1);
t224 = Ifges(5,4) * t262 - Ifges(5,2) * t261 + Ifges(5,6) * qJD(4);
t234 = -t261 * mrSges(6,2) - t262 * mrSges(6,3);
t331 = qJDD(1) * t293;
t333 = t262 * qJD(4);
t244 = -qJDD(1) * t329 + t298 * t331 + t333;
t334 = t261 * qJD(4);
t245 = qJDD(1) * t318 - t334;
t254 = t261 * mrSges(6,1) - qJD(4) * mrSges(6,3);
t232 = t261 * pkin(4) - t262 * qJ(5);
t302 = qJD(4) ^ 2;
t188 = -qJDD(4) * pkin(4) - t302 * qJ(5) + t262 * t232 + qJDD(5) - t192;
t182 = (t261 * t262 - qJDD(4)) * pkin(8) + (t245 + t334) * pkin(5) + t188;
t256 = t262 * pkin(5) - qJD(4) * pkin(8);
t260 = t261 ^ 2;
t285 = t293 ^ 2;
t247 = t296 * t270 - t294 * t271;
t321 = qJDD(3) - t247;
t209 = (-pkin(2) - t342) * qJDD(1) + (-qJ(3) + (-t285 - t286) * pkin(7)) * t303 + t321;
t344 = -2 * qJD(5);
t305 = pkin(4) * t333 + t262 * t344 + (-t245 + t334) * qJ(5) + t209;
t185 = -t260 * pkin(5) - t262 * t256 + (pkin(4) + pkin(8)) * t244 + t305;
t297 = sin(qJ(6));
t300 = cos(qJ(6));
t179 = t300 * t182 - t297 * t185;
t249 = -t297 * qJD(4) + t300 * t261;
t207 = t249 * qJD(6) + t300 * qJDD(4) + t297 * t244;
t250 = t300 * qJD(4) + t297 * t261;
t214 = -t249 * mrSges(7,1) + t250 * mrSges(7,2);
t259 = qJD(6) + t262;
t217 = -t259 * mrSges(7,2) + t249 * mrSges(7,3);
t243 = qJDD(6) + t245;
t176 = m(7) * t179 + t243 * mrSges(7,1) - t207 * mrSges(7,3) - t250 * t214 + t259 * t217;
t180 = t297 * t182 + t300 * t185;
t206 = -t250 * qJD(6) - t297 * qJDD(4) + t300 * t244;
t218 = t259 * mrSges(7,1) - t250 * mrSges(7,3);
t177 = m(7) * t180 - t243 * mrSges(7,2) + t206 * mrSges(7,3) + t249 * t214 - t259 * t218;
t164 = t300 * t176 + t297 * t177;
t315 = -t302 * pkin(4) + qJDD(4) * qJ(5) - t261 * t232 + t193;
t184 = -t244 * pkin(5) - t260 * pkin(8) + ((2 * qJD(5)) + t256) * qJD(4) + t315;
t196 = Ifges(7,5) * t250 + Ifges(7,6) * t249 + Ifges(7,3) * t259;
t198 = Ifges(7,1) * t250 + Ifges(7,4) * t249 + Ifges(7,5) * t259;
t167 = -mrSges(7,1) * t184 + mrSges(7,3) * t180 + Ifges(7,4) * t207 + Ifges(7,2) * t206 + Ifges(7,6) * t243 - t250 * t196 + t259 * t198;
t197 = Ifges(7,4) * t250 + Ifges(7,2) * t249 + Ifges(7,6) * t259;
t168 = mrSges(7,2) * t184 - mrSges(7,3) * t179 + Ifges(7,1) * t207 + Ifges(7,4) * t206 + Ifges(7,5) * t243 + t249 * t196 - t259 * t197;
t186 = qJD(4) * t344 - t315;
t221 = Ifges(6,5) * qJD(4) - Ifges(6,6) * t262 + Ifges(6,3) * t261;
t312 = -mrSges(6,2) * t188 + mrSges(6,3) * t186 - Ifges(6,1) * qJDD(4) + Ifges(6,4) * t245 - Ifges(6,5) * t244 + pkin(8) * t164 + t297 * t167 - t300 * t168 + t262 * t221;
t181 = -m(7) * t184 + t206 * mrSges(7,1) - t207 * mrSges(7,2) + t249 * t217 - t250 * t218;
t255 = t262 * mrSges(6,1) + qJD(4) * mrSges(6,2);
t313 = -m(6) * t186 + qJDD(4) * mrSges(6,3) + qJD(4) * t255 - t181;
t316 = -m(6) * t188 - t245 * mrSges(6,1) - t262 * t234 - t164;
t223 = Ifges(6,4) * qJD(4) - Ifges(6,2) * t262 + Ifges(6,6) * t261;
t337 = Ifges(5,1) * t262 - Ifges(5,4) * t261 + Ifges(5,5) * qJD(4) - t223;
t348 = t337 * t261 - mrSges(5,2) * t193 + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t254 + t316) + qJ(5) * (-t244 * mrSges(6,1) - t261 * t234 + t313) + mrSges(5,1) * t192 + t262 * t224 - Ifges(5,6) * t244 + Ifges(5,5) * t245 + Ifges(5,3) * qJDD(4) - t312;
t233 = t261 * mrSges(5,1) + t262 * mrSges(5,2);
t252 = -qJD(4) * mrSges(5,2) - t261 * mrSges(5,3);
t160 = m(5) * t192 - t245 * mrSges(5,3) - t262 * t233 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t252 - t254) * qJD(4) + t316;
t253 = qJD(4) * mrSges(5,1) - t262 * mrSges(5,3);
t171 = m(5) * t193 - qJDD(4) * mrSges(5,2) - qJD(4) * t253 + (-t233 - t234) * t261 + (-mrSges(5,3) - mrSges(6,1)) * t244 + t313;
t155 = t343 * t160 + t298 * t171;
t210 = -t293 * t227 + t336;
t323 = Ifges(4,4) * t293 + Ifges(4,2) * t295;
t324 = Ifges(4,1) * t293 + Ifges(4,4) * t295;
t346 = qJD(1) * t295;
t347 = -mrSges(4,1) * t210 + mrSges(4,2) * t211 - pkin(3) * t155 - (t323 * t335 - t324 * t346) * qJD(1) - t348;
t341 = Ifges(5,4) + Ifges(6,6);
t340 = mrSges(4,2) * t293;
t319 = mrSges(4,3) * qJDD(1) + t303 * (-mrSges(4,1) * t295 + t340);
t153 = m(4) * t210 - t293 * t319 + t155;
t325 = -t298 * t160 + t343 * t171;
t154 = m(4) * t211 + t295 * t319 + t325;
t326 = -t293 * t153 + t295 * t154;
t146 = m(3) * t248 - t303 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t326;
t220 = -qJDD(1) * pkin(2) - t303 * qJ(3) + t321;
t165 = -t297 * t176 + t300 * t177;
t190 = t244 * pkin(4) + t305;
t161 = m(6) * t190 - t244 * mrSges(6,2) - t245 * mrSges(6,3) - t261 * t254 - t262 * t255 + t165;
t309 = m(5) * t209 + t244 * mrSges(5,1) + t245 * mrSges(5,2) + t261 * t252 + t262 * t253 + t161;
t307 = -m(4) * t220 + mrSges(4,1) * t330 - t309 + (t285 * t303 + t339) * mrSges(4,3);
t157 = (mrSges(3,1) - t340) * qJDD(1) + m(3) * t247 + t307 - t303 * mrSges(3,2);
t142 = t294 * t146 + t296 * t157;
t148 = t295 * t153 + t293 * t154;
t225 = Ifges(6,1) * qJD(4) - Ifges(6,4) * t262 + Ifges(6,5) * t261;
t338 = -Ifges(5,5) * t262 + Ifges(5,6) * t261 - Ifges(5,3) * qJD(4) - t225;
t327 = t296 * t146 - t294 * t157;
t322 = Ifges(4,5) * t293 + Ifges(4,6) * t295;
t310 = -mrSges(6,1) * t186 + mrSges(6,2) * t190 - pkin(5) * t181 - pkin(8) * t165 - t300 * t167 - t297 * t168;
t143 = -mrSges(5,1) * t209 + mrSges(5,3) * t193 - pkin(4) * t161 + t338 * t262 + t341 * t245 + (-Ifges(5,2) - Ifges(6,3)) * t244 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + t337 * qJD(4) + t310;
t314 = mrSges(7,1) * t179 - mrSges(7,2) * t180 + Ifges(7,5) * t207 + Ifges(7,6) * t206 + Ifges(7,3) * t243 + t250 * t197 - t249 * t198;
t308 = mrSges(6,1) * t188 - mrSges(6,3) * t190 + pkin(5) * t164 + t314;
t149 = t308 + mrSges(5,2) * t209 - mrSges(5,3) * t192 - qJ(5) * t161 + t338 * t261 + (Ifges(5,1) + Ifges(6,2)) * t245 - t341 * t244 + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) + (-t224 + t221) * qJD(4);
t267 = t322 * qJD(1);
t135 = -mrSges(4,1) * t220 + mrSges(4,3) * t211 - pkin(3) * t309 + pkin(7) * t325 + qJDD(1) * t323 + t143 * t343 + t298 * t149 - t267 * t335;
t138 = mrSges(4,2) * t220 - mrSges(4,3) * t210 - pkin(7) * t155 + qJDD(1) * t324 - t298 * t143 + t149 * t343 + t267 * t346;
t317 = -mrSges(3,2) * t248 + qJ(3) * t326 + t295 * t135 + t293 * t138 + pkin(2) * (-mrSges(4,2) * t331 + t307) + mrSges(3,1) * t247 + Ifges(3,3) * qJDD(1);
t311 = mrSges(2,1) * t273 - mrSges(2,2) * t274 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t317;
t140 = m(2) * t274 - t303 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t327;
t139 = m(2) * t273 + qJDD(1) * mrSges(2,1) - t303 * mrSges(2,2) + t142;
t136 = -pkin(2) * t148 + (Ifges(3,6) - t322) * qJDD(1) + mrSges(3,3) * t248 - mrSges(3,1) * t292 + t303 * Ifges(3,5) + t347;
t133 = mrSges(3,2) * t292 - mrSges(3,3) * t247 + Ifges(3,5) * qJDD(1) - t303 * Ifges(3,6) - qJ(3) * t148 - t293 * t135 + t295 * t138;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t273 + Ifges(2,5) * qJDD(1) - t303 * Ifges(2,6) - qJ(2) * t142 + t296 * t133 - t294 * t136;
t131 = Ifges(2,6) * qJDD(1) + t303 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t274 + t294 * t133 + t296 * t136 - pkin(1) * (m(3) * t292 + t148) + qJ(2) * t327;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t301 * t132 - t299 * t131 - pkin(6) * (t301 * t139 + t299 * t140), t132, t133, t138, t149, -t261 * t223 - t312, t168; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t299 * t132 + t301 * t131 + pkin(6) * (-t299 * t139 + t301 * t140), t131, t136, t135, t143, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t245 + Ifges(6,6) * t244 - qJD(4) * t221 + t261 * t225 - t308, t167; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t311, t311, t317, qJDD(1) * t322 - t347, t348, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t245 + Ifges(6,3) * t244 + qJD(4) * t223 + t262 * t225 - t310, t314;];
m_new  = t1;
