% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP2
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
% Datum: 2019-05-05 14:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:47:02
% EndTime: 2019-05-05 14:47:15
% DurationCPUTime: 8.02s
% Computational Cost: add. (137330->317), mult. (302866->382), div. (0->0), fcn. (207244->10), ass. (0->129)
t298 = qJD(1) ^ 2;
t288 = sin(pkin(10));
t290 = cos(pkin(10));
t293 = sin(qJ(4));
t295 = cos(qJ(4));
t310 = t288 * t293 - t290 * t295;
t258 = t310 * qJD(1);
t294 = sin(qJ(1));
t296 = cos(qJ(1));
t270 = t294 * g(1) - g(2) * t296;
t267 = qJDD(1) * pkin(1) + t270;
t271 = -g(1) * t296 - g(2) * t294;
t268 = -pkin(1) * t298 + t271;
t289 = sin(pkin(9));
t291 = cos(pkin(9));
t251 = t289 * t267 + t291 * t268;
t235 = -pkin(2) * t298 + qJDD(1) * qJ(3) + t251;
t287 = -g(3) + qJDD(2);
t325 = qJD(1) * qJD(3);
t329 = t290 * t287 - 0.2e1 * t288 * t325;
t336 = pkin(3) * t290;
t212 = (-pkin(7) * qJDD(1) + t298 * t336 - t235) * t288 + t329;
t220 = t288 * t287 + (t235 + 0.2e1 * t325) * t290;
t324 = qJDD(1) * t290;
t283 = t290 ^ 2;
t333 = t283 * t298;
t215 = -pkin(3) * t333 + pkin(7) * t324 + t220;
t192 = t293 * t212 + t295 * t215;
t311 = t288 * t295 + t290 * t293;
t259 = t311 * qJD(1);
t246 = pkin(4) * t258 - pkin(8) * t259;
t297 = qJD(4) ^ 2;
t187 = -pkin(4) * t297 + qJDD(4) * pkin(8) - t246 * t258 + t192;
t282 = t288 ^ 2;
t250 = t291 * t267 - t289 * t268;
t313 = qJDD(3) - t250;
t218 = (-pkin(2) - t336) * qJDD(1) + (-qJ(3) + (-t282 - t283) * pkin(7)) * t298 + t313;
t326 = t259 * qJD(4);
t247 = -t310 * qJDD(1) - t326;
t327 = t258 * qJD(4);
t248 = t311 * qJDD(1) - t327;
t189 = (-t248 + t327) * pkin(8) + (-t247 + t326) * pkin(4) + t218;
t292 = sin(qJ(5));
t337 = cos(qJ(5));
t183 = -t292 * t187 + t337 * t189;
t184 = t337 * t187 + t292 * t189;
t252 = -t337 * qJD(4) + t292 * t259;
t253 = t292 * qJD(4) + t337 * t259;
t257 = qJD(5) + t258;
t196 = Ifges(7,5) * t253 + Ifges(7,6) * t257 + Ifges(7,3) * t252;
t199 = Ifges(6,4) * t253 - Ifges(6,2) * t252 + Ifges(6,6) * t257;
t201 = Ifges(6,1) * t253 - Ifges(6,4) * t252 + Ifges(6,5) * t257;
t213 = t253 * qJD(5) - t337 * qJDD(4) + t292 * t248;
t214 = -t252 * qJD(5) + t292 * qJDD(4) + t337 * t248;
t223 = mrSges(7,1) * t252 - mrSges(7,3) * t253;
t245 = qJDD(5) - t247;
t222 = pkin(5) * t252 - qJ(6) * t253;
t256 = t257 ^ 2;
t179 = -pkin(5) * t256 + qJ(6) * t245 + 0.2e1 * qJD(6) * t257 - t222 * t252 + t184;
t181 = -t245 * pkin(5) - t256 * qJ(6) + t253 * t222 + qJDD(6) - t183;
t200 = Ifges(7,1) * t253 + Ifges(7,4) * t257 + Ifges(7,5) * t252;
t308 = mrSges(7,1) * t181 - mrSges(7,3) * t179 - Ifges(7,4) * t214 - Ifges(7,2) * t245 - Ifges(7,6) * t213 - t252 * t200;
t226 = -mrSges(7,2) * t252 + mrSges(7,3) * t257;
t318 = -m(7) * t181 + t245 * mrSges(7,1) + t257 * t226;
t229 = -mrSges(7,1) * t257 + mrSges(7,2) * t253;
t323 = m(7) * t179 + t245 * mrSges(7,3) + t257 * t229;
t339 = -(-t199 + t196) * t253 + mrSges(6,1) * t183 - mrSges(6,2) * t184 + Ifges(6,5) * t214 - Ifges(6,6) * t213 + Ifges(6,3) * t245 + pkin(5) * (-t214 * mrSges(7,2) - t253 * t223 + t318) + qJ(6) * (-t213 * mrSges(7,2) - t252 * t223 + t323) + t252 * t201 - t308;
t241 = mrSges(5,1) * t258 + mrSges(5,2) * t259;
t255 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t259;
t228 = mrSges(6,1) * t257 - mrSges(6,3) * t253;
t330 = -mrSges(6,1) * t252 - mrSges(6,2) * t253 - t223;
t335 = -mrSges(6,3) - mrSges(7,2);
t169 = m(6) * t184 - t245 * mrSges(6,2) + t335 * t213 - t257 * t228 + t330 * t252 + t323;
t227 = -mrSges(6,2) * t257 - mrSges(6,3) * t252;
t171 = m(6) * t183 + t245 * mrSges(6,1) + t335 * t214 + t257 * t227 + t330 * t253 + t318;
t319 = t337 * t169 - t171 * t292;
t157 = m(5) * t192 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t247 - qJD(4) * t255 - t241 * t258 + t319;
t191 = t295 * t212 - t293 * t215;
t254 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t258;
t186 = -qJDD(4) * pkin(4) - t297 * pkin(8) + t259 * t246 - t191;
t182 = -0.2e1 * qJD(6) * t253 + (t252 * t257 - t214) * qJ(6) + (t253 * t257 + t213) * pkin(5) + t186;
t176 = m(7) * t182 + mrSges(7,1) * t213 - t214 * mrSges(7,3) + t226 * t252 - t253 * t229;
t301 = -m(6) * t186 - t213 * mrSges(6,1) - mrSges(6,2) * t214 - t252 * t227 - t228 * t253 - t176;
t166 = m(5) * t191 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t248 + qJD(4) * t254 - t241 * t259 + t301;
t151 = t293 * t157 + t295 * t166;
t219 = -t235 * t288 + t329;
t317 = -mrSges(7,1) * t182 + mrSges(7,2) * t179;
t198 = Ifges(7,4) * t253 + Ifges(7,2) * t257 + Ifges(7,6) * t252;
t332 = -Ifges(6,5) * t253 + Ifges(6,6) * t252 - Ifges(6,3) * t257 - t198;
t159 = -mrSges(6,1) * t186 + mrSges(6,3) * t184 - pkin(5) * t176 + (t200 + t201) * t257 + t332 * t253 + (Ifges(6,6) - Ifges(7,6)) * t245 + (Ifges(6,4) - Ifges(7,5)) * t214 + (-Ifges(6,2) - Ifges(7,3)) * t213 + t317;
t307 = mrSges(7,2) * t181 - mrSges(7,3) * t182 + Ifges(7,1) * t214 + Ifges(7,4) * t245 + Ifges(7,5) * t213 + t257 * t196;
t161 = mrSges(6,2) * t186 - mrSges(6,3) * t183 + Ifges(6,1) * t214 - Ifges(6,4) * t213 + Ifges(6,5) * t245 - qJ(6) * t176 - t257 * t199 + t332 * t252 + t307;
t233 = Ifges(5,4) * t259 - Ifges(5,2) * t258 + Ifges(5,6) * qJD(4);
t234 = Ifges(5,1) * t259 - Ifges(5,4) * t258 + Ifges(5,5) * qJD(4);
t303 = -mrSges(5,1) * t191 + mrSges(5,2) * t192 - Ifges(5,5) * t248 - Ifges(5,6) * t247 - Ifges(5,3) * qJDD(4) - pkin(4) * t301 - pkin(8) * t319 - t337 * t159 - t292 * t161 - t259 * t233 - t258 * t234;
t315 = Ifges(4,4) * t288 + Ifges(4,2) * t290;
t316 = Ifges(4,1) * t288 + Ifges(4,4) * t290;
t338 = -mrSges(4,1) * t219 + mrSges(4,2) * t220 - pkin(3) * t151 - (t288 * t315 - t290 * t316) * t298 + t303;
t334 = mrSges(4,2) * t288;
t309 = mrSges(4,3) * qJDD(1) + t298 * (-mrSges(4,1) * t290 + t334);
t149 = m(4) * t219 - t309 * t288 + t151;
t320 = t295 * t157 - t293 * t166;
t150 = m(4) * t220 + t309 * t290 + t320;
t321 = -t149 * t288 + t290 * t150;
t141 = m(3) * t251 - mrSges(3,1) * t298 - qJDD(1) * mrSges(3,2) + t321;
t231 = -qJDD(1) * pkin(2) - t298 * qJ(3) + t313;
t163 = t292 * t169 + t337 * t171;
t305 = m(5) * t218 - t247 * mrSges(5,1) + t248 * mrSges(5,2) + t258 * t254 + t259 * t255 + t163;
t302 = -m(4) * t231 + mrSges(4,1) * t324 - t305 + (t282 * t298 + t333) * mrSges(4,3);
t153 = t302 + (mrSges(3,1) - t334) * qJDD(1) - t298 * mrSges(3,2) + m(3) * t250;
t138 = t289 * t141 + t291 * t153;
t143 = t290 * t149 + t288 * t150;
t314 = Ifges(4,5) * t288 + Ifges(4,6) * t290;
t328 = t298 * t314;
t322 = t291 * t141 - t153 * t289;
t232 = Ifges(5,5) * t259 - Ifges(5,6) * t258 + Ifges(5,3) * qJD(4);
t144 = mrSges(5,2) * t218 - mrSges(5,3) * t191 + Ifges(5,1) * t248 + Ifges(5,4) * t247 + Ifges(5,5) * qJDD(4) - pkin(8) * t163 - qJD(4) * t233 - t292 * t159 + t337 * t161 - t258 * t232;
t145 = -mrSges(5,1) * t218 + mrSges(5,3) * t192 + Ifges(5,4) * t248 + Ifges(5,2) * t247 + Ifges(5,6) * qJDD(4) - pkin(4) * t163 + qJD(4) * t234 - t259 * t232 - t339;
t132 = -mrSges(4,1) * t231 + mrSges(4,3) * t220 - pkin(3) * t305 + pkin(7) * t320 + t315 * qJDD(1) + t293 * t144 + t295 * t145 - t288 * t328;
t134 = mrSges(4,2) * t231 - mrSges(4,3) * t219 - pkin(7) * t151 + t316 * qJDD(1) + t295 * t144 - t293 * t145 + t290 * t328;
t306 = -mrSges(3,2) * t251 + qJ(3) * t321 + t290 * t132 + t288 * t134 + pkin(2) * (-qJDD(1) * t334 + t302) + mrSges(3,1) * t250 + Ifges(3,3) * qJDD(1);
t304 = mrSges(2,1) * t270 - mrSges(2,2) * t271 + Ifges(2,3) * qJDD(1) + pkin(1) * t138 + t306;
t136 = m(2) * t271 - mrSges(2,1) * t298 - qJDD(1) * mrSges(2,2) + t322;
t135 = m(2) * t270 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t298 + t138;
t130 = -pkin(2) * t143 + (Ifges(3,6) - t314) * qJDD(1) + t298 * Ifges(3,5) - mrSges(3,1) * t287 + mrSges(3,3) * t251 + t338;
t129 = mrSges(3,2) * t287 - mrSges(3,3) * t250 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t298 - qJ(3) * t143 - t132 * t288 + t134 * t290;
t128 = -mrSges(2,2) * g(3) - mrSges(2,3) * t270 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t298 - qJ(2) * t138 + t129 * t291 - t130 * t289;
t127 = Ifges(2,6) * qJDD(1) + t298 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t271 + t289 * t129 + t291 * t130 - pkin(1) * (m(3) * t287 + t143) + qJ(2) * t322;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t296 * t128 - t294 * t127 - pkin(6) * (t135 * t296 + t136 * t294), t128, t129, t134, t144, t161, -t198 * t252 + t307; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t294 * t128 + t296 * t127 + pkin(6) * (-t135 * t294 + t136 * t296), t127, t130, t132, t145, t159, -t253 * t196 - t308; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t304, t304, t306, t314 * qJDD(1) - t338, -t303, t339, Ifges(7,5) * t214 + Ifges(7,6) * t245 + Ifges(7,3) * t213 + t253 * t198 - t257 * t200 - t317;];
m_new  = t1;
