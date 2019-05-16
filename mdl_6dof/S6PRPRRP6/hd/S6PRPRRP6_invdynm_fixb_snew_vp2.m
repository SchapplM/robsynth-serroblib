% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-05 00:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:59:19
% EndTime: 2019-05-04 23:59:33
% DurationCPUTime: 5.69s
% Computational Cost: add. (76597->292), mult. (139624->351), div. (0->0), fcn. (86797->10), ass. (0->120)
t271 = sin(pkin(10));
t273 = cos(pkin(10));
t254 = t271 * g(1) - t273 * g(2);
t255 = -t273 * g(1) - t271 * g(2);
t268 = -g(3) + qJDD(1);
t272 = sin(pkin(6));
t274 = cos(pkin(6));
t277 = sin(qJ(2));
t279 = cos(qJ(2));
t197 = -t277 * t255 + (t254 * t274 + t268 * t272) * t279;
t281 = qJD(2) ^ 2;
t288 = -t281 * qJ(3) + qJDD(3) - t197;
t319 = -pkin(2) - pkin(8);
t192 = t319 * qJDD(2) + t288;
t227 = -t272 * t254 + t274 * t268;
t276 = sin(qJ(4));
t278 = cos(qJ(4));
t187 = t276 * t192 + t278 * t227;
t249 = (mrSges(5,1) * t276 + mrSges(5,2) * t278) * qJD(2);
t303 = qJD(2) * qJD(4);
t299 = t278 * t303;
t251 = -t276 * qJDD(2) - t299;
t305 = qJD(2) * t278;
t257 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t305;
t250 = (pkin(4) * t276 - pkin(9) * t278) * qJD(2);
t280 = qJD(4) ^ 2;
t304 = t276 * qJD(2);
t182 = -t280 * pkin(4) + qJDD(4) * pkin(9) - t250 * t304 + t187;
t309 = t274 * t277;
t310 = t272 * t277;
t198 = t254 * t309 + t279 * t255 + t268 * t310;
t321 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t198;
t190 = t319 * t281 - t321;
t300 = t276 * t303;
t252 = t278 * qJDD(2) - t300;
t184 = (-t252 + t300) * pkin(9) + (-t251 + t299) * pkin(4) + t190;
t275 = sin(qJ(5));
t318 = cos(qJ(5));
t179 = t318 * t182 + t275 * t184;
t248 = t275 * qJD(4) + t318 * t305;
t213 = t248 * qJD(5) - t318 * qJDD(4) + t275 * t252;
t260 = qJD(5) + t304;
t223 = t260 * mrSges(6,1) - t248 * mrSges(6,3);
t244 = qJDD(5) - t251;
t247 = -t318 * qJD(4) + t275 * t305;
t217 = t247 * pkin(5) - t248 * qJ(6);
t259 = t260 ^ 2;
t174 = -t259 * pkin(5) + t244 * qJ(6) + 0.2e1 * qJD(6) * t260 - t247 * t217 + t179;
t224 = -t260 * mrSges(7,1) + t248 * mrSges(7,2);
t301 = m(7) * t174 + t244 * mrSges(7,3) + t260 * t224;
t218 = t247 * mrSges(7,1) - t248 * mrSges(7,3);
t306 = -t247 * mrSges(6,1) - t248 * mrSges(6,2) - t218;
t315 = -mrSges(6,3) - mrSges(7,2);
t164 = m(6) * t179 - t244 * mrSges(6,2) + t315 * t213 - t260 * t223 + t306 * t247 + t301;
t178 = -t275 * t182 + t318 * t184;
t214 = -t247 * qJD(5) + t275 * qJDD(4) + t318 * t252;
t222 = -t260 * mrSges(6,2) - t247 * mrSges(6,3);
t176 = -t244 * pkin(5) - t259 * qJ(6) + t248 * t217 + qJDD(6) - t178;
t225 = -t247 * mrSges(7,2) + t260 * mrSges(7,3);
t296 = -m(7) * t176 + t244 * mrSges(7,1) + t260 * t225;
t166 = m(6) * t178 + t244 * mrSges(6,1) + t315 * t214 + t260 * t222 + t306 * t248 + t296;
t298 = t318 * t164 - t275 * t166;
t152 = m(5) * t187 - qJDD(4) * mrSges(5,2) + t251 * mrSges(5,3) - qJD(4) * t257 - t249 * t304 + t298;
t186 = t278 * t192 - t276 * t227;
t256 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t304;
t181 = -qJDD(4) * pkin(4) - t280 * pkin(9) + t250 * t305 - t186;
t177 = -0.2e1 * qJD(6) * t248 + (t247 * t260 - t214) * qJ(6) + (t248 * t260 + t213) * pkin(5) + t181;
t171 = m(7) * t177 + t213 * mrSges(7,1) - t214 * mrSges(7,3) - t248 * t224 + t247 * t225;
t283 = -m(6) * t181 - t213 * mrSges(6,1) - t214 * mrSges(6,2) - t247 * t222 - t248 * t223 - t171;
t161 = m(5) * t186 + qJDD(4) * mrSges(5,1) - t252 * mrSges(5,3) + qJD(4) * t256 - t249 * t305 + t283;
t147 = t278 * t152 - t276 * t161;
t145 = m(4) * t227 + t147;
t203 = Ifges(7,1) * t248 + Ifges(7,4) * t260 + Ifges(7,5) * t247;
t204 = Ifges(6,1) * t248 - Ifges(6,4) * t247 + Ifges(6,5) * t260;
t295 = -mrSges(7,1) * t177 + mrSges(7,2) * t174;
t201 = Ifges(7,4) * t248 + Ifges(7,2) * t260 + Ifges(7,6) * t247;
t308 = -Ifges(6,5) * t248 + Ifges(6,6) * t247 - Ifges(6,3) * t260 - t201;
t154 = -mrSges(6,1) * t181 + mrSges(6,3) * t179 - pkin(5) * t171 + (t203 + t204) * t260 + t308 * t248 + (Ifges(6,6) - Ifges(7,6)) * t244 + (Ifges(6,4) - Ifges(7,5)) * t214 + (-Ifges(6,2) - Ifges(7,3)) * t213 + t295;
t202 = Ifges(6,4) * t248 - Ifges(6,2) * t247 + Ifges(6,6) * t260;
t199 = Ifges(7,5) * t248 + Ifges(7,6) * t260 + Ifges(7,3) * t247;
t292 = mrSges(7,2) * t176 - mrSges(7,3) * t177 + Ifges(7,1) * t214 + Ifges(7,4) * t244 + Ifges(7,5) * t213 + t260 * t199;
t157 = mrSges(6,2) * t181 - mrSges(6,3) * t178 + Ifges(6,1) * t214 - Ifges(6,4) * t213 + Ifges(6,5) * t244 - qJ(6) * t171 - t260 * t202 + t308 * t247 + t292;
t159 = t275 * t164 + t318 * t166;
t232 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t278 - Ifges(5,6) * t276) * qJD(2);
t233 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t278 - Ifges(5,2) * t276) * qJD(2);
t137 = mrSges(5,2) * t190 - mrSges(5,3) * t186 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - pkin(9) * t159 - qJD(4) * t233 - t275 * t154 + t318 * t157 - t232 * t304;
t234 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t278 - Ifges(5,4) * t276) * qJD(2);
t293 = mrSges(7,1) * t176 - mrSges(7,3) * t174 - Ifges(7,4) * t214 - Ifges(7,2) * t244 - Ifges(7,6) * t213 - t247 * t203;
t320 = -(-t202 + t199) * t248 + mrSges(6,1) * t178 - mrSges(6,2) * t179 + Ifges(6,5) * t214 - Ifges(6,6) * t213 + Ifges(6,3) * t244 + pkin(5) * (-t214 * mrSges(7,2) - t248 * t218 + t296) + qJ(6) * (-t213 * mrSges(7,2) - t247 * t218 + t301) + t247 * t204 - t293;
t141 = -mrSges(5,1) * t190 + mrSges(5,3) * t187 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * qJDD(4) - pkin(4) * t159 + qJD(4) * t234 - t232 * t305 - t320;
t155 = -m(5) * t190 + t251 * mrSges(5,1) - t252 * mrSges(5,2) - t256 * t304 - t257 * t305 - t159;
t193 = t281 * pkin(2) + t321;
t287 = -mrSges(4,1) * t193 - pkin(3) * t155 - pkin(8) * t147 - t276 * t137 - t278 * t141;
t313 = Ifges(4,5) - Ifges(3,6);
t314 = -Ifges(4,4) + Ifges(3,5);
t316 = mrSges(3,1) - mrSges(4,2);
t129 = mrSges(3,3) * t198 - pkin(2) * t145 - t313 * qJDD(2) - t316 * t227 + t314 * t281 + t287;
t146 = t276 * t152 + t278 * t161;
t195 = -qJDD(2) * pkin(2) + t288;
t291 = -m(4) * t195 + t281 * mrSges(4,3) - t146;
t143 = m(3) * t197 - t281 * mrSges(3,2) + t316 * qJDD(2) + t291;
t285 = -m(4) * t193 + t281 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t155;
t150 = m(3) * t198 - t281 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t285;
t140 = -t277 * t143 + t279 * t150;
t322 = pkin(7) * t140 + t129 * t279;
t311 = t143 * t279;
t144 = m(3) * t227 + t145;
t135 = -t272 * t144 + t150 * t309 + t274 * t311;
t289 = mrSges(4,2) * t195 - mrSges(4,3) * t193 + Ifges(4,1) * qJDD(2) - pkin(8) * t146 + t278 * t137 - t276 * t141;
t127 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t197 - mrSges(3,2) * t198 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t291) + qJ(3) * t285 + t289;
t286 = mrSges(5,1) * t186 - mrSges(5,2) * t187 + Ifges(5,5) * t252 + Ifges(5,6) * t251 + Ifges(5,3) * qJDD(4) + pkin(4) * t283 + pkin(9) * t298 + t318 * t154 + t275 * t157 + t233 * t305 + t234 * t304;
t284 = mrSges(4,1) * t195 + pkin(3) * t146 + t286;
t131 = -mrSges(3,3) * t197 - qJ(3) * t145 + t284 + t314 * qJDD(2) + (mrSges(3,2) - mrSges(4,3)) * t227 + t313 * t281;
t290 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + pkin(1) * t135 + t274 * t127 + t131 * t310 + t322 * t272;
t138 = m(2) * t255 + t140;
t134 = t274 * t144 + (t150 * t277 + t311) * t272;
t132 = m(2) * t254 + t135;
t125 = mrSges(2,2) * t268 - mrSges(2,3) * t254 - t277 * t129 + t279 * t131 + (-t134 * t272 - t135 * t274) * pkin(7);
t124 = -mrSges(2,1) * t268 + mrSges(2,3) * t255 - pkin(1) * t134 - t272 * t127 + (t131 * t277 + t322) * t274;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t273 * t125 - t271 * t124 - qJ(1) * (t273 * t132 + t271 * t138), t125, t131, t289, t137, t157, -t247 * t201 + t292; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t125 + t273 * t124 + qJ(1) * (-t271 * t132 + t273 * t138), t124, t129, mrSges(4,3) * t227 + Ifges(4,4) * qJDD(2) - t281 * Ifges(4,5) - t284, t141, t154, -t248 * t199 - t293; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t290, t290, t127, -mrSges(4,2) * t227 + t281 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t287, t286, t320, Ifges(7,5) * t214 + Ifges(7,6) * t244 + Ifges(7,3) * t213 + t248 * t201 - t260 * t203 - t295;];
m_new  = t1;
