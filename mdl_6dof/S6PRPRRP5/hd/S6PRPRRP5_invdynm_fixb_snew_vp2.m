% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRP5
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
% Datum: 2019-05-04 23:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:54:17
% EndTime: 2019-05-04 23:54:30
% DurationCPUTime: 5.63s
% Computational Cost: add. (78407->291), mult. (144144->351), div. (0->0), fcn. (89962->10), ass. (0->121)
t273 = sin(pkin(10));
t275 = cos(pkin(10));
t259 = g(1) * t273 - g(2) * t275;
t260 = -g(1) * t275 - g(2) * t273;
t270 = -g(3) + qJDD(1);
t274 = sin(pkin(6));
t276 = cos(pkin(6));
t279 = sin(qJ(2));
t282 = cos(qJ(2));
t199 = -t279 * t260 + (t259 * t276 + t270 * t274) * t282;
t284 = qJD(2) ^ 2;
t291 = -t284 * qJ(3) + qJDD(3) - t199;
t321 = -pkin(2) - pkin(8);
t193 = t321 * qJDD(2) + t291;
t231 = -t259 * t274 + t270 * t276;
t278 = sin(qJ(4));
t281 = cos(qJ(4));
t187 = t278 * t193 + t281 * t231;
t254 = (mrSges(5,1) * t278 + mrSges(5,2) * t281) * qJD(2);
t307 = qJD(2) * qJD(4);
t302 = t281 * t307;
t256 = -qJDD(2) * t278 - t302;
t308 = qJD(2) * t281;
t262 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t308;
t255 = (pkin(4) * t278 - pkin(9) * t281) * qJD(2);
t283 = qJD(4) ^ 2;
t309 = qJD(2) * t278;
t181 = -pkin(4) * t283 + qJDD(4) * pkin(9) - t255 * t309 + t187;
t312 = t276 * t279;
t313 = t274 * t279;
t200 = t259 * t312 + t282 * t260 + t270 * t313;
t323 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t200;
t190 = t321 * t284 - t323;
t303 = t278 * t307;
t257 = qJDD(2) * t281 - t303;
t184 = (-t257 + t303) * pkin(9) + (-t256 + t302) * pkin(4) + t190;
t277 = sin(qJ(5));
t280 = cos(qJ(5));
t175 = -t277 * t181 + t280 * t184;
t252 = qJD(4) * t280 - t277 * t308;
t219 = qJD(5) * t252 + qJDD(4) * t277 + t257 * t280;
t253 = qJD(4) * t277 + t280 * t308;
t221 = -mrSges(7,1) * t252 + mrSges(7,2) * t253;
t222 = -mrSges(6,1) * t252 + mrSges(6,2) * t253;
t264 = qJD(5) + t309;
t226 = -mrSges(6,2) * t264 + mrSges(6,3) * t252;
t249 = qJDD(5) - t256;
t171 = -0.2e1 * qJD(6) * t253 + (t252 * t264 - t219) * qJ(6) + (t252 * t253 + t249) * pkin(5) + t175;
t225 = -mrSges(7,2) * t264 + mrSges(7,3) * t252;
t305 = m(7) * t171 + t249 * mrSges(7,1) + t264 * t225;
t161 = m(6) * t175 + t249 * mrSges(6,1) + t264 * t226 + (-t221 - t222) * t253 + (-mrSges(6,3) - mrSges(7,3)) * t219 + t305;
t176 = t280 * t181 + t277 * t184;
t218 = -qJD(5) * t253 + qJDD(4) * t280 - t257 * t277;
t227 = pkin(5) * t264 - qJ(6) * t253;
t248 = t252 ^ 2;
t174 = -pkin(5) * t248 + qJ(6) * t218 + 0.2e1 * qJD(6) * t252 - t227 * t264 + t176;
t304 = m(7) * t174 + t218 * mrSges(7,3) + t252 * t221;
t228 = mrSges(7,1) * t264 - mrSges(7,3) * t253;
t310 = -mrSges(6,1) * t264 + mrSges(6,3) * t253 - t228;
t318 = -mrSges(6,2) - mrSges(7,2);
t164 = m(6) * t176 + t218 * mrSges(6,3) + t252 * t222 + t318 * t249 + t310 * t264 + t304;
t301 = -t161 * t277 + t280 * t164;
t153 = m(5) * t187 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t256 - qJD(4) * t262 - t254 * t309 + t301;
t186 = t193 * t281 - t278 * t231;
t261 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t309;
t180 = -qJDD(4) * pkin(4) - pkin(9) * t283 + t255 * t308 - t186;
t178 = -pkin(5) * t218 - qJ(6) * t248 + t227 * t253 + qJDD(6) + t180;
t299 = -m(7) * t178 + t218 * mrSges(7,1) + t252 * t225;
t287 = -m(6) * t180 + t218 * mrSges(6,1) + t318 * t219 + t252 * t226 + t310 * t253 + t299;
t165 = m(5) * t186 + qJDD(4) * mrSges(5,1) - t257 * mrSges(5,3) + qJD(4) * t261 - t254 * t308 + t287;
t146 = t281 * t153 - t165 * t278;
t144 = m(4) * t231 + t146;
t201 = Ifges(7,5) * t253 + Ifges(7,6) * t252 + Ifges(7,3) * t264;
t202 = Ifges(6,5) * t253 + Ifges(6,6) * t252 + Ifges(6,3) * t264;
t206 = Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t264;
t205 = Ifges(7,1) * t253 + Ifges(7,4) * t252 + Ifges(7,5) * t264;
t297 = -mrSges(7,1) * t178 + mrSges(7,3) * t174 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t249 + t264 * t205;
t148 = Ifges(6,4) * t219 + Ifges(6,2) * t218 + Ifges(6,6) * t249 + t264 * t206 - mrSges(6,1) * t180 + mrSges(6,3) * t176 - pkin(5) * (t219 * mrSges(7,2) - t299) + qJ(6) * (-t249 * mrSges(7,2) - t264 * t228 + t304) + (-pkin(5) * t228 - t201 - t202) * t253 + t297;
t168 = -t219 * mrSges(7,3) - t253 * t221 + t305;
t203 = Ifges(7,4) * t253 + Ifges(7,2) * t252 + Ifges(7,6) * t264;
t204 = Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t264;
t295 = mrSges(7,2) * t178 - mrSges(7,3) * t171 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t249 + t252 * t201;
t156 = mrSges(6,2) * t180 - mrSges(6,3) * t175 + Ifges(6,1) * t219 + Ifges(6,4) * t218 + Ifges(6,5) * t249 - qJ(6) * t168 + t252 * t202 + (-t203 - t204) * t264 + t295;
t158 = t280 * t161 + t277 * t164;
t236 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t281 - Ifges(5,6) * t278) * qJD(2);
t237 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t281 - Ifges(5,2) * t278) * qJD(2);
t136 = mrSges(5,2) * t190 - mrSges(5,3) * t186 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * qJDD(4) - pkin(9) * t158 - qJD(4) * t237 - t148 * t277 + t156 * t280 - t236 * t309;
t238 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t281 - Ifges(5,4) * t278) * qJD(2);
t296 = -mrSges(7,1) * t171 + mrSges(7,2) * t174 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t249 - t253 * t203;
t322 = mrSges(6,1) * t175 - mrSges(6,2) * t176 + Ifges(6,5) * t219 + Ifges(6,6) * t218 + Ifges(6,3) * t249 + pkin(5) * t168 + t253 * t204 - (t206 + t205) * t252 - t296;
t140 = -mrSges(5,1) * t190 + mrSges(5,3) * t187 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * qJDD(4) - pkin(4) * t158 + qJD(4) * t238 - t236 * t308 - t322;
t154 = -m(5) * t190 + mrSges(5,1) * t256 - t257 * mrSges(5,2) - t261 * t309 - t262 * t308 - t158;
t194 = t284 * pkin(2) + t323;
t290 = -mrSges(4,1) * t194 - pkin(3) * t154 - pkin(8) * t146 - t278 * t136 - t281 * t140;
t316 = Ifges(4,5) - Ifges(3,6);
t317 = -Ifges(4,4) + Ifges(3,5);
t319 = mrSges(3,1) - mrSges(4,2);
t128 = mrSges(3,3) * t200 - pkin(2) * t144 - t316 * qJDD(2) - t319 * t231 + t317 * t284 + t290;
t145 = t278 * t153 + t281 * t165;
t196 = -qJDD(2) * pkin(2) + t291;
t294 = -m(4) * t196 + t284 * mrSges(4,3) - t145;
t142 = m(3) * t199 - t284 * mrSges(3,2) + t319 * qJDD(2) + t294;
t288 = -m(4) * t194 + t284 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t154;
t151 = m(3) * t200 - mrSges(3,1) * t284 - qJDD(2) * mrSges(3,2) + t288;
t139 = -t142 * t279 + t282 * t151;
t324 = pkin(7) * t139 + t128 * t282;
t314 = t142 * t282;
t143 = m(3) * t231 + t144;
t134 = -t143 * t274 + t151 * t312 + t276 * t314;
t292 = mrSges(4,2) * t196 - mrSges(4,3) * t194 + Ifges(4,1) * qJDD(2) - pkin(8) * t145 + t281 * t136 - t278 * t140;
t126 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t199 - mrSges(3,2) * t200 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t294) + qJ(3) * t288 + t292;
t289 = mrSges(5,1) * t186 - mrSges(5,2) * t187 + Ifges(5,5) * t257 + Ifges(5,6) * t256 + Ifges(5,3) * qJDD(4) + pkin(4) * t287 + pkin(9) * t301 + t280 * t148 + t277 * t156 + t237 * t308 + t238 * t309;
t286 = mrSges(4,1) * t196 + pkin(3) * t145 + t289;
t130 = t316 * t284 + (mrSges(3,2) - mrSges(4,3)) * t231 + t317 * qJDD(2) + t286 - qJ(3) * t144 - mrSges(3,3) * t199;
t293 = mrSges(2,1) * t259 - mrSges(2,2) * t260 + pkin(1) * t134 + t276 * t126 + t130 * t313 + t324 * t274;
t137 = m(2) * t260 + t139;
t133 = t276 * t143 + (t151 * t279 + t314) * t274;
t131 = m(2) * t259 + t134;
t124 = mrSges(2,2) * t270 - mrSges(2,3) * t259 - t279 * t128 + t282 * t130 + (-t133 * t274 - t134 * t276) * pkin(7);
t123 = -mrSges(2,1) * t270 + mrSges(2,3) * t260 - pkin(1) * t133 - t274 * t126 + (t130 * t279 + t324) * t276;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t124 - t273 * t123 - qJ(1) * (t131 * t275 + t137 * t273), t124, t130, t292, t136, t156, -t203 * t264 + t295; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t273 * t124 + t275 * t123 + qJ(1) * (-t131 * t273 + t137 * t275), t123, t128, mrSges(4,3) * t231 + Ifges(4,4) * qJDD(2) - t284 * Ifges(4,5) - t286, t140, t148, -t253 * t201 + t297; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t293, t293, t126, -mrSges(4,2) * t231 + t284 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t290, t289, t322, -t252 * t205 - t296;];
m_new  = t1;
