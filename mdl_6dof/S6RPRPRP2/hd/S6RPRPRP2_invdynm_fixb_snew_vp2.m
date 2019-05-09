% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:32:39
% EndTime: 2019-05-05 17:32:55
% DurationCPUTime: 8.59s
% Computational Cost: add. (149633->337), mult. (320341->414), div. (0->0), fcn. (207588->10), ass. (0->127)
t340 = -2 * qJD(4);
t302 = sin(qJ(1));
t304 = cos(qJ(1));
t284 = t302 * g(1) - g(2) * t304;
t275 = qJDD(1) * pkin(1) + t284;
t285 = -g(1) * t304 - g(2) * t302;
t306 = qJD(1) ^ 2;
t277 = -pkin(1) * t306 + t285;
t297 = sin(pkin(9));
t299 = cos(pkin(9));
t253 = t297 * t275 + t299 * t277;
t239 = -pkin(2) * t306 + qJDD(1) * pkin(7) + t253;
t295 = -g(3) + qJDD(2);
t301 = sin(qJ(3));
t303 = cos(qJ(3));
t228 = -t301 * t239 + t303 * t295;
t329 = qJD(1) * qJD(3);
t327 = t303 * t329;
t278 = qJDD(1) * t301 + t327;
t202 = (-t278 + t327) * qJ(4) + (t301 * t303 * t306 + qJDD(3)) * pkin(3) + t228;
t229 = t303 * t239 + t301 * t295;
t279 = qJDD(1) * t303 - t301 * t329;
t332 = qJD(1) * t301;
t281 = qJD(3) * pkin(3) - qJ(4) * t332;
t294 = t303 ^ 2;
t203 = -pkin(3) * t294 * t306 + qJ(4) * t279 - qJD(3) * t281 + t229;
t296 = sin(pkin(10));
t298 = cos(pkin(10));
t265 = (t296 * t303 + t298 * t301) * qJD(1);
t193 = t298 * t202 - t296 * t203 + t265 * t340;
t264 = (t296 * t301 - t298 * t303) * qJD(1);
t194 = t296 * t202 + t298 * t203 + t264 * t340;
t242 = pkin(4) * t264 - pkin(8) * t265;
t305 = qJD(3) ^ 2;
t191 = -pkin(4) * t305 + qJDD(3) * pkin(8) - t242 * t264 + t194;
t252 = t299 * t275 - t297 * t277;
t317 = -qJDD(1) * pkin(2) - t252;
t204 = -t279 * pkin(3) + qJDD(4) + t281 * t332 + (-qJ(4) * t294 - pkin(7)) * t306 + t317;
t254 = -t278 * t296 + t279 * t298;
t255 = t278 * t298 + t279 * t296;
t196 = (qJD(3) * t264 - t255) * pkin(8) + (qJD(3) * t265 - t254) * pkin(4) + t204;
t300 = sin(qJ(5));
t337 = cos(qJ(5));
t187 = -t300 * t191 + t337 * t196;
t188 = t337 * t191 + t300 * t196;
t256 = -t337 * qJD(3) + t300 * t265;
t257 = t300 * qJD(3) + t337 * t265;
t263 = qJD(5) + t264;
t205 = Ifges(7,5) * t257 + Ifges(7,6) * t263 + Ifges(7,3) * t256;
t208 = Ifges(6,4) * t257 - Ifges(6,2) * t256 + Ifges(6,6) * t263;
t210 = Ifges(6,1) * t257 - Ifges(6,4) * t256 + Ifges(6,5) * t263;
t219 = t257 * qJD(5) - t337 * qJDD(3) + t300 * t255;
t220 = -t256 * qJD(5) + t300 * qJDD(3) + t337 * t255;
t225 = mrSges(7,1) * t256 - mrSges(7,3) * t257;
t251 = qJDD(5) - t254;
t224 = pkin(5) * t256 - qJ(6) * t257;
t262 = t263 ^ 2;
t183 = -pkin(5) * t262 + qJ(6) * t251 + 0.2e1 * qJD(6) * t263 - t224 * t256 + t188;
t185 = -t251 * pkin(5) - t262 * qJ(6) + t257 * t224 + qJDD(6) - t187;
t209 = Ifges(7,1) * t257 + Ifges(7,4) * t263 + Ifges(7,5) * t256;
t316 = mrSges(7,1) * t185 - mrSges(7,3) * t183 - Ifges(7,4) * t220 - Ifges(7,2) * t251 - Ifges(7,6) * t219 - t256 * t209;
t230 = -mrSges(7,2) * t256 + mrSges(7,3) * t263;
t321 = -m(7) * t185 + t251 * mrSges(7,1) + t263 * t230;
t233 = -mrSges(7,1) * t263 + mrSges(7,2) * t257;
t328 = m(7) * t183 + t251 * mrSges(7,3) + t263 * t233;
t339 = -(-t208 + t205) * t257 + mrSges(6,1) * t187 - mrSges(6,2) * t188 + Ifges(6,5) * t220 - Ifges(6,6) * t219 + Ifges(6,3) * t251 + pkin(5) * (-t220 * mrSges(7,2) - t257 * t225 + t321) + qJ(6) * (-t219 * mrSges(7,2) - t256 * t225 + t328) + t256 * t210 - t316;
t241 = mrSges(5,1) * t264 + mrSges(5,2) * t265;
t259 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t265;
t232 = mrSges(6,1) * t263 - mrSges(6,3) * t257;
t333 = -mrSges(6,1) * t256 - mrSges(6,2) * t257 - t225;
t336 = -mrSges(6,3) - mrSges(7,2);
t173 = m(6) * t188 - t251 * mrSges(6,2) + t336 * t219 - t263 * t232 + t333 * t256 + t328;
t231 = -mrSges(6,2) * t263 - mrSges(6,3) * t256;
t175 = m(6) * t187 + t251 * mrSges(6,1) + t336 * t220 + t263 * t231 + t333 * t257 + t321;
t323 = t337 * t173 - t175 * t300;
t161 = m(5) * t194 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t254 - qJD(3) * t259 - t241 * t264 + t323;
t258 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t264;
t190 = -qJDD(3) * pkin(4) - t305 * pkin(8) + t265 * t242 - t193;
t186 = -0.2e1 * qJD(6) * t257 + (t256 * t263 - t220) * qJ(6) + (t257 * t263 + t219) * pkin(5) + t190;
t180 = m(7) * t186 + mrSges(7,1) * t219 - t220 * mrSges(7,3) + t230 * t256 - t257 * t233;
t310 = -m(6) * t190 - t219 * mrSges(6,1) - mrSges(6,2) * t220 - t256 * t231 - t232 * t257 - t180;
t170 = m(5) * t193 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t255 + qJD(3) * t258 - t241 * t265 + t310;
t155 = t296 * t161 + t298 * t170;
t270 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t301 + Ifges(4,2) * t303) * qJD(1);
t271 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t301 + Ifges(4,4) * t303) * qJD(1);
t320 = -mrSges(7,1) * t186 + mrSges(7,2) * t183;
t207 = Ifges(7,4) * t257 + Ifges(7,2) * t263 + Ifges(7,6) * t256;
t335 = -Ifges(6,5) * t257 + Ifges(6,6) * t256 - Ifges(6,3) * t263 - t207;
t163 = -mrSges(6,1) * t190 + mrSges(6,3) * t188 - pkin(5) * t180 + (t209 + t210) * t263 + t335 * t257 + (Ifges(6,6) - Ifges(7,6)) * t251 + (Ifges(6,4) - Ifges(7,5)) * t220 + (-Ifges(6,2) - Ifges(7,3)) * t219 + t320;
t315 = mrSges(7,2) * t185 - mrSges(7,3) * t186 + Ifges(7,1) * t220 + Ifges(7,4) * t251 + Ifges(7,5) * t219 + t263 * t205;
t165 = mrSges(6,2) * t190 - mrSges(6,3) * t187 + Ifges(6,1) * t220 - Ifges(6,4) * t219 + Ifges(6,5) * t251 - qJ(6) * t180 - t263 * t208 + t335 * t256 + t315;
t236 = Ifges(5,4) * t265 - Ifges(5,2) * t264 + Ifges(5,6) * qJD(3);
t237 = Ifges(5,1) * t265 - Ifges(5,4) * t264 + Ifges(5,5) * qJD(3);
t311 = -mrSges(5,1) * t193 + mrSges(5,2) * t194 - Ifges(5,5) * t255 - Ifges(5,6) * t254 - Ifges(5,3) * qJDD(3) - pkin(4) * t310 - pkin(8) * t323 - t337 * t163 - t300 * t165 - t265 * t236 - t264 * t237;
t338 = mrSges(4,1) * t228 - mrSges(4,2) * t229 + Ifges(4,5) * t278 + Ifges(4,6) * t279 + Ifges(4,3) * qJDD(3) + pkin(3) * t155 + (t270 * t301 - t271 * t303) * qJD(1) - t311;
t276 = (-mrSges(4,1) * t303 + mrSges(4,2) * t301) * qJD(1);
t331 = qJD(1) * t303;
t283 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t331;
t153 = m(4) * t228 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t278 + qJD(3) * t283 - t276 * t332 + t155;
t282 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t332;
t324 = t298 * t161 - t170 * t296;
t154 = m(4) * t229 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t279 - qJD(3) * t282 + t276 * t331 + t324;
t325 = -t153 * t301 + t303 * t154;
t145 = m(3) * t253 - mrSges(3,1) * t306 - qJDD(1) * mrSges(3,2) + t325;
t238 = -t306 * pkin(7) + t317;
t167 = t300 * t173 + t337 * t175;
t313 = m(5) * t204 - t254 * mrSges(5,1) + mrSges(5,2) * t255 + t264 * t258 + t259 * t265 + t167;
t309 = -m(4) * t238 + t279 * mrSges(4,1) - mrSges(4,2) * t278 - t282 * t332 + t283 * t331 - t313;
t157 = m(3) * t252 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t306 + t309;
t142 = t297 * t145 + t299 * t157;
t147 = t303 * t153 + t301 * t154;
t326 = t299 * t145 - t157 * t297;
t235 = Ifges(5,5) * t265 - Ifges(5,6) * t264 + Ifges(5,3) * qJD(3);
t148 = mrSges(5,2) * t204 - mrSges(5,3) * t193 + Ifges(5,1) * t255 + Ifges(5,4) * t254 + Ifges(5,5) * qJDD(3) - pkin(8) * t167 - qJD(3) * t236 - t300 * t163 + t337 * t165 - t264 * t235;
t149 = -mrSges(5,1) * t204 + mrSges(5,3) * t194 + Ifges(5,4) * t255 + Ifges(5,2) * t254 + Ifges(5,6) * qJDD(3) - pkin(4) * t167 + qJD(3) * t237 - t265 * t235 - t339;
t269 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t301 + Ifges(4,6) * t303) * qJD(1);
t136 = -mrSges(4,1) * t238 + mrSges(4,3) * t229 + Ifges(4,4) * t278 + Ifges(4,2) * t279 + Ifges(4,6) * qJDD(3) - pkin(3) * t313 + qJ(4) * t324 + qJD(3) * t271 + t296 * t148 + t298 * t149 - t269 * t332;
t138 = mrSges(4,2) * t238 - mrSges(4,3) * t228 + Ifges(4,1) * t278 + Ifges(4,4) * t279 + Ifges(4,5) * qJDD(3) - qJ(4) * t155 - qJD(3) * t270 + t148 * t298 - t149 * t296 + t269 * t331;
t314 = mrSges(3,1) * t252 - mrSges(3,2) * t253 + Ifges(3,3) * qJDD(1) + pkin(2) * t309 + pkin(7) * t325 + t303 * t136 + t301 * t138;
t312 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t314;
t140 = m(2) * t285 - mrSges(2,1) * t306 - qJDD(1) * mrSges(2,2) + t326;
t139 = m(2) * t284 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t306 + t142;
t134 = -mrSges(3,1) * t295 + mrSges(3,3) * t253 + t306 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t147 - t338;
t133 = mrSges(3,2) * t295 - mrSges(3,3) * t252 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t306 - pkin(7) * t147 - t136 * t301 + t138 * t303;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t284 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t306 - qJ(2) * t142 + t133 * t299 - t134 * t297;
t131 = Ifges(2,6) * qJDD(1) + t306 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t285 + t297 * t133 + t299 * t134 - pkin(1) * (m(3) * t295 + t147) + qJ(2) * t326;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t132 - t302 * t131 - pkin(6) * (t139 * t304 + t140 * t302), t132, t133, t138, t148, t165, -t207 * t256 + t315; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t302 * t132 + t304 * t131 + pkin(6) * (-t139 * t302 + t140 * t304), t131, t134, t136, t149, t163, -t257 * t205 - t316; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t312, t312, t314, t338, -t311, t339, Ifges(7,5) * t220 + Ifges(7,6) * t251 + Ifges(7,3) * t219 + t257 * t207 - t263 * t209 - t320;];
m_new  = t1;
