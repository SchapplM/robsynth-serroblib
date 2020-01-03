% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:28
% EndTime: 2019-12-31 20:09:35
% DurationCPUTime: 2.80s
% Computational Cost: add. (24173->312), mult. (49222->362), div. (0->0), fcn. (25797->6), ass. (0->108)
t267 = sin(qJ(4));
t270 = cos(qJ(4));
t271 = cos(qJ(2));
t302 = qJD(1) * t271;
t238 = -t267 * qJD(2) - t270 * t302;
t268 = sin(qJ(2));
t300 = qJD(1) * qJD(2);
t295 = t268 * t300;
t244 = t271 * qJDD(1) - t295;
t199 = t238 * qJD(4) + t270 * qJDD(2) - t267 * t244;
t239 = t270 * qJD(2) - t267 * t302;
t201 = -t238 * mrSges(6,1) + t239 * mrSges(6,2);
t301 = t268 * qJD(1);
t252 = pkin(3) * t301 - qJD(2) * pkin(7);
t266 = t271 ^ 2;
t274 = qJD(1) ^ 2;
t294 = t271 * t300;
t243 = t268 * qJDD(1) + t294;
t269 = sin(qJ(1));
t272 = cos(qJ(1));
t253 = t269 * g(1) - t272 * g(2);
t291 = -qJDD(1) * pkin(1) - t253;
t315 = -2 * qJD(3);
t280 = pkin(2) * t295 + t301 * t315 + (-t243 - t294) * qJ(3) + t291;
t164 = -t252 * t301 + (-pkin(3) * t266 - pkin(6)) * t274 + (-pkin(2) - pkin(7)) * t244 + t280;
t254 = -t272 * g(1) - t269 * g(2);
t223 = -t274 * pkin(1) + qJDD(1) * pkin(6) + t254;
t208 = -t271 * g(3) - t268 * t223;
t240 = (-pkin(2) * t271 - qJ(3) * t268) * qJD(1);
t273 = qJD(2) ^ 2;
t176 = -qJDD(2) * pkin(2) - t273 * qJ(3) + t240 * t301 + qJDD(3) - t208;
t169 = (-t268 * t271 * t274 - qJDD(2)) * pkin(7) + (t243 - t294) * pkin(3) + t176;
t158 = -t267 * t164 + t270 * t169;
t237 = qJDD(4) + t243;
t256 = qJD(4) + t301;
t153 = -0.2e1 * qJD(5) * t239 + (t238 * t256 - t199) * qJ(5) + (t238 * t239 + t237) * pkin(4) + t158;
t203 = -t256 * mrSges(6,2) + t238 * mrSges(6,3);
t298 = m(6) * t153 + t237 * mrSges(6,1) + t256 * t203;
t150 = -t199 * mrSges(6,3) - t239 * t201 + t298;
t159 = t270 * t164 + t267 * t169;
t182 = Ifges(5,4) * t239 + Ifges(5,2) * t238 + Ifges(5,6) * t256;
t183 = Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t256;
t184 = Ifges(5,1) * t239 + Ifges(5,4) * t238 + Ifges(5,5) * t256;
t198 = -t239 * qJD(4) - t267 * qJDD(2) - t270 * t244;
t205 = t256 * pkin(4) - t239 * qJ(5);
t236 = t238 ^ 2;
t156 = -t236 * pkin(4) + t198 * qJ(5) + 0.2e1 * qJD(5) * t238 - t256 * t205 + t159;
t181 = Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t256;
t287 = -mrSges(6,1) * t153 + mrSges(6,2) * t156 - Ifges(6,5) * t199 - Ifges(6,6) * t198 - Ifges(6,3) * t237 - t239 * t181;
t316 = mrSges(5,1) * t158 - mrSges(5,2) * t159 + Ifges(5,5) * t199 + Ifges(5,6) * t198 + Ifges(5,3) * t237 + pkin(4) * t150 + t239 * t182 - t287 - (t184 + t183) * t238;
t202 = -t238 * mrSges(5,1) + t239 * mrSges(5,2);
t204 = -t256 * mrSges(5,2) + t238 * mrSges(5,3);
t142 = m(5) * t158 + t237 * mrSges(5,1) + t256 * t204 + (-t201 - t202) * t239 + (-mrSges(5,3) - mrSges(6,3)) * t199 + t298;
t206 = t256 * mrSges(6,1) - t239 * mrSges(6,3);
t207 = t256 * mrSges(5,1) - t239 * mrSges(5,3);
t297 = m(6) * t156 + t198 * mrSges(6,3) + t238 * t201;
t146 = m(5) * t159 + t198 * mrSges(5,3) + t238 * t202 + (-t206 - t207) * t256 + (-mrSges(5,2) - mrSges(6,2)) * t237 + t297;
t139 = t270 * t142 + t267 * t146;
t310 = t274 * pkin(6);
t170 = -t244 * pkin(2) + t280 - t310;
t314 = mrSges(4,1) * t176 - mrSges(4,3) * t170 + pkin(3) * t139 + t316;
t209 = -t268 * g(3) + t271 * t223;
t218 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t268 + Ifges(3,4) * t271) * qJD(1);
t241 = (mrSges(4,2) * t271 - mrSges(4,3) * t268) * qJD(1);
t250 = -mrSges(4,1) * t302 - qJD(2) * mrSges(4,3);
t174 = t273 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t315 - t240 * t302 - t209;
t251 = mrSges(4,1) * t301 + qJD(2) * mrSges(4,2);
t168 = -t266 * t274 * pkin(7) + t244 * pkin(3) + qJD(2) * t252 - t174;
t162 = -t198 * pkin(4) - t236 * qJ(5) + t239 * t205 + qJDD(5) + t168;
t296 = m(6) * t162 + t199 * mrSges(6,2) + t239 * t206;
t311 = -m(5) * t168 - t199 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t198 - t239 * t207 + (t203 + t204) * t238 - t296;
t278 = -m(4) * t174 + qJDD(2) * mrSges(4,3) + qJD(2) * t251 + t241 * t302 - t311;
t179 = Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t256;
t180 = Ifges(5,5) * t239 + Ifges(5,6) * t238 + Ifges(5,3) * t256;
t288 = -mrSges(6,1) * t162 + mrSges(6,3) * t156 + Ifges(6,4) * t199 + Ifges(6,2) * t198 + Ifges(6,6) * t237 + t256 * t183;
t131 = Ifges(5,4) * t199 + Ifges(5,2) * t198 + Ifges(5,6) * t237 + t256 * t184 - mrSges(5,1) * t168 + mrSges(5,3) * t159 - pkin(4) * (-t198 * mrSges(6,1) - t238 * t203 + t296) + qJ(5) * (-t237 * mrSges(6,2) - t256 * t206 + t297) + (-t180 - t179) * t239 + t288;
t286 = mrSges(6,2) * t162 - mrSges(6,3) * t153 + Ifges(6,1) * t199 + Ifges(6,4) * t198 + Ifges(6,5) * t237 + t238 * t179;
t138 = mrSges(5,2) * t168 - mrSges(5,3) * t158 + Ifges(5,1) * t199 + Ifges(5,4) * t198 + Ifges(5,5) * t237 - qJ(5) * t150 + t238 * t180 + (-t181 - t182) * t256 + t286;
t220 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t268 - Ifges(4,6) * t271) * qJD(1);
t282 = -mrSges(4,2) * t176 + mrSges(4,3) * t174 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t243 + Ifges(4,5) * t244 + pkin(7) * t139 + t267 * t131 - t270 * t138 - t220 * t302;
t284 = -m(4) * t176 - t243 * mrSges(4,1) - t139;
t219 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t268 - Ifges(4,3) * t271) * qJD(1);
t303 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t268 + Ifges(3,2) * t271) * qJD(1) - t219;
t313 = (-t271 * t218 + t303 * t268) * qJD(1) + mrSges(3,1) * t208 - mrSges(3,2) * t209 + Ifges(3,5) * t243 + Ifges(3,6) * t244 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t250 - t241 * t301 + t284) + qJ(3) * (t244 * mrSges(4,1) + t278) - t282;
t308 = Ifges(3,4) + Ifges(4,6);
t140 = -t267 * t142 + t270 * t146;
t221 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t268 - Ifges(4,5) * t271) * qJD(1);
t304 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t268 + Ifges(3,6) * t271) * qJD(1) + t221;
t242 = (-mrSges(3,1) * t271 + mrSges(3,2) * t268) * qJD(1);
t249 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t302;
t134 = m(3) * t208 - t243 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t249 - t250) * qJD(2) + (-t241 - t242) * t301 + t284;
t248 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t301;
t145 = t278 - qJDD(2) * mrSges(3,2) + t242 * t302 - qJD(2) * t248 + m(3) * t209 + (mrSges(3,3) + mrSges(4,1)) * t244;
t292 = -t268 * t134 + t271 * t145;
t289 = -m(4) * t170 - t244 * mrSges(4,2) + t251 * t301 - t140;
t135 = -t243 * mrSges(4,3) + t250 * t302 - t289;
t222 = t291 - t310;
t281 = -mrSges(4,1) * t174 + mrSges(4,2) * t170 - pkin(3) * t311 - pkin(7) * t140 - t270 * t131 - t267 * t138;
t125 = -mrSges(3,1) * t222 + mrSges(3,3) * t209 - pkin(2) * t135 + (Ifges(3,2) + Ifges(4,3)) * t244 + t308 * t243 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t218 - t220) * qJD(2) - t304 * t301 + t281;
t127 = t304 * t302 + t308 * t244 + (Ifges(3,1) + Ifges(4,2)) * t243 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t303 * qJD(2) + mrSges(3,2) * t222 - mrSges(3,3) * t208 - qJ(3) * t135 + t314;
t277 = -m(3) * t222 + t249 * t302 + t244 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t243 + (-t248 * t268 - t250 * t271) * qJD(1) + t289;
t283 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + Ifges(2,3) * qJDD(1) + pkin(1) * t277 + pkin(6) * t292 + t271 * t125 + t268 * t127;
t132 = m(2) * t253 + qJDD(1) * mrSges(2,1) - t274 * mrSges(2,2) + t277;
t130 = t271 * t134 + t268 * t145;
t128 = m(2) * t254 - t274 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t292;
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t254 + t274 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t313;
t122 = -mrSges(2,2) * g(3) - mrSges(2,3) * t253 + Ifges(2,5) * qJDD(1) - t274 * Ifges(2,6) - pkin(6) * t130 - t268 * t125 + t271 * t127;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t272 * t122 - t269 * t123 - pkin(5) * (t269 * t128 + t272 * t132), t122, t127, -t219 * t301 - t282, t138, -t256 * t181 + t286; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t269 * t122 + t272 * t123 + pkin(5) * (t272 * t128 - t269 * t132), t123, t125, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t243 - Ifges(4,6) * t244 - qJD(2) * t219 - t221 * t302 - t314, t131, -t239 * t179 + t288; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, t313, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t243 - Ifges(4,3) * t244 + qJD(2) * t220 + t221 * t301 - t281, t316, -t238 * t183 - t287;];
m_new = t1;
