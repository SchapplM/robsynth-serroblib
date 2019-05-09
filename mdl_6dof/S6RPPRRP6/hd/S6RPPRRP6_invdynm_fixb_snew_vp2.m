% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:00:12
% EndTime: 2019-05-05 15:00:18
% DurationCPUTime: 2.58s
% Computational Cost: add. (31637->295), mult. (58051->333), div. (0->0), fcn. (29439->6), ass. (0->105)
t271 = sin(qJ(1));
t273 = cos(qJ(1));
t242 = t271 * g(1) - t273 * g(2);
t276 = qJD(1) ^ 2;
t217 = -qJDD(1) * pkin(1) - t276 * qJ(2) + qJDD(2) - t242;
t207 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t217;
t243 = -t273 * g(1) - t271 * g(2);
t317 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t243;
t202 = -t276 * pkin(7) - t207;
t270 = sin(qJ(4));
t272 = cos(qJ(4));
t304 = qJD(1) * qJD(4);
t299 = t272 * t304;
t237 = -qJDD(1) * t270 - t299;
t300 = t270 * t304;
t238 = qJDD(1) * t272 - t300;
t169 = (-t238 + t300) * pkin(8) + (-t237 + t299) * pkin(4) + t202;
t208 = qJDD(3) + (-pkin(1) - qJ(3)) * t276 + t317;
t203 = -qJDD(1) * pkin(7) + t208;
t194 = -g(3) * t272 + t270 * t203;
t236 = (pkin(4) * t270 - pkin(8) * t272) * qJD(1);
t275 = qJD(4) ^ 2;
t306 = qJD(1) * t270;
t172 = -pkin(4) * t275 + qJDD(4) * pkin(8) - t236 * t306 + t194;
t269 = sin(qJ(5));
t315 = cos(qJ(5));
t166 = t315 * t169 - t269 * t172;
t167 = t269 * t169 + t315 * t172;
t305 = qJD(1) * t272;
t233 = -t315 * qJD(4) + t269 * t305;
t234 = t269 * qJD(4) + t315 * t305;
t245 = qJD(5) + t306;
t175 = Ifges(7,5) * t234 + Ifges(7,6) * t245 + Ifges(7,3) * t233;
t178 = Ifges(6,4) * t234 - Ifges(6,2) * t233 + Ifges(6,6) * t245;
t180 = Ifges(6,1) * t234 - Ifges(6,4) * t233 + Ifges(6,5) * t245;
t191 = qJD(5) * t234 - t315 * qJDD(4) + t238 * t269;
t192 = -t233 * qJD(5) + t269 * qJDD(4) + t315 * t238;
t199 = mrSges(7,1) * t233 - mrSges(7,3) * t234;
t232 = qJDD(5) - t237;
t198 = pkin(5) * t233 - qJ(6) * t234;
t244 = t245 ^ 2;
t162 = -pkin(5) * t244 + qJ(6) * t232 + 0.2e1 * qJD(6) * t245 - t198 * t233 + t167;
t164 = -t232 * pkin(5) - t244 * qJ(6) + t234 * t198 + qJDD(6) - t166;
t179 = Ifges(7,1) * t234 + Ifges(7,4) * t245 + Ifges(7,5) * t233;
t290 = mrSges(7,1) * t164 - mrSges(7,3) * t162 - Ifges(7,4) * t192 - Ifges(7,2) * t232 - Ifges(7,6) * t191 - t233 * t179;
t212 = -mrSges(7,2) * t233 + mrSges(7,3) * t245;
t295 = -m(7) * t164 + t232 * mrSges(7,1) + t245 * t212;
t211 = -mrSges(7,1) * t245 + mrSges(7,2) * t234;
t301 = m(7) * t162 + t232 * mrSges(7,3) + t245 * t211;
t316 = -(-t178 + t175) * t234 + mrSges(6,1) * t166 - mrSges(6,2) * t167 + Ifges(6,5) * t192 - Ifges(6,6) * t191 + Ifges(6,3) * t232 + pkin(5) * (-t192 * mrSges(7,2) - t234 * t199 + t295) + qJ(6) * (-t191 * mrSges(7,2) - t233 * t199 + t301) + t233 * t180 - t290;
t314 = mrSges(3,2) - mrSges(4,3);
t313 = -mrSges(6,3) - mrSges(7,2);
t312 = Ifges(3,4) - Ifges(4,5);
t311 = Ifges(2,6) - Ifges(3,5);
t310 = mrSges(4,3) * t276;
t235 = (mrSges(5,1) * t270 + mrSges(5,2) * t272) * qJD(1);
t241 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t305;
t210 = mrSges(6,1) * t245 - mrSges(6,3) * t234;
t307 = -mrSges(6,1) * t233 - mrSges(6,2) * t234 - t199;
t152 = m(6) * t167 - t232 * mrSges(6,2) + t313 * t191 - t245 * t210 + t307 * t233 + t301;
t209 = -mrSges(6,2) * t245 - mrSges(6,3) * t233;
t154 = m(6) * t166 + t232 * mrSges(6,1) + t313 * t192 + t245 * t209 + t307 * t234 + t295;
t297 = t315 * t152 - t154 * t269;
t144 = m(5) * t194 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t237 - qJD(4) * t241 - t235 * t306 + t297;
t193 = t270 * g(3) + t272 * t203;
t240 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t306;
t171 = -qJDD(4) * pkin(4) - t275 * pkin(8) + t236 * t305 - t193;
t165 = -0.2e1 * qJD(6) * t234 + (t233 * t245 - t192) * qJ(6) + (t234 * t245 + t191) * pkin(5) + t171;
t157 = m(7) * t165 + mrSges(7,1) * t191 - t192 * mrSges(7,3) - t234 * t211 + t212 * t233;
t280 = -m(6) * t171 - t191 * mrSges(6,1) - mrSges(6,2) * t192 - t233 * t209 - t210 * t234 - t157;
t149 = m(5) * t193 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t238 + qJD(4) * t240 - t235 * t305 + t280;
t131 = t270 * t144 + t272 * t149;
t146 = t269 * t152 + t315 * t154;
t177 = Ifges(7,4) * t234 + Ifges(7,2) * t245 + Ifges(7,6) * t233;
t309 = -Ifges(6,5) * t234 + Ifges(6,6) * t233 - Ifges(6,3) * t245 - t177;
t298 = t272 * t144 - t270 * t149;
t296 = m(4) * t208 + qJDD(1) * mrSges(4,2) + t131;
t294 = -mrSges(7,1) * t165 + mrSges(7,2) * t162;
t291 = -m(5) * t202 + t237 * mrSges(5,1) - t238 * mrSges(5,2) - t240 * t306 - t241 * t305 - t146;
t289 = mrSges(7,2) * t164 - mrSges(7,3) * t165 + Ifges(7,1) * t192 + Ifges(7,4) * t232 + Ifges(7,5) * t191 + t245 * t175;
t137 = -mrSges(6,1) * t171 + mrSges(6,3) * t167 - pkin(5) * t157 + (t179 + t180) * t245 + t309 * t234 + (Ifges(6,6) - Ifges(7,6)) * t232 + (Ifges(6,4) - Ifges(7,5)) * t192 + (-Ifges(6,2) - Ifges(7,3)) * t191 + t294;
t139 = mrSges(6,2) * t171 - mrSges(6,3) * t166 + Ifges(6,1) * t192 - Ifges(6,4) * t191 + Ifges(6,5) * t232 - qJ(6) * t157 - t245 * t178 + t309 * t233 + t289;
t219 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t272 - Ifges(5,6) * t270) * qJD(1);
t220 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t272 - Ifges(5,2) * t270) * qJD(1);
t122 = mrSges(5,2) * t202 - mrSges(5,3) * t193 + Ifges(5,1) * t238 + Ifges(5,4) * t237 + Ifges(5,5) * qJDD(4) - pkin(8) * t146 - qJD(4) * t220 - t269 * t137 + t315 * t139 - t219 * t306;
t221 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t272 - Ifges(5,4) * t270) * qJD(1);
t124 = -mrSges(5,1) * t202 + mrSges(5,3) * t194 + Ifges(5,4) * t238 + Ifges(5,2) * t237 + Ifges(5,6) * qJDD(4) - pkin(4) * t146 + qJD(4) * t221 - t219 * t305 - t316;
t288 = mrSges(4,1) * t207 + mrSges(4,2) * g(3) + t276 * Ifges(4,4) + Ifges(4,5) * qJDD(1) + pkin(3) * t291 + pkin(7) * t298 + t270 * t122 + t272 * t124;
t215 = pkin(1) * t276 - t317;
t287 = -m(3) * t215 + t276 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t296;
t286 = mrSges(4,2) * t208 - mrSges(4,3) * t207 + Ifges(4,1) * qJDD(1) - pkin(7) * t131 + t272 * t122 - t124 * t270;
t285 = mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t238 + Ifges(5,6) * t237 + Ifges(5,3) * qJDD(4) + pkin(4) * t280 + pkin(8) * t297 + t315 * t137 + t269 * t139 + t220 * t305 + t221 * t306;
t136 = m(4) * t207 - t276 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t291;
t284 = mrSges(3,1) * t217 + pkin(2) * t136 + t288;
t283 = -m(3) * t217 + t276 * mrSges(3,3) - t136;
t282 = mrSges(4,1) * t208 - Ifges(4,4) * qJDD(1) + pkin(3) * t131 + t285;
t281 = mrSges(3,2) * t217 - mrSges(3,3) * t215 + Ifges(3,1) * qJDD(1) - qJ(3) * t136 + t286;
t279 = -mrSges(2,2) * t243 + qJ(2) * (t287 - t310) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t283) + mrSges(2,1) * t242 + Ifges(2,3) * qJDD(1) + t281;
t278 = mrSges(3,1) * t215 + pkin(2) * (-t296 + t310) + qJ(3) * (-m(4) * g(3) + t298) - t282;
t133 = t283 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) - t276 * mrSges(2,2) + m(2) * t242;
t128 = (-m(3) - m(4)) * g(3) + t298;
t125 = m(2) * t243 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t276 + t287;
t119 = -t278 + (Ifges(2,5) - t312) * t276 + t311 * qJDD(1) + (mrSges(2,1) - t314) * g(3) + mrSges(2,3) * t243 - pkin(1) * t128;
t118 = t284 + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t242 - qJ(2) * t128 - t311 * t276;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t273 * t118 - t271 * t119 - pkin(6) * (t125 * t271 + t133 * t273), t118, t281, t286, t122, t139, -t177 * t233 + t289; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t118 + t273 * t119 + pkin(6) * (t125 * t273 - t133 * t271), t119, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t276 * Ifges(3,5) - t284, -mrSges(4,3) * g(3) - t276 * Ifges(4,5) - t282, t124, t137, -t234 * t175 - t290; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t279, t279, Ifges(3,5) * qJDD(1) + t314 * g(3) + t312 * t276 + t278, t288, t285, t316, Ifges(7,5) * t192 + Ifges(7,6) * t232 + Ifges(7,3) * t191 + t234 * t177 - t245 * t179 - t294;];
m_new  = t1;
