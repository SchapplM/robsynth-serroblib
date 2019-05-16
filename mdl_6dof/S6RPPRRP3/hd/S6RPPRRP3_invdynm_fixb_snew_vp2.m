% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:50:46
% EndTime: 2019-05-05 14:50:53
% DurationCPUTime: 3.65s
% Computational Cost: add. (50305->293), mult. (91308->343), div. (0->0), fcn. (49702->8), ass. (0->113)
t264 = sin(qJ(1));
t266 = cos(qJ(1));
t240 = t264 * g(1) - t266 * g(2);
t231 = qJDD(1) * pkin(1) + t240;
t241 = -t266 * g(1) - t264 * g(2);
t268 = qJD(1) ^ 2;
t233 = -t268 * pkin(1) + t241;
t260 = sin(pkin(9));
t261 = cos(pkin(9));
t200 = t260 * t231 + t261 * t233;
t304 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t200;
t302 = -pkin(2) - pkin(7);
t177 = t302 * t268 - t304;
t263 = sin(qJ(4));
t265 = cos(qJ(4));
t291 = qJD(1) * qJD(4);
t286 = t265 * t291;
t235 = -t263 * qJDD(1) - t286;
t287 = t263 * t291;
t236 = t265 * qJDD(1) - t287;
t166 = (-t236 + t287) * pkin(8) + (-t235 + t286) * pkin(4) + t177;
t199 = t261 * t231 - t260 * t233;
t281 = -t268 * qJ(3) + qJDD(3) - t199;
t178 = t302 * qJDD(1) + t281;
t257 = -g(3) + qJDD(2);
t172 = t263 * t178 + t265 * t257;
t234 = (pkin(4) * t263 - pkin(8) * t265) * qJD(1);
t267 = qJD(4) ^ 2;
t292 = t263 * qJD(1);
t169 = -t267 * pkin(4) + qJDD(4) * pkin(8) - t234 * t292 + t172;
t262 = sin(qJ(5));
t301 = cos(qJ(5));
t163 = t301 * t166 - t262 * t169;
t164 = t262 * t166 + t301 * t169;
t293 = qJD(1) * t265;
t229 = -t301 * qJD(4) + t262 * t293;
t230 = t262 * qJD(4) + t301 * t293;
t243 = qJD(5) + t292;
t182 = Ifges(7,5) * t230 + Ifges(7,6) * t243 + Ifges(7,3) * t229;
t185 = Ifges(6,4) * t230 - Ifges(6,2) * t229 + Ifges(6,6) * t243;
t187 = Ifges(6,1) * t230 - Ifges(6,4) * t229 + Ifges(6,5) * t243;
t196 = t230 * qJD(5) - t301 * qJDD(4) + t262 * t236;
t197 = -t229 * qJD(5) + t262 * qJDD(4) + t301 * t236;
t204 = t229 * mrSges(7,1) - t230 * mrSges(7,3);
t228 = qJDD(5) - t235;
t203 = t229 * pkin(5) - t230 * qJ(6);
t242 = t243 ^ 2;
t159 = -t242 * pkin(5) + t228 * qJ(6) + 0.2e1 * qJD(6) * t243 - t229 * t203 + t164;
t161 = -t228 * pkin(5) - t242 * qJ(6) + t230 * t203 + qJDD(6) - t163;
t186 = Ifges(7,1) * t230 + Ifges(7,4) * t243 + Ifges(7,5) * t229;
t280 = mrSges(7,1) * t161 - mrSges(7,3) * t159 - Ifges(7,4) * t197 - Ifges(7,2) * t228 - Ifges(7,6) * t196 - t229 * t186;
t209 = -t229 * mrSges(7,2) + t243 * mrSges(7,3);
t283 = -m(7) * t161 + t228 * mrSges(7,1) + t243 * t209;
t208 = -t243 * mrSges(7,1) + t230 * mrSges(7,2);
t289 = m(7) * t159 + t228 * mrSges(7,3) + t243 * t208;
t303 = -(-t185 + t182) * t230 + mrSges(6,1) * t163 - mrSges(6,2) * t164 + Ifges(6,5) * t197 - Ifges(6,6) * t196 + Ifges(6,3) * t228 + pkin(5) * (-t197 * mrSges(7,2) - t230 * t204 + t283) + qJ(6) * (-t196 * mrSges(7,2) - t229 * t204 + t289) + t229 * t187 - t280;
t300 = mrSges(3,1) - mrSges(4,2);
t299 = -mrSges(6,3) - mrSges(7,2);
t298 = -Ifges(4,4) + Ifges(3,5);
t297 = Ifges(4,5) - Ifges(3,6);
t232 = (mrSges(5,1) * t263 + mrSges(5,2) * t265) * qJD(1);
t239 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t293;
t207 = t243 * mrSges(6,1) - t230 * mrSges(6,3);
t294 = -t229 * mrSges(6,1) - t230 * mrSges(6,2) - t204;
t149 = m(6) * t164 - t228 * mrSges(6,2) + t299 * t196 - t243 * t207 + t294 * t229 + t289;
t206 = -t243 * mrSges(6,2) - t229 * mrSges(6,3);
t151 = m(6) * t163 + t228 * mrSges(6,1) + t299 * t197 + t243 * t206 + t294 * t230 + t283;
t284 = t301 * t149 - t262 * t151;
t139 = m(5) * t172 - qJDD(4) * mrSges(5,2) + t235 * mrSges(5,3) - qJD(4) * t239 - t232 * t292 + t284;
t171 = t265 * t178 - t263 * t257;
t238 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t292;
t168 = -qJDD(4) * pkin(4) - t267 * pkin(8) + t234 * t293 - t171;
t162 = -0.2e1 * qJD(6) * t230 + (t229 * t243 - t197) * qJ(6) + (t230 * t243 + t196) * pkin(5) + t168;
t156 = m(7) * t162 + t196 * mrSges(7,1) - t197 * mrSges(7,3) - t230 * t208 + t229 * t209;
t271 = -m(6) * t168 - t196 * mrSges(6,1) - t197 * mrSges(6,2) - t229 * t206 - t230 * t207 - t156;
t146 = m(5) * t171 + qJDD(4) * mrSges(5,1) - t236 * mrSges(5,3) + qJD(4) * t238 - t232 * t293 + t271;
t130 = t263 * t139 + t265 * t146;
t181 = -qJDD(1) * pkin(2) + t281;
t278 = -m(4) * t181 + t268 * mrSges(4,3) - t130;
t127 = m(3) * t199 - t268 * mrSges(3,2) + t300 * qJDD(1) + t278;
t144 = t262 * t149 + t301 * t151;
t137 = -m(5) * t177 + t235 * mrSges(5,1) - t236 * mrSges(5,2) - t238 * t292 - t239 * t293 - t144;
t179 = t268 * pkin(2) + t304;
t274 = -m(4) * t179 + t268 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t137;
t134 = m(3) * t200 - t268 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t274;
t124 = t261 * t127 + t260 * t134;
t184 = Ifges(7,4) * t230 + Ifges(7,2) * t243 + Ifges(7,6) * t229;
t296 = -Ifges(6,5) * t230 + Ifges(6,6) * t229 - Ifges(6,3) * t243 - t184;
t285 = -t260 * t127 + t261 * t134;
t131 = t265 * t139 - t263 * t146;
t129 = m(4) * t257 + t131;
t282 = -mrSges(7,1) * t162 + mrSges(7,2) * t159;
t279 = mrSges(7,2) * t161 - mrSges(7,3) * t162 + Ifges(7,1) * t197 + Ifges(7,4) * t228 + Ifges(7,5) * t196 + t243 * t182;
t140 = -mrSges(6,1) * t168 + mrSges(6,3) * t164 - pkin(5) * t156 + (t186 + t187) * t243 + t296 * t230 + (Ifges(6,6) - Ifges(7,6)) * t228 + (Ifges(6,4) - Ifges(7,5)) * t197 + (-Ifges(6,2) - Ifges(7,3)) * t196 + t282;
t142 = mrSges(6,2) * t168 - mrSges(6,3) * t163 + Ifges(6,1) * t197 - Ifges(6,4) * t196 + Ifges(6,5) * t228 - qJ(6) * t156 - t243 * t185 + t296 * t229 + t279;
t216 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t265 - Ifges(5,6) * t263) * qJD(1);
t217 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t265 - Ifges(5,2) * t263) * qJD(1);
t120 = mrSges(5,2) * t177 - mrSges(5,3) * t171 + Ifges(5,1) * t236 + Ifges(5,4) * t235 + Ifges(5,5) * qJDD(4) - pkin(8) * t144 - qJD(4) * t217 - t262 * t140 + t301 * t142 - t216 * t292;
t218 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t265 - Ifges(5,4) * t263) * qJD(1);
t125 = -mrSges(5,1) * t177 + mrSges(5,3) * t172 + Ifges(5,4) * t236 + Ifges(5,2) * t235 + Ifges(5,6) * qJDD(4) - pkin(4) * t144 + qJD(4) * t218 - t216 * t293 - t303;
t277 = mrSges(4,2) * t181 - mrSges(4,3) * t179 + Ifges(4,1) * qJDD(1) - pkin(7) * t130 + t265 * t120 - t263 * t125;
t276 = -mrSges(4,1) * t179 - pkin(3) * t137 - pkin(7) * t131 - t263 * t120 - t265 * t125;
t275 = mrSges(5,1) * t171 - mrSges(5,2) * t172 + Ifges(5,5) * t236 + Ifges(5,6) * t235 + Ifges(5,3) * qJDD(4) + pkin(4) * t271 + pkin(8) * t284 + t301 * t140 + t262 * t142 + t217 * t293 + t218 * t292;
t273 = -mrSges(3,2) * t200 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t278) + qJ(3) * t274 + mrSges(3,1) * t199 + Ifges(3,3) * qJDD(1) + t277;
t272 = mrSges(4,1) * t181 + pkin(3) * t130 + t275;
t270 = mrSges(2,1) * t240 - mrSges(2,2) * t241 + Ifges(2,3) * qJDD(1) + pkin(1) * t124 + t273;
t122 = m(2) * t241 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t121 = m(2) * t240 + qJDD(1) * mrSges(2,1) - t268 * mrSges(2,2) + t124;
t118 = -mrSges(3,3) * t199 - qJ(3) * t129 + t272 + t297 * t268 + t298 * qJDD(1) + (mrSges(3,2) - mrSges(4,3)) * t257;
t117 = mrSges(3,3) * t200 - pkin(2) * t129 - t297 * qJDD(1) - t300 * t257 + t298 * t268 + t276;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t240 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - qJ(2) * t124 - t260 * t117 + t261 * t118;
t115 = Ifges(2,6) * qJDD(1) + t268 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t241 + t260 * t118 + t261 * t117 - pkin(1) * (m(3) * t257 + t129) + qJ(2) * t285;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t266 * t116 - t264 * t115 - pkin(6) * (t266 * t121 + t264 * t122), t116, t118, t277, t120, t142, -t229 * t184 + t279; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t116 + t266 * t115 + pkin(6) * (-t264 * t121 + t266 * t122), t115, t117, mrSges(4,3) * t257 + Ifges(4,4) * qJDD(1) - t268 * Ifges(4,5) - t272, t125, t140, -t230 * t182 - t280; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t270, t270, t273, -mrSges(4,2) * t257 + t268 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t276, t275, t303, Ifges(7,5) * t197 + Ifges(7,6) * t228 + Ifges(7,3) * t196 + t230 * t184 - t243 * t186 - t282;];
m_new  = t1;
