% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:56
% EndTime: 2019-12-31 22:28:18
% DurationCPUTime: 11.42s
% Computational Cost: add. (202103->310), mult. (410635->390), div. (0->0), fcn. (288952->10), ass. (0->125)
t261 = sin(qJ(1));
t266 = cos(qJ(1));
t249 = t261 * g(1) - t266 * g(2);
t268 = qJD(1) ^ 2;
t232 = -qJDD(1) * pkin(1) - t268 * pkin(6) - t249;
t260 = sin(qJ(2));
t265 = cos(qJ(2));
t282 = qJD(1) * qJD(2);
t281 = t265 * t282;
t243 = t260 * qJDD(1) + t281;
t253 = t260 * t282;
t244 = t265 * qJDD(1) - t253;
t199 = (-t243 - t281) * pkin(7) + (-t244 + t253) * pkin(2) + t232;
t250 = -t266 * g(1) - t261 * g(2);
t233 = -t268 * pkin(1) + qJDD(1) * pkin(6) + t250;
t221 = -t260 * g(3) + t265 * t233;
t242 = (-pkin(2) * t265 - pkin(7) * t260) * qJD(1);
t267 = qJD(2) ^ 2;
t283 = t265 * qJD(1);
t202 = -t267 * pkin(2) + qJDD(2) * pkin(7) + t242 * t283 + t221;
t259 = sin(qJ(3));
t264 = cos(qJ(3));
t186 = t264 * t199 - t259 * t202;
t284 = qJD(1) * t260;
t239 = t264 * qJD(2) - t259 * t284;
t213 = t239 * qJD(3) + t259 * qJDD(2) + t264 * t243;
t238 = qJDD(3) - t244;
t240 = t259 * qJD(2) + t264 * t284;
t252 = qJD(3) - t283;
t170 = (t239 * t252 - t213) * pkin(8) + (t239 * t240 + t238) * pkin(3) + t186;
t187 = t259 * t199 + t264 * t202;
t212 = -t240 * qJD(3) + t264 * qJDD(2) - t259 * t243;
t222 = t252 * pkin(3) - t240 * pkin(8);
t237 = t239 ^ 2;
t172 = -t237 * pkin(3) + t212 * pkin(8) - t252 * t222 + t187;
t258 = sin(qJ(4));
t263 = cos(qJ(4));
t157 = t263 * t170 - t258 * t172;
t215 = t263 * t239 - t258 * t240;
t185 = t215 * qJD(4) + t258 * t212 + t263 * t213;
t216 = t258 * t239 + t263 * t240;
t234 = qJDD(4) + t238;
t251 = qJD(4) + t252;
t154 = (t215 * t251 - t185) * pkin(9) + (t215 * t216 + t234) * pkin(4) + t157;
t158 = t258 * t170 + t263 * t172;
t184 = -t216 * qJD(4) + t263 * t212 - t258 * t213;
t205 = t251 * pkin(4) - t216 * pkin(9);
t214 = t215 ^ 2;
t155 = -t214 * pkin(4) + t184 * pkin(9) - t251 * t205 + t158;
t257 = sin(qJ(5));
t262 = cos(qJ(5));
t153 = t257 * t154 + t262 * t155;
t220 = -t265 * g(3) - t260 * t233;
t201 = -qJDD(2) * pkin(2) - t267 * pkin(7) + t242 * t284 - t220;
t179 = -t212 * pkin(3) - t237 * pkin(8) + t240 * t222 + t201;
t160 = -t184 * pkin(4) - t214 * pkin(9) + t216 * t205 + t179;
t195 = t257 * t215 + t262 * t216;
t165 = -t195 * qJD(5) + t262 * t184 - t257 * t185;
t194 = t262 * t215 - t257 * t216;
t166 = t194 * qJD(5) + t257 * t184 + t262 * t185;
t246 = qJD(5) + t251;
t173 = Ifges(6,5) * t195 + Ifges(6,6) * t194 + Ifges(6,3) * t246;
t175 = Ifges(6,1) * t195 + Ifges(6,4) * t194 + Ifges(6,5) * t246;
t228 = qJDD(5) + t234;
t142 = -mrSges(6,1) * t160 + mrSges(6,3) * t153 + Ifges(6,4) * t166 + Ifges(6,2) * t165 + Ifges(6,6) * t228 - t195 * t173 + t246 * t175;
t152 = t262 * t154 - t257 * t155;
t174 = Ifges(6,4) * t195 + Ifges(6,2) * t194 + Ifges(6,6) * t246;
t143 = mrSges(6,2) * t160 - mrSges(6,3) * t152 + Ifges(6,1) * t166 + Ifges(6,4) * t165 + Ifges(6,5) * t228 + t194 * t173 - t246 * t174;
t190 = Ifges(5,5) * t216 + Ifges(5,6) * t215 + Ifges(5,3) * t251;
t192 = Ifges(5,1) * t216 + Ifges(5,4) * t215 + Ifges(5,5) * t251;
t188 = -t246 * mrSges(6,2) + t194 * mrSges(6,3);
t189 = t246 * mrSges(6,1) - t195 * mrSges(6,3);
t277 = m(6) * t160 - t165 * mrSges(6,1) + t166 * mrSges(6,2) - t194 * t188 + t195 * t189;
t178 = -t194 * mrSges(6,1) + t195 * mrSges(6,2);
t149 = m(6) * t152 + t228 * mrSges(6,1) - t166 * mrSges(6,3) - t195 * t178 + t246 * t188;
t150 = m(6) * t153 - t228 * mrSges(6,2) + t165 * mrSges(6,3) + t194 * t178 - t246 * t189;
t278 = -t257 * t149 + t262 * t150;
t129 = -mrSges(5,1) * t179 + mrSges(5,3) * t158 + Ifges(5,4) * t185 + Ifges(5,2) * t184 + Ifges(5,6) * t234 - pkin(4) * t277 + pkin(9) * t278 + t262 * t142 + t257 * t143 - t216 * t190 + t251 * t192;
t141 = t262 * t149 + t257 * t150;
t191 = Ifges(5,4) * t216 + Ifges(5,2) * t215 + Ifges(5,6) * t251;
t130 = mrSges(5,2) * t179 - mrSges(5,3) * t157 + Ifges(5,1) * t185 + Ifges(5,4) * t184 + Ifges(5,5) * t234 - pkin(9) * t141 - t257 * t142 + t262 * t143 + t215 * t190 - t251 * t191;
t206 = Ifges(4,5) * t240 + Ifges(4,6) * t239 + Ifges(4,3) * t252;
t208 = Ifges(4,1) * t240 + Ifges(4,4) * t239 + Ifges(4,5) * t252;
t203 = -t251 * mrSges(5,2) + t215 * mrSges(5,3);
t204 = t251 * mrSges(5,1) - t216 * mrSges(5,3);
t273 = m(5) * t179 - t184 * mrSges(5,1) + t185 * mrSges(5,2) - t215 * t203 + t216 * t204 + t277;
t196 = -t215 * mrSges(5,1) + t216 * mrSges(5,2);
t138 = m(5) * t157 + t234 * mrSges(5,1) - t185 * mrSges(5,3) - t216 * t196 + t251 * t203 + t141;
t139 = m(5) * t158 - t234 * mrSges(5,2) + t184 * mrSges(5,3) + t215 * t196 - t251 * t204 + t278;
t279 = -t258 * t138 + t263 * t139;
t118 = -mrSges(4,1) * t201 + mrSges(4,3) * t187 + Ifges(4,4) * t213 + Ifges(4,2) * t212 + Ifges(4,6) * t238 - pkin(3) * t273 + pkin(8) * t279 + t263 * t129 + t258 * t130 - t240 * t206 + t252 * t208;
t134 = t263 * t138 + t258 * t139;
t207 = Ifges(4,4) * t240 + Ifges(4,2) * t239 + Ifges(4,6) * t252;
t119 = mrSges(4,2) * t201 - mrSges(4,3) * t186 + Ifges(4,1) * t213 + Ifges(4,4) * t212 + Ifges(4,5) * t238 - pkin(8) * t134 - t258 * t129 + t263 * t130 + t239 * t206 - t252 * t207;
t217 = -t239 * mrSges(4,1) + t240 * mrSges(4,2);
t218 = -t252 * mrSges(4,2) + t239 * mrSges(4,3);
t132 = m(4) * t186 + t238 * mrSges(4,1) - t213 * mrSges(4,3) - t240 * t217 + t252 * t218 + t134;
t219 = t252 * mrSges(4,1) - t240 * mrSges(4,3);
t133 = m(4) * t187 - t238 * mrSges(4,2) + t212 * mrSges(4,3) + t239 * t217 - t252 * t219 + t279;
t128 = -t259 * t132 + t264 * t133;
t145 = -m(4) * t201 + t212 * mrSges(4,1) - t213 * mrSges(4,2) + t239 * t218 - t240 * t219 - t273;
t230 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t260 + Ifges(3,2) * t265) * qJD(1);
t231 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t260 + Ifges(3,4) * t265) * qJD(1);
t285 = mrSges(3,1) * t220 - mrSges(3,2) * t221 + Ifges(3,5) * t243 + Ifges(3,6) * t244 + Ifges(3,3) * qJDD(2) + pkin(2) * t145 + pkin(7) * t128 + t264 * t118 + t259 * t119 + (t260 * t230 - t265 * t231) * qJD(1);
t241 = (-mrSges(3,1) * t265 + mrSges(3,2) * t260) * qJD(1);
t247 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t284;
t126 = m(3) * t221 - qJDD(2) * mrSges(3,2) + t244 * mrSges(3,3) - qJD(2) * t247 + t241 * t283 + t128;
t248 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t283;
t144 = m(3) * t220 + qJDD(2) * mrSges(3,1) - t243 * mrSges(3,3) + qJD(2) * t248 - t241 * t284 + t145;
t280 = t265 * t126 - t260 * t144;
t127 = t264 * t132 + t259 * t133;
t229 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t260 + Ifges(3,6) * t265) * qJD(1);
t115 = mrSges(3,2) * t232 - mrSges(3,3) * t220 + Ifges(3,1) * t243 + Ifges(3,4) * t244 + Ifges(3,5) * qJDD(2) - pkin(7) * t127 - qJD(2) * t230 - t259 * t118 + t264 * t119 + t229 * t283;
t274 = -mrSges(6,1) * t152 + mrSges(6,2) * t153 - Ifges(6,5) * t166 - Ifges(6,6) * t165 - Ifges(6,3) * t228 - t195 * t174 + t194 * t175;
t271 = -mrSges(5,1) * t157 + mrSges(5,2) * t158 - Ifges(5,5) * t185 - Ifges(5,6) * t184 - Ifges(5,3) * t234 - pkin(4) * t141 - t216 * t191 + t215 * t192 + t274;
t269 = mrSges(4,1) * t186 - mrSges(4,2) * t187 + Ifges(4,5) * t213 + Ifges(4,6) * t212 + Ifges(4,3) * t238 + pkin(3) * t134 + t240 * t207 - t239 * t208 - t271;
t117 = -mrSges(3,1) * t232 + mrSges(3,3) * t221 + Ifges(3,4) * t243 + Ifges(3,2) * t244 + Ifges(3,6) * qJDD(2) - pkin(2) * t127 + qJD(2) * t231 - t229 * t284 - t269;
t272 = -m(3) * t232 + t244 * mrSges(3,1) - t243 * mrSges(3,2) - t247 * t284 + t248 * t283 - t127;
t275 = mrSges(2,1) * t249 - mrSges(2,2) * t250 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t280 + t260 * t115 + t265 * t117;
t123 = m(2) * t249 + qJDD(1) * mrSges(2,1) - t268 * mrSges(2,2) + t272;
t122 = t260 * t126 + t265 * t144;
t120 = m(2) * t250 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t280;
t113 = mrSges(2,1) * g(3) + mrSges(2,3) * t250 + t268 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t122 - t285;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - pkin(6) * t122 + t265 * t115 - t260 * t117;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t266 * t112 - t261 * t113 - pkin(5) * (t261 * t120 + t266 * t123), t112, t115, t119, t130, t143; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t261 * t112 + t266 * t113 + pkin(5) * (t266 * t120 - t261 * t123), t113, t117, t118, t129, t142; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t275, t275, t285, t269, -t271, -t274;];
m_new = t1;
