% Calculate vector of centrifugal and coriolis load on the joints for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:51:50
% EndTime: 2018-11-23 14:51:56
% DurationCPUTime: 5.90s
% Computational Cost: add. (7617->446), mult. (19972->635), div. (0->0), fcn. (16704->14), ass. (0->226)
t182 = cos(pkin(6));
t169 = qJD(1) * t182 + qJD(2);
t177 = sin(pkin(13));
t179 = sin(pkin(6));
t186 = sin(qJ(3));
t190 = cos(qJ(3));
t180 = cos(pkin(13));
t181 = cos(pkin(7));
t245 = t180 * t181;
t196 = (t177 * t190 + t186 * t245) * t179;
t178 = sin(pkin(7));
t247 = t178 * t186;
t115 = qJD(1) * t196 + t169 * t247;
t185 = sin(qJ(4));
t235 = qJD(4) * t185;
t231 = pkin(4) * t235;
t315 = -t115 + t231;
t240 = qJD(1) * t179;
t226 = t180 * t240;
t248 = t177 * t186;
t114 = t190 * (t169 * t178 + t181 * t226) - t240 * t248;
t188 = cos(qJ(5));
t189 = cos(qJ(4));
t184 = sin(qJ(5));
t244 = t184 * t185;
t159 = -t188 * t189 + t244;
t291 = -pkin(10) - pkin(9);
t227 = qJD(4) * t291;
t163 = t185 * t227;
t167 = t291 * t185;
t168 = t291 * t189;
t207 = t188 * t167 + t168 * t184;
t220 = t189 * t227;
t306 = qJD(5) * t207 + t159 * t114 + t188 * t163 + t184 * t220;
t176 = qJD(4) + qJD(5);
t133 = t176 * t159;
t160 = t184 * t189 + t185 * t188;
t134 = t176 * t160;
t314 = pkin(5) * t134 + pkin(11) * t133 + t315;
t313 = t176 * Ifges(6,6) / 0.2e1;
t261 = t176 * Ifges(6,5);
t174 = -pkin(4) * t189 - pkin(3);
t99 = qJD(3) * t174 - t114;
t312 = t99 * mrSges(6,2) + t261 / 0.2e1;
t155 = t159 * qJD(3);
t156 = t160 * qJD(3);
t262 = t156 * Ifges(6,4);
t311 = t313 + t262 / 0.2e1 - t155 * Ifges(6,2) / 0.2e1;
t183 = sin(qJ(6));
t187 = cos(qJ(6));
t128 = pkin(5) * t159 - pkin(11) * t160 + t174;
t143 = t167 * t184 - t168 * t188;
t82 = t128 * t183 + t143 * t187;
t308 = -qJD(6) * t82 - t306 * t183 + t314 * t187;
t81 = t128 * t187 - t143 * t183;
t307 = qJD(6) * t81 + t314 * t183 + t306 * t187;
t305 = qJD(5) * t143 - t160 * t114 + t163 * t184 - t188 * t220;
t139 = -t156 * t183 + t176 * t187;
t140 = t156 * t187 + t176 * t183;
t264 = t156 * mrSges(6,3);
t257 = mrSges(6,1) * t176 + mrSges(7,1) * t139 - mrSges(7,2) * t140 - t264;
t121 = t133 * qJD(3);
t79 = qJD(6) * t139 - t121 * t187;
t80 = -qJD(6) * t140 + t121 * t183;
t36 = -mrSges(7,1) * t80 + mrSges(7,2) * t79;
t104 = qJD(3) * pkin(9) + t115;
t221 = pkin(10) * qJD(3) + t104;
t146 = t169 * t181 - t178 * t226;
t243 = t185 * t146;
t68 = t189 * t221 + t243;
t195 = t68 * qJD(4);
t197 = (t190 * t245 - t248) * t179;
t246 = t178 * t190;
t101 = (qJD(1) * t197 + t169 * t246) * qJD(3);
t253 = t101 * t185;
t258 = t188 * t68;
t141 = t189 * t146;
t209 = t221 * t185;
t67 = t141 - t209;
t64 = qJD(4) * pkin(4) + t67;
t31 = t184 * t64 + t258;
t242 = qJD(4) * t141 + t189 * t101;
t40 = -qJD(4) * t209 + t242;
t9 = t184 * t40 - t188 * (-t195 - t253) + t31 * qJD(5);
t304 = m(7) * t9 + t36;
t303 = t182 * t246 + t197;
t29 = pkin(11) * t176 + t31;
t62 = pkin(5) * t155 - pkin(11) * t156 + t99;
t12 = -t183 * t29 + t187 * t62;
t122 = t134 * qJD(3);
t102 = t115 * qJD(3);
t95 = qJD(3) * t231 + t102;
t41 = pkin(5) * t122 + pkin(11) * t121 + t95;
t259 = t184 * t68;
t30 = t188 * t64 - t259;
t8 = qJD(5) * t30 - t101 * t244 - t184 * t195 + t188 * t40;
t2 = qJD(6) * t12 + t183 * t41 + t187 * t8;
t13 = t183 * t62 + t187 * t29;
t3 = -qJD(6) * t13 - t183 * t8 + t187 * t41;
t302 = -t3 * t183 + t187 * t2;
t301 = -t99 * mrSges(6,1) - t12 * mrSges(7,1) + t13 * mrSges(7,2) + t311;
t150 = qJD(6) + t155;
t217 = mrSges(7,1) * t183 + mrSges(7,2) * t187;
t28 = -pkin(5) * t176 - t30;
t204 = t28 * t217;
t212 = Ifges(7,5) * t187 - Ifges(7,6) * t183;
t271 = Ifges(7,4) * t187;
t214 = -Ifges(7,2) * t183 + t271;
t272 = Ifges(7,4) * t183;
t216 = Ifges(7,1) * t187 - t272;
t281 = t187 / 0.2e1;
t282 = -t183 / 0.2e1;
t286 = t140 / 0.2e1;
t273 = Ifges(7,4) * t140;
t72 = Ifges(7,2) * t139 + Ifges(7,6) * t150 + t273;
t137 = Ifges(7,4) * t139;
t73 = t140 * Ifges(7,1) + t150 * Ifges(7,5) + t137;
t300 = (-t12 * t187 - t13 * t183) * mrSges(7,3) + t150 * t212 / 0.2e1 + t216 * t286 + t139 * t214 / 0.2e1 + t204 + t73 * t281 + t72 * t282;
t299 = -Ifges(5,1) / 0.2e1;
t298 = -t72 / 0.2e1;
t237 = qJD(3) * t189;
t297 = -Ifges(5,4) * t237 / 0.2e1;
t296 = -t12 * t183 + t13 * t187;
t295 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t79 + Ifges(7,6) * t80;
t293 = t79 / 0.2e1;
t292 = t80 / 0.2e1;
t131 = t182 * t247 + t196;
t149 = -t178 * t179 * t180 + t181 * t182;
t89 = -t131 * t185 + t149 * t189;
t90 = t131 * t189 + t149 * t185;
t44 = t184 * t90 - t188 * t89;
t290 = t44 * t9;
t289 = t122 / 0.2e1;
t288 = -t139 / 0.2e1;
t287 = -t140 / 0.2e1;
t285 = -t150 / 0.2e1;
t284 = t155 / 0.2e1;
t283 = -t156 / 0.2e1;
t151 = t181 * t189 - t185 * t247;
t152 = t181 * t185 + t189 * t247;
t208 = t188 * t151 - t152 * t184;
t280 = t208 * t9;
t279 = t207 * t9;
t275 = mrSges(6,3) * t155;
t274 = Ifges(5,4) * t185;
t148 = Ifges(6,4) * t155;
t268 = t139 * Ifges(7,6);
t267 = t140 * Ifges(7,5);
t266 = t150 * Ifges(7,3);
t263 = t156 * Ifges(6,1);
t232 = qJD(3) * qJD(4);
t158 = (mrSges(5,1) * t185 + mrSges(5,2) * t189) * t232;
t70 = mrSges(6,1) * t122 - mrSges(6,2) * t121;
t256 = t70 + t158;
t255 = Ifges(5,5) * qJD(4);
t254 = Ifges(5,6) * qJD(4);
t252 = t102 * t303;
t251 = t102 * t190;
t250 = t155 * t183;
t249 = t155 * t187;
t126 = mrSges(6,1) * t155 + mrSges(6,2) * t156;
t241 = t126 + (-mrSges(5,1) * t189 + mrSges(5,2) * t185) * qJD(3);
t239 = qJD(3) * t185;
t238 = qJD(3) * t186;
t236 = qJD(3) * t190;
t234 = qJD(6) * t183;
t233 = qJD(6) * t187;
t225 = t178 * t238;
t224 = t178 * t236;
t223 = t255 / 0.2e1;
t222 = -t254 / 0.2e1;
t127 = pkin(5) * t156 + pkin(11) * t155;
t218 = mrSges(7,1) * t187 - mrSges(7,2) * t183;
t215 = Ifges(7,1) * t183 + t271;
t213 = Ifges(7,2) * t187 + t272;
t211 = Ifges(7,5) * t183 + Ifges(7,6) * t187;
t45 = t184 * t89 + t188 * t90;
t34 = -t183 * t45 - t187 * t303;
t35 = -t183 * t303 + t187 * t45;
t75 = t104 * t189 + t243;
t113 = t151 * t184 + t152 * t188;
t93 = -t113 * t183 - t187 * t246;
t206 = -t113 * t187 + t183 * t246;
t103 = -qJD(3) * pkin(3) - t114;
t74 = -t104 * t185 + t141;
t203 = t74 * mrSges(5,3) + t239 * t299 + t297 - t255 / 0.2e1 - t103 * mrSges(5,2);
t202 = t75 * mrSges(5,3) + t254 / 0.2e1 + (t189 * Ifges(5,2) + t274) * qJD(3) / 0.2e1 - t103 * mrSges(5,1);
t46 = mrSges(7,1) * t122 - mrSges(7,3) * t79;
t47 = -mrSges(7,2) * t122 + mrSges(7,3) * t80;
t96 = -mrSges(7,2) * t150 + mrSges(7,3) * t139;
t97 = mrSges(7,1) * t150 - mrSges(7,3) * t140;
t193 = -t97 * t233 - t96 * t234 + t187 * t47 + m(7) * (-t12 * t233 - t13 * t234 + t302) - t183 * t46;
t110 = -t148 + t261 + t263;
t24 = t79 * Ifges(7,4) + t80 * Ifges(7,2) + t122 * Ifges(7,6);
t25 = t79 * Ifges(7,1) + t80 * Ifges(7,4) + t122 * Ifges(7,5);
t71 = t266 + t267 + t268;
t192 = -Ifges(6,6) * t122 - Ifges(6,5) * t121 + t250 * t298 - t8 * mrSges(6,2) + t211 * t289 + t213 * t292 + t215 * t293 + t24 * t281 - t30 * t275 + t183 * t25 / 0.2e1 + t73 * t249 / 0.2e1 + (-mrSges(6,1) - t218) * t9 + (-t148 + t110) * t284 + (-t262 + t71) * t283 + (-Ifges(6,1) * t283 - t212 * t285 - t214 * t288 - t216 * t287 + t204 + t312) * t155 + (Ifges(7,5) * t287 - Ifges(6,2) * t284 + Ifges(7,6) * t288 + Ifges(7,3) * t285 + t301 + t313) * t156 + (-t12 * t249 - t13 * t250 + t302) * mrSges(7,3) + t300 * qJD(6);
t191 = qJD(3) ^ 2;
t166 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t237;
t165 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t239;
t144 = -mrSges(6,2) * t176 - t275;
t136 = -qJD(4) * t152 - t185 * t224;
t135 = qJD(4) * t151 + t189 * t224;
t124 = t131 * qJD(3);
t123 = t303 * qJD(3);
t118 = Ifges(7,3) * t122;
t108 = pkin(4) * t239 + t127;
t52 = -qJD(4) * t90 - t123 * t185;
t51 = qJD(4) * t89 + t123 * t189;
t49 = qJD(5) * t113 + t135 * t184 - t188 * t136;
t48 = qJD(5) * t208 + t135 * t188 + t136 * t184;
t43 = -qJD(4) * t75 - t253;
t42 = -t104 * t235 + t242;
t33 = t188 * t67 - t259;
t32 = t184 * t67 + t258;
t27 = qJD(6) * t206 - t183 * t48 + t187 * t225;
t26 = qJD(6) * t93 + t183 * t225 + t187 * t48;
t17 = t127 * t183 + t187 * t30;
t16 = t127 * t187 - t183 * t30;
t15 = t108 * t183 + t187 * t33;
t14 = t108 * t187 - t183 * t33;
t11 = qJD(5) * t45 + t184 * t51 - t188 * t52;
t10 = -qJD(5) * t44 + t184 * t52 + t188 * t51;
t5 = -qJD(6) * t35 - t10 * t183 + t124 * t187;
t4 = qJD(6) * t34 + t10 * t187 + t124 * t183;
t1 = [t10 * t144 + t52 * t165 + t51 * t166 + t34 * t46 + t35 * t47 + t44 * t36 + t4 * t96 + t5 * t97 - t256 * t303 + t241 * t124 - t257 * t11 + (-t121 * t44 - t122 * t45) * mrSges(6,3) + (-t124 * mrSges(4,1) - t123 * mrSges(4,2) + (-t185 * t90 - t189 * t89) * qJD(4) * mrSges(5,3)) * qJD(3) + m(6) * (t10 * t31 - t11 * t30 + t124 * t99 - t303 * t95 + t45 * t8 + t290) + m(7) * (t11 * t28 + t12 * t5 + t13 * t4 + t2 * t35 + t3 * t34 + t290) + m(4) * (t101 * t131 - t114 * t124 + t115 * t123 - t252) + m(5) * (t103 * t124 + t42 * t90 + t43 * t89 + t51 * t75 + t52 * t74 - t252); -t208 * t36 + t135 * t166 + t136 * t165 + t48 * t144 + t26 * t96 + t27 * t97 + t93 * t46 - t206 * t47 - t257 * t49 + (-t113 * t122 + t121 * t208) * mrSges(6,3) + (-t151 * t189 - t152 * t185) * mrSges(5,3) * t232 + m(6) * (t113 * t8 - t30 * t49 + t31 * t48 - t280) + m(7) * (t12 * t27 + t13 * t26 - t2 * t206 + t28 * t49 + t3 * t93 - t280) + m(5) * (t135 * t75 + t136 * t74 + t151 * t43 + t152 * t42) + ((-mrSges(4,2) * t191 - t256) * t190 + (-t191 * mrSges(4,1) + qJD(3) * t241) * t186 + m(6) * (-t190 * t95 + t238 * t99) + m(4) * (t101 * t186 - t114 * t238 + t115 * t236 - t251) + m(5) * (t103 * t238 - t251)) * t178; (t102 * mrSges(5,2) - t43 * mrSges(5,3) + t114 * t165 + (-pkin(9) * t166 + pkin(4) * t126 - 0.3e1 / 0.2e1 * Ifges(5,4) * t239 + t222 - t202) * qJD(4)) * t185 + (mrSges(4,1) * qJD(3) - t241) * t115 - pkin(3) * t158 - t207 * t36 - t102 * mrSges(4,1) + t81 * t46 + t82 * t47 + t308 * t97 + t307 * t96 + m(5) * (-pkin(3) * t102 + (-t43 * t185 - t75 * t235) * pkin(9)) + (t71 / 0.2e1 + t266 / 0.2e1 + t268 / 0.2e1 + t267 / 0.2e1 - t301 - t311) * t134 + t306 * t144 + (t121 * t207 - t122 * t143 + t133 * t30 - t134 * t31) * mrSges(6,3) - (-t148 / 0.2e1 + t263 / 0.2e1 + t110 / 0.2e1 + t300 + t312) * t133 + (qJD(3) * t114 - t101) * mrSges(4,2) + (-Ifges(6,1) * t121 - Ifges(6,4) * t122 + t95 * mrSges(6,2) + t25 * t281 + t24 * t282 + t212 * t289 + t216 * t293 + t214 * t292 + (mrSges(6,3) + t217) * t9 + (-t183 * t2 - t187 * t3) * mrSges(7,3) + (-mrSges(7,3) * t296 + t187 * t298 + t211 * t285 + t213 * t288 + t215 * t287 + t28 * t218 + t73 * t282) * qJD(6)) * t160 + (t118 / 0.2e1 + Ifges(6,4) * t121 + t95 * mrSges(6,1) - t8 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t122 + t295) * t159 - m(5) * (-t114 * t185 * t74 + t103 * t115) + t174 * t70 + (-t102 * mrSges(5,1) + (t223 + (-m(5) * t74 - t165) * pkin(9) + (0.3e1 / 0.2e1 * Ifges(5,4) * t189 + (-0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t185) * qJD(3) - t203) * qJD(4) + (m(5) * pkin(9) + mrSges(5,3)) * t42 + (-m(5) * t75 - t166) * t114) * t189 - t257 * t305 + (t12 * t308 + t13 * t307 + t2 * t82 + t28 * t305 + t3 * t81 - t279) * m(7) + (t143 * t8 + t174 * t95 - t305 * t30 + t306 * t31 + t315 * t99 - t279) * m(6); t75 * t165 - t74 * t166 - t33 * t144 - t14 * t97 - t15 * t96 - t42 * mrSges(5,2) + t43 * mrSges(5,1) + t192 - m(6) * (-t30 * t32 + t31 * t33) - m(7) * (t12 * t14 + t13 * t15 + t28 * t32) + t257 * t32 + t31 * t264 + t193 * (pkin(4) * t184 + pkin(11)) + ((t223 + t297 + t203) * t189 + (t222 + (t274 / 0.2e1 + (t299 + Ifges(5,2) / 0.2e1) * t189) * qJD(3) + (-m(6) * t99 - t126) * pkin(4) + t202) * t185) * qJD(3) + (m(6) * (t184 * t8 - t188 * t9) + (t121 * t188 - t122 * t184) * mrSges(6,3) + ((-m(6) * t30 + m(7) * t28 - t257) * t184 + (m(6) * t31 + m(7) * t296 - t183 * t97 + t187 * t96 + t144) * t188) * qJD(5)) * pkin(4) + t304 * (-pkin(4) * t188 - pkin(5)); -t30 * t144 - t17 * t96 - t16 * t97 + t192 + (t257 + t264) * t31 - m(7) * (t12 * t16 + t13 * t17 + t28 * t31) + t193 * pkin(11) - t304 * pkin(5); t118 - t28 * (mrSges(7,1) * t140 + mrSges(7,2) * t139) + (Ifges(7,5) * t139 - Ifges(7,6) * t140) * t285 - t12 * t96 + t13 * t97 + (Ifges(7,1) * t139 - t273) * t287 + t72 * t286 + (t12 * t139 + t13 * t140) * mrSges(7,3) + (-Ifges(7,2) * t140 + t137 + t73) * t288 + t295;];
tauc  = t1(:);
