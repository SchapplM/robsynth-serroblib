% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:22
% DurationCPUTime: 5.87s
% Computational Cost: add. (4286->423), mult. (10696->558), div. (0->0), fcn. (6968->6), ass. (0->208)
t298 = mrSges(5,2) - mrSges(4,1);
t297 = -Ifges(4,5) + Ifges(5,4);
t296 = -Ifges(4,6) + Ifges(5,5);
t173 = sin(qJ(3));
t258 = cos(qJ(3));
t259 = cos(qJ(2));
t203 = t258 * t259;
t192 = qJD(1) * t203;
t174 = sin(qJ(2));
t225 = qJD(1) * t174;
t136 = t173 * t225 - t192;
t171 = qJD(2) + qJD(3);
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t118 = t136 * t172 + t171 * t175;
t131 = Ifges(4,4) * t136;
t146 = t173 * t259 + t174 * t258;
t137 = t146 * qJD(1);
t133 = qJD(5) + t137;
t287 = t133 * Ifges(6,3);
t117 = t136 * t175 - t171 * t172;
t288 = t117 * Ifges(6,6);
t295 = t137 * Ifges(4,1) + t171 * Ifges(4,5) + t118 * Ifges(6,5) - t131 + t287 + t288;
t273 = pkin(3) + pkin(8);
t154 = (-pkin(7) - pkin(6)) * t174;
t148 = qJD(1) * t154;
t142 = qJD(2) * pkin(2) + t148;
t217 = t259 * pkin(6);
t155 = pkin(7) * t259 + t217;
t149 = t155 * qJD(1);
t227 = t173 * t149;
t109 = -t258 * t142 + t227;
t254 = t137 * pkin(4);
t191 = t109 + t254;
t292 = qJD(4) + t191;
t53 = -t171 * t273 + t292;
t167 = -pkin(2) * t259 - pkin(1);
t153 = qJD(1) * t167;
t180 = -t137 * qJ(4) + t153;
t54 = t136 * t273 + t180;
t12 = -t172 * t54 + t175 * t53;
t13 = t172 * t53 + t175 * t54;
t194 = t12 * t172 - t13 * t175;
t115 = t171 * t146;
t102 = t115 * qJD(1);
t226 = t173 * t174;
t199 = t171 * t226;
t101 = qJD(1) * t199 - t171 * t192;
t220 = qJD(1) * qJD(2);
t208 = t174 * t220;
t162 = pkin(2) * t208;
t187 = qJ(4) * t101 - qJD(4) * t137 + t162;
t10 = t102 * t273 + t187;
t139 = t258 * t149;
t110 = t173 * t142 + t139;
t181 = qJD(2) * t155;
t141 = t258 * t181;
t150 = qJD(2) * t154;
t143 = qJD(1) * t150;
t47 = qJD(1) * t141 + qJD(3) * t110 + t173 * t143;
t21 = -t101 * pkin(4) + t47;
t235 = qJD(5) * t13;
t2 = -t10 * t172 + t175 * t21 - t235;
t1 = qJD(5) * t12 + t10 * t175 + t172 * t21;
t256 = t1 * t172;
t178 = m(6) * (-qJD(5) * t194 + t2 * t175 + t256);
t51 = qJD(5) * t117 + t102 * t172;
t22 = -mrSges(6,1) * t101 - mrSges(6,3) * t51;
t52 = -qJD(5) * t118 + t102 * t175;
t23 = mrSges(6,2) * t101 + mrSges(6,3) * t52;
t294 = t172 * t23 + t175 * t22 + t178;
t198 = mrSges(6,1) * t175 - mrSges(6,2) * t172;
t250 = Ifges(6,4) * t118;
t42 = Ifges(6,2) * t117 + Ifges(6,6) * t133 + t250;
t255 = t136 * pkin(4);
t93 = -t171 * qJ(4) - t110;
t64 = -t93 - t255;
t293 = -t175 * t42 / 0.2e1 + t198 * t64;
t79 = t136 * pkin(3) + t180;
t291 = t153 * mrSges(4,1) - t79 * mrSges(5,2);
t290 = t153 * mrSges(4,2) - t13 * mrSges(6,2) + t109 * mrSges(4,3) - t79 * mrSges(5,3);
t275 = t51 / 0.2e1;
t274 = t52 / 0.2e1;
t272 = -t101 / 0.2e1;
t252 = mrSges(5,1) * t136;
t123 = -mrSges(5,3) * t171 + t252;
t63 = -mrSges(6,1) * t117 + mrSges(6,2) * t118;
t236 = t63 - t123;
t121 = -mrSges(4,2) * t171 - mrSges(4,3) * t136;
t285 = t121 - t123;
t244 = t137 * mrSges(4,3);
t284 = -t137 * mrSges(5,1) - t298 * t171 - t244;
t211 = qJD(1) * t259;
t283 = t174 * pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t211) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t225) * t217;
t281 = -qJD(4) - t109;
t119 = -t258 * t154 + t155 * t173;
t120 = t173 * t154 + t155 * t258;
t280 = -t101 * t119 - t102 * t120 + t146 * t47;
t80 = -pkin(3) * t171 - t281;
t279 = m(5) * t80 - t284;
t278 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t276 = Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t272;
t270 = -t117 / 0.2e1;
t269 = -t118 / 0.2e1;
t268 = t118 / 0.2e1;
t267 = -t133 / 0.2e1;
t266 = -t136 / 0.2e1;
t265 = t136 / 0.2e1;
t264 = -t137 / 0.2e1;
t263 = t137 / 0.2e1;
t261 = -t171 / 0.2e1;
t260 = t171 / 0.2e1;
t257 = pkin(2) * t173;
t251 = Ifges(3,4) * t174;
t249 = Ifges(6,4) * t172;
t248 = Ifges(6,4) * t175;
t247 = t101 * mrSges(5,1);
t245 = t119 * t47;
t243 = t137 * Ifges(4,4);
t242 = t137 * Ifges(5,6);
t239 = t175 * mrSges(6,3);
t237 = t93 * t137;
t232 = t115 * t172;
t231 = t115 * t175;
t230 = t137 * t172;
t145 = -t203 + t226;
t229 = t145 * t172;
t228 = t145 * t175;
t224 = qJD(2) * t174;
t223 = qJD(3) * t173;
t222 = qJD(5) * t172;
t221 = qJD(5) * t175;
t219 = Ifges(6,5) * t51 + Ifges(6,6) * t52 - Ifges(6,3) * t101;
t216 = t258 * pkin(2);
t170 = pkin(2) * t224;
t169 = pkin(2) * t225;
t214 = Ifges(3,4) * t259;
t210 = qJD(2) * t259;
t209 = qJD(3) * t258;
t206 = -t222 / 0.2e1;
t111 = t148 * t173 + t139;
t204 = pkin(2) * t209;
t166 = -t216 - pkin(3);
t202 = qJD(1) * t210;
t197 = Ifges(6,1) * t172 + t248;
t196 = Ifges(6,2) * t175 + t249;
t195 = Ifges(6,5) * t172 + Ifges(6,6) * t175;
t104 = pkin(3) * t137 + qJ(4) * t136;
t185 = -t146 * qJ(4) + t167;
t73 = t145 * t273 + t185;
t91 = pkin(4) * t146 + t119;
t26 = t172 * t91 + t175 * t73;
t25 = -t172 * t73 + t175 * t91;
t76 = -mrSges(6,2) * t133 + mrSges(6,3) * t117;
t77 = mrSges(6,1) * t133 - mrSges(6,3) * t118;
t193 = -t172 * t77 + t175 * t76;
t81 = t104 + t169;
t190 = Ifges(3,5) * t259 - Ifges(3,6) * t174;
t112 = t148 * t258 - t227;
t189 = t145 * t221 + t232;
t188 = t145 * t222 - t231;
t114 = -t171 * t203 + t199;
t186 = qJ(4) * t114 - qJD(4) * t146 + t170;
t138 = t173 * t181;
t46 = -qJD(1) * t138 + t142 * t209 + t258 * t143 - t149 * t223;
t55 = -t258 * t150 - t154 * t209 + t155 * t223 + t138;
t184 = pkin(1) * (mrSges(3,1) * t174 + mrSges(3,2) * t259);
t183 = t174 * (Ifges(3,1) * t259 - t251);
t182 = (Ifges(3,2) * t259 + t251) * qJD(1);
t39 = -qJD(4) * t171 - t46;
t56 = qJD(3) * t120 + t173 * t150 + t141;
t130 = Ifges(5,6) * t136;
t16 = -pkin(4) * t102 - t39;
t116 = Ifges(6,4) * t117;
t43 = Ifges(6,1) * t118 + Ifges(6,5) * t133 + t116;
t7 = t51 * Ifges(6,4) + t52 * Ifges(6,2) - t101 * Ifges(6,6);
t86 = t171 * Ifges(5,5) + t136 * Ifges(5,3) - t242;
t87 = t171 * Ifges(5,4) - t137 * Ifges(5,2) + t130;
t88 = -t136 * Ifges(4,2) + t171 * Ifges(4,6) + t243;
t177 = -t39 * mrSges(5,3) - t46 * mrSges(4,2) + t80 * t252 + (Ifges(6,5) * t175 - Ifges(6,6) * t172) * t272 + (-Ifges(6,2) * t172 + t248) * t274 + (Ifges(6,1) * t175 - t249) * t275 + t175 * t276 - t172 * t7 / 0.2e1 + t16 * (mrSges(6,1) * t172 + mrSges(6,2) * t175) + t298 * t47 + (t206 - t230 / 0.2e1) * t43 + (t130 + t87) * t266 + (-t243 + t86) * t264 + (t242 + t88) * t263 + (-Ifges(4,1) * t264 - Ifges(6,5) * t269 + Ifges(5,2) * t263 - Ifges(6,6) * t270 - Ifges(6,3) * t267 + t297 * t261 + t290) * t136 + (Ifges(5,3) * t266 - t13 * t239 + t195 * t267 + t196 * t270 + t197 * t269 + t296 * t261 - t291 + t293) * t137 + (mrSges(6,1) * t136 + (t222 + t230) * mrSges(6,3)) * t12 + t296 * t102 + t297 * t101 + t293 * qJD(5) + (-Ifges(4,2) * t137 - t131 + t295) * t265 - (t117 * t196 + t118 * t197 + t133 * t195) * qJD(5) / 0.2e1;
t168 = Ifges(3,4) * t211;
t164 = qJ(4) + t257;
t158 = t204 + qJD(4);
t135 = Ifges(3,1) * t225 + Ifges(3,5) * qJD(2) + t168;
t134 = Ifges(3,6) * qJD(2) + t182;
t132 = t137 * pkin(8);
t108 = t145 * pkin(3) + t185;
t107 = -mrSges(5,2) * t136 - mrSges(5,3) * t137;
t106 = mrSges(4,1) * t136 + mrSges(4,2) * t137;
t92 = -t145 * pkin(4) + t120;
t75 = t112 - t254;
t74 = t111 - t255;
t72 = t110 - t255;
t69 = t104 + t132;
t62 = t132 + t81;
t34 = pkin(3) * t115 + t186;
t31 = -t114 * pkin(4) + t56;
t30 = -pkin(4) * t115 - t55;
t24 = pkin(3) * t102 + t187;
t20 = t172 * t72 + t175 * t69;
t19 = -t172 * t69 + t175 * t72;
t18 = t172 * t74 + t175 * t62;
t17 = -t172 * t62 + t175 * t74;
t15 = t115 * t273 + t186;
t9 = -mrSges(6,1) * t52 + mrSges(6,2) * t51;
t4 = -qJD(5) * t26 - t15 * t172 + t175 * t31;
t3 = qJD(5) * t25 + t15 * t175 + t172 * t31;
t5 = [(-0.2e1 * t184 + t183) * t220 + (t115 * t93 + t280) * mrSges(5,1) + m(4) * (t120 * t46 + 0.2e1 * t153 * t170 + t245) + m(5) * (t108 * t24 - t120 * t39 + t34 * t79 + t245) + (t135 + qJD(1) * (Ifges(3,1) * t174 + t214)) * t210 / 0.2e1 - (t182 + t134) * t224 / 0.2e1 + (qJD(5) * t43 + t7) * t228 / 0.2e1 + (t265 * Ifges(5,6) - t295 / 0.2e1 + Ifges(5,2) * t264 - t287 / 0.2e1 - t288 / 0.2e1 - t12 * mrSges(6,1) + t87 / 0.2e1 - Ifges(4,1) * t263 - Ifges(4,4) * t266 - Ifges(6,5) * t268 + t297 * t260 - t290 - t80 * mrSges(5,1)) * t114 + (t145 * t101 + t115 * t264) * Ifges(5,6) + t30 * t63 + t25 * t22 + t26 * t23 + t117 * (Ifges(6,4) * t189 - Ifges(6,2) * t188) / 0.2e1 + t133 * (Ifges(6,5) * t189 - Ifges(6,6) * t188) / 0.2e1 + m(6) * (t1 * t26 + t12 * t4 + t13 * t3 + t16 * t92 + t2 * t25 + t30 * t64) + t106 * t170 + (Ifges(5,3) * t265 + t86 / 0.2e1 - t88 / 0.2e1 - Ifges(4,4) * t263 - Ifges(4,2) * t266 + t296 * t260 + t291) * t115 + (t1 * t228 - t12 * t189 - t13 * t188 - t2 * t229) * mrSges(6,3) + (mrSges(4,1) * t162 + t39 * mrSges(5,1) - t24 * mrSges(5,2) - t46 * mrSges(4,3) - t16 * t198 + t195 * t272 + t196 * t274 + t197 * t275 + t42 * t206 + (Ifges(5,3) + Ifges(4,2)) * t102 + (-t272 + t101 / 0.2e1) * Ifges(4,4)) * t145 + (Ifges(6,1) * t189 - Ifges(6,4) * t188) * t268 + t43 * t232 / 0.2e1 + t42 * t231 / 0.2e1 + (-Ifges(3,2) * t174 + t214) * t202 + t3 * t76 + t4 * t77 + t92 * t9 + t34 * t107 + t108 * (-mrSges(5,2) * t102 + mrSges(5,3) * t101) + t64 * (mrSges(6,1) * t188 + mrSges(6,2) * t189) + t229 * t276 + (Ifges(6,6) * t274 + Ifges(6,5) * t275 - t24 * mrSges(5,3) + mrSges(4,2) * t162 + (Ifges(6,3) + Ifges(4,1)) * t272 + t278 + t219 / 0.2e1 + (-Ifges(5,2) - Ifges(4,1) / 0.2e1) * t101 + (-Ifges(5,6) - Ifges(4,4)) * t102) * t146 + t167 * (mrSges(4,1) * t102 - mrSges(4,2) * t101) + (m(4) * t109 + t279) * t56 + (-t110 * t115 + t280) * mrSges(4,3) + (t190 * qJD(2) / 0.2e1 - t283) * qJD(2) + (-m(4) * t110 + m(5) * t93 - t285) * t55; (t221 * t76 - t222 * t77 + t294) * (-pkin(8) + t166) - (-Ifges(3,2) * t225 + t135 + t168) * t211 / 0.2e1 + (-t13 * t221 - t256) * mrSges(6,3) + (-t102 * t164 - t237) * mrSges(5,1) + t177 + (-t111 * t80 - t164 * t39 + t166 * t47 - t79 * t81 + (t112 - t158) * t93) * m(5) + (-t12 * t17 - t13 * t18 + t16 * t164 + (-t75 + t158) * t64) * m(6) + (-t109 * t111 - t110 * t112 - t153 * t169 + (-t258 * t47 + t173 * t46 + (t109 * t173 + t110 * t258) * qJD(3)) * pkin(2)) * m(4) + (t283 + (-t183 / 0.2e1 + t184) * qJD(1)) * qJD(1) + (-mrSges(3,1) * t202 + mrSges(3,2) * t208) * pkin(6) + (t101 * t216 - t102 * t257) * mrSges(4,3) - t166 * t247 - t2 * t239 + t134 * t225 / 0.2e1 - t190 * t220 / 0.2e1 - t106 * t169 - Ifges(3,6) * t208 + (m(6) * (t12 * t175 + t13 * t172) + t175 * t77 + t172 * t76 + t279) * pkin(2) * t223 - t75 * t63 - t18 * t76 - t17 * t77 - t81 * t107 + Ifges(3,5) * t202 + t110 * t244 + t164 * t9 + t284 * t111 - t285 * t112 + t236 * t158 + t121 * t204; qJ(4) * t9 + (-(qJD(5) * t76 + t22) * t273 + (-t2 - t235) * mrSges(6,3)) * t175 + (-t1 * mrSges(6,3) - (-qJD(5) * t77 + t23) * t273) * t172 + (t244 + t284) * t110 + t285 * t109 + t236 * qJD(4) + t177 + (pkin(3) * t101 - qJ(4) * t102 - t237) * mrSges(5,1) - t273 * t178 + t191 * t63 - t20 * t76 - t19 * t77 - t104 * t107 + (t16 * qJ(4) - t12 * t19 - t13 * t20 + t292 * t64) * m(6) + (-t47 * pkin(3) - t39 * qJ(4) - t79 * t104 - t80 * t110 + t281 * t93) * m(5); -t247 - t236 * t171 + t193 * qJD(5) + (t107 + t193) * t137 - m(6) * (t137 * t194 + t64 * t171) + (t79 * t137 + t93 * t171 + t47) * m(5) + t294; -t64 * (mrSges(6,1) * t118 + mrSges(6,2) * t117) + (Ifges(6,1) * t117 - t250) * t269 + t42 * t268 + (Ifges(6,5) * t117 - Ifges(6,6) * t118) * t267 - t12 * t76 + t13 * t77 + (t117 * t12 + t118 * t13) * mrSges(6,3) + t219 + (-Ifges(6,2) * t118 + t116 + t43) * t270 + t278;];
tauc = t5(:);
