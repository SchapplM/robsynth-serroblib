% Calculate time derivative of joint inertia matrix for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:09
% EndTime: 2019-03-09 16:39:19
% DurationCPUTime: 4.38s
% Computational Cost: add. (7039->419), mult. (15887->590), div. (0->0), fcn. (14824->8), ass. (0->186)
t197 = sin(qJ(3));
t198 = sin(qJ(2));
t200 = cos(qJ(3));
t201 = cos(qJ(2));
t176 = t197 * t201 + t200 * t198;
t285 = qJD(2) + qJD(3);
t152 = t285 * t176;
t298 = Ifges(7,4) + Ifges(6,5);
t302 = t298 * t152;
t301 = Ifges(6,1) + Ifges(7,1);
t300 = -Ifges(6,4) + Ifges(7,5);
t296 = -Ifges(6,6) + Ifges(7,6);
t289 = (mrSges(6,3) + mrSges(7,2));
t299 = 2 * t289;
t297 = Ifges(7,2) + Ifges(6,3);
t175 = t197 * t198 - t200 * t201;
t151 = t285 * t175;
t194 = sin(pkin(10));
t195 = cos(pkin(10));
t196 = sin(qJ(5));
t199 = cos(qJ(5));
t174 = t194 * t199 + t195 * t196;
t169 = t174 * qJD(5);
t237 = t199 * t195;
t173 = t194 * t196 - t237;
t50 = t151 * t173 - t169 * t176;
t230 = qJD(5) * t199;
t231 = qJD(5) * t196;
t241 = t176 * t195;
t242 = t176 * t194;
t51 = -t151 * t174 + t230 * t241 - t231 * t242;
t295 = t300 * t51 + t301 * t50 + t302;
t121 = t174 * t176;
t122 = t173 * t176;
t294 = t300 * t121 - t301 * t122 + t298 * t175;
t293 = -mrSges(5,1) * t195 + mrSges(5,2) * t194 - mrSges(4,1);
t168 = t173 * qJD(5);
t292 = -t301 * t168 + t300 * t169;
t291 = t300 * t173 + t301 * t174;
t290 = -t298 * t168 + t296 * t169;
t288 = m(7) * qJD(6);
t269 = -pkin(8) - pkin(7);
t182 = t269 * t198;
t183 = t269 * t201;
t287 = t200 * t182 + t183 * t197;
t190 = -pkin(2) * t201 - pkin(1);
t138 = t175 * pkin(3) - t176 * qJ(4) + t190;
t157 = t182 * t197 - t183 * t200;
t88 = t195 * t138 - t157 * t194;
t89 = t194 * t138 + t195 * t157;
t286 = -t194 * t88 + t195 * t89;
t284 = m(7) * qJ(6) + mrSges(7,3);
t146 = mrSges(6,1) * t173 + mrSges(6,2) * t174;
t252 = pkin(2) * qJD(3);
t283 = ((t146 + t293) * t197 - mrSges(4,2) * t200) * t252;
t244 = t151 * t195;
t74 = pkin(2) * qJD(2) * t198 + pkin(3) * t152 + qJ(4) * t151 - qJD(4) * t176;
t225 = qJD(2) * t269;
t178 = t198 * t225;
t216 = t201 * t225;
t96 = qJD(3) * t287 + t200 * t178 + t197 * t216;
t31 = -t194 * t96 + t195 * t74;
t21 = pkin(4) * t152 + pkin(9) * t244 + t31;
t63 = pkin(4) * t175 - pkin(9) * t241 + t88;
t73 = -pkin(9) * t242 + t89;
t257 = t196 * t63 + t199 * t73;
t245 = t151 * t194;
t32 = t194 * t74 + t195 * t96;
t28 = pkin(9) * t245 + t32;
t6 = -qJD(5) * t257 - t196 * t28 + t199 * t21;
t282 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t281 = -mrSges(6,2) + t284;
t280 = 2 * m(5);
t279 = 2 * m(6);
t278 = 0.2e1 * m(7);
t277 = 0.2e1 * pkin(2);
t97 = qJD(3) * t157 + t178 * t197 - t200 * t216;
t276 = 0.2e1 * t97;
t193 = t195 ^ 2;
t127 = t169 * mrSges(7,1) + t168 * mrSges(7,3);
t275 = 0.2e1 * t127;
t128 = t169 * mrSges(6,1) - t168 * mrSges(6,2);
t274 = 0.2e1 * t128;
t145 = mrSges(7,1) * t173 - mrSges(7,3) * t174;
t273 = 0.2e1 * t145;
t272 = 0.2e1 * t190;
t271 = m(4) / 0.2e1;
t263 = pkin(2) * t200;
t262 = t31 * mrSges(5,3);
t261 = -pkin(9) - qJ(4);
t186 = pkin(2) * t197 + qJ(4);
t260 = -pkin(9) - t186;
t33 = mrSges(6,1) * t152 - mrSges(6,3) * t50;
t34 = -t152 * mrSges(7,1) + t50 * mrSges(7,2);
t259 = t34 - t33;
t35 = -mrSges(6,2) * t152 - mrSges(6,3) * t51;
t36 = -mrSges(7,2) * t51 + mrSges(7,3) * t152;
t258 = t35 + t36;
t255 = Ifges(5,4) * t194;
t254 = Ifges(5,4) * t195;
t253 = Ifges(5,2) * t194;
t251 = t287 * t97;
t160 = t168 * mrSges(7,2);
t250 = t194 * t31;
t248 = t195 * t32;
t102 = -mrSges(7,2) * t121 + mrSges(7,3) * t175;
t99 = -mrSges(6,2) * t175 - mrSges(6,3) * t121;
t246 = t99 + t102;
t243 = t287 * t197;
t239 = t186 * t195;
t184 = t200 * t252 + qJD(4);
t238 = t196 * t184;
t100 = mrSges(6,1) * t175 + mrSges(6,3) * t122;
t101 = -mrSges(7,1) * t175 - mrSges(7,2) * t122;
t236 = t101 - t100;
t84 = -mrSges(5,1) * t245 - mrSges(5,2) * t244;
t233 = t194 ^ 2 + t193;
t232 = qJD(3) * t197;
t229 = t196 * qJD(4);
t228 = t199 * qJD(4);
t227 = 0.2e1 * t201;
t226 = pkin(2) * t232;
t187 = -pkin(4) * t195 - pkin(3);
t18 = t51 * mrSges(6,1) + t50 * mrSges(6,2);
t17 = t51 * mrSges(7,1) - t50 * mrSges(7,3);
t224 = t261 * t194;
t223 = t260 * t194;
t191 = t195 * pkin(9);
t171 = t191 + t239;
t125 = t171 * t196 - t199 * t223;
t126 = t199 * t171 + t196 * t223;
t220 = qJD(5) * t260;
t86 = -t171 * t231 + t184 * t237 + (t199 * t220 - t238) * t194;
t87 = t171 * t230 + t195 * t238 + (t199 * t184 + t196 * t220) * t194;
t222 = t125 * t87 + t126 * t86;
t221 = qJD(5) * t261;
t181 = qJ(4) * t195 + t191;
t114 = t195 * t228 - t181 * t231 + (t199 * t221 - t229) * t194;
t115 = t195 * t229 + t181 * t230 + (t196 * t221 + t228) * t194;
t154 = t181 * t196 - t199 * t224;
t155 = t199 * t181 + t196 * t224;
t219 = t155 * t114 + t115 * t154;
t218 = t233 * t184;
t217 = t233 * qJD(4);
t118 = pkin(4) * t242 - t287;
t214 = Ifges(5,5) * t195 - Ifges(5,6) * t194;
t25 = -t196 * t73 + t199 * t63;
t211 = 0.2e1 * t233 * mrSges(5,3);
t210 = t297 * t152 + t296 * t51 + t298 * t50;
t5 = t196 * t21 + t199 * t28 + t63 * t230 - t231 * t73;
t209 = t114 * t126 + t115 * t125 + t154 * t87 + t155 * t86;
t208 = t127 + t128;
t103 = pkin(5) * t169 + qJ(6) * t168 - qJD(6) * t174;
t133 = pkin(5) * t173 - qJ(6) * t174 + t187;
t129 = -Ifges(7,5) * t168 + Ifges(7,3) * t169;
t130 = -Ifges(6,4) * t168 - Ifges(6,2) * t169;
t147 = Ifges(7,5) * t174 + Ifges(7,3) * t173;
t148 = Ifges(6,4) * t174 - Ifges(6,2) * t173;
t207 = t292 * t174 + (t129 - t130) * t173 + (t147 - t148) * t169 - t291 * t168;
t69 = -pkin(4) * t245 + t97;
t204 = (pkin(5) * t168 - qJ(6) * t169 - qJD(6) * t173) * mrSges(7,2) + t290;
t1 = qJ(6) * t152 + qJD(6) * t175 + t5;
t13 = Ifges(7,5) * t50 + Ifges(7,6) * t152 + Ifges(7,3) * t51;
t14 = Ifges(6,4) * t50 - Ifges(6,2) * t51 + Ifges(6,6) * t152;
t23 = qJ(6) * t175 + t257;
t24 = -pkin(5) * t175 - t25;
t3 = -pkin(5) * t152 - t6;
t42 = pkin(5) * t121 + qJ(6) * t122 + t118;
t61 = t152 * Ifges(5,6) - (-t253 + t254) * t151;
t62 = t152 * Ifges(5,5) - (Ifges(5,1) * t195 - t255) * t151;
t65 = -Ifges(7,5) * t122 + Ifges(7,6) * t175 + Ifges(7,3) * t121;
t66 = -Ifges(6,4) * t122 - Ifges(6,2) * t121 + Ifges(6,6) * t175;
t8 = pkin(5) * t51 - qJ(6) * t50 + qJD(6) * t122 + t69;
t202 = (Ifges(5,5) * t194 + Ifges(5,6) * t195 + t296 * t173 + t298 * t174) * t152 / 0.2e1 - t294 * t168 / 0.2e1 + t295 * t174 / 0.2e1 + t293 * t97 + t290 * t175 / 0.2e1 + t291 * t50 / 0.2e1 - t292 * t122 / 0.2e1 + (Ifges(5,2) * t195 + t255) * t245 / 0.2e1 - (Ifges(5,1) * t194 + t254) * t244 / 0.2e1 + (t13 / 0.2e1 - t14 / 0.2e1) * t173 + (t65 / 0.2e1 - t66 / 0.2e1) * t169 + (t147 / 0.2e1 - t148 / 0.2e1) * t51 + (t129 / 0.2e1 - t130 / 0.2e1) * t121 + (t168 * t25 - t169 * t257 - t173 * t5 - t174 * t6) * mrSges(6,3) + (-t1 * t173 - t169 * t23 + t174 * t3) * mrSges(7,2) + mrSges(5,3) * t248 - t24 * t160 - t96 * mrSges(4,2) + t42 * t127 + t118 * t128 + t8 * t145 + t69 * t146 - Ifges(4,5) * t151 - Ifges(4,6) * t152 + t194 * t62 / 0.2e1 + t195 * t61 / 0.2e1;
t189 = -pkin(3) - t263;
t179 = t187 - t263;
t140 = mrSges(5,1) * t175 - mrSges(5,3) * t241;
t139 = -mrSges(5,2) * t175 - mrSges(5,3) * t242;
t134 = (mrSges(5,1) * t194 + mrSges(5,2) * t195) * t176;
t120 = t133 - t263;
t98 = t103 + t226;
t91 = mrSges(5,1) * t152 + mrSges(5,3) * t244;
t90 = -mrSges(5,2) * t152 + mrSges(5,3) * t245;
t77 = mrSges(6,1) * t121 - mrSges(6,2) * t122;
t76 = mrSges(7,1) * t121 + mrSges(7,3) * t122;
t2 = [(t65 - t66) * t51 + (-0.2e1 * mrSges(4,3) * t96 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + t297) * t152 + t210) * t175 + t294 * t50 + (t296 * t152 + t13 - t14) * t121 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t201) * t227 + (m(4) * pkin(2) * t272 + (mrSges(4,1) * t175 + mrSges(4,2) * t176) * t277 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t198 + (Ifges(3,1) - Ifges(3,2)) * t227) * t198) * qJD(2) + (t118 * t69 + t25 * t6 + t257 * t5) * t279 + 0.2e1 * t257 * t35 + (t1 * t23 + t24 * t3 + t42 * t8) * t278 - (t295 + t302) * t122 + 0.2e1 * (t151 * t287 - t152 * t157) * mrSges(4,3) - 0.2e1 * t287 * t84 + (mrSges(4,1) * t152 - mrSges(4,2) * t151) * t272 + 0.2e1 * m(4) * (t157 * t96 - t251) + (t31 * t88 + t32 * t89 - t251) * t280 + (mrSges(4,3) * t276 - t194 * t61 + t195 * t62 + (-(2 * Ifges(4,4)) + t214) * t152 - (Ifges(5,1) * t193 + (2 * Ifges(4,1)) + (t253 - 0.2e1 * t254) * t194) * t151) * t176 + t134 * t276 - 0.2e1 * (-Ifges(4,4) + t214) * t151 * t175 + 0.2e1 * t25 * t33 + 0.2e1 * t24 * t34 + 0.2e1 * t23 * t36 + 0.2e1 * t42 * t17 + 0.2e1 * t8 * t76 + 0.2e1 * t69 * t77 + 0.2e1 * t89 * t90 + 0.2e1 * t88 * t91 + 0.2e1 * t5 * t99 + 0.2e1 * t6 * t100 + 0.2e1 * t3 * t101 + 0.2e1 * t1 * t102 + 0.2e1 * t118 * t18 + 0.2e1 * t32 * t139 + 0.2e1 * t31 * t140; ((t197 * t96 - t200 * t97) * t271 + ((t157 * t200 - t243) * t271 - m(5) * t243 / 0.2e1 + m(6) * t118 * t197 / 0.2e1) * qJD(3)) * t277 + t202 + t236 * t87 + t246 * t86 + (t184 * t139 + t186 * t90) * t195 + (-t184 * t140 - t186 * t91 - t262) * t194 + t258 * t126 + t259 * t125 + m(7) * (t1 * t126 + t120 * t8 + t125 * t3 + t23 * t86 + t24 * t87 + t42 * t98) + (Ifges(3,5) * t201 - Ifges(3,6) * t198 + (-mrSges(3,1) * t201 + mrSges(3,2) * t198) * pkin(7)) * qJD(2) + m(5) * (t286 * t184 - t186 * t250 + t189 * t97 + t32 * t239) + m(6) * (-t125 * t6 + t126 * t5 + t179 * t69 - t25 * t87 + t257 * t86) + ((t134 + t77) * t232 + (t151 * t200 - t152 * t197 + (-t175 * t200 + t176 * t197) * qJD(3)) * mrSges(4,3)) * pkin(2) + t98 * t76 + t120 * t17 + t179 * t18 + t189 * t84; t207 + 0.2e1 * t283 + t184 * t211 + (t186 * t218 + t189 * t226) * t280 + (t179 * t226 + t222) * t279 + (t120 * t98 + t222) * t278 + t120 * t275 + t98 * t273 + t179 * t274 + (-t125 * t168 - t126 * t169 - t86 * t173 + t87 * t174) * t299; m(5) * (-pkin(3) * t97 + t286 * qJD(4) + (t248 - t250) * qJ(4)) + t202 + (qJ(4) * t90 + qJD(4) * t139) * t195 + (-qJ(4) * t91 - qJD(4) * t140 - t262) * t194 + t258 * t155 + t259 * t154 + t236 * t115 + t246 * t114 + m(6) * (t114 * t257 - t115 * t25 - t154 * t6 + t155 * t5 + t187 * t69) + m(7) * (t1 * t155 + t103 * t42 + t114 * t23 + t115 * t24 + t133 * t8 + t154 * t3) - pkin(3) * t84 + t103 * t76 + t133 * t17 + t187 * t18; t207 + m(5) * (-pkin(3) * t226 + qJ(4) * t218 + t186 * t217) + m(6) * (t187 * t226 + t209) + m(7) * (t103 * t120 + t133 * t98 + t209) + t283 + (t103 + t98) * t145 + (t179 + t187) * t128 + (t120 + t133) * t127 + (t218 + t217) * mrSges(5,3) + t289 * ((t115 + t87) * t174 + (-t114 - t86) * t173 + (-t126 - t155) * t169 - (t125 + t154) * t168); t207 + (t103 * t133 + t219) * t278 + t219 * t279 + (qJ(4) * t233 * t280 + t211) * qJD(4) + t133 * t275 + t103 * t273 + t187 * t274 + (-t114 * t173 + t115 * t174 - t154 * t168 - t155 * t169) * t299; m(5) * t97 + m(6) * t69 + m(7) * t8 + t17 + t18 + t84; m(7) * t98 + (m(5) + m(6)) * t226 + t208; m(7) * t103 + t208; 0; -pkin(5) * t34 + m(7) * (-pkin(5) * t3 + qJ(6) * t1 + qJD(6) * t23) + qJD(6) * t102 + qJ(6) * t36 + t1 * mrSges(7,3) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t6 * mrSges(6,1) + t210; t126 * t288 + t281 * t86 + t282 * t87 + t204; t114 * t281 + t115 * t282 + t155 * t288 + t204; 0; 0.2e1 * t284 * qJD(6); m(7) * t3 + t34; m(7) * t87 - t160; m(7) * t115 - t160; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
