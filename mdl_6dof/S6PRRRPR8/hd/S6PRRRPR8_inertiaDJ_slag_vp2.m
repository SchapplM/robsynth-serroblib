% Calculate time derivative of joint inertia matrix for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:18
% EndTime: 2019-03-08 23:48:29
% DurationCPUTime: 4.68s
% Computational Cost: add. (4474->559), mult. (12813->801), div. (0->0), fcn. (12066->12), ass. (0->225)
t283 = -Ifges(6,1) - Ifges(5,3);
t282 = Ifges(5,5) - Ifges(6,4);
t281 = Ifges(5,6) - Ifges(6,5);
t173 = sin(qJ(6));
t177 = cos(qJ(6));
t178 = cos(qJ(4));
t227 = qJD(6) * t178;
t174 = sin(qJ(4));
t229 = qJD(4) * t174;
t190 = t173 * t227 + t177 * t229;
t189 = t173 * t229 - t177 * t227;
t140 = -t178 * mrSges(5,1) + t174 * mrSges(5,2);
t280 = -m(5) * pkin(3) - mrSges(4,1) + t140;
t246 = qJ(5) * t174;
t203 = -pkin(4) * t178 - t246;
t137 = -pkin(3) + t203;
t139 = t178 * mrSges(6,2) - t174 * mrSges(6,3);
t279 = m(6) * t137 + t139;
t169 = sin(pkin(7));
t175 = sin(qJ(3));
t230 = qJD(3) * t175;
t219 = t169 * t230;
t171 = cos(pkin(7));
t245 = t169 * t175;
t114 = t171 * t174 + t178 * t245;
t179 = cos(qJ(3));
t231 = qJD(3) * t169;
t218 = t179 * t231;
t81 = qJD(4) * t114 + t174 * t218;
t61 = mrSges(6,1) * t81 - mrSges(6,3) * t219;
t64 = -mrSges(5,2) * t219 - mrSges(5,3) * t81;
t277 = -t61 + t64;
t223 = t174 * t245;
t113 = -t178 * t171 + t223;
t244 = t169 * t179;
t84 = mrSges(6,1) * t113 + mrSges(6,3) * t244;
t86 = mrSges(5,2) * t244 - mrSges(5,3) * t113;
t276 = t84 - t86;
t85 = mrSges(6,1) * t114 - mrSges(6,2) * t244;
t87 = -mrSges(5,1) * t244 - mrSges(5,3) * t114;
t275 = t85 - t87;
t141 = mrSges(7,1) * t173 + mrSges(7,2) * t177;
t274 = mrSges(6,3) + t141;
t273 = -mrSges(4,1) * t171 + mrSges(5,1) * t113 + mrSges(5,2) * t114 + mrSges(4,3) * t245;
t160 = pkin(9) * t245;
t256 = pkin(2) * t179;
t115 = t171 * t256 - t160;
t180 = cos(qJ(2));
t233 = t179 * t180;
t176 = sin(qJ(2));
t239 = t175 * t176;
t272 = t171 * t233 - t239;
t116 = t171 * t175 * pkin(2) + pkin(9) * t244;
t102 = pkin(10) * t171 + t116;
t103 = (-pkin(3) * t179 - pkin(10) * t175 - pkin(2)) * t169;
t108 = (pkin(3) * t175 - pkin(10) * t179) * t231;
t109 = t115 * qJD(3);
t228 = qJD(4) * t178;
t25 = -t102 * t228 - t103 * t229 + t108 * t178 - t174 * t109;
t20 = -pkin(4) * t219 - t25;
t80 = qJD(4) * t223 - t171 * t228 - t178 * t218;
t62 = -t80 * mrSges(6,1) + mrSges(6,2) * t219;
t271 = m(6) * t20 + t62;
t192 = t102 * t229 - t103 * t228 - t174 * t108 - t178 * t109;
t17 = -t169 * (qJ(5) * t230 - qJD(5) * t179) + t192;
t172 = cos(pkin(6));
t170 = sin(pkin(6));
t243 = t170 * t180;
t112 = -t169 * t243 + t171 * t172;
t237 = t176 * t179;
t238 = t175 * t180;
t193 = t171 * t238 + t237;
t224 = t172 * t245;
t187 = t170 * t193 + t224;
t183 = t187 * t174;
t49 = -t112 * t178 + t183;
t74 = -t170 * t272 - t172 * t244;
t22 = -t173 * t74 + t177 * t49;
t23 = t173 * t49 + t177 * t74;
t232 = qJD(2) * t176;
t221 = t170 * t232;
t210 = t169 * t221;
t47 = t172 * t218 + (t272 * qJD(3) + (-t171 * t239 + t233) * qJD(2)) * t170;
t50 = t112 * t174 + t178 * t187;
t14 = qJD(4) * t50 + t47 * t174 - t178 * t210;
t46 = t172 * t219 + (t193 * qJD(3) + (t171 * t237 + t238) * qJD(2)) * t170;
t3 = -qJD(6) * t23 + t14 * t177 - t173 * t46;
t4 = qJD(6) * t22 + t14 * t173 + t177 * t46;
t270 = qJD(6) * (t173 * t22 - t177 * t23) - t173 * t4 - t177 * t3;
t255 = mrSges(7,3) * t178;
t132 = mrSges(7,1) * t174 + t173 * t255;
t133 = -mrSges(7,2) * t174 - t177 * t255;
t262 = pkin(4) + pkin(11);
t118 = -t178 * t262 - pkin(3) - t246;
t261 = pkin(5) + pkin(10);
t150 = t261 * t174;
t72 = -t118 * t173 + t150 * t177;
t73 = t118 * t177 + t150 * t173;
t269 = m(7) * (-t173 * t72 + t177 * t73) + t177 * t133 - t173 * t132;
t268 = 0.2e1 * m(7);
t267 = -0.2e1 * mrSges(4,3);
t82 = t113 * t177 + t173 * t244;
t35 = qJD(6) * t82 + t173 * t81 + t177 * t219;
t265 = t35 / 0.2e1;
t264 = t82 / 0.2e1;
t194 = -t113 * t173 + t177 * t244;
t263 = -t194 / 0.2e1;
t204 = -Ifges(7,5) * t173 - Ifges(7,6) * t177;
t260 = t204 * qJD(6) / 0.2e1;
t251 = Ifges(7,4) * t177;
t145 = -Ifges(7,2) * t173 + t251;
t259 = t145 / 0.2e1;
t258 = -t173 / 0.2e1;
t257 = -t177 / 0.2e1;
t15 = t174 * t210 - qJD(4) * t183 + (qJD(4) * t112 + t47) * t178;
t5 = t50 * t15;
t26 = t74 * t46;
t56 = t178 * t102 + t174 * t103;
t254 = Ifges(5,4) * t174;
t253 = Ifges(5,4) * t178;
t252 = Ifges(7,4) * t173;
t250 = Ifges(6,6) * t174;
t249 = Ifges(6,6) * t178;
t110 = t116 * qJD(3);
t248 = t110 * t74;
t147 = Ifges(7,1) * t177 - t252;
t241 = t173 * t147;
t240 = t173 * t262;
t235 = t177 * t145;
t234 = t177 * t262;
t36 = qJD(6) * t194 - t173 * t219 + t177 * t81;
t6 = Ifges(7,5) * t35 + Ifges(7,6) * t36 - Ifges(7,3) * t80;
t225 = m(6) * pkin(10) + mrSges(6,1);
t151 = t261 * t178;
t220 = t169 ^ 2 * t232;
t213 = t244 / 0.2e1;
t55 = -t174 * t102 + t103 * t178;
t211 = pkin(4) * t229 - qJD(5) * t174;
t52 = pkin(4) * t244 - t55;
t30 = pkin(5) * t114 + pkin(11) * t244 + t52;
t101 = t160 + (-pkin(3) - t256) * t171;
t188 = -qJ(5) * t114 + t101;
t37 = t113 * t262 + t188;
t10 = t173 * t30 + t177 * t37;
t9 = -t173 * t37 + t177 * t30;
t208 = t10 * t177 - t173 * t9;
t207 = mrSges(7,1) * t177 - mrSges(7,2) * t173;
t206 = Ifges(7,1) * t173 + t251;
t205 = Ifges(7,2) * t177 + t252;
t202 = t14 * t174 + t15 * t178;
t135 = qJD(4) * t151;
t96 = (pkin(11) * t174 - qJ(5) * t178) * qJD(4) + t211;
t33 = qJD(6) * t72 + t135 * t173 + t177 * t96;
t34 = -qJD(6) * t73 + t135 * t177 - t173 * t96;
t200 = -t173 * t33 - t177 * t34;
t53 = -mrSges(7,2) * t114 + mrSges(7,3) * t82;
t54 = mrSges(7,1) * t114 + mrSges(7,3) * t194;
t199 = -t173 * t54 + t177 * t53;
t197 = t15 * qJ(5) + t50 * qJD(5);
t51 = qJ(5) * t244 - t56;
t196 = t283 * t219 + t281 * t81 + t282 * t80;
t28 = -Ifges(7,4) * t194 + Ifges(7,2) * t82 + Ifges(7,6) * t114;
t29 = -Ifges(7,1) * t194 + Ifges(7,4) * t82 + Ifges(7,5) * t114;
t195 = t257 * t28 + t258 * t29;
t65 = t189 * Ifges(7,5) + t190 * Ifges(7,6) + Ifges(7,3) * t228;
t184 = m(7) * t270;
t182 = qJ(5) * t80 - qJD(5) * t114 + t110;
t167 = Ifges(5,5) * t228;
t166 = Ifges(6,5) * t229;
t154 = Ifges(4,5) * t218;
t148 = Ifges(5,1) * t174 + t253;
t146 = Ifges(5,2) * t178 + t254;
t144 = Ifges(7,5) * t177 - Ifges(7,6) * t173;
t143 = -Ifges(6,2) * t174 - t249;
t142 = -Ifges(6,3) * t178 - t250;
t134 = t261 * t229;
t131 = (Ifges(5,1) * t178 - t254) * qJD(4);
t130 = t206 * qJD(6);
t129 = (-Ifges(5,2) * t174 + t253) * qJD(4);
t128 = t205 * qJD(6);
t126 = (-Ifges(6,2) * t178 + t250) * qJD(4);
t125 = (Ifges(6,3) * t174 - t249) * qJD(4);
t124 = (mrSges(5,1) * t174 + mrSges(5,2) * t178) * qJD(4);
t123 = (-mrSges(6,2) * t174 - mrSges(6,3) * t178) * qJD(4);
t122 = t207 * qJD(6);
t121 = -mrSges(4,2) * t171 + mrSges(4,3) * t244;
t117 = t207 * t178;
t111 = -qJ(5) * t228 + t211;
t107 = (mrSges(4,1) * t175 + mrSges(4,2) * t179) * t231;
t106 = Ifges(7,5) * t174 - t178 * t206;
t105 = Ifges(7,6) * t174 - t178 * t205;
t104 = Ifges(7,3) * t174 + t178 * t204;
t92 = mrSges(7,1) * t228 - mrSges(7,3) * t189;
t91 = -mrSges(7,2) * t228 + mrSges(7,3) * t190;
t71 = -mrSges(7,1) * t190 + mrSges(7,2) * t189;
t69 = -mrSges(6,2) * t113 - mrSges(6,3) * t114;
t67 = -t147 * t227 + (Ifges(7,5) * t178 + t174 * t206) * qJD(4);
t66 = -t145 * t227 + (Ifges(7,6) * t178 + t174 * t205) * qJD(4);
t63 = mrSges(5,1) * t219 + mrSges(5,3) * t80;
t60 = Ifges(5,1) * t114 - Ifges(5,4) * t113 - Ifges(5,5) * t244;
t59 = Ifges(5,4) * t114 - Ifges(5,2) * t113 - Ifges(5,6) * t244;
t58 = -Ifges(6,4) * t244 - Ifges(6,2) * t114 + Ifges(6,6) * t113;
t57 = -Ifges(6,5) * t244 - Ifges(6,6) * t114 + Ifges(6,3) * t113;
t48 = pkin(4) * t113 + t188;
t45 = -mrSges(7,1) * t82 - mrSges(7,2) * t194;
t44 = mrSges(5,1) * t81 - mrSges(5,2) * t80;
t43 = -mrSges(6,2) * t81 + mrSges(6,3) * t80;
t42 = -Ifges(5,1) * t80 - Ifges(5,4) * t81 + Ifges(5,5) * t219;
t41 = -Ifges(5,4) * t80 - Ifges(5,2) * t81 + Ifges(5,6) * t219;
t40 = Ifges(6,4) * t219 + Ifges(6,2) * t80 + Ifges(6,6) * t81;
t39 = Ifges(6,5) * t219 + Ifges(6,6) * t80 + Ifges(6,3) * t81;
t38 = -pkin(5) * t113 - t51;
t27 = -Ifges(7,5) * t194 + Ifges(7,6) * t82 + Ifges(7,3) * t114;
t21 = pkin(4) * t81 + t182;
t19 = mrSges(7,2) * t80 + mrSges(7,3) * t36;
t18 = -mrSges(7,1) * t80 - mrSges(7,3) * t35;
t16 = t262 * t81 + t182;
t13 = -pkin(5) * t81 - t17;
t12 = -pkin(5) * t80 - t219 * t262 - t25;
t11 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t8 = Ifges(7,1) * t35 + Ifges(7,4) * t36 - Ifges(7,5) * t80;
t7 = Ifges(7,4) * t35 + Ifges(7,2) * t36 - Ifges(7,6) * t80;
t2 = -qJD(6) * t10 + t12 * t177 - t16 * t173;
t1 = qJD(6) * t9 + t12 * t173 + t16 * t177;
t24 = [0.2e1 * m(7) * (t22 * t3 + t23 * t4 + t5) + 0.2e1 * m(4) * (t47 * t224 + t26 + (t112 * t169 * t232 + t193 * t47) * t170) + 0.2e1 * (m(6) + m(5)) * (t14 * t49 + t26 + t5); -qJD(2) * mrSges(3,2) * t243 - t187 * mrSges(4,3) * t219 + t47 * t121 + t112 * t107 + t4 * t53 + t3 * t54 + t22 * t18 + t23 * t19 + m(7) * (t1 * t23 + t10 * t4 + t2 * t22 + t3 * t9) + m(5) * t248 + t170 * (-mrSges(4,1) * t179 + mrSges(4,2) * t175) * t220 - mrSges(3,1) * t221 + m(4) * (t109 * t224 + t248 + t116 * t47 + (-pkin(2) * t220 + t109 * t193) * t170) + (-m(5) * t25 + t271 - t63) * t49 + (-m(5) * t55 + m(6) * t52 + t275) * t14 + (m(6) * t21 + mrSges(4,3) * t218 + t43 + t44) * t74 + (-m(5) * t192 - m(6) * t17 + m(7) * t13 + t11 + t277) * t50 + (-m(4) * t115 + m(5) * t101 + m(6) * t48 + t273 + t69) * t46 + (m(5) * t56 - m(6) * t51 + m(7) * t38 - t276 + t45) * t15; 0.2e1 * m(5) * (t101 * t110 - t192 * t56 + t25 * t55) - 0.2e1 * t192 * t86 - t194 * t8 + (t1 * t10 + t13 * t38 + t2 * t9) * t268 + 0.2e1 * t109 * t121 + 0.2e1 * t101 * t44 + t82 * t7 + 0.2e1 * t17 * t84 + 0.2e1 * t20 * t85 + 0.2e1 * t25 * t87 + 0.2e1 * t51 * t61 + 0.2e1 * t52 * t62 + 0.2e1 * t55 * t63 + 0.2e1 * t56 * t64 + 0.2e1 * t21 * t69 + 0.2e1 * t1 * t53 + 0.2e1 * t2 * t54 + 0.2e1 * t13 * t45 + 0.2e1 * t48 * t43 + t35 * t29 + t36 * t28 + 0.2e1 * t38 * t11 + 0.2e1 * t9 * t18 + 0.2e1 * t10 * t19 + (t42 + t6 - t40) * t114 + (t39 - t41) * t113 + 0.2e1 * m(4) * (t109 * t116 - t110 * t115) + 0.2e1 * m(6) * (t17 * t51 + t20 * t52 + t21 * t48) + (t57 - t59) * t81 + (t58 - t60 - t27) * t80 + t171 * t154 + 0.2e1 * t273 * t110 + (-0.2e1 * pkin(2) * t107 + t196 * t179 + ((0.2e1 * Ifges(4,4) * t244 + Ifges(4,5) * t171 + t115 * t267) * t179 + (-0.2e1 * Ifges(4,4) * t245 + t116 * t267 - 0.2e1 * Ifges(4,6) * t171 + t282 * t114 - t281 * t113 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t283) * t244) * t175) * qJD(3)) * t169; -t47 * mrSges(4,2) + t15 * t117 + t3 * t132 + t4 * t133 + t22 * t92 + t23 * t91 + t50 * t71 + m(7) * (-t134 * t50 + t15 * t151 + t22 * t34 + t23 * t33 + t3 * t72 + t4 * t73) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (t228 * t49 - t229 * t50 + t202) * pkin(10) + (m(6) * t111 + t123 + t124) * t74 + (t279 + t280) * t46 + (mrSges(5,3) + mrSges(6,1)) * ((-t174 * t50 + t178 * t49) * qJD(4) + t202); t280 * t110 + (t142 / 0.2e1 - t146 / 0.2e1) * t81 + (t143 / 0.2e1 - t148 / 0.2e1 - t104 / 0.2e1) * t80 + m(7) * (t1 * t73 + t10 * t33 + t13 * t151 - t134 * t38 + t2 * t72 + t34 * t9) + (-t39 / 0.2e1 + t41 / 0.2e1 + t7 * t257 + t8 * t258 - t192 * mrSges(5,3) - t17 * mrSges(6,1) + (t29 * t257 + t173 * t28 / 0.2e1) * qJD(6)) * t178 + (t277 * t178 + (t62 - t63) * t174 + (t174 * t276 + t178 * t275) * qJD(4) + m(6) * (-t17 * t178 + t20 * t174 + t228 * t52 + t229 * t51) + m(5) * (-t25 * t174 - t178 * t192 - t228 * t55 - t229 * t56)) * pkin(10) + (-t126 / 0.2e1 + t131 / 0.2e1 + t65 / 0.2e1) * t114 + (t125 / 0.2e1 - t129 / 0.2e1) * t113 + t67 * t263 + t66 * t264 + t106 * t265 - t134 * t45 + t137 * t43 + t21 * t139 + t151 * t11 + t48 * t123 + t101 * t124 + t2 * t132 + t1 * t133 + t13 * t117 - t109 * mrSges(4,2) + t111 * t69 + t10 * t91 + t9 * t92 + t36 * t105 / 0.2e1 + t72 * t18 + t73 * t19 + t38 * t71 + t33 * t53 + t34 * t54 - pkin(3) * t44 + t154 + (-t40 / 0.2e1 + t42 / 0.2e1 + t6 / 0.2e1 + t20 * mrSges(6,1) - t25 * mrSges(5,3)) * t174 + m(6) * (t111 * t48 + t137 * t21) + ((t52 * mrSges(6,1) - t55 * mrSges(5,3) - t58 / 0.2e1 + t60 / 0.2e1 + t27 / 0.2e1 + Ifges(6,4) * t213) * t178 + (t51 * mrSges(6,1) - t56 * mrSges(5,3) + Ifges(5,6) * t213 + t57 / 0.2e1 - t59 / 0.2e1 - t195) * t174) * qJD(4) + ((-t166 / 0.2e1 - t167 / 0.2e1) * t179 + (-Ifges(4,6) + (Ifges(5,6) / 0.2e1 - Ifges(6,5) / 0.2e1) * t178 + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t174) * t230) * t169; 0.2e1 * t137 * t123 + 0.2e1 * t151 * t71 - 0.2e1 * pkin(3) * t124 + 0.2e1 * t34 * t132 + 0.2e1 * t33 * t133 - 0.2e1 * t134 * t117 + 0.2e1 * t73 * t91 + 0.2e1 * t72 * t92 + (-t134 * t151 + t33 * t73 + t34 * t72) * t268 + 0.2e1 * t279 * t111 + (-t126 + t131 + t65 + (t105 * t177 + t106 * t173 + t142 - t146) * qJD(4)) * t174 + (-t173 * t67 - t177 * t66 - t125 + t129 + (t105 * t173 - t106 * t177) * qJD(6) + (t104 - t143 + t148) * qJD(4)) * t178; t50 * t122 + (-mrSges(5,1) + mrSges(6,2)) * t14 + (-mrSges(5,2) + t274) * t15 + m(6) * (-pkin(4) * t14 + t197) + m(7) * t197 + t262 * t184 + t270 * mrSges(7,3); (-t208 * mrSges(7,3) - (m(7) * t208 + t199) * t262 + t195) * qJD(6) - t196 + t13 * t141 - t80 * t144 / 0.2e1 + t36 * t259 + t147 * t265 + t38 * t122 + t114 * t260 - t128 * t264 - t130 * t263 - pkin(4) * t62 + t192 * mrSges(5,2) + t25 * mrSges(5,1) - t17 * mrSges(6,3) + t20 * mrSges(6,2) + m(6) * (-pkin(4) * t20 - qJ(5) * t17 - qJD(5) * t51) + m(7) * (qJ(5) * t13 + qJD(5) * t38 - t1 * t240 - t2 * t234) + (-t61 + t11) * qJ(5) + (-t84 + t45) * qJD(5) + (-t7 / 0.2e1 - t262 * t19 - t1 * mrSges(7,3)) * t173 + (t8 / 0.2e1 - t262 * t18 - t2 * mrSges(7,3)) * t177; t166 + t167 + t177 * t67 / 0.2e1 + t66 * t258 + t174 * t260 - t134 * t141 + t151 * t122 + qJD(5) * t117 + qJ(5) * t71 + m(7) * (-qJ(5) * t134 + qJD(5) * t151 - t234 * t34 - t240 * t33) - t92 * t234 - t91 * t240 + t200 * mrSges(7,3) + (t225 * qJD(5) - t128 * t257 - t130 * t258) * t178 + ((-t73 * mrSges(7,3) - t178 * t147 / 0.2e1 - t105 / 0.2e1) * t177 + (t178 * t259 + t72 * mrSges(7,3) - t106 / 0.2e1) * t173 - t269 * t262) * qJD(6) + ((-pkin(4) * mrSges(6,1) + t144 / 0.2e1 - Ifges(6,4)) * t178 + (t235 / 0.2e1 + t241 / 0.2e1 - qJ(5) * mrSges(6,1) - Ifges(5,6)) * t174 + (m(6) * t203 + t139 + t140) * pkin(10)) * qJD(4); 0.2e1 * qJ(5) * t122 + t128 * t173 - t130 * t177 + (-t235 - t241) * qJD(6) + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t274) * qJD(5); m(6) * t14 - t184; t173 * t19 + t177 * t18 + t199 * qJD(6) + m(7) * (qJD(6) * t208 + t1 * t173 + t177 * t2) + t271; -m(7) * t200 + qJD(6) * t269 + t173 * t91 + t177 * t92 + t225 * t228; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t34 - mrSges(7,2) * t33 + t65; ((mrSges(7,2) * t262 - Ifges(7,6)) * t177 + (mrSges(7,1) * t262 - Ifges(7,5)) * t173) * qJD(6); -t141 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t24(1) t24(2) t24(4) t24(7) t24(11) t24(16); t24(2) t24(3) t24(5) t24(8) t24(12) t24(17); t24(4) t24(5) t24(6) t24(9) t24(13) t24(18); t24(7) t24(8) t24(9) t24(10) t24(14) t24(19); t24(11) t24(12) t24(13) t24(14) t24(15) t24(20); t24(16) t24(17) t24(18) t24(19) t24(20) t24(21);];
Mq  = res;
