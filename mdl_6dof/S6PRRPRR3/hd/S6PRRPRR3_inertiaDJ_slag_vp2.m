% Calculate time derivative of joint inertia matrix for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:24
% EndTime: 2019-03-08 22:02:33
% DurationCPUTime: 4.33s
% Computational Cost: add. (6678->522), mult. (19791->790), div. (0->0), fcn. (19942->14), ass. (0->238)
t163 = sin(pkin(7));
t174 = cos(qJ(3));
t232 = t163 * t174;
t286 = 0.2e1 * t232;
t168 = sin(qJ(6));
t172 = cos(qJ(6));
t169 = sin(qJ(5));
t215 = qJD(6) * t169;
t173 = cos(qJ(5));
t217 = qJD(5) * t173;
t282 = t168 * t215 - t172 * t217;
t164 = sin(pkin(6));
t166 = cos(pkin(7));
t167 = cos(pkin(6));
t175 = cos(qJ(2));
t118 = -t163 * t164 * t175 + t167 * t166;
t162 = sin(pkin(13));
t165 = cos(pkin(13));
t225 = t174 * t175;
t170 = sin(qJ(3));
t171 = sin(qJ(2));
t229 = t170 * t171;
t182 = t166 * t225 - t229;
t87 = t164 * t182 + t167 * t232;
t227 = t171 * t174;
t228 = t170 * t175;
t181 = t166 * t228 + t227;
t233 = t163 * t170;
t88 = t164 * t181 + t167 * t233;
t52 = t162 * t87 + t165 * t88;
t185 = t118 * t173 - t169 * t52;
t237 = qJD(5) * t185;
t285 = t168 / 0.2e1;
t263 = t172 / 0.2e1;
t260 = pkin(2) * t166;
t152 = t170 * t260;
t120 = pkin(9) * t232 + t152;
t106 = qJ(4) * t232 + t120;
t153 = t174 * t260;
t254 = -pkin(9) - qJ(4);
t200 = t254 * t170;
t95 = pkin(3) * t166 + t163 * t200 + t153;
t70 = t165 * t106 + t162 * t95;
t56 = pkin(10) * t166 + t70;
t231 = t165 * t174;
t110 = t162 * t233 - t163 * t231;
t111 = (t162 * t174 + t165 * t170) * t163;
t126 = (-pkin(3) * t174 - pkin(2)) * t163;
t75 = pkin(4) * t110 - pkin(10) * t111 + t126;
t252 = t169 * t75 + t173 * t56;
t281 = qJD(5) * t252;
t283 = (t232 * t254 - t152) * qJD(3) - qJD(4) * t233;
t146 = qJD(3) * t153;
t84 = t146 + (qJD(3) * t200 + qJD(4) * t174) * t163;
t50 = t283 * t162 + t165 * t84;
t107 = qJD(3) * t111;
t222 = qJD(3) * t163;
t108 = (-t162 * t170 + t231) * t222;
t221 = qJD(3) * t170;
t208 = t163 * t221;
t196 = pkin(3) * t208;
t72 = pkin(4) * t107 - pkin(10) * t108 + t196;
t14 = -t169 * t50 + t173 * t72 - t281;
t284 = m(6) * (t14 + t281);
t224 = t168 ^ 2 + t172 ^ 2;
t214 = qJD(6) * t172;
t180 = -t168 * t217 - t169 * t214;
t137 = -mrSges(7,1) * t172 + mrSges(7,2) * t168;
t280 = -m(7) * pkin(5) - mrSges(6,1) + t137;
t279 = 0.2e1 * m(7);
t278 = -2 * mrSges(4,3);
t277 = -2 * mrSges(5,3);
t276 = -2 * Ifges(5,4);
t155 = pkin(3) * t162 + pkin(10);
t275 = 0.2e1 * t155;
t274 = 0.2e1 * t173;
t273 = m(6) / 0.2e1;
t272 = m(7) / 0.2e1;
t92 = t111 * t169 - t173 * t166;
t67 = -qJD(5) * t92 + t108 * t173;
t93 = t111 * t173 + t166 * t169;
t74 = t110 * t168 + t172 * t93;
t28 = -qJD(6) * t74 + t107 * t172 - t168 * t67;
t271 = t28 / 0.2e1;
t247 = Ifges(7,4) * t168;
t140 = Ifges(7,2) * t172 + t247;
t246 = Ifges(7,4) * t172;
t189 = -Ifges(7,2) * t168 + t246;
t81 = -t140 * t215 + (Ifges(7,6) * t169 + t173 * t189) * qJD(5);
t270 = t81 / 0.2e1;
t142 = Ifges(7,1) * t168 + t246;
t190 = Ifges(7,1) * t172 - t247;
t82 = -t142 * t215 + (Ifges(7,5) * t169 + t173 * t190) * qJD(5);
t269 = t82 / 0.2e1;
t114 = -Ifges(7,5) * t173 + t169 * t190;
t268 = t114 / 0.2e1;
t267 = Ifges(7,5) * t285 + Ifges(7,6) * t263;
t266 = t142 / 0.2e1;
t265 = -t168 / 0.2e1;
t264 = -t172 / 0.2e1;
t46 = mrSges(6,1) * t107 - mrSges(6,3) * t67;
t73 = t110 * t172 - t168 * t93;
t27 = qJD(6) * t73 + t107 * t168 + t172 * t67;
t8 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t262 = t8 - t46;
t156 = -pkin(3) * t165 - pkin(4);
t261 = m(6) * t156;
t259 = pkin(5) * t169;
t258 = pkin(11) * t173;
t223 = qJD(2) * t164;
t209 = t171 * t223;
t195 = t163 * t209;
t41 = t118 * t169 + t173 * t52;
t236 = qJD(5) * t41;
t65 = -t167 * t208 + (-t181 * qJD(3) + (-t166 * t227 - t228) * qJD(2)) * t164;
t207 = t174 * t222;
t66 = t167 * t207 + (t182 * qJD(3) + (-t166 * t229 + t225) * qJD(2)) * t164;
t34 = t162 * t65 + t165 * t66;
t16 = t169 * t34 - t173 * t195 + t236;
t257 = t16 * t185;
t33 = t162 * t66 - t165 * t65;
t51 = t162 * t88 - t165 * t87;
t256 = t33 * t51;
t49 = t162 * t84 - t165 * t283;
t255 = t49 * t51;
t39 = -mrSges(7,1) * t73 + mrSges(7,2) * t74;
t77 = mrSges(6,1) * t110 - mrSges(6,3) * t93;
t253 = t39 - t77;
t251 = -mrSges(5,1) * t166 + mrSges(6,1) * t92 + mrSges(6,2) * t93 + mrSges(5,3) * t111;
t250 = mrSges(7,3) * t169;
t249 = Ifges(6,4) * t169;
t248 = Ifges(6,4) * t173;
t245 = Ifges(7,6) * t168;
t244 = t107 * Ifges(6,5);
t243 = t107 * Ifges(6,6);
t242 = t110 * Ifges(6,6);
t15 = t169 * t195 + t173 * t34 + t237;
t241 = t15 * t173;
t240 = t16 * t169;
t239 = -mrSges(6,1) * t173 + mrSges(6,2) * t169 - mrSges(5,1);
t36 = -t169 * t56 + t173 * t75;
t29 = -pkin(5) * t110 - t36;
t238 = qJD(5) * t29;
t122 = -pkin(5) * t173 - pkin(11) * t169 + t156;
t230 = t168 * t173;
t90 = t122 * t172 - t155 * t230;
t235 = qJD(6) * t90;
t226 = t172 * t173;
t91 = t122 * t168 + t155 * t226;
t234 = qJD(6) * t91;
t220 = qJD(5) * t168;
t219 = qJD(5) * t169;
t218 = qJD(5) * t172;
t216 = qJD(6) * t168;
t213 = m(7) * t274;
t68 = qJD(5) * t93 + t108 * t169;
t5 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t68;
t212 = Ifges(6,5) * t67 - Ifges(6,6) * t68 + Ifges(6,3) * t107;
t10 = -pkin(5) * t107 - t14;
t211 = -m(7) * t10 - t8;
t210 = Ifges(4,5) * t207 + Ifges(5,5) * t108 - Ifges(5,6) * t107;
t205 = t155 * t219;
t69 = -t162 * t106 + t165 * t95;
t78 = t107 * mrSges(5,1) + t108 * mrSges(5,2);
t136 = (-t258 + t259) * qJD(5);
t58 = t136 * t168 - t172 * t205 + t235;
t198 = t58 - t235;
t59 = t136 * t172 + t168 * t205 - t234;
t197 = -t59 - t234;
t30 = pkin(11) * t110 + t252;
t55 = -pkin(4) * t166 - t69;
t35 = pkin(5) * t92 - pkin(11) * t93 + t55;
t11 = -t168 * t30 + t172 * t35;
t19 = pkin(5) * t68 - pkin(11) * t67 + t49;
t13 = t169 * t72 + t173 * t50 + t75 * t217 - t219 * t56;
t9 = pkin(11) * t107 + t13;
t1 = qJD(6) * t11 + t168 * t19 + t172 * t9;
t12 = t168 * t35 + t172 * t30;
t2 = -qJD(6) * t12 - t168 * t9 + t172 * t19;
t194 = t1 * t172 - t168 * t2;
t20 = -t168 * t41 + t172 * t51;
t3 = qJD(6) * t20 + t15 * t172 + t168 * t33;
t21 = t168 * t51 + t172 * t41;
t4 = -qJD(6) * t21 - t15 * t168 + t172 * t33;
t193 = -t168 * t4 + t172 * t3;
t192 = t169 * mrSges(6,1) + t173 * mrSges(6,2);
t191 = mrSges(7,1) * t168 + mrSges(7,2) * t172;
t17 = mrSges(7,1) * t68 - mrSges(7,3) * t27;
t18 = -mrSges(7,2) * t68 + mrSges(7,3) * t28;
t188 = -t168 * t17 + t172 * t18;
t23 = Ifges(7,4) * t74 + Ifges(7,2) * t73 + Ifges(7,6) * t92;
t24 = Ifges(7,1) * t74 + Ifges(7,4) * t73 + Ifges(7,5) * t92;
t184 = t23 * t265 + t24 * t263;
t183 = -t185 * t217 + t240;
t47 = -mrSges(6,2) * t107 - mrSges(6,3) * t68;
t178 = t47 + m(6) * (-qJD(5) * t36 + t13) + t253 * qJD(5);
t176 = -t11 * t214 - t12 * t216 + t194;
t80 = -t282 * Ifges(7,5) + Ifges(7,6) * t180 + Ifges(7,3) * t219;
t159 = Ifges(6,5) * t217;
t158 = Ifges(7,5) * t214;
t143 = Ifges(6,1) * t169 + t248;
t141 = Ifges(6,2) * t173 + t249;
t135 = -mrSges(7,1) * t173 - t172 * t250;
t134 = mrSges(7,2) * t173 - t168 * t250;
t133 = (Ifges(6,1) * t173 - t249) * qJD(5);
t132 = t190 * qJD(6);
t131 = (-Ifges(6,2) * t169 + t248) * qJD(5);
t130 = t189 * qJD(6);
t129 = -Ifges(7,6) * t216 + t158;
t128 = t192 * qJD(5);
t127 = t191 * qJD(6);
t125 = -mrSges(4,2) * t166 + mrSges(4,3) * t232;
t124 = mrSges(4,1) * t166 - mrSges(4,3) * t233;
t121 = t191 * t169;
t119 = -pkin(9) * t233 + t153;
t117 = t120 * qJD(3);
t116 = -pkin(9) * t208 + t146;
t115 = (mrSges(4,1) * t170 + mrSges(4,2) * t174) * t222;
t113 = -Ifges(7,6) * t173 + t169 * t189;
t112 = -Ifges(7,3) * t173 + (Ifges(7,5) * t172 - t245) * t169;
t105 = -mrSges(7,2) * t219 + mrSges(7,3) * t180;
t104 = mrSges(7,1) * t219 + t282 * mrSges(7,3);
t96 = -mrSges(5,2) * t166 - mrSges(5,3) * t110;
t94 = t118 * t195;
t86 = t180 * mrSges(7,1) + t282 * mrSges(7,2);
t79 = mrSges(5,1) * t110 + mrSges(5,2) * t111;
t76 = -mrSges(6,2) * t110 - mrSges(6,3) * t92;
t45 = Ifges(6,1) * t93 - Ifges(6,4) * t92 + Ifges(6,5) * t110;
t44 = Ifges(6,4) * t93 - Ifges(6,2) * t92 + t242;
t43 = mrSges(7,1) * t92 - mrSges(7,3) * t74;
t42 = -mrSges(7,2) * t92 + mrSges(7,3) * t73;
t38 = mrSges(6,1) * t68 + mrSges(6,2) * t67;
t32 = Ifges(6,1) * t67 - Ifges(6,4) * t68 + t244;
t31 = Ifges(6,4) * t67 - Ifges(6,2) * t68 + t243;
t22 = Ifges(7,5) * t74 + Ifges(7,6) * t73 + Ifges(7,3) * t92;
t7 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t68;
t6 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t68;
t25 = [0.2e1 * m(7) * (t20 * t4 + t21 * t3 - t257) + 0.2e1 * m(6) * (t15 * t41 + t256 - t257) + 0.2e1 * m(5) * (t34 * t52 + t256 + t94) + 0.2e1 * m(4) * (t65 * t87 + t66 * t88 + t94); t65 * t124 + t66 * t125 + t15 * t76 + t20 * t17 + t21 * t18 + t3 * t42 + t34 * t96 + t51 * t38 + t4 * t43 + t41 * t47 - t262 * t185 + t251 * t33 + t253 * t16 + (t78 + t115) * t118 + (-mrSges(3,1) * t171 - mrSges(3,2) * t175) * t223 + (-t107 * t52 + t108 * t51) * mrSges(5,3) + ((t79 + (-mrSges(4,1) * t174 + mrSges(4,2) * t170) * t163) * t209 + (-t170 * t88 - t174 * t87) * qJD(3) * mrSges(4,3)) * t163 + m(5) * (-t33 * t69 + t34 * t70 + t255 + t50 * t52 + (pkin(3) * t118 * t221 + t126 * t209) * t163) + m(4) * (-pkin(2) * t163 ^ 2 * t209 + t116 * t88 - t117 * t87 + t119 * t65 + t120 * t66) + m(6) * (t13 * t41 + t14 * t185 + t15 * t252 - t16 * t36 + t33 * t55 + t255) + m(7) * (t1 * t21 - t10 * t185 + t11 * t4 + t12 * t3 + t16 * t29 + t2 * t20); 0.2e1 * m(6) * (t13 * t252 + t14 * t36 + t49 * t55) + 0.2e1 * t252 * t47 + 0.2e1 * m(4) * (t116 * t120 - t117 * t119) + 0.2e1 * t251 * t49 + (t22 - t44) * t68 + (t70 * t277 + t111 * t276 + Ifges(6,5) * t93 - Ifges(5,6) * t166 - Ifges(6,6) * t92 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t110) * t107 + (0.2e1 * Ifges(5,1) * t111 + Ifges(5,5) * t166 + t110 * t276 + t277 * t69) * t108 + t110 * t212 + t210 * t166 + (t5 - t31) * t92 + (-0.2e1 * pkin(2) * t115 + ((Ifges(4,4) * t286 + Ifges(4,5) * t166 + t119 * t278) * t174 + (-0.2e1 * Ifges(4,4) * t233 - 0.2e1 * Ifges(4,6) * t166 + t120 * t278 + (Ifges(4,1) - Ifges(4,2)) * t286 + 0.2e1 * (m(5) * t126 + t79) * pkin(3)) * t170) * qJD(3)) * t163 + (t1 * t12 + t10 * t29 + t11 * t2) * t279 + 0.2e1 * t11 * t17 + 0.2e1 * t12 * t18 + t27 * t24 + t28 * t23 + 0.2e1 * t29 * t8 + 0.2e1 * t10 * t39 + 0.2e1 * t1 * t42 + 0.2e1 * t2 * t43 + 0.2e1 * t36 * t46 + 0.2e1 * t55 * t38 + t67 * t45 + t73 * t6 + t74 * t7 + 0.2e1 * t13 * t76 + 0.2e1 * t14 * t77 + 0.2e1 * m(5) * (-t49 * t69 + t50 * t70) + t93 * t32 + 0.2e1 * t50 * t96 - 0.2e1 * t117 * t124 + 0.2e1 * t116 * t125 + 0.2e1 * t126 * t78; t65 * mrSges(4,1) - t66 * mrSges(4,2) - t34 * mrSges(5,2) + t20 * t104 + t21 * t105 + t16 * t121 + t51 * t128 + t3 * t134 + t4 * t135 + t185 * t86 + t239 * t33 + m(7) * (t20 * t59 + t21 * t58 + t3 * t91 + t4 * t90) + t33 * t261 + (t183 * t272 + (-t219 * t41 + t183 + t241) * t273) * t275 + m(5) * (t162 * t34 - t165 * t33) * pkin(3) + (t241 + t240 + (-t169 * t41 - t173 * t185) * qJD(5)) * mrSges(6,3); (t7 * t263 + t6 * t265 - t14 * mrSges(6,3) + t32 / 0.2e1 + t244 / 0.2e1 + (t23 * t264 + t24 * t265) * qJD(6) + (-t252 * mrSges(6,3) - t44 / 0.2e1 + t22 / 0.2e1 - t242 / 0.2e1) * qJD(5) + (-qJD(5) * t76 - t211 - t284 - t46) * t155) * t169 + m(7) * (t1 * t91 + t11 * t59 + t12 * t58 + t2 * t90) + (t80 / 0.2e1 - t131 / 0.2e1) * t92 + (t112 / 0.2e1 - t141 / 0.2e1) * t68 + t210 + t110 * t159 / 0.2e1 + (t239 + t261) * t49 + (t13 * mrSges(6,3) - t5 / 0.2e1 + t31 / 0.2e1 + t243 / 0.2e1 + (t45 / 0.2e1 - t36 * mrSges(6,3) + t184) * qJD(5) + (m(7) * t238 + t178) * t155) * t173 - Ifges(4,6) * t208 + t27 * t268 + t74 * t269 + t73 * t270 + t113 * t271 + (m(5) * (t162 * t50 - t165 * t49) + (-t107 * t162 - t108 * t165) * mrSges(5,3)) * pkin(3) - t50 * mrSges(5,2) + t58 * t42 + t59 * t43 - t29 * t86 + t90 * t17 + t91 * t18 + t11 * t104 + t12 * t105 - t116 * mrSges(4,2) - t117 * mrSges(4,1) + t10 * t121 + t55 * t128 + t93 * t133 / 0.2e1 + t1 * t134 + t2 * t135 + t67 * t143 / 0.2e1 + t156 * t38; 0.2e1 * t58 * t134 + 0.2e1 * t91 * t105 + 0.2e1 * t59 * t135 + 0.2e1 * t90 * t104 + (t58 * t91 + t59 * t90) * t279 + 0.2e1 * t156 * t128 + (t131 - t80 + (-t113 * t168 + t114 * t172 + t121 * t275 + t143) * qJD(5)) * t173 + (-0.2e1 * t155 * t86 - t168 * t81 + t172 * t82 + t133 + (-t113 * t172 - t114 * t168) * qJD(6) + (t155 ^ 2 * t213 + t112 - t141) * qJD(5)) * t169; m(5) * t195 + ((-t20 * t220 + t21 * t218 - t16) * t272 + (-t16 + t236) * t273) * t274 + 0.2e1 * ((-t20 * t214 - t21 * t216 + t193 - t237) * t272 + (t15 - t237) * t273) * t169; m(5) * t196 + ((-t168 * t43 + t172 * t42 + t76) * qJD(5) + m(7) * (-t11 * t220 + t12 * t218 - t10) + t284 - t262) * t173 + ((-t168 * t42 - t172 * t43) * qJD(6) + m(7) * (t176 + t238) + t178 + t188) * t169 + t78; t173 * t86 + (m(7) * (-t168 * t59 + t172 * t58 - t214 * t90 - t216 * t91) - t134 * t216 + t172 * t105 - t135 * t214 - t168 * t104) * t169 + (m(7) * (t226 * t91 - t230 * t90 + (t169 ^ 2 - t173 ^ 2) * t155) + t134 * t226 - t135 * t230 + t169 * t121) * qJD(5); (-0.1e1 + t224) * t213 * t219; -t15 * mrSges(6,2) - t185 * t127 + (m(7) * pkin(11) + mrSges(7,3)) * ((-t168 * t21 - t172 * t20) * qJD(6) + t193) + t280 * t16; -t13 * mrSges(6,2) + t14 * mrSges(6,1) + t29 * t127 + t92 * t129 / 0.2e1 + t73 * t130 / 0.2e1 + t74 * t132 / 0.2e1 + t10 * t137 + t68 * t267 + t140 * t271 + t27 * t266 + t7 * t285 + t6 * t263 + t184 * qJD(6) + t211 * pkin(5) + ((-t11 * t172 - t12 * t168) * qJD(6) + t194) * mrSges(7,3) + (m(7) * t176 - t214 * t43 - t216 * t42 + t188) * pkin(11) + t212; pkin(5) * t86 + t159 + (-t129 / 0.2e1 + t280 * t155 * qJD(5)) * t173 + (qJD(6) * t268 + t217 * t266 + t270 + t198 * mrSges(7,3) + (m(7) * t198 - qJD(6) * t135 + t105) * pkin(11)) * t172 + (-qJD(6) * t113 / 0.2e1 - t140 * t217 / 0.2e1 + t269 + t197 * mrSges(7,3) + (m(7) * t197 - qJD(6) * t134 - t104) * pkin(11)) * t168 + (t132 * t263 + t130 * t265 + t155 * t127 + (t140 * t264 + t142 * t265) * qJD(6) + (t155 * mrSges(6,2) - Ifges(6,6) + t267) * qJD(5)) * t169; -t173 * t127 + (t169 * t137 + m(7) * (t224 * t258 - t259) + t224 * t173 * mrSges(7,3) - t192) * qJD(5); -0.2e1 * pkin(5) * t127 + t130 * t172 + t132 * t168 + (-t140 * t168 + t142 * t172) * qJD(6); mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t59 - mrSges(7,2) * t58 + t80; t86; t158 + (pkin(11) * t137 - t245) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
