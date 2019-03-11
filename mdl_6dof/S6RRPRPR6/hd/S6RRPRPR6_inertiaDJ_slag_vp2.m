% Calculate time derivative of joint inertia matrix for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:14
% EndTime: 2019-03-09 10:38:26
% DurationCPUTime: 4.81s
% Computational Cost: add. (6496->529), mult. (17911->748), div. (0->0), fcn. (17286->10), ass. (0->227)
t170 = sin(pkin(6));
t278 = 0.2e1 * t170;
t277 = Ifges(6,1) + Ifges(5,3);
t276 = -Ifges(6,4) + Ifges(5,5);
t275 = Ifges(6,5) - Ifges(5,6);
t173 = sin(qJ(6));
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t214 = qJD(6) * t177;
t176 = cos(qJ(6));
t218 = qJD(4) * t176;
t185 = t173 * t214 + t174 * t218;
t219 = qJD(4) * t174;
t274 = t173 * t219 - t176 * t214;
t254 = -t173 / 0.2e1;
t142 = mrSges(7,1) * t173 + mrSges(7,2) * t176;
t273 = mrSges(6,3) + t142;
t175 = sin(qJ(2));
t172 = cos(pkin(6));
t252 = pkin(1) * t172;
t158 = t175 * t252;
t178 = cos(qJ(2));
t228 = t170 * t178;
t229 = t170 * t175;
t245 = -pkin(8) - qJ(3);
t272 = (t245 * t228 - t158) * qJD(2) - qJD(3) * t229;
t169 = sin(pkin(11));
t171 = cos(pkin(11));
t112 = (t169 * t178 + t171 * t175) * t170;
t109 = qJD(2) * t112;
t221 = qJD(2) * t170;
t227 = t171 * t178;
t110 = (-t169 * t175 + t227) * t221;
t94 = t112 * t177 + t172 * t174;
t71 = qJD(4) * t94 + t110 * t174;
t111 = t169 * t229 - t170 * t227;
t93 = t112 * t174 - t172 * t177;
t78 = t111 * t176 + t173 * t93;
t27 = -qJD(6) * t78 - t109 * t173 + t176 * t71;
t72 = -t112 * t219 + (qJD(4) * t172 + t110) * t177;
t17 = -mrSges(7,2) * t72 + mrSges(7,3) * t27;
t77 = -t111 * t173 + t176 * t93;
t28 = qJD(6) * t77 + t109 * t176 + t173 * t71;
t18 = mrSges(7,1) * t72 - mrSges(7,3) * t28;
t43 = -mrSges(7,2) * t94 + mrSges(7,3) * t77;
t44 = mrSges(7,1) * t94 - mrSges(7,3) * t78;
t190 = t173 * t44 - t176 * t43;
t270 = -t190 * qJD(6) + t173 * t17 + t176 * t18;
t217 = qJD(4) * t177;
t159 = t178 * t252;
t152 = qJD(2) * t159;
t203 = t245 * t175;
t90 = t152 + (qJD(2) * t203 + qJD(3) * t178) * t170;
t56 = t169 * t272 + t171 * t90;
t123 = pkin(8) * t228 + t158;
t108 = qJ(3) * t228 + t123;
t95 = pkin(2) * t172 + t170 * t203 + t159;
t74 = t171 * t108 + t169 * t95;
t62 = pkin(9) * t172 + t74;
t211 = t175 * t221;
t199 = pkin(2) * t211;
t76 = pkin(3) * t109 - pkin(9) * t110 + t199;
t128 = (-pkin(2) * t178 - pkin(1)) * t170;
t79 = pkin(3) * t111 - pkin(9) * t112 + t128;
t15 = -t174 * t56 + t177 * t76 - t62 * t217 - t79 * t219;
t12 = -pkin(4) * t109 - t15;
t37 = t174 * t79 + t177 * t62;
t29 = -qJ(5) * t111 - t37;
t50 = mrSges(5,1) * t109 - mrSges(5,3) * t72;
t52 = t72 * mrSges(6,1) + t109 * mrSges(6,2);
t269 = m(5) * (t37 * qJD(4) + t15) - m(6) * (t29 * qJD(4) + t12) - t52 + t50;
t268 = 2 * m(7);
t267 = -2 * mrSges(3,3);
t266 = -2 * mrSges(4,3);
t265 = -2 * Ifges(4,4);
t264 = 0.2e1 * t112;
t263 = 0.2e1 * t128;
t262 = m(4) * pkin(2);
t261 = t28 / 0.2e1;
t260 = t77 / 0.2e1;
t238 = Ifges(7,4) * t173;
t148 = Ifges(7,1) * t176 - t238;
t237 = Ifges(7,4) * t176;
t195 = Ifges(7,1) * t173 + t237;
t88 = -t148 * t214 + (Ifges(7,5) * t177 + t174 * t195) * qJD(4);
t259 = t88 / 0.2e1;
t258 = pkin(4) + pkin(10);
t193 = -Ifges(7,5) * t173 - Ifges(7,6) * t176;
t257 = t193 * qJD(6) / 0.2e1;
t256 = Ifges(7,5) * t176 / 0.2e1 + Ifges(7,6) * t254;
t146 = -Ifges(7,2) * t173 + t237;
t255 = t146 / 0.2e1;
t253 = -t176 / 0.2e1;
t251 = pkin(4) * t177;
t55 = t169 * t90 - t171 * t272;
t250 = t55 * mrSges(4,1);
t249 = t55 * mrSges(5,1);
t248 = t55 * mrSges(5,2);
t247 = t56 * mrSges(4,2);
t161 = pkin(2) * t169 + pkin(9);
t246 = pkin(5) + t161;
t80 = mrSges(6,1) * t93 - mrSges(6,3) * t111;
t82 = -mrSges(5,2) * t111 - mrSges(5,3) * t93;
t243 = t80 - t82;
t81 = mrSges(6,1) * t94 + mrSges(6,2) * t111;
t83 = mrSges(5,1) * t111 - mrSges(5,3) * t94;
t242 = t81 - t83;
t241 = mrSges(7,3) * t177;
t240 = Ifges(5,4) * t174;
t239 = Ifges(5,4) * t177;
t236 = Ifges(6,6) * t174;
t235 = Ifges(6,6) * t177;
t234 = t111 * Ifges(6,4);
t233 = t111 * Ifges(5,6);
t117 = -pkin(8) * t211 + t152;
t232 = t117 * mrSges(3,2);
t118 = t123 * qJD(2);
t231 = t118 * mrSges(3,1);
t230 = qJ(5) * t174;
t226 = t173 * t148;
t225 = t173 * t258;
t224 = t176 * t146;
t223 = t176 * t258;
t222 = t173 ^ 2 + t176 ^ 2;
t127 = t246 * t177;
t120 = qJD(4) * t127;
t220 = qJD(4) * t173;
t216 = qJD(6) * t173;
t215 = qJD(6) * t176;
t213 = t174 * qJD(5);
t7 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t72;
t212 = Ifges(3,5) * t178 * t221 + Ifges(4,5) * t110 - Ifges(4,6) * t109;
t162 = -pkin(2) * t171 - pkin(3);
t205 = m(6) * t161 + mrSges(6,1);
t204 = m(7) * t222;
t73 = -t169 * t108 + t171 * t95;
t36 = -t174 * t62 + t177 * t79;
t200 = pkin(4) * t219 - t213;
t183 = -qJ(5) * t72 - qJD(5) * t94 + t55;
t13 = t258 * t71 + t183;
t3 = pkin(5) * t72 - t258 * t109 - t15;
t19 = pkin(5) * t94 - t258 * t111 - t36;
t61 = -pkin(3) * t172 - t73;
t184 = -qJ(5) * t94 + t61;
t21 = t258 * t93 + t184;
t5 = -t173 * t21 + t176 * t19;
t1 = qJD(6) * t5 + t13 * t176 + t173 * t3;
t6 = t173 * t19 + t176 * t21;
t2 = -qJD(6) * t6 - t13 * t173 + t176 * t3;
t198 = t1 * t173 + t176 * t2;
t197 = t173 * t5 - t176 * t6;
t196 = mrSges(7,1) * t176 - mrSges(7,2) * t173;
t141 = t177 * mrSges(6,2) - t174 * mrSges(6,3);
t194 = Ifges(7,2) * t176 + t238;
t107 = (pkin(10) * t174 - qJ(5) * t177) * qJD(4) + t200;
t187 = t162 - t230;
t113 = -t258 * t177 + t187;
t126 = t246 * t174;
t84 = -t113 * t173 + t126 * t176;
t41 = qJD(6) * t84 + t107 * t176 + t120 * t173;
t85 = t113 * t176 + t126 * t173;
t42 = -qJD(6) * t85 - t107 * t173 + t120 * t176;
t191 = t173 * t41 + t176 * t42;
t189 = -t173 * t84 + t176 * t85;
t188 = t109 * t277 + t275 * t71 + t276 * t72;
t23 = Ifges(7,4) * t78 + Ifges(7,2) * t77 + Ifges(7,6) * t94;
t24 = Ifges(7,1) * t78 + Ifges(7,4) * t77 + Ifges(7,5) * t94;
t186 = t23 * t253 + t24 * t254;
t14 = t174 * t76 + t177 * t56 + t79 * t217 - t62 * t219;
t86 = Ifges(7,5) * t274 + Ifges(7,6) * t185 + Ifges(7,3) * t217;
t10 = -qJ(5) * t109 - qJD(5) * t111 - t14;
t30 = -pkin(4) * t111 - t36;
t49 = -mrSges(5,2) * t109 - mrSges(5,3) * t71;
t51 = mrSges(6,1) * t71 - mrSges(6,3) * t109;
t181 = t49 + m(5) * (-t36 * qJD(4) + t14) + m(6) * (t30 * qJD(4) - t10) - t51;
t105 = -mrSges(7,2) * t217 + mrSges(7,3) * t185;
t106 = mrSges(7,1) * t217 - mrSges(7,3) * t274;
t139 = mrSges(7,1) * t174 + t173 * t241;
t140 = -mrSges(7,2) * t174 - t176 * t241;
t180 = t173 * t105 + t176 * t106 - t139 * t216 + t140 * t215;
t165 = Ifges(5,5) * t217;
t164 = Ifges(6,5) * t219;
t149 = Ifges(5,1) * t174 + t239;
t147 = Ifges(5,2) * t177 + t240;
t144 = -Ifges(6,2) * t174 - t235;
t143 = -Ifges(6,3) * t177 - t236;
t138 = (Ifges(5,1) * t177 - t240) * qJD(4);
t137 = t195 * qJD(6);
t136 = (-Ifges(5,2) * t174 + t239) * qJD(4);
t135 = t194 * qJD(6);
t133 = (-Ifges(6,2) * t177 + t236) * qJD(4);
t132 = (Ifges(6,3) * t174 - t235) * qJD(4);
t131 = (t174 * mrSges(5,1) + t177 * mrSges(5,2)) * qJD(4);
t130 = (-t174 * mrSges(6,2) - t177 * mrSges(6,3)) * qJD(4);
t129 = t196 * qJD(6);
t125 = t196 * t177;
t124 = t187 - t251;
t122 = -pkin(8) * t229 + t159;
t121 = qJ(5) * t217 - t200;
t119 = t246 * t219;
t116 = Ifges(7,5) * t174 - t177 * t195;
t115 = Ifges(7,6) * t174 - t177 * t194;
t114 = Ifges(7,3) * t174 + t177 * t193;
t100 = t110 * mrSges(4,2);
t91 = mrSges(7,1) * t185 - mrSges(7,2) * t274;
t87 = -t146 * t214 + (Ifges(7,6) * t177 + t174 * t194) * qJD(4);
t60 = -mrSges(6,2) * t93 - mrSges(6,3) * t94;
t48 = Ifges(5,1) * t94 - Ifges(5,4) * t93 + Ifges(5,5) * t111;
t47 = Ifges(5,4) * t94 - Ifges(5,2) * t93 + t233;
t46 = -Ifges(6,2) * t94 + Ifges(6,6) * t93 + t234;
t45 = Ifges(6,5) * t111 - Ifges(6,6) * t94 + Ifges(6,3) * t93;
t40 = -mrSges(7,1) * t77 + mrSges(7,2) * t78;
t39 = -mrSges(6,2) * t71 - mrSges(6,3) * t72;
t38 = mrSges(5,1) * t71 + mrSges(5,2) * t72;
t35 = pkin(4) * t93 + t184;
t34 = Ifges(5,1) * t72 - Ifges(5,4) * t71 + t109 * Ifges(5,5);
t33 = Ifges(5,4) * t72 - Ifges(5,2) * t71 + t109 * Ifges(5,6);
t32 = t109 * Ifges(6,4) - Ifges(6,2) * t72 + Ifges(6,6) * t71;
t31 = t109 * Ifges(6,5) - Ifges(6,6) * t72 + Ifges(6,3) * t71;
t22 = Ifges(7,5) * t78 + Ifges(7,6) * t77 + Ifges(7,3) * t94;
t20 = -pkin(5) * t93 - t29;
t16 = pkin(4) * t71 + t183;
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t28;
t9 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t72;
t8 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t72;
t4 = -pkin(5) * t71 - t10;
t25 = [(t212 - 0.2e1 * t231 - 0.2e1 * t232 - 0.2e1 * t247 - 0.2e1 * t250) * t172 + 0.2e1 * m(3) * (t117 * t123 - t118 * t122) + 0.2e1 * m(5) * (t14 * t37 + t15 * t36 + t55 * t61) + 0.2e1 * m(6) * (t10 * t29 + t12 * t30 + t16 * t35) + t55 * mrSges(4,3) * t264 + (mrSges(4,1) * t263 + t74 * t266 + t112 * t265 - Ifges(4,6) * t172 + t276 * t94 + t275 * t93 + ((2 * Ifges(4,2)) + t277) * t111) * t109 + (0.2e1 * (t117 * t178 + t118 * t175) * mrSges(3,3) + ((t122 * t267 + Ifges(3,5) * t172 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t178) * t278) * t178 + (t262 * t263 + 0.2e1 * pkin(2) * (mrSges(4,1) * t111 + mrSges(4,2) * t112) + t123 * t267 - 0.2e1 * Ifges(3,6) * t172 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t175 + (Ifges(3,1) - Ifges(3,2)) * t178) * t278) * t175) * qJD(2)) * t170 + (t22 + t48 - t46) * t72 + (t45 - t47) * t71 + 0.2e1 * m(4) * (-t55 * t73 + t56 * t74) + (-t32 + t34 + t7 + 0.2e1 * t248) * t94 + (t31 - t33 + 0.2e1 * t249) * t93 + t100 * t263 + (t1 * t6 + t2 * t5 + t20 * t4) * t268 + (Ifges(4,1) * t264 + Ifges(4,5) * t172 + t111 * t265 + t73 * t266) * t110 + (t56 * t266 + t188) * t111 + 0.2e1 * t6 * t17 + 0.2e1 * t5 * t18 + 0.2e1 * t20 * t11 + t27 * t23 + t28 * t24 + 0.2e1 * t35 * t39 + 0.2e1 * t4 * t40 + 0.2e1 * t1 * t43 + 0.2e1 * t2 * t44 + 0.2e1 * t37 * t49 + 0.2e1 * t36 * t50 + 0.2e1 * t29 * t51 + 0.2e1 * t30 * t52 + 0.2e1 * t16 * t60 + 0.2e1 * t61 * t38 + t77 * t8 + t78 * t9 + 0.2e1 * t10 * t80 + 0.2e1 * t12 * t81 + 0.2e1 * t14 * t82 + 0.2e1 * t15 * t83; -Ifges(3,6) * t211 + m(7) * (t1 * t85 - t119 * t20 + t127 * t4 + t2 * t84 + t41 * t6 + t42 * t5) + m(6) * (-t121 * t35 + t124 * t16) - t247 + t212 - t250 + (t86 / 0.2e1 - t133 / 0.2e1 + t138 / 0.2e1) * t94 + (t132 / 0.2e1 - t136 / 0.2e1) * t93 + (t114 / 0.2e1 - t144 / 0.2e1 + t149 / 0.2e1) * t72 + (t143 / 0.2e1 - t147 / 0.2e1) * t71 - t232 - t231 + (t165 / 0.2e1 + t164 / 0.2e1) * t111 + (t12 * mrSges(6,1) - t15 * mrSges(5,3) + t7 / 0.2e1 - t32 / 0.2e1 + t34 / 0.2e1 + t248 + (-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t109 + (t45 / 0.2e1 - t47 / 0.2e1 - t37 * mrSges(5,3) + t29 * mrSges(6,1) - t233 / 0.2e1 - t186) * qJD(4) + (t243 * qJD(4) - t269) * t161) * t174 + (-t109 * t169 - t110 * t171) * pkin(2) * mrSges(4,3) + (m(5) * t55 + t38) * t162 + t78 * t259 + t87 * t260 + t116 * t261 + (t169 * t56 - t171 * t55) * t262 + (t8 * t253 + t9 * t254 - t10 * mrSges(6,1) + t14 * mrSges(5,3) - t31 / 0.2e1 + t33 / 0.2e1 - t249 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t109 + (t173 * t23 / 0.2e1 + t24 * t253) * qJD(6) + (-t46 / 0.2e1 + t48 / 0.2e1 + t22 / 0.2e1 - t36 * mrSges(5,3) + t30 * mrSges(6,1) - t234 / 0.2e1) * qJD(4) + (t242 * qJD(4) + t181) * t161) * t177 + t41 * t43 + t42 * t44 + t84 * t18 + t85 * t17 - t20 * t91 + t6 * t105 + t5 * t106 + t27 * t115 / 0.2e1 - t119 * t40 - t121 * t60 + t124 * t39 + t4 * t125 + t127 * t11 + t35 * t130 + t61 * t131 + t2 * t139 + t1 * t140 + t16 * t141; (-t119 * t127 + t41 * t85 + t42 * t84) * t268 + 0.2e1 * t85 * t105 + 0.2e1 * t84 * t106 - 0.2e1 * t119 * t125 - 0.2e1 * t127 * t91 + 0.2e1 * t124 * t130 + 0.2e1 * t42 * t139 + 0.2e1 * t41 * t140 + 0.2e1 * t162 * t131 + 0.2e1 * (-m(6) * t124 - t141) * t121 + (-t133 + t138 + t86 + (t176 * t115 + t173 * t116 + t143 - t147) * qJD(4)) * t174 + (-t173 * t88 - t176 * t87 - t132 + t136 + (t115 * t173 - t116 * t176) * qJD(6) + (t114 - t144 + t149) * qJD(4)) * t177; m(4) * t199 + t109 * mrSges(4,1) + t100 + (t11 + (t173 * t43 + t176 * t44 + t242) * qJD(4) + m(7) * (t5 * t218 + t6 * t220 + t4) + t181) * t174 + ((t40 - t243) * qJD(4) + m(7) * (qJD(4) * t20 - t6 * t215 + t5 * t216 - t198) + t269 - t270) * t177; (m(7) * (t84 * t218 + t85 * t220 - t119) + t140 * t220 + t139 * t218 - t91) * t174 + (m(7) * (-t85 * t215 + t84 * t216 + t120 - t191) + qJD(4) * t125 - t180) * t177; (0.1e1 - t222) * t174 * t217 * t268; t188 + (t197 * mrSges(7,3) - (-m(7) * t197 - t190) * t258 + t186) * qJD(6) + (-t258 * t18 - t2 * mrSges(7,3) + t9 / 0.2e1) * t176 + (-t258 * t17 - t1 * mrSges(7,3) - t8 / 0.2e1) * t173 + (t40 - t80) * qJD(5) + (t11 - t51) * qJ(5) + m(7) * (qJ(5) * t4 + qJD(5) * t20 - t1 * t225 - t2 * t223) + m(6) * (-pkin(4) * t12 - qJ(5) * t10 - qJD(5) * t29) - t10 * mrSges(6,3) + t12 * mrSges(6,2) - t14 * mrSges(5,2) + t15 * mrSges(5,1) - pkin(4) * t52 + t20 * t129 + t94 * t257 - t135 * t260 - t78 * t137 / 0.2e1 + t4 * t142 + t72 * t256 + t27 * t255 + t148 * t261; m(7) * (-qJ(5) * t119 + qJD(5) * t127 - t42 * t223 - t41 * t225) - t106 * t223 - t105 * t225 + t164 + t165 - qJ(5) * t91 + qJD(5) * t125 + t127 * t129 - t119 * t142 + t87 * t254 + t174 * t257 + t176 * t259 - t191 * mrSges(7,3) + (t205 * qJD(5) - t135 * t253 - t137 * t254) * t177 + ((-t115 / 0.2e1 - t85 * mrSges(7,3) - t177 * t148 / 0.2e1) * t176 + (-t116 / 0.2e1 + t84 * mrSges(7,3) + t177 * t255) * t173 - (m(7) * t189 - t173 * t139 + t176 * t140) * t258) * qJD(6) + ((-pkin(4) * mrSges(6,1) - Ifges(6,4) + t256) * t177 + (-Ifges(5,6) - qJ(5) * mrSges(6,1) + t226 / 0.2e1 + t224 / 0.2e1) * t174 + (m(6) * (-t230 - t251) + t174 * mrSges(5,2) - t177 * mrSges(5,1) + t141) * t161) * qJD(4); t174 * t129 + m(7) * t213 + m(6) * t121 + ((m(7) * qJ(5) - mrSges(5,2) + t273) * t177 + (-t222 * mrSges(7,3) - t204 * t258 - mrSges(5,1) + mrSges(6,2)) * t174) * qJD(4); 0.2e1 * qJ(5) * t129 + t135 * t173 - t137 * t176 + (-t224 - t226) * qJD(6) + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t273) * qJD(5); m(7) * (-qJD(6) * t197 + t198) + m(6) * t12 + t52 + t270; m(7) * (qJD(6) * t189 + t191) + t205 * t217 + t180; (m(6) + t204) * t219; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t42 - mrSges(7,2) * t41 + t86; t91; ((mrSges(7,2) * t258 - Ifges(7,6)) * t176 + (mrSges(7,1) * t258 - Ifges(7,5)) * t173) * qJD(6); -t142 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
