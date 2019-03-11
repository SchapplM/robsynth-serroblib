% Calculate time derivative of joint inertia matrix for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:18
% EndTime: 2019-03-08 22:30:24
% DurationCPUTime: 2.62s
% Computational Cost: add. (3000->381), mult. (7091->573), div. (0->0), fcn. (6109->10), ass. (0->175)
t204 = pkin(4) + pkin(8);
t127 = cos(qJ(5));
t123 = sin(qJ(5));
t128 = cos(qJ(3));
t177 = qJD(5) * t128;
t167 = t123 * t177;
t124 = sin(qJ(3));
t181 = qJD(3) * t124;
t138 = t127 * t181 + t167;
t99 = mrSges(6,1) * t123 + mrSges(6,2) * t127;
t218 = mrSges(5,3) + t99;
t160 = pkin(3) * t181 - qJD(4) * t124;
t63 = (pkin(9) * t124 - qJ(4) * t128) * qJD(3) + t160;
t180 = qJD(3) * t128;
t93 = t204 * t180;
t161 = -t123 * t63 + t127 * t93;
t130 = -pkin(3) - pkin(9);
t190 = qJ(4) * t124;
t76 = t128 * t130 - pkin(2) - t190;
t162 = pkin(10) * t128 - t76;
t102 = t204 * t124;
t83 = t123 * t102;
t10 = (-pkin(10) * t123 * t124 + pkin(5) * t128) * qJD(3) + (t127 * t162 - t83) * qJD(5) + t161;
t178 = qJD(5) * t127;
t179 = qJD(5) * t123;
t14 = t102 * t178 + t123 * t93 + t127 * t63 - t179 * t76;
t11 = pkin(10) * t138 + t14;
t122 = sin(qJ(6));
t126 = cos(qJ(6));
t84 = t127 * t102;
t35 = pkin(5) * t124 + t123 * t162 + t84;
t184 = t127 * t128;
t45 = t127 * t76 + t83;
t39 = -pkin(10) * t184 + t45;
t12 = -t122 * t39 + t126 * t35;
t2 = qJD(6) * t12 + t10 * t122 + t11 * t126;
t13 = t122 * t35 + t126 * t39;
t3 = -qJD(6) * t13 + t10 * t126 - t11 * t122;
t217 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t98 = t128 * mrSges(5,2) - t124 * mrSges(5,3);
t215 = -t128 * mrSges(4,1) + t124 * mrSges(4,2) + t98;
t15 = -qJD(5) * t45 + t161;
t214 = t123 * t14 + t127 * t15;
t213 = qJD(5) + qJD(6);
t120 = sin(pkin(6));
t129 = cos(qJ(2));
t188 = t120 * t129;
t121 = cos(pkin(6));
t125 = sin(qJ(2));
t189 = t120 * t125;
t170 = t124 * t189;
t71 = -t121 * t128 + t170;
t141 = -t123 * t71 + t127 * t188;
t183 = qJD(2) * t125;
t169 = t120 * t183;
t182 = qJD(2) * t129;
t168 = t120 * t182;
t72 = t121 * t124 + t128 * t189;
t51 = qJD(3) * t72 + t124 * t168;
t16 = qJD(5) * t141 - t123 * t169 + t127 * t51;
t53 = t123 * t188 + t127 * t71;
t17 = qJD(5) * t53 + t123 * t51 + t127 * t169;
t212 = qJD(5) * (t123 * t53 + t127 * t141) - t123 * t17 - t127 * t16;
t186 = t126 * t127;
t142 = t122 * t123 - t186;
t143 = t122 * t127 + t126 * t123;
t176 = qJD(6) * t122;
t46 = -t122 * t179 - t123 * t176 + t186 * t213;
t47 = t213 * t143;
t211 = qJD(6) * (-t122 * t142 - t126 * t143) - t122 * t46 + t126 * t47;
t210 = t213 * t128;
t209 = 2 * m(7);
t208 = m(5) / 0.2e1;
t207 = m(7) * pkin(5);
t206 = -t143 / 0.2e1;
t205 = -t142 / 0.2e1;
t202 = -t123 / 0.2e1;
t201 = t124 / 0.2e1;
t200 = -t127 / 0.2e1;
t52 = -qJD(3) * t170 + (qJD(3) * t121 + t168) * t128;
t32 = t72 * t52;
t199 = pkin(10) - t130;
t198 = -Ifges(7,5) * t47 - Ifges(7,6) * t46;
t197 = Ifges(6,4) * t123;
t196 = Ifges(6,4) * t127;
t195 = Ifges(6,5) * t124;
t156 = Ifges(6,1) * t123 + t196;
t65 = -t128 * t156 + t195;
t193 = t123 * t65;
t155 = Ifges(6,2) * t127 + t197;
t64 = t124 * Ifges(6,6) - t128 * t155;
t191 = t127 * t64;
t101 = Ifges(6,1) * t127 - t197;
t187 = t123 * t101;
t100 = -Ifges(6,2) * t123 + t196;
t185 = t127 * t100;
t103 = t204 * t128;
t175 = qJD(6) * t126;
t174 = 2 * mrSges(7,3);
t30 = t142 * t210 + t143 * t181;
t31 = -t142 * t181 + t143 * t210;
t173 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t180;
t172 = m(5) * pkin(8) + mrSges(5,1);
t23 = t122 * t141 + t126 * t53;
t5 = qJD(6) * t23 + t122 * t16 + t126 * t17;
t24 = t122 * t53 - t126 * t141;
t6 = -qJD(6) * t24 - t122 * t17 + t126 * t16;
t171 = t6 * mrSges(7,1) - t5 * mrSges(7,2);
t166 = t127 * t177;
t165 = t123 * t181;
t163 = -t47 * mrSges(7,1) - t46 * mrSges(7,2);
t95 = t199 * t127;
t159 = -t142 * t47 - t143 * t46;
t157 = mrSges(6,1) * t127 - mrSges(6,2) * t123;
t154 = -Ifges(6,5) * t123 - Ifges(6,6) * t127;
t153 = -pkin(3) * t128 - t190;
t94 = t199 * t123;
t56 = -t122 * t95 - t126 * t94;
t55 = t122 * t94 - t126 * t95;
t44 = -t123 * t76 + t84;
t149 = t123 * t44 - t127 * t45;
t90 = mrSges(6,3) * t123 * t128 + mrSges(6,1) * t124;
t91 = -mrSges(6,2) * t124 - mrSges(6,3) * t184;
t147 = -t123 * t90 + t127 * t91;
t146 = t51 * t124 + t52 * t128;
t79 = t199 * t179;
t80 = qJD(5) * t95;
t21 = qJD(6) * t55 + t122 * t79 - t126 * t80;
t22 = -qJD(6) * t56 + t122 * t80 + t126 * t79;
t145 = t22 * mrSges(7,1) - t21 * mrSges(7,2) + t198;
t144 = t52 * qJ(4) + t72 * qJD(4);
t140 = Ifges(6,5) * t165 + t138 * Ifges(6,6) + Ifges(6,3) * t180 + t173;
t137 = t165 - t166;
t136 = -t47 * t12 + t46 * t13 - t142 * t3 + t143 * t2;
t135 = -t142 * t6 + t143 * t5 - t47 * t23 + t46 * t24;
t134 = t142 * t22 - t143 * t21 - t46 * t56 + t47 * t55;
t131 = m(6) * t212;
t111 = pkin(5) * t123 + qJ(4);
t107 = pkin(5) * t178 + qJD(4);
t96 = -pkin(2) + t153;
t92 = t204 * t181;
t89 = t156 * qJD(5);
t88 = t155 * qJD(5);
t87 = (mrSges(4,1) * t124 + mrSges(4,2) * t128) * qJD(3);
t86 = (-mrSges(5,2) * t124 - mrSges(5,3) * t128) * qJD(3);
t85 = t157 * qJD(5);
t77 = (-mrSges(7,1) * t122 - mrSges(7,2) * t126) * qJD(6) * pkin(5);
t75 = t157 * t128;
t70 = pkin(5) * t184 + t103;
t69 = -qJ(4) * t180 + t160;
t67 = t143 * t128;
t66 = t142 * t128;
t62 = mrSges(6,1) * t180 - mrSges(6,3) * t137;
t61 = -mrSges(6,2) * t180 + mrSges(6,3) * t138;
t60 = mrSges(7,1) * t124 + t67 * mrSges(7,3);
t59 = -mrSges(7,2) * t124 + t66 * mrSges(7,3);
t57 = -pkin(5) * t167 + (-pkin(5) * t127 - t204) * t181;
t50 = -Ifges(7,1) * t142 - Ifges(7,4) * t143;
t49 = -Ifges(7,4) * t142 - Ifges(7,2) * t143;
t48 = mrSges(7,1) * t143 - mrSges(7,2) * t142;
t40 = -mrSges(6,1) * t138 + mrSges(6,2) * t137;
t38 = -mrSges(7,1) * t66 - mrSges(7,2) * t67;
t37 = -t101 * t177 + (Ifges(6,5) * t128 + t124 * t156) * qJD(3);
t36 = -t100 * t177 + (Ifges(6,6) * t128 + t124 * t155) * qJD(3);
t34 = -Ifges(7,1) * t67 + Ifges(7,4) * t66 + Ifges(7,5) * t124;
t33 = -Ifges(7,4) * t67 + Ifges(7,2) * t66 + Ifges(7,6) * t124;
t27 = -Ifges(7,1) * t47 - Ifges(7,4) * t46;
t26 = -Ifges(7,4) * t47 - Ifges(7,2) * t46;
t25 = mrSges(7,1) * t46 - mrSges(7,2) * t47;
t19 = -mrSges(7,2) * t180 + t31 * mrSges(7,3);
t18 = mrSges(7,1) * t180 - t30 * mrSges(7,3);
t9 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t180;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t31 + Ifges(7,6) * t180;
t1 = [0.2e1 * m(7) * (t23 * t6 + t24 * t5 + t32) + 0.2e1 * m(6) * (-t141 * t17 + t16 * t53 + t32) + 0.2e1 * (m(5) + m(4)) * (-t120 ^ 2 * t125 * t182 + t71 * t51 + t32); t16 * t90 + t17 * t91 + t23 * t18 + t24 * t19 + t5 * t59 + t53 * t62 - t141 * t61 + t6 * t60 + (t40 + t9) * t72 + (t38 + t75) * t52 + m(7) * (t12 * t6 + t13 * t5 + t2 * t24 + t23 * t3 + t52 * t70 + t57 * t72) + m(6) * (t103 * t52 - t14 * t141 + t15 * t53 + t44 * t16 + t45 * t17 - t72 * t92) + 0.2e1 * (m(4) / 0.2e1 + t208) * (t180 * t71 - t181 * t72 + t146) * pkin(8) + (mrSges(4,3) + mrSges(5,1)) * ((-t124 * t72 + t128 * t71) * qJD(3) + t146) + ((-t86 - t87) * t129 + (-t129 * mrSges(3,2) + (-mrSges(3,1) + t215) * t125) * qJD(2) - m(4) * pkin(2) * t183 + 0.2e1 * (-t129 * t69 + t183 * t96) * t208) * t120; 0.2e1 * t12 * t18 + 0.2e1 * t13 * t19 + t31 * t33 + t30 * t34 + 0.2e1 * t57 * t38 + 0.2e1 * t2 * t59 + 0.2e1 * t3 * t60 + 0.2e1 * t45 * t61 + 0.2e1 * t44 * t62 + t66 * t7 - t67 * t8 + 0.2e1 * t70 * t9 - 0.2e1 * pkin(2) * t87 + 0.2e1 * t15 * t90 + 0.2e1 * t14 * t91 - 0.2e1 * t92 * t75 + 0.2e1 * t96 * t86 + 0.2e1 * t103 * t40 + 0.2e1 * (m(5) * t96 + t98) * t69 + (t12 * t3 + t13 * t2 + t57 * t70) * t209 + 0.2e1 * m(6) * (-t103 * t92 + t45 * t14 + t44 * t15) + ((t193 + t191 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t124) * qJD(3) + t140) * t124 + (-t123 * t37 - t127 * t36 + (t123 * t64 + (-t65 - t195) * t127) * qJD(5) + (-Ifges(7,5) * t67 + Ifges(7,6) * t66 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t154) * t128 + (-(2 * Ifges(4,2)) - (2 * Ifges(5,3)) + (2 * Ifges(4,1)) + (2 * Ifges(5,2)) + Ifges(7,3) + Ifges(6,3)) * t124) * qJD(3)) * t128; (t25 + t85) * t72 + (-mrSges(4,1) + mrSges(5,2)) * t51 + (-mrSges(4,2) + t48 + t218) * t52 + m(7) * (t107 * t72 + t111 * t52 + t21 * t24 + t22 * t23 + t56 * t5 + t55 * t6) + m(6) * t144 + m(5) * (-pkin(3) * t51 + t144) - t130 * t131 - t135 * mrSges(7,3) + t212 * mrSges(6,3); (-t15 * mrSges(6,3) + t130 * t62 + t37 / 0.2e1) * t127 + (-t14 * mrSges(6,3) + t130 * t61 - t36 / 0.2e1) * t123 + t7 * t206 + t198 * t201 + t8 * t205 + (qJD(4) * t172 - t200 * t88 - t202 * t89) * t128 + m(7) * (t107 * t70 + t111 * t57 + t22 * t12 + t21 * t13 + t56 * t2 + t55 * t3) + ((-pkin(3) * mrSges(5,1) + Ifges(7,5) * t205 + Ifges(7,6) * t206 + Ifges(6,5) * t127 / 0.2e1 + Ifges(6,6) * t202 + Ifges(4,5) - Ifges(5,4)) * t128 + (t185 / 0.2e1 - qJ(4) * mrSges(5,1) + t187 / 0.2e1 - Ifges(4,6) + Ifges(5,5)) * t124 + (m(5) * t153 + t215) * pkin(8)) * qJD(3) + m(6) * (-qJ(4) * t92 + qJD(4) * t103 + t130 * t214) + (-t193 / 0.2e1 - t191 / 0.2e1 + t154 * t201 + (t101 * t200 + t123 * t100 / 0.2e1) * t128 + t149 * mrSges(6,3) + (-m(6) * t149 + t147) * t130) * qJD(5) - t136 * mrSges(7,3) + qJ(4) * t40 - t46 * t33 / 0.2e1 - t47 * t34 / 0.2e1 + t31 * t49 / 0.2e1 + t30 * t50 / 0.2e1 + t55 * t18 + t56 * t19 + t57 * t48 + t21 * t59 + t22 * t60 + t66 * t26 / 0.2e1 - t67 * t27 / 0.2e1 + t70 * t25 + qJD(4) * t75 - t92 * t99 + t103 * t85 + t107 * t38 + t111 * t9; 0.2e1 * t107 * t48 + 0.2e1 * t111 * t25 - t47 * t50 - t142 * t27 - t46 * t49 - t143 * t26 + (t107 * t111 + t56 * t21 + t55 * t22) * t209 + 0.2e1 * qJ(4) * t85 + t123 * t88 - t127 * t89 + (-t185 - t187) * qJD(5) + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t218) * qJD(4) + t134 * t174; m(5) * t51 + m(7) * t135 - t131; t123 * t61 + t127 * t62 - t142 * t18 + t143 * t19 + t46 * t59 - t47 * t60 + t147 * qJD(5) + t172 * t180 + m(7) * t136 + m(6) * (-qJD(5) * t149 + t214); -m(7) * t134 + t159 * t174; -0.2e1 * m(7) * t159; t16 * mrSges(6,1) - t17 * mrSges(6,2) + (t122 * t5 + t126 * t6 + (-t122 * t23 + t126 * t24) * qJD(6)) * t207 + t171; -Ifges(6,5) * t166 + t15 * mrSges(6,1) - t14 * mrSges(6,2) + (-t60 * t176 + t126 * t18 + m(7) * (-t12 * t176 + t122 * t2 + t126 * t3 + t13 * t175) + t59 * t175 + t122 * t19) * pkin(5) + t140 + t217; ((-mrSges(6,2) * t130 - Ifges(6,6)) * t127 + (-mrSges(6,1) * t130 - Ifges(6,5)) * t123) * qJD(5) + (m(7) * (t122 * t21 + t126 * t22 + (-t122 * t55 + t126 * t56) * qJD(6)) + t211 * mrSges(7,3)) * pkin(5) + t145; -t99 * qJD(5) - t207 * t211 + t163; 0.2e1 * t77; t171; t173 + t217; t145; t163; t77; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
