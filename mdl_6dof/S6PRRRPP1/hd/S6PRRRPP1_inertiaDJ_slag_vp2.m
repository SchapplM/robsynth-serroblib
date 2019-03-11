% Calculate time derivative of joint inertia matrix for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:54
% EndTime: 2019-03-08 22:43:02
% DurationCPUTime: 3.24s
% Computational Cost: add. (2869->437), mult. (7787->633), div. (0->0), fcn. (6701->10), ass. (0->191)
t237 = Ifges(7,4) + Ifges(6,5);
t236 = Ifges(7,6) - Ifges(6,6);
t154 = sin(qJ(4));
t155 = sin(qJ(3));
t157 = cos(qJ(4));
t192 = qJD(4) * t157;
t158 = cos(qJ(3));
t195 = qJD(3) * t158;
t164 = t154 * t195 + t155 * t192;
t239 = -Ifges(7,2) - Ifges(5,3) - Ifges(6,3);
t238 = mrSges(6,3) + mrSges(7,2);
t151 = sin(pkin(11));
t206 = cos(pkin(11));
t177 = t206 * t157;
t178 = t206 * t154;
t193 = qJD(4) * t155;
t183 = t154 * t193;
t184 = t157 * t195;
t232 = t183 - t184;
t51 = t151 * t232 - t177 * t193 - t178 * t195;
t114 = t151 * t157 + t178;
t205 = t151 * t154;
t167 = t177 - t205;
t52 = -t114 * t193 + t167 * t195;
t15 = -mrSges(7,1) * t51 - t52 * mrSges(7,3);
t16 = -t51 * mrSges(6,1) + mrSges(6,2) * t52;
t235 = t15 + t16;
t105 = t114 * qJD(4);
t106 = t167 * qJD(4);
t59 = mrSges(7,1) * t105 - mrSges(7,3) * t106;
t60 = mrSges(6,1) * t105 + mrSges(6,2) * t106;
t234 = t59 + t60;
t233 = m(7) * qJD(6);
t129 = -pkin(3) * t158 - pkin(9) * t155 - pkin(2);
t200 = t157 * t158;
t139 = pkin(8) * t200;
t92 = t129 * t154 + t139;
t221 = pkin(4) * t151;
t140 = qJ(6) + t221;
t231 = m(7) * t140 + mrSges(7,3);
t188 = t206 * pkin(4);
t143 = -t188 - pkin(5);
t224 = m(6) * pkin(4);
t230 = m(7) * t143 - t206 * t224 - mrSges(6,1) - mrSges(7,1);
t229 = t151 * t224 - mrSges(6,2) + t231;
t228 = 2 * m(5);
t227 = 0.2e1 * m(6);
t226 = 0.2e1 * m(7);
t225 = 0.2e1 * pkin(8);
t213 = Ifges(5,4) * t154;
t132 = Ifges(5,2) * t157 + t213;
t223 = -t132 / 0.2e1;
t222 = -t154 / 0.2e1;
t220 = pkin(8) * t154;
t218 = -qJ(5) - pkin(9);
t191 = qJD(5) * t157;
t124 = (pkin(3) * t155 - pkin(9) * t158) * qJD(3);
t196 = qJD(3) * t155;
t198 = t124 * t157 + t196 * t220;
t19 = -t155 * t191 + (pkin(4) * t155 - qJ(5) * t200) * qJD(3) + (-t139 + (qJ(5) * t155 - t129) * t154) * qJD(4) + t198;
t199 = t124 * t154 + t129 * t192;
t201 = t155 * t157;
t25 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t201 + (-qJD(5) * t155 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t158) * t154 + t199;
t7 = t151 * t19 + t206 * t25;
t33 = mrSges(7,2) * t51 + mrSges(7,3) * t196;
t34 = -mrSges(6,2) * t196 + mrSges(6,3) * t51;
t217 = t33 + t34;
t35 = mrSges(6,1) * t196 - mrSges(6,3) * t52;
t36 = -mrSges(7,1) * t196 + mrSges(7,2) * t52;
t216 = -t35 + t36;
t116 = t157 * t129;
t69 = -qJ(5) * t201 + t116 + (-pkin(4) - t220) * t158;
t202 = t154 * t155;
t81 = -qJ(5) * t202 + t92;
t28 = t151 * t69 + t206 * t81;
t96 = t114 * t155;
t84 = -mrSges(7,2) * t96 - mrSges(7,3) * t158;
t85 = mrSges(6,2) * t158 - mrSges(6,3) * t96;
t215 = t84 + t85;
t97 = t167 * t155;
t86 = -mrSges(6,1) * t158 - mrSges(6,3) * t97;
t87 = mrSges(7,1) * t158 + mrSges(7,2) * t97;
t214 = -t86 + t87;
t212 = Ifges(5,4) * t157;
t211 = t106 * mrSges(7,2);
t153 = cos(pkin(6));
t152 = sin(pkin(6));
t156 = sin(qJ(2));
t204 = t152 * t156;
t107 = -t153 * t158 + t155 * t204;
t108 = t153 * t155 + t158 * t204;
t159 = cos(qJ(2));
t197 = qJD(2) * t159;
t186 = t152 * t197;
t77 = qJD(3) * t108 + t155 * t186;
t43 = t107 * t77;
t210 = t158 * Ifges(5,6);
t209 = t77 * t155;
t78 = -qJD(3) * t107 + t158 * t186;
t208 = t78 * t158;
t207 = -mrSges(5,1) * t157 + mrSges(5,2) * t154 - mrSges(4,1);
t203 = t152 * t159;
t125 = pkin(4) * t202 + pkin(8) * t155;
t194 = qJD(4) * t154;
t190 = pkin(4) * t194;
t189 = pkin(9) * t194;
t149 = pkin(8) * t195;
t90 = pkin(4) * t164 + t149;
t144 = -pkin(4) * t157 - pkin(3);
t187 = qJD(2) * t204;
t179 = qJD(4) * t218;
t104 = t154 * t179 + t191;
t163 = -qJD(5) * t154 + t157 * t179;
t53 = t104 * t151 - t163 * t206;
t54 = t104 * t206 + t151 * t163;
t131 = t218 * t157;
t82 = -t131 * t151 - t178 * t218;
t83 = -t131 * t206 + t205 * t218;
t181 = t53 * t82 + t54 * t83;
t180 = Ifges(5,6) * t154 + (2 * Ifges(4,4));
t37 = (-t157 * t196 - t158 * t194) * pkin(8) + t199;
t91 = -t158 * t220 + t116;
t176 = -qJD(4) * t91 + t37;
t174 = mrSges(5,1) * t154 + mrSges(5,2) * t157;
t173 = Ifges(5,1) * t157 - t213;
t133 = Ifges(5,1) * t154 + t212;
t172 = -Ifges(5,2) * t154 + t212;
t171 = Ifges(5,5) * t154 + Ifges(5,6) * t157;
t168 = -t108 * t157 + t154 * t203;
t79 = -t108 * t154 - t157 * t203;
t29 = -t151 * t168 - t206 * t79;
t30 = t151 * t79 - t168 * t206;
t21 = qJD(4) * t168 - t78 * t154 + t157 * t187;
t22 = qJD(4) * t79 + t154 * t187 + t78 * t157;
t6 = t151 * t22 - t206 * t21;
t8 = t151 * t21 + t206 * t22;
t170 = t29 * t53 + t30 * t54 + t6 * t82 + t8 * t83;
t169 = t107 * t195 + t209;
t5 = -t151 * t25 + t19 * t206;
t27 = -t151 * t81 + t206 * t69;
t166 = -Ifges(5,5) * t184 + t196 * t239 + t236 * t51 - t237 * t52;
t160 = -t21 * t154 + t22 * t157 + (t154 * t168 - t157 * t79) * qJD(4);
t148 = Ifges(5,5) * t192;
t123 = -mrSges(5,1) * t158 - mrSges(5,3) * t201;
t122 = mrSges(5,2) * t158 - mrSges(5,3) * t202;
t121 = t173 * qJD(4);
t120 = t172 * qJD(4);
t119 = (mrSges(4,1) * t155 + mrSges(4,2) * t158) * qJD(3);
t118 = t174 * qJD(4);
t111 = t174 * t155;
t103 = Ifges(7,4) * t106;
t102 = Ifges(6,5) * t106;
t101 = Ifges(6,6) * t105;
t100 = Ifges(7,6) * t105;
t95 = -Ifges(5,5) * t158 + t155 * t173;
t94 = t155 * t172 - t210;
t89 = -mrSges(5,2) * t196 - mrSges(5,3) * t164;
t88 = mrSges(5,1) * t196 + mrSges(5,3) * t232;
t76 = Ifges(6,1) * t114 + Ifges(6,4) * t167;
t75 = Ifges(7,1) * t114 - Ifges(7,5) * t167;
t74 = Ifges(6,4) * t114 + Ifges(6,2) * t167;
t73 = Ifges(7,5) * t114 - Ifges(7,3) * t167;
t72 = -mrSges(6,1) * t167 + mrSges(6,2) * t114;
t71 = -mrSges(7,1) * t167 - mrSges(7,3) * t114;
t67 = mrSges(5,1) * t164 - mrSges(5,2) * t232;
t65 = -pkin(5) * t167 - qJ(6) * t114 + t144;
t64 = Ifges(6,1) * t106 - Ifges(6,4) * t105;
t63 = Ifges(7,1) * t106 + Ifges(7,5) * t105;
t62 = Ifges(6,4) * t106 - Ifges(6,2) * t105;
t61 = Ifges(7,5) * t106 + Ifges(7,3) * t105;
t58 = -t133 * t193 + (Ifges(5,5) * t155 + t158 * t173) * qJD(3);
t57 = -t132 * t193 + (Ifges(5,6) * t155 + t158 * t172) * qJD(3);
t56 = mrSges(6,1) * t96 + mrSges(6,2) * t97;
t55 = mrSges(7,1) * t96 - mrSges(7,3) * t97;
t42 = Ifges(6,1) * t97 - Ifges(6,4) * t96 - Ifges(6,5) * t158;
t41 = Ifges(7,1) * t97 - Ifges(7,4) * t158 + Ifges(7,5) * t96;
t40 = Ifges(6,4) * t97 - Ifges(6,2) * t96 - Ifges(6,6) * t158;
t39 = Ifges(7,5) * t97 - Ifges(7,6) * t158 + Ifges(7,3) * t96;
t38 = -qJD(4) * t92 + t198;
t32 = pkin(5) * t96 - qJ(6) * t97 + t125;
t31 = pkin(5) * t105 - qJ(6) * t106 - qJD(6) * t114 + t190;
t24 = pkin(5) * t158 - t27;
t23 = -qJ(6) * t158 + t28;
t14 = Ifges(6,1) * t52 + Ifges(6,4) * t51 + Ifges(6,5) * t196;
t13 = Ifges(7,1) * t52 + Ifges(7,4) * t196 - Ifges(7,5) * t51;
t12 = Ifges(6,4) * t52 + Ifges(6,2) * t51 + Ifges(6,6) * t196;
t11 = Ifges(7,5) * t52 + Ifges(7,6) * t196 - Ifges(7,3) * t51;
t9 = -pkin(5) * t51 - qJ(6) * t52 - qJD(6) * t97 + t90;
t4 = -pkin(5) * t196 - t5;
t3 = qJ(6) * t196 - qJD(6) * t158 + t7;
t1 = [0.2e1 * m(5) * (-t168 * t22 + t21 * t79 + t43) + 0.2e1 * m(4) * (-t152 ^ 2 * t156 * t197 + t108 * t78 + t43) + 0.2e1 * (m(6) + m(7)) * (t29 * t6 + t30 * t8 + t43); t22 * t122 + t21 * t123 + t79 * t88 - t168 * t89 + t215 * t8 + t214 * t6 + t217 * t30 + t216 * t29 + (t55 + t56 + t111) * t77 + (t67 + t235) * t107 + (-t159 * t119 + (-t159 * mrSges(3,2) + (-mrSges(4,1) * t158 + mrSges(4,2) * t155 - mrSges(3,1)) * t156) * qJD(2)) * t152 + (t209 + t208 + (t107 * t158 - t108 * t155) * qJD(3)) * mrSges(4,3) + m(7) * (t107 * t9 + t23 * t8 + t24 * t6 + t29 * t4 + t3 * t30 + t32 * t77) + m(6) * (t107 * t90 + t125 * t77 - t27 * t6 + t28 * t8 - t29 * t5 + t30 * t7) + m(5) * (-t168 * t37 + t21 * t91 + t22 * t92 + t38 * t79) - m(4) * pkin(2) * t187 + (m(5) * t169 / 0.2e1 + m(4) * (-t108 * t196 + t169 + t208) / 0.2e1) * t225; (t13 + t14) * t97 + (t11 - t12) * t96 + (t41 + t42) * t52 + (t40 - t39) * t51 + (t67 * t225 - t154 * t57 + t157 * t58 + (-t154 * t95 - t157 * t94 + t158 * t171) * qJD(4) + ((Ifges(5,5) * t157 - t180) * t155 + t237 * t97 + t236 * t96 + (pkin(8) ^ 2 * t228 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t239) * t158) * qJD(3)) * t155 + (t23 * t3 + t24 * t4 + t32 * t9) * t226 + (t125 * t90 + t27 * t5 + t28 * t7) * t227 + (t37 * t92 + t38 * t91) * t228 + ((t111 * t225 - t154 * t94 + t157 * t95 + t158 * t180) * qJD(3) + t166) * t158 + 0.2e1 * t32 * t15 + 0.2e1 * t23 * t33 + 0.2e1 * t28 * t34 + 0.2e1 * t27 * t35 + 0.2e1 * t24 * t36 + 0.2e1 * t9 * t55 + 0.2e1 * t3 * t84 + 0.2e1 * t7 * t85 + 0.2e1 * t5 * t86 + 0.2e1 * t4 * t87 + 0.2e1 * t90 * t56 + 0.2e1 * t91 * t88 + 0.2e1 * t92 * t89 - 0.2e1 * pkin(2) * t119 + 0.2e1 * t37 * t122 + 0.2e1 * t38 * t123 + 0.2e1 * t125 * t16; -t78 * mrSges(4,2) + (t118 + t234) * t107 + (t71 + t72 + t207) * t77 + m(6) * (t107 * t190 + t144 * t77 + t170) + m(7) * (t107 * t31 + t65 * t77 + t170) + t160 * mrSges(5,3) + t238 * (-t105 * t30 + t106 * t29 + t114 * t6 + t167 * t8) + (-pkin(3) * t77 + pkin(9) * t160) * m(5); (t75 / 0.2e1 + t76 / 0.2e1) * t52 + (-t73 / 0.2e1 + t74 / 0.2e1) * t51 + (t13 / 0.2e1 + t14 / 0.2e1) * t114 + (t41 / 0.2e1 + t42 / 0.2e1) * t106 + (t39 / 0.2e1 - t40 / 0.2e1) * t105 + m(7) * (t23 * t54 + t24 * t53 + t3 * t83 + t31 * t32 + t4 * t82 + t65 * t9) + (t63 / 0.2e1 + t64 / 0.2e1) * t97 + (t61 / 0.2e1 - t62 / 0.2e1) * t96 + m(5) * (-pkin(9) * t154 * t38 - pkin(3) * t149 - t189 * t92) + m(6) * (t125 * t190 + t144 * t90 - t27 * t53 + t28 * t54 - t5 * t82 + t7 * t83) + (t133 * t195 / 0.2e1 + qJD(4) * t95 / 0.2e1 + t57 / 0.2e1 + t176 * mrSges(5,3) + (m(5) * t176 - qJD(4) * t123 + t89) * pkin(9)) * t157 + (-t148 / 0.2e1 - t102 / 0.2e1 + t101 / 0.2e1 - t103 / 0.2e1 - t100 / 0.2e1 + (pkin(8) * t207 + Ifges(4,5)) * qJD(3)) * t158 + t214 * t53 + t215 * t54 + t216 * t82 + t217 * t83 + (t195 * t223 - t38 * mrSges(5,3) - pkin(9) * t88 + t58 / 0.2e1 + (pkin(4) * t56 - t92 * mrSges(5,3) - pkin(9) * t122 - t94 / 0.2e1 + t210 / 0.2e1) * qJD(4)) * t154 + t31 * t55 + t32 * t59 + t65 * t15 - pkin(3) * t67 + t9 * t71 + t90 * t72 + t125 * t60 + t144 * t16 - (t11 / 0.2e1 - t12 / 0.2e1) * t167 + (-t105 * t28 - t106 * t27 - t114 * t5 + t167 * t7) * mrSges(6,3) + (-t105 * t23 + t106 * t24 + t114 * t4 + t167 * t3) * mrSges(7,2) + (t120 * t222 - Ifges(4,6) * qJD(3) + t157 * t121 / 0.2e1 + (t133 * t222 + t157 * t223) * qJD(4) + (mrSges(4,2) * qJD(3) + t118) * pkin(8) + (t237 * t114 - t167 * t236 + t171) * qJD(3) / 0.2e1) * t155; -0.2e1 * pkin(3) * t118 + t157 * t120 + t154 * t121 + 0.2e1 * t144 * t60 + 0.2e1 * t31 * t71 + 0.2e1 * t65 * t59 + (t63 + t64) * t114 - (t61 - t62) * t167 + (t75 + t76) * t106 + (t73 - t74) * t105 + (t157 * t133 + (0.2e1 * pkin(4) * t72 - t132) * t154) * qJD(4) + (t31 * t65 + t181) * t226 + (t144 * t190 + t181) * t227 + 0.2e1 * t238 * (-t105 * t83 + t106 * t82 + t114 * t53 + t167 * t54); t21 * mrSges(5,1) - t22 * mrSges(5,2) + t229 * t8 + t230 * t6 + t233 * t30; -t166 + (t151 * t7 + t206 * t5) * t224 + m(7) * (qJD(6) * t23 + t140 * t3 + t143 * t4) - Ifges(5,5) * t183 + t35 * t188 + t34 * t221 + t3 * mrSges(7,3) - t4 * mrSges(7,1) + t5 * mrSges(6,1) - t7 * mrSges(6,2) - t37 * mrSges(5,2) + t38 * mrSges(5,1) + qJD(6) * t84 + t140 * t33 + t143 * t36 - t164 * Ifges(5,6); t83 * t233 - pkin(9) * mrSges(5,1) * t192 + mrSges(5,2) * t189 - Ifges(5,6) * t194 + t143 * t211 + t100 - t101 + t102 + t103 + t148 + t229 * t54 + t230 * t53 + (-t105 * t221 - t106 * t188) * mrSges(6,3) + (qJD(6) * t167 - t105 * t140) * mrSges(7,2); 0.2e1 * t231 * qJD(6); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t77; m(6) * t90 + m(7) * t9 + t235; m(6) * t190 + m(7) * t31 + t234; 0; 0; m(7) * t6; m(7) * t4 + t36; m(7) * t53 + t211; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
