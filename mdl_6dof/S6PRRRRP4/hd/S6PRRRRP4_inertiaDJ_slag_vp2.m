% Calculate time derivative of joint inertia matrix for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:12
% EndTime: 2019-03-09 00:13:19
% DurationCPUTime: 3.72s
% Computational Cost: add. (4139->479), mult. (10774->687), div. (0->0), fcn. (9513->10), ass. (0->196)
t245 = -mrSges(6,1) - mrSges(7,1);
t244 = -mrSges(6,2) + mrSges(7,3);
t241 = Ifges(7,4) + Ifges(6,5);
t243 = -Ifges(7,2) - Ifges(6,3);
t240 = Ifges(7,6) - Ifges(6,6);
t159 = sin(qJ(3));
t158 = sin(qJ(4));
t163 = cos(qJ(3));
t207 = qJD(3) * t163;
t190 = t158 * t207;
t162 = cos(qJ(4));
t204 = qJD(4) * t162;
t171 = t159 * t204 + t190;
t242 = mrSges(6,3) + mrSges(7,2);
t157 = sin(qJ(5));
t161 = cos(qJ(5));
t118 = t157 * t162 + t158 * t161;
t106 = t118 * t159;
t189 = t162 * t207;
t208 = qJD(3) * t159;
t239 = -Ifges(5,5) * t189 - Ifges(5,3) * t208;
t129 = -pkin(3) * t163 - pkin(9) * t159 - pkin(2);
t211 = t162 * t163;
t142 = pkin(8) * t211;
t103 = t158 * t129 + t142;
t205 = qJD(4) * t159;
t192 = t158 * t205;
t172 = t189 - t192;
t238 = qJD(4) + qJD(5);
t116 = t162 * t129;
t212 = t159 * t162;
t229 = pkin(8) * t158;
t74 = -pkin(10) * t212 + t116 + (-pkin(4) - t229) * t163;
t213 = t158 * t159;
t92 = -pkin(10) * t213 + t103;
t225 = t157 * t74 + t161 * t92;
t127 = (pkin(3) * t159 - pkin(9) * t163) * qJD(3);
t210 = t162 * t127 + t208 * t229;
t29 = (pkin(4) * t159 - pkin(10) * t211) * qJD(3) + (-t142 + (pkin(10) * t159 - t129) * t158) * qJD(4) + t210;
t206 = qJD(4) * t158;
t58 = t158 * t127 + t129 * t204 + (-t162 * t208 - t163 * t206) * pkin(8);
t45 = -t171 * pkin(10) + t58;
t12 = -qJD(5) * t225 - t157 * t45 + t161 * t29;
t237 = 2 * m(5);
t236 = 2 * m(6);
t235 = 2 * m(7);
t234 = 0.2e1 * pkin(4);
t233 = 0.2e1 * pkin(8);
t232 = -pkin(10) - pkin(9);
t222 = Ifges(5,4) * t158;
t133 = Ifges(5,2) * t162 + t222;
t231 = -t133 / 0.2e1;
t230 = -t158 / 0.2e1;
t203 = qJD(5) * t157;
t52 = -t203 * t213 + (t238 * t212 + t190) * t161 + t172 * t157;
t30 = -mrSges(7,2) * t52 + mrSges(7,3) * t208;
t33 = -mrSges(6,2) * t208 - mrSges(6,3) * t52;
t227 = t30 + t33;
t117 = t157 * t158 - t161 * t162;
t51 = -t238 * t106 - t117 * t207;
t31 = mrSges(6,1) * t208 - mrSges(6,3) * t51;
t32 = -mrSges(7,1) * t208 + t51 * mrSges(7,2);
t226 = t32 - t31;
t95 = -mrSges(7,2) * t106 - mrSges(7,3) * t163;
t96 = mrSges(6,2) * t163 - mrSges(6,3) * t106;
t224 = t95 + t96;
t107 = t117 * t159;
t97 = -mrSges(6,1) * t163 + mrSges(6,3) * t107;
t98 = mrSges(7,1) * t163 - mrSges(7,2) * t107;
t223 = -t97 + t98;
t221 = Ifges(5,4) * t162;
t220 = Ifges(5,6) * t158;
t156 = cos(pkin(6));
t155 = sin(pkin(6));
t160 = sin(qJ(2));
t215 = t155 * t160;
t108 = -t156 * t163 + t159 * t215;
t109 = t156 * t159 + t163 * t215;
t164 = cos(qJ(2));
t209 = qJD(2) * t164;
t193 = t155 * t209;
t88 = t109 * qJD(3) + t159 * t193;
t60 = t108 * t88;
t219 = t163 * Ifges(5,6);
t218 = t88 * t159;
t89 = -t108 * qJD(3) + t163 * t193;
t217 = t89 * t163;
t132 = -mrSges(5,1) * t162 + mrSges(5,2) * t158;
t216 = -mrSges(4,1) + t132;
t214 = t155 * t164;
t128 = pkin(4) * t213 + t159 * pkin(8);
t202 = qJD(5) * t161;
t201 = pkin(4) * t206;
t200 = pkin(4) * t203;
t199 = pkin(4) * t202;
t152 = pkin(8) * t207;
t198 = t232 * t158;
t173 = -t109 * t162 + t158 * t214;
t90 = -t109 * t158 - t162 * t214;
t37 = -t157 * t173 - t161 * t90;
t197 = t37 * t203;
t135 = t232 * t162;
t93 = -t157 * t135 - t161 * t198;
t196 = t93 * t203;
t101 = t171 * pkin(4) + t152;
t147 = -pkin(4) * t162 - pkin(3);
t195 = qJD(4) * t232;
t194 = qJD(2) * t215;
t126 = t158 * t195;
t184 = t162 * t195;
t55 = -t93 * qJD(5) + t161 * t126 + t157 * t184;
t94 = -t161 * t135 + t157 * t198;
t56 = t94 * qJD(5) + t157 * t126 - t161 * t184;
t188 = t94 * t55 + t56 * t93;
t187 = (2 * Ifges(4,4)) + t220;
t102 = -t163 * t229 + t116;
t186 = -qJD(4) * t102 + t58;
t183 = mrSges(5,1) * t158 + mrSges(5,2) * t162;
t182 = Ifges(5,1) * t162 - t222;
t134 = Ifges(5,1) * t158 + t221;
t181 = -Ifges(5,2) * t158 + t221;
t180 = Ifges(5,5) * t158 + Ifges(5,6) * t162;
t34 = -t157 * t92 + t161 * t74;
t38 = t157 * t90 - t161 * t173;
t177 = t243 * t208 - t240 * t52 - t241 * t51;
t24 = t173 * qJD(4) - t89 * t158 + t162 * t194;
t25 = t90 * qJD(4) + t158 * t194 + t89 * t162;
t5 = -t37 * qJD(5) + t157 * t24 + t161 * t25;
t6 = t38 * qJD(5) + t157 * t25 - t161 * t24;
t176 = t244 * t5 + t245 * t6;
t175 = t37 * t56 + t55 * t38 + t94 * t5 + t6 * t93;
t174 = t108 * t207 + t218;
t11 = t157 * t29 + t161 * t45 + t74 * t202 - t92 * t203;
t81 = t238 * t118;
t75 = Ifges(7,6) * t81;
t76 = Ifges(6,6) * t81;
t80 = t238 * t117;
t77 = Ifges(6,5) * t80;
t78 = Ifges(7,4) * t80;
t169 = t244 * t55 + t245 * t56 + t75 - t76 - t77 - t78;
t8 = qJ(6) * t208 - qJD(6) * t163 + t11;
t9 = -pkin(5) * t208 - t12;
t168 = t12 * mrSges(6,1) - t9 * mrSges(7,1) - t11 * mrSges(6,2) + t8 * mrSges(7,3) - t177;
t166 = -t24 * t158 + t25 * t162 + (t158 * t173 - t162 * t90) * qJD(4);
t140 = qJD(6) + t199;
t165 = -mrSges(6,2) * t199 + t140 * mrSges(7,3) + t245 * t200;
t154 = qJD(6) * mrSges(7,3);
t151 = Ifges(5,5) * t204;
t146 = -pkin(4) * t161 - pkin(5);
t144 = pkin(4) * t157 + qJ(6);
t125 = -mrSges(5,1) * t163 - mrSges(5,3) * t212;
t124 = mrSges(5,2) * t163 - mrSges(5,3) * t213;
t123 = t182 * qJD(4);
t122 = t181 * qJD(4);
t121 = (mrSges(4,1) * t159 + mrSges(4,2) * t163) * qJD(3);
t120 = t183 * qJD(4);
t113 = t183 * t159;
t105 = -Ifges(5,5) * t163 + t182 * t159;
t104 = t181 * t159 - t219;
t100 = -mrSges(5,2) * t208 - t171 * mrSges(5,3);
t99 = mrSges(5,1) * t208 - t172 * mrSges(5,3);
t87 = Ifges(6,1) * t118 - Ifges(6,4) * t117;
t86 = Ifges(7,1) * t118 + Ifges(7,5) * t117;
t85 = Ifges(6,4) * t118 - Ifges(6,2) * t117;
t84 = Ifges(7,5) * t118 + Ifges(7,3) * t117;
t83 = mrSges(6,1) * t117 + mrSges(6,2) * t118;
t82 = mrSges(7,1) * t117 - mrSges(7,3) * t118;
t72 = t171 * mrSges(5,1) + t172 * mrSges(5,2);
t71 = pkin(5) * t117 - qJ(6) * t118 + t147;
t68 = mrSges(6,1) * t106 - mrSges(6,2) * t107;
t67 = mrSges(7,1) * t106 + mrSges(7,3) * t107;
t66 = -t134 * t205 + (Ifges(5,5) * t159 + t182 * t163) * qJD(3);
t65 = -t133 * t205 + (Ifges(5,6) * t159 + t181 * t163) * qJD(3);
t64 = -Ifges(6,1) * t107 - Ifges(6,4) * t106 - Ifges(6,5) * t163;
t63 = -Ifges(7,1) * t107 - Ifges(7,4) * t163 + Ifges(7,5) * t106;
t62 = -Ifges(6,4) * t107 - Ifges(6,2) * t106 - Ifges(6,6) * t163;
t61 = -Ifges(7,5) * t107 - Ifges(7,6) * t163 + Ifges(7,3) * t106;
t59 = -t103 * qJD(4) + t210;
t57 = pkin(5) * t106 + qJ(6) * t107 + t128;
t44 = -Ifges(6,1) * t80 - Ifges(6,4) * t81;
t43 = -Ifges(7,1) * t80 + Ifges(7,5) * t81;
t42 = -Ifges(6,4) * t80 - Ifges(6,2) * t81;
t41 = -Ifges(7,5) * t80 + Ifges(7,3) * t81;
t40 = mrSges(6,1) * t81 - mrSges(6,2) * t80;
t39 = mrSges(7,1) * t81 + mrSges(7,3) * t80;
t28 = pkin(5) * t163 - t34;
t27 = -qJ(6) * t163 + t225;
t22 = pkin(5) * t81 + qJ(6) * t80 - qJD(6) * t118 + t201;
t19 = mrSges(6,1) * t52 + mrSges(6,2) * t51;
t18 = mrSges(7,1) * t52 - mrSges(7,3) * t51;
t17 = Ifges(6,1) * t51 - Ifges(6,4) * t52 + Ifges(6,5) * t208;
t16 = Ifges(7,1) * t51 + Ifges(7,4) * t208 + Ifges(7,5) * t52;
t15 = Ifges(6,4) * t51 - Ifges(6,2) * t52 + Ifges(6,6) * t208;
t14 = Ifges(7,5) * t51 + Ifges(7,6) * t208 + Ifges(7,3) * t52;
t13 = pkin(5) * t52 - qJ(6) * t51 + qJD(6) * t107 + t101;
t1 = [0.2e1 * m(5) * (-t173 * t25 + t24 * t90 + t60) + 0.2e1 * m(4) * (-t155 ^ 2 * t160 * t209 + t109 * t89 + t60) + 0.2e1 * (m(7) + m(6)) * (t37 * t6 + t38 * t5 + t60); -t173 * t100 + t25 * t124 + t24 * t125 + t90 * t99 + t223 * t6 + t224 * t5 + t227 * t38 + t226 * t37 + (t113 + t67 + t68) * t88 + (t72 + t18 + t19) * t108 + (-t164 * t121 + (-t164 * mrSges(3,2) + (-mrSges(4,1) * t163 + mrSges(4,2) * t159 - mrSges(3,1)) * t160) * qJD(2)) * t155 + (t218 + t217 + (t108 * t163 - t109 * t159) * qJD(3)) * mrSges(4,3) + m(5) * (t102 * t24 + t103 * t25 - t173 * t58 + t59 * t90) - m(4) * pkin(2) * t194 + m(7) * (t108 * t13 + t27 * t5 + t28 * t6 + t37 * t9 + t38 * t8 + t57 * t88) + m(6) * (t101 * t108 + t11 * t38 - t12 * t37 + t128 * t88 + t225 * t5 - t34 * t6) + (m(5) * t174 / 0.2e1 + m(4) * (-t109 * t208 + t174 + t217) / 0.2e1) * t233; (t101 * t128 + t11 * t225 + t12 * t34) * t236 + 0.2e1 * t225 * t33 - (t17 + t16) * t107 + (t14 - t15) * t106 + ((-t158 * t104 + t162 * t105 + t113 * t233 + t187 * t163) * qJD(3) + t177 + t239) * t163 + (t13 * t57 + t27 * t8 + t28 * t9) * t235 + (t102 * t59 + t103 * t58) * t237 + (t61 - t62) * t52 + 0.2e1 * t128 * t19 - 0.2e1 * pkin(2) * t121 + 0.2e1 * t58 * t124 + 0.2e1 * t59 * t125 + (t63 + t64) * t51 + 0.2e1 * t102 * t99 + 0.2e1 * t103 * t100 + 0.2e1 * t8 * t95 + 0.2e1 * t11 * t96 + 0.2e1 * t12 * t97 + 0.2e1 * t9 * t98 + 0.2e1 * t101 * t68 + 0.2e1 * t13 * t67 + 0.2e1 * t57 * t18 + 0.2e1 * t27 * t30 + 0.2e1 * t28 * t32 + 0.2e1 * t34 * t31 + (t72 * t233 - t158 * t65 + t162 * t66 + (-t162 * t104 - t158 * t105 + t163 * t180) * qJD(4) + ((Ifges(5,5) * t162 - t187) * t159 - t241 * t107 + t240 * t106 + (pkin(8) ^ 2 * t237 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) + t243) * t163) * qJD(3)) * t159; -t89 * mrSges(4,2) + (t120 + t39 + t40) * t108 + (t82 + t83 + t216) * t88 + m(7) * (t108 * t22 + t71 * t88 + t175) + m(6) * (t108 * t201 + t147 * t88 + t175) + t166 * mrSges(5,3) + t242 * (-t117 * t5 + t118 * t6 - t37 * t80 - t38 * t81) + (-pkin(3) * t88 + t166 * pkin(9)) * m(5); m(6) * (t101 * t147 + t11 * t94 - t12 * t93 + t128 * t201 + t225 * t55 - t34 * t56) + m(5) * (-pkin(3) * t152 + (-t103 * t206 - t59 * t158) * pkin(9)) + (-t11 * t117 - t118 * t12 - t225 * t81 + t34 * t80) * mrSges(6,3) + (-t117 * t8 + t118 * t9 - t27 * t81 - t28 * t80) * mrSges(7,2) + (t77 / 0.2e1 + t76 / 0.2e1 + t78 / 0.2e1 - t75 / 0.2e1 - t151 / 0.2e1 + (t216 * pkin(8) + Ifges(4,5)) * qJD(3)) * t163 + (t162 * t123 / 0.2e1 + t122 * t230 - Ifges(4,6) * qJD(3) + (t134 * t230 + t162 * t231) * qJD(4) + (mrSges(4,2) * qJD(3) + t120) * pkin(8) + (t240 * t117 + t241 * t118 + t180) * qJD(3) / 0.2e1) * t159 + (t207 * t231 + t66 / 0.2e1 - pkin(9) * t99 - t59 * mrSges(5,3) + (-pkin(9) * t124 - t103 * mrSges(5,3) + pkin(4) * t68 + t219 / 0.2e1 - t104 / 0.2e1) * qJD(4)) * t158 + t224 * t55 + t226 * t93 + t227 * t94 + t223 * t56 + (t134 * t207 / 0.2e1 + t65 / 0.2e1 + qJD(4) * t105 / 0.2e1 + t186 * mrSges(5,3) + (m(5) * t186 - qJD(4) * t125 + t100) * pkin(9)) * t162 + m(7) * (t13 * t71 + t22 * t57 + t27 * t55 + t28 * t56 + t8 * t94 + t9 * t93) + (t84 / 0.2e1 - t85 / 0.2e1) * t52 + (t86 / 0.2e1 + t87 / 0.2e1) * t51 + (t61 / 0.2e1 - t62 / 0.2e1) * t81 - (t63 / 0.2e1 + t64 / 0.2e1) * t80 + t147 * t19 + t128 * t40 + t101 * t83 + t13 * t82 - pkin(3) * t72 + t22 * t67 + t71 * t18 + t57 * t39 + (t16 / 0.2e1 + t17 / 0.2e1) * t118 + (t14 / 0.2e1 - t15 / 0.2e1) * t117 - (t43 / 0.2e1 + t44 / 0.2e1) * t107 + (t41 / 0.2e1 - t42 / 0.2e1) * t106; -0.2e1 * pkin(3) * t120 + t162 * t122 + t158 * t123 + 0.2e1 * t147 * t40 + 0.2e1 * t22 * t82 + 0.2e1 * t71 * t39 + (t84 - t85) * t81 - (t86 + t87) * t80 + (t43 + t44) * t118 + (t41 - t42) * t117 + (t162 * t134 + (t83 * t234 - t133) * t158) * qJD(4) + (t147 * t201 + t188) * t236 + (t22 * t71 + t188) * t235 + 0.2e1 * t242 * (-t117 * t55 + t118 * t56 - t80 * t93 - t81 * t94); t24 * mrSges(5,1) - t25 * mrSges(5,2) + m(7) * (t140 * t38 + t144 * t5 + t146 * t6) + (m(7) * t197 / 0.2e1 + m(6) * (t157 * t5 - t161 * t6 + t38 * t202 + t197) / 0.2e1) * t234 + t176; m(7) * (t140 * t27 + t144 * t8 + t146 * t9) + ((m(6) * t12 + t31 + (m(6) * t225 + t96) * qJD(5)) * t161 + (m(6) * t11 + t33 + (-m(6) * t34 + m(7) * t28 + t223) * qJD(5)) * t157) * pkin(4) - t171 * Ifges(5,6) - Ifges(5,5) * t192 + t140 * t95 + t144 * t30 + t146 * t32 - t58 * mrSges(5,2) + t59 * mrSges(5,1) + t168 - t239; m(7) * (t140 * t94 + t144 * t55 + t146 * t56) + t151 + (t132 * pkin(9) - t220) * qJD(4) + (-t117 * t140 - t144 * t81 - t146 * t80) * mrSges(7,2) + (t118 * mrSges(7,2) * t203 + m(7) * t196 + m(6) * (t157 * t55 - t161 * t56 + t94 * t202 + t196) + (-t157 * t81 + t161 * t80 + (-t117 * t161 + t118 * t157) * qJD(5)) * mrSges(6,3)) * pkin(4) + t169; 0.2e1 * m(7) * (t140 * t144 + t146 * t200) + 0.2e1 * t165; m(7) * (-pkin(5) * t6 + qJ(6) * t5 + qJD(6) * t38) + t176; m(7) * (-pkin(5) * t9 + qJ(6) * t8 + qJD(6) * t27) + qJD(6) * t95 + qJ(6) * t30 - pkin(5) * t32 + t168; m(7) * (-pkin(5) * t56 + qJ(6) * t55 + qJD(6) * t94) + (pkin(5) * t80 - qJ(6) * t81 - qJD(6) * t117) * mrSges(7,2) + t169; m(7) * (-pkin(5) * t200 + qJ(6) * t140 + qJD(6) * t144) + t154 + t165; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t154; m(7) * t6; m(7) * t9 + t32; m(7) * t56 - t80 * mrSges(7,2); m(7) * t200; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
