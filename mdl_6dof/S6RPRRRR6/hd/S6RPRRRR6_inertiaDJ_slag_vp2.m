% Calculate time derivative of joint inertia matrix for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:18
% EndTime: 2019-03-09 07:12:27
% DurationCPUTime: 4.60s
% Computational Cost: add. (12062->483), mult. (26580->723), div. (0->0), fcn. (27073->10), ass. (0->190)
t173 = sin(pkin(11));
t174 = cos(pkin(11));
t178 = sin(qJ(3));
t182 = cos(qJ(3));
t149 = t173 * t178 - t182 * t174;
t141 = t149 * qJD(3);
t150 = t173 * t182 + t178 * t174;
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t210 = qJD(4) * t181;
t188 = -t177 * t141 + t150 * t210;
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t190 = t176 * t177 - t180 * t181;
t104 = t190 * t150;
t203 = -pkin(2) * t174 - pkin(1);
t112 = pkin(3) * t149 - pkin(8) * t150 + t203;
t229 = pkin(7) + qJ(2);
t161 = t229 * t173;
t162 = t229 * t174;
t127 = -t178 * t161 + t162 * t182;
t120 = t181 * t127;
t78 = t177 * t112 + t120;
t142 = t150 * qJD(3);
t213 = t181 * t141;
t244 = -Ifges(5,5) * t213 + Ifges(5,3) * t142;
t243 = -t182 * t161 - t162 * t178;
t234 = -pkin(9) - pkin(8);
t165 = t234 * t177;
t166 = t234 * t181;
t129 = t176 * t165 - t180 * t166;
t194 = mrSges(5,1) * t177 + mrSges(5,2) * t181;
t156 = t194 * qJD(4);
t242 = qJD(4) + qJD(5);
t152 = t176 * t181 + t177 * t180;
t122 = t242 * t152;
t48 = -t122 * t150 + t141 * t190;
t49 = t104 * t242 + t152 * t141;
t241 = Ifges(6,5) * t48 + Ifges(6,6) * t49 + Ifges(6,3) * t142;
t240 = 2 * m(6);
t239 = 2 * m(7);
t238 = -2 * mrSges(4,3);
t236 = -0.2e1 * t243;
t235 = m(6) * pkin(4);
t231 = -t150 / 0.2e1;
t227 = Ifges(5,4) * t177;
t163 = Ifges(5,2) * t181 + t227;
t230 = -t163 / 0.2e1;
t175 = sin(qJ(6));
t179 = cos(qJ(6));
t115 = -t152 * t175 - t179 * t190;
t121 = t242 * t190;
t59 = qJD(6) * t115 - t121 * t179 - t122 * t175;
t116 = t152 * t179 - t175 * t190;
t60 = -qJD(6) * t116 + t121 * t175 - t122 * t179;
t228 = Ifges(7,5) * t59 + Ifges(7,6) * t60;
t218 = t150 * t181;
t77 = t181 * t112 - t127 * t177;
t56 = pkin(4) * t149 - pkin(9) * t218 + t77;
t219 = t150 * t177;
t68 = -pkin(9) * t219 + t78;
t38 = t176 * t56 + t180 * t68;
t226 = Ifges(5,4) * t181;
t225 = Ifges(5,6) * t177;
t168 = pkin(4) * t180 + pkin(5);
t206 = qJD(6) * t179;
t207 = qJD(6) * t175;
t216 = t175 * t176;
t110 = t168 * t206 + (-t176 * t207 + (t179 * t180 - t216) * qJD(5)) * pkin(4);
t224 = t110 * mrSges(7,2);
t96 = qJD(2) * t150 + qJD(3) * t127;
t223 = t243 * t96;
t222 = t142 * Ifges(5,5);
t221 = t142 * Ifges(5,6);
t220 = t149 * Ifges(5,6);
t215 = t176 * t179;
t212 = -Ifges(6,5) * t121 - Ifges(6,6) * t122;
t211 = qJD(4) * t177;
t209 = qJD(5) * t176;
t208 = qJD(5) * t180;
t103 = t152 * t150;
t72 = -t103 * t179 + t104 * t175;
t18 = qJD(6) * t72 + t175 * t49 + t179 * t48;
t73 = -t103 * t175 - t104 * t179;
t19 = -qJD(6) * t73 - t175 * t48 + t179 * t49;
t205 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t142;
t204 = pkin(4) * t211;
t169 = -pkin(4) * t181 - pkin(3);
t202 = qJD(4) * t234;
t201 = t150 * t211;
t32 = -t60 * mrSges(7,1) + t59 * mrSges(7,2);
t199 = -(2 * Ifges(4,4)) - t225;
t37 = -t176 * t68 + t180 * t56;
t198 = t142 * mrSges(4,1) - t141 * mrSges(4,2);
t82 = t122 * mrSges(6,1) - t121 * mrSges(6,2);
t111 = -t168 * t207 + (-t176 * t206 + (-t175 * t180 - t215) * qJD(5)) * pkin(4);
t108 = t111 * mrSges(7,1);
t197 = t108 - t224;
t109 = pkin(3) * t142 + pkin(8) * t141;
t95 = -t149 * qJD(2) + qJD(3) * t243;
t196 = t181 * t109 - t177 * t95;
t128 = t180 * t165 + t166 * t176;
t97 = pkin(4) * t219 - t243;
t193 = Ifges(5,1) * t181 - t227;
t192 = -Ifges(5,2) * t177 + t226;
t26 = pkin(5) * t149 + pkin(10) * t104 + t37;
t27 = -pkin(10) * t103 + t38;
t12 = -t175 * t27 + t179 * t26;
t13 = t175 * t26 + t179 * t27;
t159 = t177 * t202;
t160 = t181 * t202;
t86 = t180 * t159 + t176 * t160 + t165 * t208 + t166 * t209;
t66 = -pkin(10) * t122 + t86;
t87 = -qJD(5) * t129 - t159 * t176 + t180 * t160;
t67 = pkin(10) * t121 + t87;
t100 = -pkin(10) * t190 + t129;
t99 = -pkin(10) * t152 + t128;
t69 = -t100 * t175 + t179 * t99;
t21 = qJD(6) * t69 + t175 * t67 + t179 * t66;
t70 = t100 * t179 + t175 * t99;
t22 = -qJD(6) * t70 - t175 * t66 + t179 * t67;
t191 = t22 * mrSges(7,1) - t21 * mrSges(7,2) + t228;
t30 = pkin(9) * t213 + pkin(4) * t142 + (-t120 + (pkin(9) * t150 - t112) * t177) * qJD(4) + t196;
t42 = t177 * t109 + t112 * t210 - t127 * t211 + t181 * t95;
t36 = -pkin(9) * t188 + t42;
t11 = -qJD(5) * t38 - t176 * t36 + t180 * t30;
t4 = pkin(5) * t142 - pkin(10) * t48 + t11;
t10 = t176 * t30 + t180 * t36 + t56 * t208 - t209 * t68;
t5 = pkin(10) * t49 + t10;
t2 = qJD(6) * t12 + t175 * t4 + t179 * t5;
t3 = -qJD(6) * t13 - t175 * t5 + t179 * t4;
t189 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t205;
t187 = t201 + t213;
t186 = -t32 - t82;
t185 = (-mrSges(6,1) * t176 - mrSges(6,2) * t180) * qJD(5) * pkin(4);
t184 = t87 * mrSges(6,1) - t86 * mrSges(6,2) + t191 + t212;
t183 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + t189 + t241;
t74 = pkin(4) * t188 + t96;
t170 = Ifges(5,5) * t210;
t164 = Ifges(5,1) * t177 + t226;
t158 = t193 * qJD(4);
t157 = t192 * qJD(4);
t146 = (-mrSges(7,1) * t175 - mrSges(7,2) * t179) * qJD(6) * pkin(5);
t140 = pkin(4) * t215 + t168 * t175;
t139 = -pkin(4) * t216 + t168 * t179;
t132 = pkin(5) * t190 + t169;
t125 = Ifges(6,1) * t152 - Ifges(6,4) * t190;
t124 = Ifges(6,4) * t152 - Ifges(6,2) * t190;
t123 = mrSges(6,1) * t190 + mrSges(6,2) * t152;
t114 = mrSges(5,1) * t149 - mrSges(5,3) * t218;
t113 = -mrSges(5,2) * t149 - mrSges(5,3) * t219;
t107 = pkin(5) * t122 + t204;
t93 = Ifges(5,5) * t149 + t150 * t193;
t92 = t150 * t192 + t220;
t91 = mrSges(6,1) * t149 + mrSges(6,3) * t104;
t90 = -mrSges(6,2) * t149 - mrSges(6,3) * t103;
t89 = -mrSges(5,2) * t142 - mrSges(5,3) * t188;
t88 = mrSges(5,1) * t142 + mrSges(5,3) * t187;
t84 = -Ifges(6,1) * t121 - Ifges(6,4) * t122;
t83 = -Ifges(6,4) * t121 - Ifges(6,2) * t122;
t81 = Ifges(7,1) * t116 + Ifges(7,4) * t115;
t80 = Ifges(7,4) * t116 + Ifges(7,2) * t115;
t79 = -mrSges(7,1) * t115 + mrSges(7,2) * t116;
t76 = mrSges(5,1) * t188 - mrSges(5,2) * t187;
t75 = mrSges(6,1) * t103 - mrSges(6,2) * t104;
t71 = pkin(5) * t103 + t97;
t64 = -Ifges(6,1) * t104 - Ifges(6,4) * t103 + Ifges(6,5) * t149;
t63 = -Ifges(6,4) * t104 - Ifges(6,2) * t103 + Ifges(6,6) * t149;
t62 = -Ifges(5,1) * t187 - Ifges(5,4) * t188 + t222;
t61 = -Ifges(5,4) * t187 - Ifges(5,2) * t188 + t221;
t58 = mrSges(7,1) * t149 - mrSges(7,3) * t73;
t57 = -mrSges(7,2) * t149 + mrSges(7,3) * t72;
t45 = -mrSges(6,2) * t142 + mrSges(6,3) * t49;
t44 = mrSges(6,1) * t142 - mrSges(6,3) * t48;
t43 = -qJD(4) * t78 + t196;
t41 = -mrSges(7,1) * t72 + mrSges(7,2) * t73;
t40 = Ifges(7,1) * t73 + Ifges(7,4) * t72 + Ifges(7,5) * t149;
t39 = Ifges(7,4) * t73 + Ifges(7,2) * t72 + Ifges(7,6) * t149;
t34 = Ifges(7,1) * t59 + Ifges(7,4) * t60;
t33 = Ifges(7,4) * t59 + Ifges(7,2) * t60;
t31 = -t49 * pkin(5) + t74;
t25 = -mrSges(6,1) * t49 + mrSges(6,2) * t48;
t24 = Ifges(6,1) * t48 + Ifges(6,4) * t49 + t142 * Ifges(6,5);
t23 = Ifges(6,4) * t48 + Ifges(6,2) * t49 + Ifges(6,6) * t142;
t15 = -mrSges(7,2) * t142 + mrSges(7,3) * t19;
t14 = mrSges(7,1) * t142 - mrSges(7,3) * t18;
t8 = -mrSges(7,1) * t19 + mrSges(7,2) * t18;
t7 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t142 * Ifges(7,5);
t6 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t142 * Ifges(7,6);
t1 = [0.2e1 * t203 * t198 + 0.2e1 * m(4) * (t127 * t95 - t223) + 0.2e1 * m(5) * (t42 * t78 + t43 * t77 - t223) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t173 ^ 2 + t174 ^ 2) + 0.2e1 * t42 * t113 + 0.2e1 * t43 * t114 - t103 * t23 - t104 * t24 + 0.2e1 * t10 * t90 + 0.2e1 * t11 * t91 + 0.2e1 * t97 * t25 + 0.2e1 * t77 * t88 + 0.2e1 * t78 * t89 + 0.2e1 * t71 * t8 + t72 * t6 + t73 * t7 + 0.2e1 * t74 * t75 + t49 * t63 + t48 * t64 + 0.2e1 * t2 * t57 + 0.2e1 * t3 * t58 + 0.2e1 * t37 * t44 + 0.2e1 * t38 * t45 + t19 * t39 + t18 * t40 + 0.2e1 * t31 * t41 + 0.2e1 * t12 * t14 + 0.2e1 * t13 * t15 + (-t199 * t141 + t95 * t238 + t205 + t241 + t244) * t149 + (-t177 * t61 - 0.2e1 * Ifges(4,1) * t141 + t181 * t62 + (t149 * (-Ifges(5,5) * t177 - Ifges(5,6) * t181) - t177 * t93 - t181 * t92) * qJD(4) + 0.2e1 * (t194 + mrSges(4,3)) * t96) * t150 + (Ifges(7,5) * t73 + Ifges(7,6) * t72 + t127 * t238 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t149 + (Ifges(5,5) * t181 + t199) * t150 - Ifges(6,5) * t104 - Ifges(6,6) * t103) * t142 + (t10 * t38 + t11 * t37 + t74 * t97) * t240 + t76 * t236 + (t12 * t3 + t13 * t2 + t31 * t71) * t239 - (mrSges(4,3) * t236 - t177 * t92 + t181 * t93) * t141; t115 * t14 + t116 * t15 - t121 * t90 - t122 * t91 - t190 * t44 + t152 * t45 + t177 * t89 + t181 * t88 + t59 * t57 + t60 * t58 + (t113 * t181 - t114 * t177) * qJD(4) + m(7) * (t115 * t3 + t116 * t2 + t12 * t60 + t13 * t59) + m(6) * (t10 * t152 - t11 * t190 - t121 * t38 - t122 * t37) + m(5) * (t177 * t42 + t181 * t43 + (-t177 * t77 + t181 * t78) * qJD(4)) + t198; 0.2e1 * m(6) * (-t121 * t152 + t122 * t190) + 0.2e1 * m(7) * (t115 * t60 + t116 * t59); (-m(5) * t96 - t76) * pkin(3) + (-t10 * t190 - t11 * t152 + t121 * t37 - t122 * t38) * mrSges(6,3) + (t228 + t212 + t170) * t149 / 0.2e1 + (Ifges(6,5) * t152 + Ifges(7,5) * t116 - Ifges(6,6) * t190 + Ifges(7,6) * t115) * t142 / 0.2e1 - t190 * t23 / 0.2e1 + m(6) * (t10 * t129 + t11 * t128 + t169 * t74 + t37 * t87 + t38 * t86) + t169 * t25 + t152 * t24 / 0.2e1 - Ifges(4,6) * t142 - Ifges(4,5) * t141 - t122 * t63 / 0.2e1 + t74 * t123 + t49 * t124 / 0.2e1 + t48 * t125 / 0.2e1 + t128 * t44 + t129 * t45 + t132 * t8 + t116 * t7 / 0.2e1 - t121 * t64 / 0.2e1 + t115 * t6 / 0.2e1 - t103 * t83 / 0.2e1 - t104 * t84 / 0.2e1 + t107 * t41 + t87 * t91 - t95 * mrSges(4,2) - t96 * mrSges(4,1) + t97 * t82 + t86 * t90 + t31 * t79 + t19 * t80 / 0.2e1 + t18 * t81 / 0.2e1 + t69 * t14 + t70 * t15 + t71 * t32 + t72 * t33 / 0.2e1 + t73 * t34 / 0.2e1 + t60 * t39 / 0.2e1 + t21 * t57 + t22 * t58 + t59 * t40 / 0.2e1 - t243 * t156 + m(7) * (t107 * t71 + t12 * t22 + t13 * t21 + t132 * t31 + t2 * t70 + t3 * t69) + (t157 * t231 - t141 * t230 - t43 * mrSges(5,3) + t222 / 0.2e1 + t96 * mrSges(5,2) + t62 / 0.2e1 + (-m(5) * t43 - t88) * pkin(8) + (-t92 / 0.2e1 - t220 / 0.2e1 + t164 * t231 - t78 * mrSges(5,3) + pkin(4) * t75 + t97 * t235 + (-m(5) * t78 - t113) * pkin(8)) * qJD(4)) * t177 + (t221 / 0.2e1 - t96 * mrSges(5,1) + t61 / 0.2e1 + t150 * t158 / 0.2e1 - t141 * t164 / 0.2e1 + t42 * mrSges(5,3) + (t93 / 0.2e1 + t150 * t230 - t77 * mrSges(5,3)) * qJD(4) + (m(5) * (-t77 * qJD(4) + t42) + t89 - qJD(4) * t114) * pkin(8)) * t181 + (t115 * t2 - t116 * t3 - t12 * t59 + t13 * t60) * mrSges(7,3); m(6) * (-t121 * t129 - t122 * t128 + t152 * t86 - t190 * t87) + m(7) * (t115 * t22 + t116 * t21 + t59 * t70 + t60 * t69); -0.2e1 * pkin(3) * t156 + 0.2e1 * t107 * t79 + t115 * t33 + t116 * t34 - t121 * t125 - t122 * t124 + 0.2e1 * t132 * t32 - t190 * t83 + t152 * t84 + t181 * t157 + t177 * t158 + 0.2e1 * t169 * t82 + t59 * t81 + t60 * t80 + (t181 * t164 + (0.2e1 * pkin(4) * t123 - t163) * t177) * qJD(4) + (t128 * t87 + t129 * t86 + t169 * t204) * t240 + (t107 * t132 + t21 * t70 + t22 * t69) * t239 + 0.2e1 * (t115 * t21 - t116 * t22 - t59 * t69 + t60 * t70) * mrSges(7,3) + 0.2e1 * (t121 * t128 - t122 * t129 - t152 * t87 - t190 * t86) * mrSges(6,3); t139 * t14 + t140 * t15 + t110 * t57 + t111 * t58 - t42 * mrSges(5,2) + t43 * mrSges(5,1) + t183 + (m(6) * (t10 * t176 + t11 * t180 + t208 * t38 - t209 * t37) + t180 * t44 + t176 * t45 - t91 * t209 + t90 * t208) * pkin(4) - Ifges(5,5) * t201 + m(7) * (t110 * t13 + t111 * t12 + t139 * t3 + t140 * t2) - t188 * Ifges(5,6) + t244; m(7) * (t110 * t116 + t111 * t115 + t139 * t60 + t140 * t59) - t156 + (-t121 * t176 - t122 * t180 + (t152 * t180 + t176 * t190) * qJD(5)) * t235 + t186; m(7) * (t110 * t70 + t111 * t69 + t139 * t22 + t140 * t21) + t170 + (-t225 + (-mrSges(5,1) * t181 + mrSges(5,2) * t177) * pkin(8)) * qJD(4) + (t110 * t115 - t111 * t116 - t139 * t59 + t140 * t60) * mrSges(7,3) + (m(6) * (t176 * t86 + t180 * t87 + (-t128 * t176 + t129 * t180) * qJD(5)) + (t180 * t121 - t176 * t122 + (t152 * t176 - t180 * t190) * qJD(5)) * mrSges(6,3)) * pkin(4) + t184; (t110 * t140 + t111 * t139) * t239 - 0.2e1 * t224 + 0.2e1 * t108 + 0.2e1 * t185; (m(7) * (-t12 * t207 + t13 * t206 + t175 * t2 + t179 * t3) + t57 * t206 + t175 * t15 - t58 * t207 + t179 * t14) * pkin(5) + t183; m(7) * (t175 * t59 + t179 * t60 + (-t115 * t175 + t116 * t179) * qJD(6)) * pkin(5) + t186; (m(7) * (t175 * t21 + t179 * t22 + (-t175 * t69 + t179 * t70) * qJD(6)) + (t175 * t60 - t179 * t59 + (t115 * t179 + t116 * t175) * qJD(6)) * mrSges(7,3)) * pkin(5) + t184; t185 + (m(7) * (t110 * t175 + t111 * t179 - t139 * t207 + t140 * t206) - mrSges(7,2) * t206 - mrSges(7,1) * t207) * pkin(5) + t197; 0.2e1 * t146; t189; -t32; t191; t197; t146; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
