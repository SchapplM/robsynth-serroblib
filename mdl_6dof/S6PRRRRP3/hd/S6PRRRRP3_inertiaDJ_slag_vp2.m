% Calculate time derivative of joint inertia matrix for
% S6PRRRRP3
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:26
% EndTime: 2019-03-09 00:06:35
% DurationCPUTime: 4.12s
% Computational Cost: add. (4162->472), mult. (10896->674), div. (0->0), fcn. (9722->10), ass. (0->195)
t244 = mrSges(6,1) + mrSges(7,1);
t240 = Ifges(6,5) + Ifges(7,5);
t239 = Ifges(6,6) + Ifges(7,6);
t243 = -Ifges(6,3) - Ifges(7,3);
t163 = sin(qJ(4));
t164 = sin(qJ(3));
t167 = cos(qJ(4));
t202 = qJD(4) * t167;
t168 = cos(qJ(3));
t205 = qJD(3) * t168;
t175 = t163 * t205 + t164 * t202;
t166 = cos(qJ(5));
t242 = pkin(4) * t166;
t238 = pkin(4) * qJD(5);
t162 = sin(qJ(5));
t180 = t162 * t163 - t166 * t167;
t111 = t180 * t164;
t230 = -pkin(10) - pkin(9);
t143 = t230 * t163;
t144 = t230 * t167;
t97 = t162 * t143 - t166 * t144;
t190 = t167 * t205;
t206 = qJD(3) * t164;
t237 = -Ifges(5,5) * t190 - Ifges(5,3) * t206;
t139 = -pkin(3) * t168 - pkin(9) * t164 - pkin(2);
t210 = t167 * t168;
t150 = pkin(8) * t210;
t106 = t163 * t139 + t150;
t236 = qJD(4) + qJD(5);
t235 = 2 * m(5);
t234 = 2 * m(6);
t233 = 2 * m(7);
t232 = 0.2e1 * pkin(8);
t231 = m(7) * pkin(5);
t221 = Ifges(5,4) * t163;
t141 = Ifges(5,2) * t167 + t221;
t229 = -t141 / 0.2e1;
t228 = -t163 / 0.2e1;
t200 = qJD(5) * t166;
t161 = cos(pkin(6));
t160 = sin(pkin(6));
t165 = sin(qJ(2));
t214 = t160 * t165;
t113 = t161 * t164 + t168 * t214;
t169 = cos(qJ(2));
t213 = t160 * t169;
t177 = -t113 * t167 + t163 * t213;
t93 = -t113 * t163 - t167 * t213;
t39 = t162 * t93 - t166 * t177;
t195 = qJD(2) * t214;
t112 = -t161 * t168 + t164 * t214;
t207 = qJD(2) * t169;
t194 = t160 * t207;
t92 = -qJD(3) * t112 + t168 * t194;
t25 = qJD(4) * t177 - t92 * t163 + t167 * t195;
t26 = qJD(4) * t93 + t163 * t195 + t92 * t167;
t38 = t162 * t177 + t166 * t93;
t8 = qJD(5) * t38 + t162 * t25 + t166 * t26;
t227 = (t162 * t8 + t200 * t39) * pkin(4);
t226 = pkin(8) * t163;
t224 = -mrSges(7,2) - mrSges(6,2);
t125 = t162 * t167 + t163 * t166;
t53 = t236 * t111 - t125 * t205;
t33 = -mrSges(7,2) * t206 + mrSges(7,3) * t53;
t34 = -mrSges(6,2) * t206 + mrSges(6,3) * t53;
t223 = t33 + t34;
t123 = t167 * t139;
t211 = t164 * t167;
t75 = -pkin(10) * t211 + t123 + (-pkin(4) - t226) * t168;
t212 = t163 * t164;
t95 = -pkin(10) * t212 + t106;
t36 = t162 * t75 + t166 * t95;
t110 = t125 * t164;
t98 = mrSges(7,2) * t168 - mrSges(7,3) * t110;
t99 = mrSges(6,2) * t168 - mrSges(6,3) * t110;
t222 = t98 + t99;
t220 = Ifges(5,4) * t167;
t219 = Ifges(5,6) * t163;
t91 = qJD(3) * t113 + t164 * t194;
t59 = t112 * t91;
t218 = t168 * Ifges(5,6);
t217 = t91 * t164;
t216 = t92 * t168;
t140 = -mrSges(5,1) * t167 + mrSges(5,2) * t163;
t215 = -mrSges(4,1) + t140;
t100 = -mrSges(7,1) * t168 + mrSges(7,3) * t111;
t101 = -mrSges(6,1) * t168 + mrSges(6,3) * t111;
t209 = t100 + t101;
t137 = (pkin(3) * t164 - pkin(9) * t168) * qJD(3);
t208 = t167 * t137 + t206 * t226;
t138 = pkin(4) * t212 + t164 * pkin(8);
t204 = qJD(4) * t163;
t203 = qJD(4) * t164;
t201 = qJD(5) * t162;
t199 = pkin(4) * t204;
t158 = pkin(8) * t205;
t30 = (pkin(4) * t164 - pkin(10) * t210) * qJD(3) + (-t150 + (pkin(10) * t164 - t139) * t163) * qJD(4) + t208;
t57 = t163 * t137 + t139 * t202 + (-t167 * t206 - t168 * t204) * pkin(8);
t46 = -pkin(10) * t175 + t57;
t12 = -qJD(5) * t36 - t162 * t46 + t166 * t30;
t84 = t236 * t125;
t52 = -t164 * t84 - t180 * t205;
t3 = pkin(5) * t206 - qJ(6) * t52 + qJD(6) * t111 + t12;
t31 = mrSges(7,1) * t206 - mrSges(7,3) * t52;
t198 = m(7) * t3 + t31;
t197 = t38 * t201;
t104 = t175 * pkin(4) + t158;
t153 = -pkin(4) * t167 - pkin(3);
t196 = qJD(4) * t230;
t193 = t163 * t203;
t17 = -t53 * mrSges(7,1) + t52 * mrSges(7,2);
t83 = t236 * t180;
t40 = t84 * mrSges(7,1) - t83 * mrSges(7,2);
t189 = t224 * t166;
t188 = (2 * Ifges(4,4)) + t219;
t35 = -t162 * t95 + t166 * t75;
t105 = -t168 * t226 + t123;
t187 = -t105 * qJD(4) + t57;
t96 = t166 * t143 + t144 * t162;
t135 = t163 * t196;
t136 = t167 * t196;
t56 = -t97 * qJD(5) - t135 * t162 + t166 * t136;
t21 = qJ(6) * t83 - qJD(6) * t125 + t56;
t186 = m(7) * t21 + t83 * mrSges(7,3);
t185 = mrSges(5,1) * t163 + mrSges(5,2) * t167;
t184 = Ifges(5,1) * t167 - t221;
t142 = Ifges(5,1) * t163 + t220;
t183 = -Ifges(5,2) * t163 + t220;
t182 = Ifges(5,5) * t163 + Ifges(5,6) * t167;
t9 = -qJD(5) * t39 - t162 * t26 + t166 * t25;
t181 = t224 * t8 + t244 * t9;
t179 = t243 * t206 - t239 * t53 - t240 * t52;
t178 = t112 * t205 + t217;
t11 = t162 * t30 + t166 * t46 + t75 * t200 - t201 * t95;
t55 = t166 * t135 + t162 * t136 + t143 * t200 + t144 * t201;
t176 = t190 - t193;
t20 = -qJ(6) * t84 - qJD(6) * t180 + t55;
t77 = Ifges(7,6) * t84;
t78 = Ifges(6,6) * t84;
t79 = Ifges(7,5) * t83;
t80 = Ifges(6,5) * t83;
t173 = t56 * mrSges(6,1) + t21 * mrSges(7,1) - t55 * mrSges(6,2) - t20 * mrSges(7,2) - t77 - t78 - t79 - t80;
t172 = -t162 * t84 + (t125 * t162 - t166 * t180) * qJD(5);
t4 = qJ(6) * t53 - qJD(6) * t110 + t11;
t171 = t12 * mrSges(6,1) + t3 * mrSges(7,1) - t11 * mrSges(6,2) - t4 * mrSges(7,2) - t179;
t170 = -t163 * t25 + t167 * t26 + (t163 * t177 - t167 * t93) * qJD(4);
t157 = Ifges(5,5) * t202;
t152 = pkin(5) + t242;
t134 = -mrSges(5,1) * t168 - mrSges(5,3) * t211;
t133 = mrSges(5,2) * t168 - mrSges(5,3) * t212;
t132 = t184 * qJD(4);
t131 = t183 * qJD(4);
t130 = (mrSges(4,1) * t164 + mrSges(4,2) * t168) * qJD(3);
t129 = t185 * qJD(4);
t119 = t185 * t164;
t109 = -Ifges(5,5) * t168 + t164 * t184;
t108 = t164 * t183 - t218;
t107 = pkin(5) * t180 + t153;
t103 = -mrSges(5,2) * t206 - mrSges(5,3) * t175;
t102 = mrSges(5,1) * t206 - mrSges(5,3) * t176;
t90 = Ifges(6,1) * t125 - Ifges(6,4) * t180;
t89 = Ifges(7,1) * t125 - Ifges(7,4) * t180;
t88 = Ifges(6,4) * t125 - Ifges(6,2) * t180;
t87 = Ifges(7,4) * t125 - Ifges(7,2) * t180;
t86 = mrSges(6,1) * t180 + mrSges(6,2) * t125;
t85 = mrSges(7,1) * t180 + mrSges(7,2) * t125;
t81 = pkin(5) * t110 + t138;
t74 = mrSges(5,1) * t175 + mrSges(5,2) * t176;
t70 = pkin(5) * t84 + t199;
t69 = mrSges(6,1) * t110 - mrSges(6,2) * t111;
t68 = mrSges(7,1) * t110 - mrSges(7,2) * t111;
t67 = -qJ(6) * t180 + t97;
t66 = -qJ(6) * t125 + t96;
t65 = -t142 * t203 + (Ifges(5,5) * t164 + t168 * t184) * qJD(3);
t64 = -t141 * t203 + (Ifges(5,6) * t164 + t168 * t183) * qJD(3);
t63 = -Ifges(6,1) * t111 - Ifges(6,4) * t110 - Ifges(6,5) * t168;
t62 = -Ifges(7,1) * t111 - Ifges(7,4) * t110 - Ifges(7,5) * t168;
t61 = -Ifges(6,4) * t111 - Ifges(6,2) * t110 - Ifges(6,6) * t168;
t60 = -Ifges(7,4) * t111 - Ifges(7,2) * t110 - Ifges(7,6) * t168;
t58 = -t106 * qJD(4) + t208;
t45 = -Ifges(6,1) * t83 - Ifges(6,4) * t84;
t44 = -Ifges(7,1) * t83 - Ifges(7,4) * t84;
t43 = -Ifges(6,4) * t83 - Ifges(6,2) * t84;
t42 = -Ifges(7,4) * t83 - Ifges(7,2) * t84;
t41 = mrSges(6,1) * t84 - mrSges(6,2) * t83;
t32 = mrSges(6,1) * t206 - mrSges(6,3) * t52;
t24 = -pkin(5) * t53 + t104;
t23 = -qJ(6) * t110 + t36;
t22 = -pkin(5) * t168 + qJ(6) * t111 + t35;
t18 = -mrSges(6,1) * t53 + mrSges(6,2) * t52;
t16 = Ifges(6,1) * t52 + Ifges(6,4) * t53 + Ifges(6,5) * t206;
t15 = Ifges(7,1) * t52 + Ifges(7,4) * t53 + Ifges(7,5) * t206;
t14 = Ifges(6,4) * t52 + Ifges(6,2) * t53 + Ifges(6,6) * t206;
t13 = Ifges(7,4) * t52 + Ifges(7,2) * t53 + Ifges(7,6) * t206;
t1 = [0.2e1 * m(5) * (-t177 * t26 + t25 * t93 + t59) + 0.2e1 * m(4) * (-t160 ^ 2 * t165 * t207 + t113 * t92 + t59) + 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t38 * t9 + t39 * t8 + t59); t93 * t102 - t177 * t103 + t26 * t133 + t25 * t134 + t209 * t9 + t222 * t8 + t223 * t39 + (t31 + t32) * t38 + (t119 + t68 + t69) * t91 + (t74 + t17 + t18) * t112 + (-t169 * t130 + (-t169 * mrSges(3,2) + (-mrSges(4,1) * t168 + mrSges(4,2) * t164 - mrSges(3,1)) * t165) * qJD(2)) * t160 + (t217 + t216 + (t112 * t168 - t113 * t164) * qJD(3)) * mrSges(4,3) + m(5) * (t105 * t25 + t106 * t26 - t177 * t57 + t58 * t93) - m(4) * pkin(2) * t195 + m(7) * (t112 * t24 + t22 * t9 + t23 * t8 + t3 * t38 + t39 * t4 + t81 * t91) + m(6) * (t104 * t112 + t11 * t39 + t12 * t38 + t138 * t91 + t35 * t9 + t36 * t8) + (m(5) * t178 / 0.2e1 + m(4) * (-t113 * t206 + t178 + t216) / 0.2e1) * t232; -(t15 + t16) * t111 - (t14 + t13) * t110 + (t60 + t61) * t53 + (t62 + t63) * t52 + ((-t163 * t108 + t167 * t109 + t119 * t232 + t168 * t188) * qJD(3) + t179 + t237) * t168 + (t22 * t3 + t23 * t4 + t24 * t81) * t233 + (t104 * t138 + t11 * t36 + t12 * t35) * t234 + (t105 * t58 + t106 * t57) * t235 + (t74 * t232 - t163 * t64 + t167 * t65 + (-t167 * t108 - t163 * t109 + t168 * t182) * qJD(4) + ((Ifges(5,5) * t167 - t188) * t164 - t240 * t111 - t239 * t110 + (pkin(8) ^ 2 * t235 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) + t243) * t168) * qJD(3)) * t164 - 0.2e1 * pkin(2) * t130 + 0.2e1 * t57 * t133 + 0.2e1 * t58 * t134 + 0.2e1 * t138 * t18 + 0.2e1 * t4 * t98 + 0.2e1 * t11 * t99 + 0.2e1 * t3 * t100 + 0.2e1 * t12 * t101 + 0.2e1 * t104 * t69 + 0.2e1 * t105 * t102 + 0.2e1 * t106 * t103 + 0.2e1 * t81 * t17 + 0.2e1 * t24 * t68 + 0.2e1 * t35 * t32 + 0.2e1 * t36 * t34 + 0.2e1 * t22 * t31 + 0.2e1 * t23 * t33; -t92 * mrSges(4,2) + (t129 + t40 + t41) * t112 + (t85 + t86 + t215) * t91 + m(6) * (t112 * t199 + t153 * t91 + t38 * t56 + t39 * t55 + t8 * t97 + t9 * t96) + m(7) * (t107 * t91 + t112 * t70 + t20 * t39 + t21 * t38 + t66 * t9 + t67 * t8) + t170 * mrSges(5,3) + (mrSges(7,3) + mrSges(6,3)) * (-t125 * t9 - t180 * t8 + t38 * t83 - t39 * t84) + (-pkin(3) * t91 + pkin(9) * t170) * m(5); (t79 / 0.2e1 + t77 / 0.2e1 + t80 / 0.2e1 + t78 / 0.2e1 - t157 / 0.2e1 + (pkin(8) * t215 + Ifges(4,5)) * qJD(3)) * t168 + (-t11 * t180 - t12 * t125 + t35 * t83 - t36 * t84) * mrSges(6,3) + (-t125 * t3 - t180 * t4 + t22 * t83 - t23 * t84) * mrSges(7,3) + (t167 * t132 / 0.2e1 + t131 * t228 - Ifges(4,6) * qJD(3) + (t142 * t228 + t167 * t229) * qJD(4) + (qJD(3) * mrSges(4,2) + t129) * pkin(8) + (t240 * t125 - t180 * t239 + t182) * qJD(3) / 0.2e1) * t164 - (t14 / 0.2e1 + t13 / 0.2e1) * t180 + m(5) * (-pkin(3) * t158 + (-t106 * t204 - t58 * t163) * pkin(9)) + m(7) * (t107 * t24 + t20 * t23 + t21 * t22 + t3 * t66 + t4 * t67 + t70 * t81) - (t60 / 0.2e1 + t61 / 0.2e1) * t84 - (t62 / 0.2e1 + t63 / 0.2e1) * t83 + t153 * t18 + (t87 / 0.2e1 + t88 / 0.2e1) * t53 + t138 * t41 + (t89 / 0.2e1 + t90 / 0.2e1) * t52 + t20 * t98 + t55 * t99 + t21 * t100 + t56 * t101 + t104 * t86 + t107 * t17 + t24 * t85 + t96 * t32 + t97 * t34 + t81 * t40 + t70 * t68 - pkin(3) * t74 + t66 * t31 + t67 * t33 + (t205 * t229 + t65 / 0.2e1 - t58 * mrSges(5,3) - pkin(9) * t102 + (-pkin(9) * t133 - t106 * mrSges(5,3) + pkin(4) * t69 + t218 / 0.2e1 - t108 / 0.2e1) * qJD(4)) * t163 + (t142 * t205 / 0.2e1 + t64 / 0.2e1 + qJD(4) * t109 / 0.2e1 + t187 * mrSges(5,3) + (m(5) * t187 - qJD(4) * t134 + t103) * pkin(9)) * t167 + m(6) * (t104 * t153 + t11 * t97 + t12 * t96 + t138 * t199 + t35 * t56 + t36 * t55) + (t15 / 0.2e1 + t16 / 0.2e1) * t125 - (t43 / 0.2e1 + t42 / 0.2e1) * t110 - (t44 / 0.2e1 + t45 / 0.2e1) * t111; -0.2e1 * pkin(3) * t129 + 0.2e1 * t107 * t40 + t167 * t131 + t163 * t132 + 0.2e1 * t153 * t41 + 0.2e1 * t70 * t85 - (t87 + t88) * t84 - (t89 + t90) * t83 + (t44 + t45) * t125 - (t43 + t42) * t180 + (t167 * t142 + (0.2e1 * pkin(4) * t86 - t141) * t163) * qJD(4) + (t153 * t199 + t55 * t97 + t56 * t96) * t234 + (t107 * t70 + t20 * t67 + t21 * t66) * t233 + 0.2e1 * (-t125 * t21 - t180 * t20 + t66 * t83 - t67 * t84) * mrSges(7,3) + 0.2e1 * (-t125 * t56 - t180 * t55 + t83 * t96 - t84 * t97) * mrSges(6,3); t25 * mrSges(5,1) - t26 * mrSges(5,2) + m(6) * ((t166 * t9 - t197) * pkin(4) + t227) + m(7) * (-pkin(4) * t197 + t152 * t9 + t227) + t181; t198 * t152 + (t166 * t32 + t223 * t162 + (-t162 * t209 + t166 * t222) * qJD(5) + m(6) * (t11 * t162 + t12 * t166 + t200 * t36 - t201 * t35) + m(7) * (t162 * t4 + t200 * t23 - t201 * t22)) * pkin(4) + t171 - t175 * Ifges(5,6) - Ifges(5,5) * t193 - t57 * mrSges(5,2) + t58 * mrSges(5,1) - t237; t157 + t186 * t152 + (pkin(9) * t140 - t219) * qJD(4) + (m(6) * (t162 * t55 + t166 * t56 + t200 * t97 - t201 * t96) + m(7) * (t162 * t20 + t200 * t67 - t201 * t66) + t172 * mrSges(7,3) + (t166 * t83 + t172) * mrSges(6,3)) * pkin(4) + t173; 0.2e1 * (t189 + ((-t152 + t242) * m(7) - t244) * t162) * t238; t231 * t9 + t181; pkin(5) * t198 + t171; pkin(5) * t186 + t173; (t189 + (-t231 - t244) * t162) * t238; 0; m(7) * t91; m(7) * t24 + t17; m(7) * t70 + t40; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
