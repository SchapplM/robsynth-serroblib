% Calculate time derivative of joint inertia matrix for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:29
% EndTime: 2019-03-08 22:54:34
% DurationCPUTime: 3.31s
% Computational Cost: add. (1680->430), mult. (4594->571), div. (0->0), fcn. (3480->8), ass. (0->188)
t237 = -Ifges(7,4) - Ifges(6,5);
t210 = -mrSges(6,1) - mrSges(7,1);
t209 = Ifges(5,5) + Ifges(7,5);
t245 = -Ifges(6,4) + t209;
t128 = sin(qJ(4));
t129 = sin(qJ(3));
t131 = cos(qJ(4));
t181 = qJD(4) * t131;
t132 = cos(qJ(3));
t184 = qJD(3) * t132;
t138 = t128 * t184 + t129 * t181;
t84 = -t131 * mrSges(5,1) + t128 * mrSges(5,2);
t244 = -m(5) * pkin(3) - mrSges(4,1) + t84;
t126 = cos(pkin(6));
t125 = sin(pkin(6));
t130 = sin(qJ(2));
t192 = t125 * t130;
t141 = t126 * t129 + t132 * t192;
t136 = t141 * qJD(3);
t213 = cos(qJ(2));
t169 = qJD(2) * t213;
t166 = t125 * t169;
t25 = t129 * t166 + t136;
t197 = t129 * t25;
t50 = -t126 * t132 + t129 * t192;
t26 = -t50 * qJD(3) + t132 * t166;
t243 = -t129 * t136 + t132 * t26 + t197;
t242 = -Ifges(6,1) - Ifges(7,1) - Ifges(5,3);
t194 = qJ(5) * t128;
t146 = -pkin(4) * t131 - t194;
t78 = -pkin(3) + t146;
t83 = t131 * mrSges(6,2) - t128 * mrSges(6,3);
t241 = m(6) * t78 + t83;
t183 = qJD(4) * t128;
t185 = qJD(3) * t129;
t240 = pkin(8) * (t131 * t185 + t132 * t183);
t239 = pkin(8) * t184;
t238 = mrSges(6,3) + mrSges(7,2);
t34 = -mrSges(5,2) * t185 - mrSges(5,3) * t138;
t36 = mrSges(6,1) * t138 - mrSges(6,3) * t185;
t235 = t34 - t36;
t170 = t131 * t184;
t182 = qJD(4) * t129;
t173 = t128 * t182;
t140 = t170 - t173;
t33 = mrSges(5,1) * t185 - mrSges(5,3) * t140;
t195 = mrSges(6,1) * t170 + mrSges(6,2) * t185;
t38 = -mrSges(6,1) * t173 + t195;
t234 = t38 - t33;
t199 = Ifges(7,6) * t128;
t147 = Ifges(7,3) * t131 + t199;
t206 = Ifges(5,4) * t128;
t156 = Ifges(5,1) * t131 - t206;
t233 = (t147 + t156) * qJD(4);
t200 = Ifges(6,6) * t131;
t150 = Ifges(6,3) * t128 - t200;
t198 = Ifges(7,6) * t131;
t151 = Ifges(7,2) * t128 + t198;
t232 = (t150 + t151) * qJD(4);
t189 = t129 * t131;
t70 = -mrSges(5,1) * t132 - mrSges(5,3) * t189;
t74 = mrSges(6,1) * t189 - mrSges(6,2) * t132;
t231 = -t70 + t74;
t191 = t128 * t129;
t72 = mrSges(6,1) * t191 + mrSges(6,3) * t132;
t73 = -mrSges(7,1) * t191 - mrSges(7,2) * t132;
t230 = -t72 + t73;
t175 = t125 * t213;
t229 = -qJD(4) * t175 + t26;
t174 = qJD(2) * t192;
t227 = qJD(4) * t141 - t174;
t37 = -mrSges(7,1) * t138 + mrSges(7,2) * t185;
t190 = t128 * t132;
t109 = pkin(8) * t190;
t77 = (pkin(3) * t129 - pkin(9) * t132) * qJD(3);
t82 = -pkin(3) * t132 - t129 * pkin(9) - pkin(2);
t208 = t128 * t77 + t82 * t181;
t161 = -qJD(5) * t132 + t208;
t6 = (-pkin(5) * t189 - t109) * qJD(4) + (-pkin(5) * t190 + (-pkin(8) * t131 + qJ(5)) * t129) * qJD(3) + t161;
t225 = m(7) * t6 + t37;
t188 = t131 * t132;
t110 = pkin(8) * t188;
t162 = qJD(4) * t110 - t131 * t77 + t82 * t183;
t212 = pkin(8) * t128;
t176 = -pkin(4) - t212;
t3 = -pkin(5) * t173 + qJD(6) * t132 + (pkin(5) * t188 + (-qJ(6) + t176) * t129) * qJD(3) + t162;
t9 = t176 * t185 + t162;
t224 = m(6) * t9 + m(7) * t3;
t223 = 0.2e1 * m(5);
t222 = 0.2e1 * m(7);
t221 = 0.2e1 * pkin(8);
t220 = 0.2e1 * mrSges(7,1);
t219 = m(6) / 0.2e1;
t217 = pkin(5) + pkin(9);
t13 = t50 * t25;
t127 = -pkin(4) - qJ(6);
t205 = Ifges(5,4) * t131;
t204 = Ifges(6,4) * t132;
t203 = Ifges(5,6) * t131;
t202 = Ifges(5,6) * t132;
t201 = Ifges(6,6) * t128;
t40 = t128 * t82 + t110;
t193 = qJD(3) * mrSges(7,3);
t187 = pkin(4) * t191 + t129 * pkin(8);
t186 = qJ(5) * qJD(5);
t180 = qJD(6) * t128;
t43 = -Ifges(5,5) * t132 + t129 * t156;
t44 = -Ifges(7,5) * t132 + t129 * t147;
t154 = -Ifges(6,2) * t131 + t201;
t47 = t129 * t154 - t204;
t179 = t47 - t43 - t44;
t93 = t217 * t131;
t178 = t50 * t184;
t177 = mrSges(7,1) * t183;
t39 = t131 * t82 - t109;
t168 = pkin(4) * t183 - qJD(5) * t128;
t167 = pkin(4) * t138 + qJ(5) * t173 + t239;
t148 = -Ifges(7,3) * t128 + t198;
t153 = Ifges(6,2) * t128 + t200;
t91 = Ifges(5,1) * t128 + t205;
t164 = -t148 / 0.2e1 + t153 / 0.2e1 + t91 / 0.2e1;
t149 = Ifges(6,3) * t131 + t201;
t152 = Ifges(7,2) * t131 - t199;
t90 = Ifges(5,2) * t131 + t206;
t163 = -t149 / 0.2e1 - t152 / 0.2e1 - t90 / 0.2e1;
t31 = qJ(5) * t132 - t40;
t4 = t128 * t229 + t131 * t227;
t5 = -t128 * t227 + t131 * t229;
t160 = t4 * t128 + t5 * t131;
t159 = mrSges(5,1) * t128 + mrSges(5,2) * t131;
t158 = -mrSges(6,2) * t128 - mrSges(6,3) * t131;
t157 = -mrSges(7,2) * t131 + mrSges(7,3) * t128;
t155 = -Ifges(5,2) * t128 + t205;
t28 = -t128 * t175 + t131 * t141;
t145 = qJ(5) * t5 + qJD(5) * t28;
t42 = t129 * t155 - t202;
t45 = -Ifges(6,5) * t132 + t129 * t150;
t46 = -Ifges(7,4) * t132 + t129 * t151;
t144 = -t42 + t45 + t46 + t202;
t143 = -qJ(5) * t131 + qJ(6) * t128;
t142 = t131 * (m(6) * pkin(9) - t210);
t135 = -Ifges(6,4) * t173 + t237 * t138 - t209 * t170 + t242 * t185;
t124 = t132 * pkin(4);
t118 = Ifges(7,4) * t183;
t117 = Ifges(5,5) * t181;
t116 = Ifges(6,5) * t183;
t115 = Ifges(7,5) * t181;
t96 = mrSges(7,1) * t170;
t92 = t217 * t128;
t85 = -mrSges(7,2) * t128 - mrSges(7,3) * t131;
t76 = qJD(4) * t93;
t75 = t217 * t183;
t71 = mrSges(7,1) * t189 + mrSges(7,3) * t132;
t69 = mrSges(5,2) * t132 - mrSges(5,3) * t191;
t67 = t155 * qJD(4);
t66 = t154 * qJD(4);
t62 = (mrSges(4,1) * t129 + mrSges(4,2) * t132) * qJD(3);
t61 = t159 * qJD(4);
t60 = t158 * qJD(4);
t59 = t157 * qJD(4);
t55 = t159 * t129;
t54 = t158 * t129;
t53 = t157 * t129;
t52 = t127 * t131 - pkin(3) - t194;
t49 = -qJ(5) * t181 + t168;
t48 = -qJ(5) * t189 + t187;
t35 = t96 + (-t177 - t193) * t129;
t32 = t124 - t39;
t30 = t129 * t143 + t187;
t29 = qJD(4) * t143 - qJD(6) * t131 + t168;
t27 = t128 * t141 + t131 * t175;
t24 = -pkin(5) * t191 - t31;
t23 = -mrSges(7,2) * t140 + mrSges(7,3) * t138;
t22 = mrSges(5,1) * t138 + mrSges(5,2) * t140;
t21 = -mrSges(6,2) * t138 - mrSges(6,3) * t140;
t20 = qJ(6) * t132 + t109 + t124 + (pkin(5) * t129 - t82) * t131;
t19 = t153 * t182 + (Ifges(6,4) * t129 + t132 * t154) * qJD(3);
t18 = t152 * t182 + (Ifges(7,4) * t129 + t132 * t151) * qJD(3);
t17 = t149 * t182 + (Ifges(6,5) * t129 + t132 * t150) * qJD(3);
t16 = t148 * t182 + (Ifges(7,5) * t129 + t132 * t147) * qJD(3);
t15 = -t91 * t182 + (Ifges(5,5) * t129 + t132 * t156) * qJD(3);
t14 = -t90 * t182 + (Ifges(5,6) * t129 + t132 * t155) * qJD(3);
t12 = t185 * t212 - t162;
t11 = t208 - t240;
t10 = (-qJ(5) * t184 - qJD(5) * t129) * t131 + t167;
t8 = -qJ(5) * t185 - t161 + t240;
t7 = t143 * t184 + (t180 + (qJ(6) * qJD(4) - qJD(5)) * t131) * t129 + t167;
t1 = [0.2e1 * m(4) * (-t125 ^ 2 * t130 * t169 + t141 * t26 + t13) + 0.2e1 * (m(6) + m(7) + m(5)) * (t27 * t4 + t28 * t5 + t13); -mrSges(3,2) * t166 - t62 * t175 + (m(4) * t243 + m(5) * (t178 + t197)) * pkin(8) + (-m(4) * pkin(2) - mrSges(4,1) * t132 + t129 * mrSges(4,2) - mrSges(3,1)) * t174 + (m(4) * t239 + m(6) * t10 + m(7) * t7 + t21 + t22 + t23) * t50 + (m(5) * t40 - m(6) * t31 + m(7) * t24 + t230 + t69) * t5 + (-m(5) * t39 + m(6) * t32 + m(7) * t20 + t231 + t71) * t4 + (m(5) * t11 - m(6) * t8 + t225 + t235) * t28 + (-m(5) * t12 + t224 + t234 + t35) * t27 + (m(6) * t48 + m(7) * t30 + t53 + t54 + t55) * t25 + (t178 + t243) * mrSges(4,3); -0.2e1 * pkin(2) * t62 + 0.2e1 * t11 * t69 + 0.2e1 * t12 * t70 + 0.2e1 * t3 * t71 + 0.2e1 * t8 * t72 + 0.2e1 * t6 * t73 + 0.2e1 * t9 * t74 + 0.2e1 * t7 * t53 + 0.2e1 * t10 * t54 + 0.2e1 * t20 * t35 + 0.2e1 * t31 * t36 + 0.2e1 * t24 * t37 + 0.2e1 * t32 * t38 + 0.2e1 * t39 * t33 + 0.2e1 * t40 * t34 + 0.2e1 * t48 * t21 + 0.2e1 * t30 * t23 + 0.2e1 * m(6) * (t10 * t48 + t31 * t8 + t32 * t9) + (t20 * t3 + t24 * t6 + t30 * t7) * t222 + (t40 * t11 + t39 * t12) * t223 + ((0.2e1 * Ifges(4,4) * t132 + t55 * t221 + (-t179 + t204) * t131 + t144 * t128) * qJD(3) + t135) * t132 + (t22 * t221 + (t15 + t16 - t19) * t131 + (-t14 + t17 + t18) * t128 + (t144 * t131 + (t132 * t209 + t179) * t128) * qJD(4) + ((-0.2e1 * Ifges(4,4) + t245 * t131 + (-Ifges(5,6) - t237) * t128) * t129 + (pkin(8) ^ 2 * t223 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t242) * t132) * qJD(3)) * t129; -t26 * mrSges(4,2) + m(7) * (t27 * t76 - t28 * t75 + t4 * t92 + t5 * t93) + 0.2e1 * (t219 + m(5) / 0.2e1) * (t181 * t27 - t183 * t28 + t160) * pkin(9) + (m(6) * t49 + m(7) * t29 + t59 + t60 + t61) * t50 + (m(7) * t52 + t241 + t244 + t85) * t25 + (mrSges(5,3) - t210) * ((-t128 * t28 + t131 * t27) * qJD(4) + t160); -pkin(3) * t22 + t10 * t83 + t78 * t21 + t52 * t23 + t29 * t53 + t30 * t59 + t92 * t35 + t93 * t37 + t48 * t60 + t49 * t54 + t7 * t85 + t76 * t71 - t75 * t73 + m(7) * (t20 * t76 - t24 * t75 + t29 * t30 + t3 * t92 + t52 * t7 + t6 * t93) + m(6) * (t10 * t78 + t48 * t49) + (-t116 / 0.2e1 - t118 / 0.2e1 - t115 / 0.2e1 - t117 / 0.2e1 + (pkin(8) * t244 + Ifges(4,5)) * qJD(3)) * t132 + (t14 / 0.2e1 - t17 / 0.2e1 - t18 / 0.2e1 + t6 * mrSges(7,1) - t8 * mrSges(6,1) + t11 * mrSges(5,3) + t164 * t184 + (t204 / 0.2e1 + t43 / 0.2e1 + t44 / 0.2e1 - t47 / 0.2e1 + t32 * mrSges(6,1) - t39 * mrSges(5,3) + t20 * mrSges(7,1)) * qJD(4) + (t231 * qJD(4) + m(6) * (qJD(4) * t32 - t8) + m(5) * (-qJD(4) * t39 + t11) + t235) * pkin(9)) * t131 + (t15 / 0.2e1 + t16 / 0.2e1 - t19 / 0.2e1 + t9 * mrSges(6,1) - t12 * mrSges(5,3) + t3 * mrSges(7,1) + t163 * t184 + (t202 / 0.2e1 - t42 / 0.2e1 + t45 / 0.2e1 + t46 / 0.2e1 - t40 * mrSges(5,3) + t31 * mrSges(6,1) - t24 * mrSges(7,1)) * qJD(4) + ((-t69 + t72) * qJD(4) + m(6) * (qJD(4) * t31 + t9) + m(5) * (-qJD(4) * t40 - t12) + t234) * pkin(9)) * t128 + (-t131 * t66 / 0.2e1 - t128 * t67 / 0.2e1 - Ifges(4,6) * qJD(3) + (qJD(3) * mrSges(4,2) + t61) * pkin(8) + (-t128 * t164 + t131 * t163) * qJD(4) + t232 * t128 / 0.2e1 + t233 * t131 / 0.2e1 + (t245 * t128 + t237 * t131 + t203) * qJD(3) / 0.2e1) * t129; -0.2e1 * pkin(3) * t61 + 0.2e1 * t78 * t60 + 0.2e1 * t29 * t85 + 0.2e1 * t52 * t59 + (t29 * t52 - t75 * t93 + t76 * t92) * t222 + 0.2e1 * t241 * t49 + (-t75 * t220 - t232 + t67) * t131 + (t76 * t220 + t233 - t66) * t128 + ((t220 * t92 - t148 + t153 + t91) * t131 + (-0.2e1 * mrSges(7,1) * t93 - t149 - t152 - t90) * t128) * qJD(4); (-mrSges(5,2) + t238) * t5 + (-mrSges(5,1) + mrSges(6,2) - mrSges(7,3)) * t4 + m(6) * (-pkin(4) * t4 + t145) + m(7) * (-qJD(6) * t27 + t127 * t4 + t145); t127 * t35 - qJD(6) * t71 - pkin(4) * t38 - t8 * mrSges(6,3) + t9 * mrSges(6,2) - t11 * mrSges(5,2) + t12 * mrSges(5,1) + t6 * mrSges(7,2) - t3 * mrSges(7,3) - t135 + t230 * qJD(5) + (-t36 + t37) * qJ(5) + m(7) * (qJ(5) * t6 + qJD(5) * t24 - qJD(6) * t20 + t127 * t3) + m(6) * (-pkin(4) * t9 - qJ(5) * t8 - qJD(5) * t31) + (-Ifges(6,4) * t131 - Ifges(5,6) * t128) * t184 + (-t209 * t128 - t203) * t182; -mrSges(7,1) * t180 + m(7) * (-qJ(5) * t75 + qJD(5) * t93 - qJD(6) * t92 + t127 * t76) + t117 + t118 + t115 + t116 - t75 * mrSges(7,2) - t76 * mrSges(7,3) + qJD(5) * t142 + ((-pkin(4) * mrSges(6,1) + t127 * mrSges(7,1) - Ifges(6,4)) * t131 + (qJ(5) * t210 - Ifges(5,6)) * t128 + (m(6) * t146 + t83 + t84) * pkin(9)) * qJD(4); 0.2e1 * m(6) * t186 + 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (-qJD(6) * t127 + t186) + 0.2e1 * t238 * qJD(5); 0.2e1 * (t219 + m(7) / 0.2e1) * t4; t96 + (t183 * t210 - t193) * t129 + t195 + t224; m(7) * t76 + qJD(4) * t142; -m(7) * qJD(6); 0; m(7) * t5; t225; -m(7) * t75 - t177; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
