% Calculate time derivative of joint inertia matrix for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:30
% EndTime: 2019-03-08 19:03:37
% DurationCPUTime: 3.38s
% Computational Cost: add. (4170->423), mult. (12065->650), div. (0->0), fcn. (12331->14), ass. (0->187)
t131 = sin(qJ(5));
t135 = cos(qJ(5));
t109 = -mrSges(6,1) * t135 + mrSges(6,2) * t131;
t221 = -m(6) * pkin(4) - mrSges(5,1) + t109;
t124 = sin(pkin(13));
t126 = sin(pkin(6));
t129 = cos(pkin(6));
t133 = sin(qJ(3));
t137 = cos(qJ(3));
t127 = cos(pkin(13));
t128 = cos(pkin(7));
t180 = t127 * t128;
t125 = sin(pkin(7));
t181 = t125 * t137;
t220 = (-t124 * t133 + t137 * t180) * t126 + t129 * t181;
t132 = sin(qJ(4));
t136 = cos(qJ(4));
t142 = -t125 * t126 * t127 + t128 * t129;
t182 = t125 * t133;
t61 = t129 * t182 + (t124 * t137 + t133 * t180) * t126;
t45 = t142 * t132 + t61 * t136;
t54 = t220 * qJD(3);
t19 = t45 * qJD(4) + t54 * t132;
t44 = t61 * t132 - t142 * t136;
t175 = qJD(3) * t137;
t158 = t125 * t175;
t87 = t128 * t132 + t136 * t182;
t68 = t87 * qJD(4) + t132 * t158;
t86 = -t136 * t128 + t132 * t182;
t197 = t86 * t19 + t68 * t44;
t219 = m(6) * pkin(10) + mrSges(6,3);
t130 = sin(qJ(6));
t134 = cos(qJ(6));
t148 = t130 * t131 - t134 * t135;
t85 = t148 * t132;
t108 = -pkin(4) * t136 - pkin(10) * t132 - pkin(3);
t177 = t135 * t136;
t118 = pkin(9) * t177;
t80 = t131 * t108 + t118;
t218 = t80 * qJD(5);
t173 = qJD(4) * t136;
t156 = t135 * t173;
t174 = qJD(4) * t132;
t217 = -Ifges(6,5) * t156 - Ifges(6,3) * t174;
t215 = qJD(5) + qJD(6);
t120 = -pkin(5) * t135 - pkin(4);
t96 = t130 * t135 + t131 * t134;
t64 = mrSges(7,1) * t148 + mrSges(7,2) * t96;
t214 = m(7) * t120 + t221 + t64;
t213 = 0.2e1 * m(6);
t212 = 0.2e1 * m(7);
t211 = 0.2e1 * pkin(9);
t210 = m(5) / 0.2e1;
t209 = m(6) / 0.2e1;
t208 = m(5) * pkin(3);
t205 = m(7) * pkin(5);
t204 = -pkin(11) - pkin(10);
t194 = Ifges(6,4) * t131;
t111 = Ifges(6,2) * t135 + t194;
t203 = -t111 / 0.2e1;
t202 = -t131 / 0.2e1;
t201 = pkin(9) * t131;
t14 = t44 * t19;
t55 = t61 * qJD(3);
t200 = t55 * t220;
t48 = t86 * t68;
t63 = t215 * t96;
t39 = -t63 * t132 - t148 * t173;
t40 = -t96 * t173 + t215 * t85;
t17 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t170 = qJD(5) * t135;
t140 = t131 * t173 + t132 * t170;
t171 = qJD(5) * t132;
t157 = t131 * t171;
t141 = t156 - t157;
t56 = t140 * mrSges(6,1) + t141 * mrSges(6,2);
t198 = t17 + t56;
t62 = t215 * t148;
t33 = mrSges(7,1) * t63 - mrSges(7,2) * t62;
t153 = mrSges(6,1) * t131 + mrSges(6,2) * t135;
t97 = t153 * qJD(5);
t196 = t33 + t97;
t84 = t96 * t132;
t53 = mrSges(7,1) * t84 - mrSges(7,2) * t85;
t90 = t153 * t132;
t195 = t53 + t90;
t193 = Ifges(6,4) * t135;
t192 = Ifges(6,6) * t131;
t191 = t136 * Ifges(6,6);
t190 = t19 * t132;
t20 = -t44 * qJD(4) + t54 * t136;
t189 = t20 * t136;
t67 = -t86 * qJD(4) + t136 * t158;
t188 = t67 * t136;
t187 = t68 * t132;
t186 = -mrSges(5,1) * t136 + mrSges(5,2) * t132 - mrSges(4,1);
t106 = (pkin(4) * t132 - pkin(10) * t136) * qJD(4);
t184 = t135 * t106 + t174 * t201;
t179 = t131 * t132;
t178 = t132 * t135;
t176 = qJD(3) * t133;
t172 = qJD(5) * t131;
t169 = qJD(6) * t130;
t168 = qJD(6) * t134;
t166 = -Ifges(7,5) * t39 - Ifges(7,6) * t40 - Ifges(7,3) * t174;
t165 = pkin(5) * t172;
t21 = -t131 * t45 - t135 * t220;
t22 = -t131 * t220 + t135 * t45;
t12 = -t130 * t22 + t134 * t21;
t4 = -t22 * qJD(5) - t131 * t20 + t135 * t55;
t5 = t21 * qJD(5) + t131 * t55 + t135 * t20;
t2 = t12 * qJD(6) + t130 * t4 + t134 * t5;
t13 = t130 * t21 + t134 * t22;
t3 = -t13 * qJD(6) - t130 * t5 + t134 * t4;
t164 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t159 = t125 * t176;
t69 = -t131 * t87 - t135 * t181;
t24 = t69 * qJD(5) + t131 * t159 + t135 * t67;
t143 = t131 * t181 - t135 * t87;
t25 = t143 * qJD(5) - t131 * t67 + t135 * t159;
t31 = t130 * t143 + t134 * t69;
t7 = t31 * qJD(6) + t130 * t25 + t134 * t24;
t32 = t130 * t69 - t134 * t143;
t8 = -t32 * qJD(6) - t130 * t24 + t134 * t25;
t163 = t8 * mrSges(7,1) - t7 * mrSges(7,2);
t160 = qJD(5) * t204;
t155 = (2 * Ifges(5,4)) + t192;
t46 = t131 * t106 + t108 * t170 + (-t135 * t174 - t136 * t172) * pkin(9);
t94 = t135 * t108;
t79 = -t136 * t201 + t94;
t154 = -t79 * qJD(5) + t46;
t152 = Ifges(6,1) * t135 - t194;
t112 = Ifges(6,1) * t131 + t193;
t151 = -Ifges(6,2) * t131 + t193;
t150 = Ifges(6,5) * t131 + Ifges(6,6) * t135;
t57 = -pkin(11) * t178 + t94 + (-pkin(5) - t201) * t136;
t71 = -pkin(11) * t179 + t80;
t29 = -t130 * t71 + t134 * t57;
t30 = t130 * t57 + t134 * t71;
t104 = t131 * t160;
t105 = t135 * t160;
t114 = t204 * t131;
t115 = t204 * t135;
t72 = t114 * t134 + t115 * t130;
t42 = t72 * qJD(6) + t104 * t134 + t105 * t130;
t73 = t114 * t130 - t115 * t134;
t43 = -t73 * qJD(6) - t104 * t130 + t105 * t134;
t58 = Ifges(7,6) * t63;
t59 = Ifges(7,5) * t62;
t149 = t43 * mrSges(7,1) - t42 * mrSges(7,2) - t58 - t59;
t26 = (pkin(5) * t132 - pkin(11) * t177) * qJD(4) + (-t118 + (pkin(11) * t132 - t108) * t131) * qJD(5) + t184;
t36 = -t140 * pkin(11) + t46;
t10 = t29 * qJD(6) + t130 * t26 + t134 * t36;
t11 = -t30 * qJD(6) - t130 * t36 + t134 * t26;
t147 = t11 * mrSges(7,1) - t10 * mrSges(7,2) - t166;
t146 = -t137 * t55 - t176 * t220;
t145 = t44 * t173 + t190;
t144 = t86 * t173 + t187;
t123 = Ifges(6,5) * t170;
t107 = (pkin(5) * t131 + pkin(9)) * t132;
t103 = -mrSges(6,1) * t136 - mrSges(6,3) * t178;
t102 = mrSges(6,2) * t136 - mrSges(6,3) * t179;
t100 = t152 * qJD(5);
t99 = t151 * qJD(5);
t98 = (mrSges(5,1) * t132 + mrSges(5,2) * t136) * qJD(4);
t92 = (-mrSges(7,1) * t130 - mrSges(7,2) * t134) * qJD(6) * pkin(5);
t83 = -Ifges(6,5) * t136 + t152 * t132;
t82 = t151 * t132 - t191;
t78 = t140 * pkin(5) + pkin(9) * t173;
t77 = -mrSges(6,2) * t174 - t140 * mrSges(6,3);
t76 = mrSges(6,1) * t174 - t141 * mrSges(6,3);
t75 = -mrSges(7,1) * t136 + mrSges(7,3) * t85;
t74 = mrSges(7,2) * t136 - mrSges(7,3) * t84;
t66 = Ifges(7,1) * t96 - Ifges(7,4) * t148;
t65 = Ifges(7,4) * t96 - Ifges(7,2) * t148;
t52 = -t112 * t171 + (Ifges(6,5) * t132 + t152 * t136) * qJD(4);
t51 = -t111 * t171 + (Ifges(6,6) * t132 + t151 * t136) * qJD(4);
t50 = -Ifges(7,1) * t85 - Ifges(7,4) * t84 - Ifges(7,5) * t136;
t49 = -Ifges(7,4) * t85 - Ifges(7,2) * t84 - Ifges(7,6) * t136;
t47 = t184 - t218;
t35 = -Ifges(7,1) * t62 - Ifges(7,4) * t63;
t34 = -Ifges(7,4) * t62 - Ifges(7,2) * t63;
t28 = -mrSges(7,2) * t174 + mrSges(7,3) * t40;
t27 = mrSges(7,1) * t174 - mrSges(7,3) * t39;
t16 = Ifges(7,1) * t39 + Ifges(7,4) * t40 + Ifges(7,5) * t174;
t15 = Ifges(7,4) * t39 + Ifges(7,2) * t40 + Ifges(7,6) * t174;
t1 = [0.2e1 * m(7) * (t12 * t3 + t13 * t2 + t14) + 0.2e1 * m(6) * (t21 * t4 + t22 * t5 + t14) + 0.2e1 * m(5) * (t20 * t45 + t14 - t200) + 0.2e1 * m(4) * (t54 * t61 - t200); m(7) * (t12 * t8 + t13 * t7 + t2 * t32 + t3 * t31 + t197) + m(6) * (-t143 * t5 + t21 * t25 + t22 * t24 + t4 * t69 + t197) + m(5) * (t87 * t20 + t67 * t45 + t197) + 0.2e1 * (t146 * t210 + m(4) * (t133 * t54 + t61 * t175 + t146) / 0.2e1) * t125; 0.2e1 * m(7) * (t31 * t8 + t32 * t7 + t48) + 0.2e1 * m(6) * (-t143 * t24 + t25 * t69 + t48) + 0.2e1 * m(5) * (-t125 ^ 2 * t133 * t175 + t87 * t67 + t48); -t54 * mrSges(4,2) + t5 * t102 + t4 * t103 + t12 * t27 + t13 * t28 + t2 * t74 + t21 * t76 + t22 * t77 + t3 * t75 - t220 * t98 + t198 * t44 + t195 * t19 + m(7) * (t10 * t13 + t107 * t19 + t11 * t12 + t2 * t30 + t29 * t3 + t44 * t78) + m(6) * (t21 * t47 + t22 * t46 + t4 * t79 + t5 * t80) + (t145 * t209 + (-t45 * t174 + t145 + t189) * t210) * t211 + (t190 + t189 + (-t132 * t45 + t136 * t44) * qJD(4)) * mrSges(5,3) + (t186 - t208) * t55; t24 * t102 + t25 * t103 + t31 * t27 + t32 * t28 + t69 * t76 + t7 * t74 - t143 * t77 + t8 * t75 + t198 * t86 + t195 * t68 + (-t137 * t98 + (-mrSges(4,2) * t137 + t186 * t133) * qJD(3)) * t125 + m(7) * (t10 * t32 + t107 * t68 + t11 * t31 + t29 * t8 + t30 * t7 + t78 * t86) + m(6) * (-t143 * t46 + t24 * t80 + t25 * t79 + t47 * t69) - t159 * t208 + (t144 * t209 + (-t87 * t174 + t144 + t188) * t210) * t211 + (t187 + t188 + (-t132 * t87 + t136 * t86) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t98 + 0.2e1 * t10 * t74 + 0.2e1 * t46 * t102 + 0.2e1 * t47 * t103 + 0.2e1 * t107 * t17 + 0.2e1 * t11 * t75 - t84 * t15 - t85 * t16 + 0.2e1 * t29 * t27 + 0.2e1 * t30 * t28 + t39 * t50 + t40 * t49 + 0.2e1 * t78 * t53 + 0.2e1 * t79 * t76 + 0.2e1 * t80 * t77 + (t46 * t80 + t47 * t79) * t213 + (t10 * t30 + t107 * t78 + t11 * t29) * t212 + ((-t131 * t82 + t135 * t83 + t155 * t136 + t90 * t211) * qJD(4) + t166 + t217) * t136 + (t56 * t211 - t131 * t51 + t135 * t52 + (-t131 * t83 - t135 * t82 + t136 * t150) * qJD(5) + (-Ifges(7,5) * t85 - Ifges(7,6) * t84 + (Ifges(6,5) * t135 - t155) * t132 + (pkin(9) ^ 2 * t213 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) - Ifges(6,3) - Ifges(7,3)) * t136) * qJD(4)) * t132; -t20 * mrSges(5,2) + t196 * t44 + m(7) * (t12 * t43 + t13 * t42 + t44 * t165 + t2 * t73 + t3 * t72) + (t12 * t62 - t13 * t63 - t148 * t2 - t3 * t96) * mrSges(7,3) + t214 * t19 + t219 * (-t4 * t131 + t5 * t135 + (-t131 * t22 - t135 * t21) * qJD(5)); -t67 * mrSges(5,2) + t196 * t86 + m(7) * (t86 * t165 + t31 * t43 + t32 * t42 + t7 * t73 + t72 * t8) + (-t148 * t7 + t31 * t62 - t32 * t63 - t8 * t96) * mrSges(7,3) + t214 * t68 + t219 * (-t25 * t131 + t24 * t135 + (t131 * t143 - t135 * t69) * qJD(5)); t120 * t17 + t107 * t33 - t148 * t15 / 0.2e1 + t96 * t16 / 0.2e1 - t84 * t34 / 0.2e1 - t85 * t35 / 0.2e1 + t73 * t28 + t42 * t74 + t43 * t75 + t78 * t64 - t63 * t49 / 0.2e1 + t40 * t65 / 0.2e1 + t39 * t66 / 0.2e1 + t72 * t27 - pkin(4) * t56 - t62 * t50 / 0.2e1 + m(7) * (t10 * t73 + t11 * t72 + t120 * t78 + t29 * t43 + t30 * t42) + (-t123 / 0.2e1 + t59 / 0.2e1 + t58 / 0.2e1 + (t221 * pkin(9) + Ifges(5,5)) * qJD(4)) * t136 + (t173 * t203 + t52 / 0.2e1 - t47 * mrSges(6,3) + (-t80 * mrSges(6,3) + t191 / 0.2e1 - t82 / 0.2e1 + (m(7) * t107 + t53) * pkin(5)) * qJD(5) + (-qJD(5) * t102 + m(6) * (-t47 - t218) - t76) * pkin(10)) * t131 + (-t10 * t148 - t11 * t96 + t29 * t62 - t30 * t63) * mrSges(7,3) + (t112 * t173 / 0.2e1 + t51 / 0.2e1 + qJD(5) * t83 / 0.2e1 + t154 * mrSges(6,3) + (m(6) * t154 - qJD(5) * t103 + t77) * pkin(10)) * t135 + (t135 * t100 / 0.2e1 + t99 * t202 - Ifges(5,6) * qJD(4) + (t112 * t202 + t135 * t203) * qJD(5) + (qJD(4) * mrSges(5,2) + t97) * pkin(9) + (Ifges(7,5) * t96 - Ifges(7,6) * t148 + t150) * qJD(4) / 0.2e1) * t132; 0.2e1 * t64 * t165 + 0.2e1 * t120 * t33 - t62 * t66 + t96 * t35 - t63 * t65 - t148 * t34 + (t120 * t165 + t42 * t73 + t43 * t72) * t212 + t131 * t100 - 0.2e1 * pkin(4) * t97 - t111 * t172 + (qJD(5) * t112 + t99) * t135 + 0.2e1 * (-t148 * t42 - t43 * t96 + t62 * t72 - t63 * t73) * mrSges(7,3); t4 * mrSges(6,1) - t5 * mrSges(6,2) + (t130 * t2 + t134 * t3 + (-t12 * t130 + t13 * t134) * qJD(6)) * t205 + t164; t25 * mrSges(6,1) - t24 * mrSges(6,2) + (t130 * t7 + t134 * t8 + (-t130 * t31 + t134 * t32) * qJD(6)) * t205 + t163; -Ifges(6,5) * t157 + t47 * mrSges(6,1) - t46 * mrSges(6,2) - t140 * Ifges(6,6) + (t74 * t168 + t130 * t28 - t75 * t169 + t134 * t27 + m(7) * (t10 * t130 + t11 * t134 + t30 * t168 - t29 * t169)) * pkin(5) + t147 - t217; t123 + (t109 * pkin(10) - t192) * qJD(5) + (m(7) * (t130 * t42 + t134 * t43 + (-t130 * t72 + t134 * t73) * qJD(6)) + (-t130 * t63 + t134 * t62 + (t130 * t96 - t134 * t148) * qJD(6)) * mrSges(7,3)) * pkin(5) + t149; 0.2e1 * t92; t164; t163; t147; t149; t92; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
