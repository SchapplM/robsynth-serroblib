% Calculate time derivative of joint inertia matrix for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:40
% EndTime: 2019-03-09 05:34:47
% DurationCPUTime: 3.77s
% Computational Cost: add. (2787->459), mult. (6114->653), div. (0->0), fcn. (4656->6), ass. (0->178)
t211 = Ifges(6,2) + Ifges(5,3);
t128 = sin(qJ(4));
t131 = cos(qJ(4));
t207 = Ifges(6,4) + Ifges(5,5);
t210 = -Ifges(6,6) * t128 - t131 * t207;
t199 = 2 * qJD(2);
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t134 = -pkin(1) - pkin(7);
t174 = t129 * t134;
t209 = -qJD(4) * t174 + qJD(2) + (pkin(3) * t132 + pkin(8) * t129) * qJD(3);
t127 = sin(qJ(6));
t130 = cos(qJ(6));
t145 = t127 * t131 - t128 * t130;
t204 = qJD(4) - qJD(6);
t208 = t204 * t145;
t171 = t128 ^ 2 + t131 ^ 2;
t176 = qJ(5) * t131;
t198 = pkin(4) + pkin(5);
t140 = -t128 * t198 + t176;
t206 = t134 + t140;
t205 = -t127 * mrSges(7,1) - t130 * mrSges(7,2);
t177 = qJ(5) * t128;
t203 = -t131 * t198 - t177;
t147 = pkin(4) * t131 + t177;
t164 = qJD(5) * t131;
t202 = qJD(4) * t147 - t164;
t170 = qJD(3) * t129;
t157 = t128 * t170;
t166 = qJD(4) * t132;
t158 = t131 * t166;
t169 = qJD(3) * t132;
t201 = Ifges(5,6) * t157 + Ifges(6,6) * t158 + t211 * t169;
t200 = 2 * m(7);
t197 = pkin(8) - pkin(9);
t96 = -qJ(5) * t127 - t130 * t198;
t67 = qJD(5) * t130 + qJD(6) * t96;
t196 = t67 * mrSges(7,2);
t97 = qJ(5) * t130 - t127 * t198;
t68 = -qJD(5) * t127 - qJD(6) * t97;
t195 = t68 * mrSges(7,1);
t144 = t127 * t128 + t130 * t131;
t42 = t204 * t144;
t21 = t132 * t42 + t145 * t170;
t22 = t132 * t208 - t144 * t170;
t193 = -Ifges(7,5) * t22 - Ifges(7,6) * t21;
t160 = t128 * t166;
t138 = t131 * t170 + t160;
t55 = mrSges(5,1) * t169 + mrSges(5,3) * t138;
t56 = -mrSges(6,1) * t169 - mrSges(6,2) * t138;
t192 = -t55 + t56;
t137 = t157 - t158;
t57 = -mrSges(5,2) * t169 + mrSges(5,3) * t137;
t58 = mrSges(6,2) * t137 + mrSges(6,3) * t169;
t191 = t57 + t58;
t184 = Ifges(6,5) * t131;
t148 = Ifges(6,3) * t128 + t184;
t182 = Ifges(6,6) * t129;
t63 = t132 * t148 + t182;
t186 = Ifges(5,4) * t131;
t149 = -Ifges(5,2) * t128 + t186;
t180 = t129 * Ifges(5,6);
t64 = t132 * t149 + t180;
t190 = t63 - t64;
t175 = t128 * t132;
t90 = -mrSges(5,2) * t129 - mrSges(5,3) * t175;
t93 = -mrSges(6,2) * t175 + mrSges(6,3) * t129;
t189 = t90 + t93;
t173 = t131 * t132;
t91 = mrSges(5,1) * t129 - mrSges(5,3) * t173;
t92 = -mrSges(6,1) * t129 + mrSges(6,2) * t173;
t188 = -t91 + t92;
t187 = Ifges(5,4) * t128;
t185 = Ifges(6,5) * t128;
t104 = -t131 * mrSges(5,1) + t128 * mrSges(5,2);
t178 = t104 - mrSges(4,1);
t98 = pkin(3) * t129 - pkin(8) * t132 + qJ(2);
t60 = t128 * t98 + t131 * t174;
t172 = t171 * pkin(8) * t169;
t168 = qJD(4) * t128;
t167 = qJD(4) * t131;
t165 = qJD(5) * t128;
t49 = t129 * qJ(5) + t60;
t110 = t197 * t131;
t162 = t129 * t169;
t161 = t134 * t169;
t105 = -Ifges(6,3) * t131 + t185;
t106 = Ifges(5,2) * t131 + t187;
t156 = t105 / 0.2e1 - t106 / 0.2e1;
t107 = Ifges(6,1) * t128 - t184;
t108 = Ifges(5,1) * t128 + t186;
t155 = -t107 / 0.2e1 - t108 / 0.2e1;
t113 = t128 * t174;
t59 = t131 * t98 - t113;
t72 = t144 * t132;
t19 = qJD(3) * t72 + t129 * t208;
t71 = t145 * t132;
t20 = -qJD(3) * t71 + t129 * t42;
t154 = t20 * mrSges(7,1) - t19 * mrSges(7,2);
t153 = t128 * mrSges(5,1) + t131 * mrSges(5,2);
t103 = -t131 * mrSges(6,1) - t128 * mrSges(6,3);
t152 = t128 * mrSges(6,1) - t131 * mrSges(6,3);
t151 = Ifges(5,1) * t131 - t187;
t150 = Ifges(6,1) * t131 + t185;
t146 = pkin(4) * t128 - t176;
t36 = t113 + (-pkin(9) * t132 - t98) * t131 - t198 * t129;
t39 = pkin(9) * t175 + t49;
t7 = -t127 * t39 + t130 * t36;
t8 = t127 * t36 + t130 * t39;
t109 = t197 * t128;
t51 = t109 * t130 - t110 * t127;
t52 = t109 * t127 + t110 * t130;
t27 = -t128 * t161 + t209 * t131 - t98 * t168;
t65 = Ifges(6,4) * t129 + t132 * t150;
t66 = Ifges(5,5) * t129 + t132 * t151;
t143 = -t129 * t207 - t65 - t66;
t142 = -t134 + t146;
t6 = pkin(9) * t160 + (pkin(9) * t129 * t131 - t132 * t198) * qJD(3) - t27;
t26 = t209 * t128 + t131 * t161 + t98 * t167;
t16 = qJ(5) * t169 + t129 * qJD(5) + t26;
t9 = -pkin(9) * t137 + t16;
t1 = qJD(6) * t7 + t127 * t6 + t130 * t9;
t2 = -qJD(6) * t8 - t127 * t9 + t130 * t6;
t141 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) + t193;
t94 = t197 * t168;
t95 = qJD(4) * t110;
t23 = qJD(6) * t51 + t127 * t95 - t130 * t94;
t24 = -qJD(6) * t52 + t127 * t94 + t130 * t95;
t40 = Ifges(7,6) * t208;
t41 = Ifges(7,5) * t42;
t139 = t24 * mrSges(7,1) - t23 * mrSges(7,2) - t40 + t41;
t135 = m(6) * t164 + (-m(6) * t147 + t103 + t104) * qJD(4);
t122 = Ifges(6,4) * t167;
t121 = Ifges(5,5) * t167;
t119 = Ifges(6,6) * t168;
t99 = -pkin(3) - t147;
t89 = t151 * qJD(4);
t88 = t150 * qJD(4);
t87 = t149 * qJD(4);
t86 = t148 * qJD(4);
t85 = t153 * qJD(4);
t84 = t152 * qJD(4);
t79 = pkin(3) - t203;
t78 = t153 * t132;
t77 = t152 * t132;
t73 = qJD(4) * t146 - t165;
t70 = t144 * t129;
t69 = t145 * t129;
t62 = t142 * t132;
t61 = qJD(4) * t140 + t165;
t54 = -mrSges(7,1) * t129 - mrSges(7,3) * t72;
t53 = mrSges(7,2) * t129 - mrSges(7,3) * t71;
t50 = -pkin(4) * t129 - t59;
t48 = t206 * t132;
t46 = -Ifges(7,1) * t145 - Ifges(7,4) * t144;
t45 = -Ifges(7,4) * t145 - Ifges(7,2) * t144;
t44 = mrSges(7,1) * t144 - mrSges(7,2) * t145;
t38 = -mrSges(5,1) * t137 - mrSges(5,2) * t138;
t37 = -mrSges(6,1) * t137 + mrSges(6,3) * t138;
t35 = mrSges(7,1) * t71 + mrSges(7,2) * t72;
t34 = -t108 * t166 + (Ifges(5,5) * t132 - t129 * t151) * qJD(3);
t33 = -t107 * t166 + (Ifges(6,4) * t132 - t129 * t150) * qJD(3);
t32 = -t106 * t166 + (Ifges(5,6) * t132 - t129 * t149) * qJD(3);
t31 = -t105 * t166 + (Ifges(6,6) * t132 - t129 * t148) * qJD(3);
t30 = Ifges(7,1) * t72 - Ifges(7,4) * t71 - Ifges(7,5) * t129;
t29 = Ifges(7,4) * t72 - Ifges(7,2) * t71 - Ifges(7,6) * t129;
t28 = t132 * t202 - t142 * t170;
t25 = -pkin(4) * t169 - t27;
t15 = Ifges(7,1) * t42 - Ifges(7,4) * t208;
t14 = Ifges(7,4) * t42 - Ifges(7,2) * t208;
t13 = mrSges(7,1) * t208 + mrSges(7,2) * t42;
t12 = (qJD(4) * t203 + t164) * t132 - t206 * t170;
t11 = -mrSges(7,1) * t169 - mrSges(7,3) * t22;
t10 = mrSges(7,2) * t169 + mrSges(7,3) * t21;
t5 = -mrSges(7,1) * t21 + mrSges(7,2) * t22;
t4 = Ifges(7,1) * t22 + Ifges(7,4) * t21 - Ifges(7,5) * t169;
t3 = Ifges(7,4) * t22 + Ifges(7,2) * t21 - Ifges(7,6) * t169;
t17 = [(mrSges(4,1) * t199 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * Ifges(4,4) * t129 + 0.2e1 * t134 * t78 + (-t182 - t190) * t128 + t143 * t131) * qJD(3) + t193 + t201) * t129 + 0.2e1 * t8 * t10 + 0.2e1 * t7 * t11 + (t1 * t8 + t12 * t48 + t2 * t7) * t200 + 0.2e1 * m(5) * (t60 * t26 + t59 * t27) + (mrSges(4,2) * t199 - 0.2e1 * t134 * t38 + (t33 + t34) * t131 + (t31 - t32) * t128 + ((-t180 + t190) * t131 + t143 * t128) * qJD(4) + (-Ifges(7,5) * t72 + Ifges(7,6) * t71 + 0.2e1 * qJ(2) * mrSges(4,1) + (-Ifges(5,6) * t128 - 0.2e1 * Ifges(4,4) - t210) * t132 + (-0.2e1 * m(5) * t134 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + (2 * Ifges(7,3)) + t211) * t129) * qJD(3)) * t132 + 0.2e1 * m(6) * (t16 * t49 + t25 * t50 + t28 * t62) + t21 * t29 + t22 * t30 + 0.2e1 * t12 * t35 + 0.2e1 * t48 * t5 + 0.2e1 * t1 * t53 + 0.2e1 * t2 * t54 + 0.2e1 * t50 * t56 + 0.2e1 * t49 * t58 + 0.2e1 * t59 * t55 + 0.2e1 * t60 * t57 + 0.2e1 * t62 * t37 - t71 * t3 + t72 * t4 + 0.2e1 * t28 * t77 + 0.2e1 * t26 * t90 + 0.2e1 * t27 * t91 + 0.2e1 * t25 * t92 + 0.2e1 * t16 * t93 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t199; m(7) * (t1 * t70 + t19 * t8 - t2 * t69 + t20 * t7) + t19 * t53 + t20 * t54 - t69 * t11 + t70 * t10 + (t5 - t37 - t38 - m(6) * t28 + m(7) * t12 + (t189 * t131 + t188 * t128 + m(6) * (t128 * t50 + t131 * t49) + m(5) * (-t128 * t59 + t131 * t60)) * qJD(3)) * t132 + (t191 * t131 + t192 * t128 + (-t128 * t189 + t131 * t188) * qJD(4) + m(6) * (t128 * t25 + t131 * t16 + t167 * t50 - t168 * t49) + m(5) * (-t27 * t128 + t131 * t26 - t167 * t59 - t168 * t60 - 0.2e1 * t161) + (m(6) * t62 - m(7) * t48 - t35 + t77 + t78) * qJD(3)) * t129; (t19 * t70 - t20 * t69 - t162) * t200 + 0.4e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t171) * t162; -t144 * t3 / 0.2e1 + m(7) * (t1 * t52 + t12 * t79 + t2 * t51 + t23 * t8 + t24 * t7 + t48 * t61) + (-t134 * t85 + (t88 / 0.2e1 + t89 / 0.2e1) * t131 + (t86 / 0.2e1 - t87 / 0.2e1) * t128 + (Ifges(7,5) * t145 / 0.2e1 + Ifges(7,6) * t144 / 0.2e1 - Ifges(4,6) - t134 * mrSges(4,2) + (Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t131 + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t128) * qJD(3) + (t128 * t155 + t131 * t156) * qJD(4)) * t132 - t145 * t4 / 0.2e1 + (t26 * mrSges(5,3) + t16 * mrSges(6,2) - t31 / 0.2e1 + t32 / 0.2e1 + t155 * t170 + (t50 * mrSges(6,2) - t59 * mrSges(5,3) + t65 / 0.2e1 + t66 / 0.2e1) * qJD(4) + (t188 * qJD(4) + m(6) * (t50 * qJD(4) + t16) + m(5) * (-t59 * qJD(4) + t26) + t191) * pkin(8)) * t131 + (t25 * mrSges(6,2) - t27 * mrSges(5,3) + t33 / 0.2e1 + t34 / 0.2e1 - t156 * t170 + (-t49 * mrSges(6,2) - t60 * mrSges(5,3) - t64 / 0.2e1 + t63 / 0.2e1 - t180 / 0.2e1) * qJD(4) + (-t189 * qJD(4) + m(6) * (-t49 * qJD(4) + t25) + m(5) * (-t60 * qJD(4) - t27) + t192) * pkin(8)) * t128 + (-t1 * t144 + t145 * t2 - t208 * t8 - t7 * t42) * mrSges(7,3) - t208 * t29 / 0.2e1 + (t121 / 0.2e1 - t41 / 0.2e1 + t40 / 0.2e1 + t122 / 0.2e1 + t119 / 0.2e1 + (-Ifges(4,5) + (-m(5) * pkin(3) + t178) * t134) * qJD(3)) * t129 + m(6) * (t28 * t99 + t62 * t73) - pkin(3) * t38 + t42 * t30 / 0.2e1 + t12 * t44 + t21 * t45 / 0.2e1 + t22 * t46 / 0.2e1 + t48 * t13 + t51 * t11 + t52 * t10 + t23 * t53 + t24 * t54 + t61 * t35 - t71 * t14 / 0.2e1 + t72 * t15 / 0.2e1 + t73 * t77 + t79 * t5 + t62 * t84 + t99 * t37 + t28 * t103; (t13 - t84 - t85) * t132 + ((t103 - t44 + t178) * t129 + (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t171) * t132) * qJD(3) + m(6) * (-t132 * t73 + t170 * t99 + t172) + m(7) * (t132 * t61 - t170 * t79 + t19 * t52 + t20 * t51 + t23 * t70 - t24 * t69) + m(5) * (-pkin(3) * t170 + t172) + (-t144 * t19 + t145 * t20 - t208 * t70 + t42 * t69) * mrSges(7,3); (t23 * t52 + t24 * t51 + t61 * t79) * t200 + 0.2e1 * t61 * t44 + 0.2e1 * t79 * t13 - t208 * t45 - t144 * t14 + t42 * t46 - t145 * t15 + 0.2e1 * t99 * t84 - 0.2e1 * pkin(3) * t85 + (t87 - t86) * t131 + (t88 + t89) * t128 + ((t107 + t108) * t131 + (t105 - t106) * t128) * qJD(4) + 0.2e1 * (m(6) * t99 + t103) * t73 + 0.2e1 * (-t144 * t23 + t145 * t24 - t208 * t52 - t42 * t51) * mrSges(7,3); t16 * mrSges(6,3) + t141 + t201 + m(7) * (t1 * t97 + t2 * t96 + t67 * t8 + t68 * t7) + m(6) * (-pkin(4) * t25 + qJ(5) * t16 + qJD(5) * t49) + (Ifges(7,3) * t132 + t210 * t129) * qJD(3) + (-Ifges(5,6) * t131 - t128 * t207) * t166 - t25 * mrSges(6,1) - t26 * mrSges(5,2) + t27 * mrSges(5,1) - pkin(4) * t56 + qJ(5) * t58 + t67 * t53 + t68 * t54 + qJD(5) * t93 + t96 * t11 + t97 * t10; m(7) * (t19 * t97 + t20 * t96 + t67 * t70 - t68 * t69) + (-m(6) * t146 - t152 - t153) * t169 + t135 * t129 - t154; m(7) * (t23 * t97 + t24 * t96 + t51 * t68 + t52 * t67) - Ifges(5,6) * t168 + t122 + t119 + t121 - t202 * mrSges(6,2) + (-t144 * t67 + t145 * t68 - t208 * t97 - t42 * t96) * mrSges(7,3) + t135 * pkin(8) - t139; 0.2e1 * t196 - 0.2e1 * t195 + (t67 * t97 + t68 * t96) * t200 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t127 * t10 + t130 * t11 + (-t127 * t54 + t130 * t53) * qJD(6) + m(7) * (t1 * t127 + t130 * t2 + (-t127 * t7 + t130 * t8) * qJD(6)) + m(6) * t25 + t56; m(7) * (t127 * t19 + t130 * t20 + (t127 * t69 + t130 * t70) * qJD(6)) + (t128 * t169 + t129 * t167) * m(6); m(7) * (t127 * t23 + t130 * t24 + (-t127 * t51 + t130 * t52) * qJD(6)) + (m(6) * pkin(8) + mrSges(6,2)) * t167 + (-t127 * t208 - t130 * t42 + (-t127 * t145 - t130 * t144) * qJD(6)) * mrSges(7,3); m(7) * (t127 * t67 + t130 * t68) + (m(7) * (-t127 * t96 + t130 * t97) - t205) * qJD(6); 0; -Ifges(7,3) * t169 - t141; t154; t139; t195 - t196; t205 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
