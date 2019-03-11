% Calculate time derivative of joint inertia matrix for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:30
% EndTime: 2019-03-09 13:12:36
% DurationCPUTime: 3.31s
% Computational Cost: add. (13283->332), mult. (27998->502), div. (0->0), fcn. (30008->10), ass. (0->165)
t135 = sin(qJ(6));
t131 = t135 ^ 2;
t139 = cos(qJ(6));
t132 = t139 ^ 2;
t221 = t131 + t132;
t137 = sin(qJ(4));
t141 = cos(qJ(4));
t133 = sin(pkin(11));
t134 = cos(pkin(11));
t138 = sin(qJ(2));
t142 = cos(qJ(2));
t110 = t133 * t142 + t134 * t138;
t196 = -qJ(3) - pkin(7);
t118 = t196 * t138;
t119 = t196 * t142;
t87 = t134 * t118 + t119 * t133;
t82 = -pkin(8) * t110 + t87;
t109 = -t133 * t138 + t134 * t142;
t88 = t133 * t118 - t134 * t119;
t83 = pkin(8) * t109 + t88;
t55 = -t137 * t83 + t141 * t82;
t175 = qJD(6) * t139;
t136 = sin(qJ(5));
t140 = cos(qJ(5));
t103 = t110 * qJD(2);
t104 = t109 * qJD(2);
t85 = t109 * t141 - t110 * t137;
t64 = qJD(4) * t85 - t103 * t137 + t104 * t141;
t86 = t109 * t137 + t110 * t141;
t65 = -qJD(4) * t86 - t103 * t141 - t104 * t137;
t66 = t136 * t86 - t140 * t85;
t35 = -qJD(5) * t66 + t136 * t65 + t140 * t64;
t67 = t136 * t85 + t140 * t86;
t153 = t135 * t35 + t67 * t175;
t125 = pkin(2) * t134 + pkin(3);
t206 = pkin(2) * t133;
t99 = t141 * t125 - t137 * t206;
t220 = t221 * t140;
t154 = -pkin(9) * t86 + t55;
t56 = t137 * t82 + t141 * t83;
t45 = pkin(9) * t85 + t56;
t22 = t136 * t154 + t140 * t45;
t128 = -pkin(2) * t142 - pkin(1);
t92 = -t109 * pkin(3) + t128;
t71 = -t85 * pkin(4) + t92;
t40 = t66 * pkin(5) - t67 * pkin(10) + t71;
t16 = -t135 * t22 + t139 * t40;
t17 = t135 * t40 + t139 * t22;
t219 = -t135 * t16 + t139 * t17;
t163 = qJD(2) * t196;
t101 = qJD(3) * t142 + t138 * t163;
t102 = -t138 * qJD(3) + t142 * t163;
t80 = -t101 * t133 + t102 * t134;
t149 = -pkin(8) * t104 + t80;
t81 = t101 * t134 + t102 * t133;
t150 = -pkin(8) * t103 + t81;
t38 = t55 * qJD(4) + t137 * t149 + t141 * t150;
t39 = -qJD(4) * t56 - t137 * t150 + t141 * t149;
t218 = t39 * mrSges(5,1) - t38 * mrSges(5,2) + Ifges(5,5) * t64 + Ifges(5,6) * t65;
t194 = mrSges(7,1) * t139;
t117 = mrSges(7,2) * t135 - t194;
t190 = pkin(4) * qJD(5);
t107 = t136 * t117 * t190;
t172 = t140 * t190;
t161 = mrSges(7,3) * t172;
t123 = t131 * t161;
t124 = t132 * t161;
t159 = mrSges(7,1) * t135 + mrSges(7,2) * t139;
t114 = t159 * qJD(6);
t127 = -pkin(4) * t140 - pkin(5);
t96 = t127 * t114;
t217 = t107 + t123 + t124 + t96;
t216 = 2 * m(5);
t215 = 2 * m(6);
t214 = 2 * m(7);
t213 = -2 * mrSges(6,3);
t144 = -t64 * pkin(9) + t39;
t19 = pkin(9) * t65 + t38;
t11 = qJD(5) * t22 + t136 * t19 - t140 * t144;
t212 = 0.2e1 * t11;
t21 = t136 * t45 - t140 * t154;
t211 = 0.2e1 * t21;
t130 = qJD(2) * t138 * pkin(2);
t89 = pkin(3) * t103 + t130;
t47 = -pkin(4) * t65 + t89;
t210 = 0.2e1 * t47;
t209 = 0.2e1 * t128;
t208 = m(4) * pkin(2);
t205 = pkin(5) * t114;
t204 = t11 * t21;
t10 = -qJD(5) * t21 + t136 * t144 + t140 * t19;
t36 = qJD(5) * t67 + t136 * t64 - t140 * t65;
t13 = pkin(5) * t36 - pkin(10) * t35 + t47;
t2 = qJD(6) * t16 + t10 * t139 + t13 * t135;
t203 = t139 * t2;
t100 = t125 * t137 + t141 * t206;
t97 = pkin(4) + t99;
t79 = t140 * t100 + t136 * t97;
t94 = t99 * qJD(4);
t95 = t100 * qJD(4);
t52 = qJD(5) * t79 + t136 * t94 + t140 * t95;
t202 = t21 * t52;
t179 = qJD(6) * t17;
t3 = -t10 * t135 + t13 * t139 - t179;
t201 = t3 * t135;
t78 = -t100 * t136 + t140 * t97;
t51 = qJD(5) * t78 - t136 * t95 + t140 * t94;
t199 = t51 * mrSges(6,2);
t198 = t94 * mrSges(5,2);
t197 = t95 * mrSges(5,1);
t181 = t139 * t35;
t195 = Ifges(7,5) * t181 + Ifges(7,3) * t36;
t193 = Ifges(7,4) * t135;
t192 = Ifges(7,4) * t139;
t191 = Ifges(7,6) * t135;
t189 = t131 * t51;
t188 = t132 * t51;
t185 = t135 * t67;
t180 = t139 * t67;
t176 = qJD(6) * t135;
t174 = 0.2e1 * t142;
t171 = t67 * t176;
t152 = t171 - t181;
t12 = mrSges(7,1) * t153 - mrSges(7,2) * t152;
t169 = m(7) * t11 + t12;
t168 = -t65 * mrSges(5,1) + t64 * mrSges(5,2);
t167 = t36 * mrSges(6,1) + t35 * mrSges(6,2);
t166 = -t176 / 0.2e1;
t165 = t103 * mrSges(4,1) + t104 * mrSges(4,2);
t164 = -(2 * Ifges(6,4)) - t191;
t162 = t221 * t51;
t160 = -t136 * mrSges(6,1) - t140 * mrSges(6,2);
t158 = Ifges(7,1) * t139 - t193;
t157 = -Ifges(7,2) * t135 + t192;
t156 = Ifges(7,5) * t135 + Ifges(7,6) * t139;
t42 = -mrSges(7,2) * t66 - mrSges(7,3) * t185;
t43 = mrSges(7,1) * t66 - mrSges(7,3) * t180;
t155 = -t135 * t43 + t139 * t42;
t115 = t157 * qJD(6);
t116 = t158 * qJD(6);
t120 = Ifges(7,2) * t139 + t193;
t121 = Ifges(7,1) * t135 + t192;
t151 = t139 * t115 + t135 * t116 - t120 * t176 + t121 * t175;
t148 = -t201 + (-t135 * t17 - t139 * t16) * qJD(6);
t46 = t52 * t117;
t48 = mrSges(7,3) * t189;
t49 = mrSges(7,3) * t188;
t50 = t52 * mrSges(6,1);
t76 = -pkin(5) - t78;
t70 = t76 * t114;
t147 = t151 + t46 + t48 + t49 - t50 + t70 - t199;
t129 = Ifges(7,5) * t175;
t33 = t66 * Ifges(7,6) + t157 * t67;
t34 = t66 * Ifges(7,5) + t158 * t67;
t6 = -Ifges(7,4) * t152 - Ifges(7,2) * t153 + t36 * Ifges(7,6);
t7 = -Ifges(7,1) * t152 - Ifges(7,4) * t153 + t36 * Ifges(7,5);
t146 = -t10 * mrSges(6,2) + mrSges(7,3) * t203 + t21 * t114 + t33 * t166 + t34 * t175 / 0.2e1 + Ifges(6,5) * t35 + t135 * t7 / 0.2e1 + t139 * t6 / 0.2e1 - t115 * t185 / 0.2e1 + t116 * t180 / 0.2e1 + t66 * (-Ifges(7,6) * t176 + t129) / 0.2e1 + (t156 / 0.2e1 - Ifges(6,6)) * t36 - t153 * t120 / 0.2e1 + (t181 / 0.2e1 + t67 * t166) * t121 + (t117 - mrSges(6,1)) * t11;
t14 = mrSges(7,1) * t36 + mrSges(7,3) * t152;
t15 = -mrSges(7,2) * t36 - mrSges(7,3) * t153;
t145 = -t42 * t176 - t43 * t175 - t135 * t14 + t139 * t15 + m(7) * (-t16 * t175 - t17 * t176 - t201 + t203);
t143 = mrSges(7,3) * t148 + t146;
t126 = pkin(4) * t136 + pkin(10);
t77 = pkin(10) + t79;
t41 = t159 * t67;
t1 = [0.2e1 * (t64 * t85 + t65 * t86) * Ifges(5,4) + 0.2e1 * (t38 * t85 - t39 * t86 - t55 * t64 + t56 * t65) * mrSges(5,3) + 0.2e1 * (-t103 * t110 + t104 * t109) * Ifges(4,4) + 0.2e1 * (-t103 * t88 - t104 * t87 + t109 * t81 - t110 * t80) * mrSges(4,3) + t22 * t36 * t213 + t165 * t209 + t12 * t211 + t41 * t212 + (t16 * t3 + t17 * t2 + t204) * t214 + (t10 * t22 + t47 * t71 + t204) * t215 + (t38 * t56 + t39 * t55 + t89 * t92) * t216 + 0.2e1 * t85 * Ifges(5,2) * t65 - 0.2e1 * t109 * Ifges(4,2) * t103 + 0.2e1 * t64 * t86 * Ifges(5,1) + 0.2e1 * t110 * t104 * Ifges(4,1) + 0.2e1 * t89 * (-mrSges(5,1) * t85 + mrSges(5,2) * t86) + 0.2e1 * t71 * t167 + 0.2e1 * t92 * t168 + 0.2e1 * m(4) * (t80 * t87 + t81 * t88) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t142) * t174 + (-0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (-mrSges(4,1) * t109 + mrSges(4,2) * t110) + t208 * t209 - 0.2e1 * Ifges(3,4) * t138 + (Ifges(3,1) - Ifges(3,2)) * t174) * t138) * qJD(2) + (mrSges(6,3) * t211 - t135 * t33 + t139 * t34) * t35 + (mrSges(6,2) * t210 + mrSges(6,3) * t212 + 0.2e1 * Ifges(6,1) * t35 - t135 * t6 + t139 * t7 + (Ifges(7,5) * t139 + t164) * t36 + (-t135 * t34 - t139 * t33 - t156 * t66) * qJD(6)) * t67 + (mrSges(6,1) * t210 + t10 * t213 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t36 + t164 * t35 + t195) * t66 + 0.2e1 * t16 * t14 + 0.2e1 * t17 * t15 + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43; (Ifges(3,5) * t142 - Ifges(3,6) * t138 + (-mrSges(3,1) * t142 + mrSges(3,2) * t138) * pkin(7)) * qJD(2) + (t77 * t15 + t51 * t42 + (-mrSges(7,3) * t16 - t43 * t77) * qJD(6)) * t139 + m(5) * (t100 * t38 + t39 * t99 - t55 * t95 + t56 * t94) + (t100 * t65 - t64 * t99 + t85 * t94 + t86 * t95) * mrSges(5,3) + t146 + t218 + (t133 * t81 + t134 * t80) * t208 + (-t35 * t78 - t36 * t79 - t51 * t66 + t52 * t67) * mrSges(6,3) + (t11 * t76 + t202 + (t148 + t203) * t77 + t219 * t51) * m(7) + (-t103 * t133 - t104 * t134) * pkin(2) * mrSges(4,3) + Ifges(4,5) * t104 - Ifges(4,6) * t103 + t76 * t12 + t80 * mrSges(4,1) - t81 * mrSges(4,2) + (-t51 * t43 + (-qJD(6) * t42 - t14) * t77 + (-t3 - t179) * mrSges(7,3)) * t135 + m(6) * (t10 * t79 - t11 * t78 + t22 * t51 + t202) + t52 * t41; -0.2e1 * t197 - 0.2e1 * t198 - 0.2e1 * t199 + 0.2e1 * t46 + 0.2e1 * t48 + 0.2e1 * t49 - 0.2e1 * t50 + 0.2e1 * t70 + (t162 * t77 + t52 * t76) * t214 + (t51 * t79 - t52 * t78) * t215 + (t100 * t94 - t95 * t99) * t216 + t151; m(4) * t130 + t135 * t15 + t139 * t14 + t155 * qJD(6) + m(7) * (t219 * qJD(6) + t135 * t2 + t139 * t3) + m(6) * t47 + m(5) * t89 + t165 + t167 + t168; 0; 0; t143 + (m(6) * (t10 * t136 - t11 * t140) + (-t136 * t36 - t140 * t35) * mrSges(6,3) + ((m(6) * t22 + m(7) * t219 - t66 * mrSges(6,3) + t155) * t140 + (mrSges(6,3) * t67 + t41 + (m(6) + m(7)) * t21) * t136) * qJD(5)) * pkin(4) + t169 * t127 + t145 * t126 + t218; m(7) * (t127 * t52 + (t188 + t189) * t126) + (m(6) * (t136 * t51 - t140 * t52) + (m(7) * (t136 * t76 + t220 * t77) + m(6) * (-t136 * t78 + t140 * t79) + t160) * qJD(5)) * pkin(4) + t147 - t198 - t197 + t217; 0; 0.2e1 * t107 + 0.2e1 * t123 + 0.2e1 * t124 + 0.2e1 * t96 + 0.2e1 * (m(7) * (t220 * t126 + t127 * t136) + t160) * t190 + t151; -pkin(5) * t169 + pkin(10) * t145 + t143; m(7) * (-pkin(5) * t52 + pkin(10) * t162) - t205 + t147; 0; -t205 + (m(7) * (-pkin(5) * t136 + t220 * pkin(10)) + t160) * t190 + t151 + t217; t151 - 0.2e1 * t205; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t171 - Ifges(7,6) * t153 + t195; t129 - t159 * t51 + (-t77 * t194 + (mrSges(7,2) * t77 - Ifges(7,6)) * t135) * qJD(6); -t114; t129 - t159 * t172 + (t117 * t126 - t191) * qJD(6); t129 + (pkin(10) * t117 - t191) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
