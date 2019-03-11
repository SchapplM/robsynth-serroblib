% Calculate time derivative of joint inertia matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:28
% EndTime: 2019-03-09 12:45:35
% DurationCPUTime: 3.84s
% Computational Cost: add. (4139->422), mult. (8955->593), div. (0->0), fcn. (7299->6), ass. (0->187)
t248 = Ifges(6,5) + Ifges(7,5);
t247 = Ifges(6,6) + Ifges(7,6);
t162 = cos(qJ(4));
t159 = sin(qJ(4));
t163 = cos(qJ(2));
t200 = qJD(4) * t163;
t196 = t159 * t200;
t160 = sin(qJ(2));
t204 = qJD(2) * t160;
t173 = t162 * t204 + t196;
t249 = mrSges(6,1) + mrSges(7,1);
t246 = Ifges(6,3) + Ifges(7,3);
t158 = sin(qJ(5));
t199 = qJD(5) * t158;
t202 = qJD(4) * t159;
t161 = cos(qJ(5));
t207 = t161 * t162;
t238 = qJD(4) + qJD(5);
t78 = -t158 * t202 - t159 * t199 + t207 * t238;
t178 = t158 * t162 + t161 * t159;
t79 = t238 * t178;
t245 = -t247 * t78 - t248 * t79;
t244 = 2 * qJ(3);
t243 = m(4) * pkin(7);
t242 = pkin(4) * t161;
t240 = pkin(4) * qJD(5);
t164 = -pkin(2) - pkin(8);
t209 = qJ(3) * t160;
t113 = t164 * t163 - pkin(1) - t209;
t229 = pkin(3) + pkin(7);
t137 = t229 * t160;
t122 = t159 * t137;
t77 = t162 * t113 + t122;
t225 = pkin(9) - t164;
t131 = t225 * t159;
t132 = t225 * t162;
t87 = -t161 * t131 - t158 * t132;
t203 = qJD(2) * t163;
t130 = t229 * t203;
t201 = qJD(4) * t162;
t188 = pkin(2) * t204 - qJD(3) * t160;
t96 = (pkin(8) * t160 - qJ(3) * t163) * qJD(2) + t188;
t27 = -t113 * t202 + t159 * t130 + t137 * t201 + t162 * t96;
t189 = t162 * t130 - t159 * t96;
t28 = -qJD(4) * t77 + t189;
t239 = t159 * t27 + t162 * t28;
t237 = t163 * t238;
t193 = t159 * t204;
t236 = Ifges(5,5) * t193 + t173 * Ifges(5,6) + Ifges(5,3) * t203;
t235 = 2 * m(6);
t234 = 2 * m(7);
t233 = -0.2e1 * pkin(1);
t104 = -qJ(3) * t203 + t188;
t232 = 0.2e1 * t104;
t182 = -pkin(2) * t163 - t209;
t133 = -pkin(1) + t182;
t231 = -0.2e1 * t133;
t230 = m(7) * pkin(5);
t227 = -t163 / 0.2e1;
t156 = t163 * pkin(7);
t226 = -mrSges(6,2) - mrSges(7,2);
t177 = t158 * t159 - t207;
t48 = -t177 * t204 + t178 * t237;
t31 = -mrSges(7,2) * t203 + mrSges(7,3) * t48;
t32 = -mrSges(6,2) * t203 + mrSges(6,3) * t48;
t224 = t31 + t32;
t123 = t162 * t137;
t190 = pkin(9) * t163 - t113;
t58 = pkin(4) * t160 + t190 * t159 + t123;
t205 = t162 * t163;
t65 = -pkin(9) * t205 + t77;
t26 = t158 * t58 + t161 * t65;
t101 = t177 * t163;
t90 = -mrSges(7,2) * t160 + mrSges(7,3) * t101;
t91 = -mrSges(6,2) * t160 + mrSges(6,3) * t101;
t221 = t90 + t91;
t102 = t178 * t163;
t92 = mrSges(7,1) * t160 + mrSges(7,3) * t102;
t93 = mrSges(6,1) * t160 + mrSges(6,3) * t102;
t220 = t92 + t93;
t219 = Ifges(5,4) * t159;
t218 = Ifges(5,4) * t162;
t217 = t158 * t78;
t185 = Ifges(5,1) * t159 + t218;
t214 = t160 * Ifges(5,5);
t99 = -t185 * t163 + t214;
t215 = t159 * t99;
t213 = t161 * t79;
t184 = Ifges(5,2) * t162 + t219;
t98 = t160 * Ifges(5,6) - t184 * t163;
t211 = t162 * t98;
t198 = qJD(5) * t161;
t210 = (t178 * t198 + t217) * pkin(4);
t136 = Ifges(5,1) * t162 - t219;
t208 = t159 * t136;
t135 = -Ifges(5,2) * t159 + t218;
t206 = t162 * t135;
t146 = t159 * pkin(4) + qJ(3);
t138 = t163 * pkin(3) + t156;
t142 = pkin(4) * t201 + qJD(3);
t47 = t177 * t237 + t178 * t204;
t20 = (-pkin(9) * t159 * t160 + pkin(4) * t163) * qJD(2) + (t190 * t162 - t122) * qJD(4) + t189;
t22 = t173 * pkin(9) + t27;
t6 = -t26 * qJD(5) - t158 * t22 + t161 * t20;
t2 = pkin(5) * t203 - qJ(6) * t47 + qJD(6) * t102 + t6;
t29 = mrSges(7,1) * t203 - mrSges(7,3) * t47;
t197 = m(7) * t2 + t29;
t108 = pkin(4) * t205 + t138;
t195 = t162 * t200;
t194 = t177 * t199;
t11 = -t48 * mrSges(7,1) + t47 * mrSges(7,2);
t36 = t78 * mrSges(7,1) - t79 * mrSges(7,2);
t191 = t226 * t161;
t25 = -t158 * t65 + t161 * t58;
t86 = t131 * t158 - t161 * t132;
t118 = t225 * t202;
t119 = qJD(4) * t132;
t35 = -qJD(5) * t87 + t161 * t118 + t119 * t158;
t15 = qJ(6) * t79 + qJD(6) * t177 + t35;
t187 = m(7) * t15 + mrSges(7,3) * t79;
t186 = mrSges(5,1) * t162 - mrSges(5,2) * t159;
t134 = mrSges(5,1) * t159 + mrSges(5,2) * t162;
t183 = -Ifges(5,5) * t159 - Ifges(5,6) * t162;
t181 = -t177 * t79 - t178 * t78;
t76 = -t113 * t159 + t123;
t180 = t76 * t159 - t77 * t162;
t127 = mrSges(5,3) * t159 * t163 + mrSges(5,1) * t160;
t128 = -mrSges(5,2) * t160 - mrSges(5,3) * t205;
t179 = -t159 * t127 + t162 * t128;
t176 = t226 * t78 - t249 * t79;
t175 = t203 * t246 + t247 * t48 + t248 * t47;
t5 = t158 * t20 + t161 * t22 + t58 * t198 - t65 * t199;
t34 = t158 * t118 - t161 * t119 + t131 * t199 - t132 * t198;
t172 = t193 - t195;
t16 = pkin(5) * t160 + qJ(6) * t102 + t25;
t17 = qJ(6) * t101 + t26;
t3 = qJ(6) * t48 + qJD(6) * t101 + t5;
t171 = t16 * t79 - t17 * t78 + t177 * t2 - t178 * t3;
t170 = t177 * t6 - t178 * t5 + t25 * t79 - t26 * t78;
t14 = -qJ(6) * t78 - qJD(6) * t178 + t34;
t56 = qJ(6) * t177 + t86;
t57 = -qJ(6) * t178 + t87;
t169 = -t14 * t178 + t15 * t177 + t56 * t79 - t57 * t78;
t168 = t177 * t35 - t178 * t34 - t78 * t87 + t79 * t86;
t167 = t35 * mrSges(6,1) + t15 * mrSges(7,1) - t34 * mrSges(6,2) - t14 * mrSges(7,2) + t245;
t166 = -t217 + (-t158 * t177 - t161 * t178) * qJD(5);
t165 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t175;
t88 = -pkin(4) * t196 + (-pkin(4) * t162 - t229) * t204;
t149 = pkin(5) + t242;
t129 = t229 * t204;
t126 = t185 * qJD(4);
t125 = t184 * qJD(4);
t124 = t186 * qJD(4);
t112 = t186 * t163;
t97 = pkin(5) * t178 + t146;
t95 = mrSges(5,1) * t203 - t172 * mrSges(5,3);
t94 = -mrSges(5,2) * t203 + t173 * mrSges(5,3);
t85 = -Ifges(6,1) * t177 - Ifges(6,4) * t178;
t84 = -Ifges(7,1) * t177 - Ifges(7,4) * t178;
t83 = -Ifges(6,4) * t177 - Ifges(6,2) * t178;
t82 = -Ifges(7,4) * t177 - Ifges(7,2) * t178;
t81 = mrSges(6,1) * t178 - mrSges(6,2) * t177;
t80 = mrSges(7,1) * t178 - mrSges(7,2) * t177;
t67 = -t173 * mrSges(5,1) + t172 * mrSges(5,2);
t66 = -pkin(5) * t101 + t108;
t64 = pkin(5) * t78 + t142;
t63 = -mrSges(6,1) * t101 - mrSges(6,2) * t102;
t62 = -mrSges(7,1) * t101 - mrSges(7,2) * t102;
t61 = -t136 * t200 + (t163 * Ifges(5,5) + t185 * t160) * qJD(2);
t60 = -t135 * t200 + (t163 * Ifges(5,6) + t184 * t160) * qJD(2);
t55 = -Ifges(6,1) * t102 + Ifges(6,4) * t101 + Ifges(6,5) * t160;
t54 = -Ifges(7,1) * t102 + Ifges(7,4) * t101 + Ifges(7,5) * t160;
t53 = -Ifges(6,4) * t102 + Ifges(6,2) * t101 + Ifges(6,6) * t160;
t52 = -Ifges(7,4) * t102 + Ifges(7,2) * t101 + Ifges(7,6) * t160;
t41 = -Ifges(6,1) * t79 - Ifges(6,4) * t78;
t40 = -Ifges(7,1) * t79 - Ifges(7,4) * t78;
t39 = -Ifges(6,4) * t79 - Ifges(6,2) * t78;
t38 = -Ifges(7,4) * t79 - Ifges(7,2) * t78;
t37 = mrSges(6,1) * t78 - mrSges(6,2) * t79;
t30 = mrSges(6,1) * t203 - mrSges(6,3) * t47;
t23 = -pkin(5) * t48 + t88;
t12 = -mrSges(6,1) * t48 + mrSges(6,2) * t47;
t10 = Ifges(6,1) * t47 + Ifges(6,4) * t48 + Ifges(6,5) * t203;
t9 = Ifges(7,1) * t47 + Ifges(7,4) * t48 + Ifges(7,5) * t203;
t8 = Ifges(6,4) * t47 + Ifges(6,2) * t48 + Ifges(6,6) * t203;
t7 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t203;
t1 = [0.2e1 * m(5) * (-t129 * t138 + t27 * t77 + t28 * t76) + ((mrSges(3,1) * t233 + mrSges(4,2) * t231 + t215 + t211 + 0.2e1 * (-Ifges(4,6) - Ifges(3,4)) * t160) * t160 + (mrSges(3,2) * t233 + mrSges(4,3) * t231 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t183) * t163 - t248 * t102 + t247 * t101 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3) + t246) * t160) * t163) * qJD(2) + m(4) * t133 * t232 + (t16 * t2 + t17 * t3 + t23 * t66) * t234 + (t108 * t88 + t25 * t6 + t26 * t5) * t235 - (t9 + t10) * t102 + (t8 + t7) * t101 + (-0.2e1 * t104 * mrSges(4,3) + t175 + t236) * t160 + (t52 + t53) * t48 + (t54 + t55) * t47 + 0.2e1 * t138 * t67 + 0.2e1 * t28 * t127 + 0.2e1 * t27 * t128 - 0.2e1 * t129 * t112 + 0.2e1 * t108 * t12 + 0.2e1 * t3 * t90 + 0.2e1 * t5 * t91 + 0.2e1 * t2 * t92 + 0.2e1 * t6 * t93 + 0.2e1 * t77 * t94 + 0.2e1 * t76 * t95 + 0.2e1 * t88 * t63 + 0.2e1 * t23 * t62 + 0.2e1 * t66 * t11 + 0.2e1 * t16 * t29 + 0.2e1 * t25 * t30 + 0.2e1 * t17 * t31 + 0.2e1 * t26 * t32 + (mrSges(4,2) * t232 - t159 * t61 - t162 * t60 + (t159 * t98 + (-t99 - t214) * t162) * qJD(4)) * t163; (t61 / 0.2e1 - t28 * mrSges(5,3) - t125 * t227 + t164 * t95) * t162 + (-t60 / 0.2e1 - t27 * mrSges(5,3) - t126 * t227 + t164 * t94) * t159 + (qJD(4) * t183 + t245) * t160 / 0.2e1 + m(5) * (-qJ(3) * t129 + t164 * t239) + t170 * mrSges(6,3) + t171 * mrSges(7,3) - (t40 / 0.2e1 + t41 / 0.2e1) * t102 + (t82 / 0.2e1 + t83 / 0.2e1) * t48 + (t38 / 0.2e1 + t39 / 0.2e1) * t101 + (t84 / 0.2e1 + t85 / 0.2e1) * t47 + (-t215 / 0.2e1 - t211 / 0.2e1 + (t159 * t135 / 0.2e1 - t162 * t136 / 0.2e1) * t163 + t180 * mrSges(5,3) + (-m(5) * t180 + t179) * t164) * qJD(4) + (m(4) * t156 + m(5) * t138 + t163 * mrSges(4,1) + t112) * qJD(3) - (t8 / 0.2e1 + t7 / 0.2e1) * t178 - (t9 / 0.2e1 + t10 / 0.2e1) * t177 + (t182 * t243 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5) * t162 / 0.2e1 - Ifges(5,6) * t159 / 0.2e1 - pkin(2) * mrSges(4,1) - (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t177 - (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t178 + (-mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t163 + (Ifges(4,5) - Ifges(3,6) - qJ(3) * mrSges(4,1) + t208 / 0.2e1 + t206 / 0.2e1 + (mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t160) * qJD(2) - t129 * t134 + t138 * t124 + t142 * t63 + t146 * t12 + t108 * t37 + t14 * t90 + t34 * t91 + t15 * t92 + t35 * t93 + t97 * t11 + t23 * t80 + t86 * t30 + t87 * t32 + t88 * t81 + t64 * t62 + t66 * t36 + qJ(3) * t67 + t56 * t29 + t57 * t31 - (t54 / 0.2e1 + t55 / 0.2e1) * t79 + (-t52 / 0.2e1 - t53 / 0.2e1) * t78 + m(7) * (t14 * t17 + t15 * t16 + t2 * t56 + t23 * t97 + t3 * t57 + t64 * t66) + m(6) * (t108 * t142 + t146 * t88 + t25 * t35 + t26 * t34 + t5 * t87 + t6 * t86); t124 * t244 + t159 * t125 - t162 * t126 + 0.2e1 * t142 * t81 + 0.2e1 * t146 * t37 + 0.2e1 * t97 * t36 + 0.2e1 * t64 * t80 - (t84 + t85) * t79 + (-t82 - t83) * t78 - (t40 + t41) * t177 - (t38 + t39) * t178 + (-t206 - t208) * qJD(4) + (t14 * t57 + t15 * t56 + t64 * t97) * t234 + (t142 * t146 + t34 * t87 + t35 * t86) * t235 + 0.2e1 * t169 * mrSges(7,3) + 0.2e1 * t168 * mrSges(6,3) + (0.2e1 * mrSges(4,3) + 0.2e1 * t134 + (m(4) + m(5)) * t244) * qJD(3); t159 * t94 + t162 * t95 - t220 * t79 + t221 * t78 - (t29 + t30) * t177 + t224 * t178 + t179 * qJD(4) + (mrSges(4,1) + t243) * t203 - m(7) * t171 - m(6) * t170 + m(5) * (-t180 * qJD(4) + t239); -m(6) * t168 - m(7) * t169 + 0.2e1 * (mrSges(7,3) + mrSges(6,3)) * t181; -0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t181; t197 * t149 + (t161 * t30 + t224 * t158 + (-t220 * t158 + t221 * t161) * qJD(5) + m(6) * (t158 * t5 + t161 * t6 + t26 * t198 - t25 * t199) + m(7) * (t158 * t3 - t16 * t199 + t17 * t198)) * pkin(4) + t165 - t27 * mrSges(5,2) + t28 * mrSges(5,1) - Ifges(5,5) * t195 + t236; t187 * t149 + ((-mrSges(5,2) * t164 - Ifges(5,6)) * t162 + (-mrSges(5,1) * t164 - Ifges(5,5)) * t159) * qJD(4) + (m(7) * (t14 * t158 + t57 * t198 - t56 * t199) + m(6) * (t158 * t34 + t161 * t35 + t87 * t198 - t86 * t199) + t166 * mrSges(7,3) + (t166 + t213) * mrSges(6,3)) * pkin(4) + t167; -t134 * qJD(4) + m(6) * ((t194 - t213) * pkin(4) + t210) + m(7) * (pkin(4) * t194 - t149 * t79 + t210) + t176; 0.2e1 * (t191 + ((-t149 + t242) * m(7) - t249) * t158) * t240; t197 * pkin(5) + t165; t187 * pkin(5) + t167; -t79 * t230 + t176; (t191 + (-t230 - t249) * t158) * t240; 0; m(7) * t23 + t11; m(7) * t64 + t36; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
