% Calculate time derivative of joint inertia matrix for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:14
% EndTime: 2019-03-08 23:31:22
% DurationCPUTime: 3.81s
% Computational Cost: add. (3047->493), mult. (7798->703), div. (0->0), fcn. (6590->10), ass. (0->206)
t218 = Ifges(6,4) + Ifges(5,5);
t239 = -Ifges(6,2) - Ifges(5,3);
t141 = sin(qJ(3));
t140 = sin(qJ(4));
t145 = cos(qJ(3));
t187 = qJD(3) * t145;
t177 = t140 * t187;
t144 = cos(qJ(4));
t184 = qJD(4) * t144;
t151 = t141 * t184 + t177;
t208 = Ifges(6,5) * t140;
t164 = Ifges(6,1) * t144 + t208;
t210 = Ifges(5,4) * t140;
t165 = Ifges(5,1) * t144 - t210;
t238 = (t164 + t165) * qJD(4);
t139 = sin(qJ(6));
t143 = cos(qJ(6));
t159 = t139 * t144 - t140 * t143;
t80 = t159 * t141;
t237 = -mrSges(7,1) * t139 - mrSges(7,2) * t143;
t113 = -pkin(3) * t145 - pkin(9) * t141 - pkin(2);
t223 = pkin(8) * t140;
t128 = t145 * t223;
t137 = t145 * pkin(4);
t44 = pkin(5) * t145 + t128 + t137 + (-pkin(10) * t141 - t113) * t144;
t194 = t140 * t141;
t192 = t144 * t145;
t129 = pkin(8) * t192;
t72 = t140 * t113 + t129;
t62 = -qJ(5) * t145 + t72;
t47 = pkin(10) * t194 + t62;
t10 = -t139 * t47 + t143 * t44;
t107 = (pkin(3) * t141 - pkin(9) * t145) * qJD(3);
t186 = qJD(4) * t140;
t168 = qJD(4) * t129 - t107 * t144 + t113 * t186;
t185 = qJD(4) * t141;
t176 = t140 * t185;
t180 = -pkin(4) - t223;
t12 = pkin(10) * t176 + (-pkin(10) * t192 + (-pkin(5) + t180) * t141) * qJD(3) + t168;
t188 = qJD(3) * t141;
t131 = qJ(5) * t188;
t193 = t141 * t144;
t213 = t140 * t107 + t113 * t184;
t15 = t131 + (-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t193 + (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t140) * t145 + t213;
t1 = qJD(6) * t10 + t12 * t139 + t143 * t15;
t11 = t139 * t44 + t143 * t47;
t2 = -qJD(6) * t11 + t12 * t143 - t139 * t15;
t236 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t235 = qJD(4) - qJD(6);
t198 = qJ(5) * t140;
t228 = pkin(4) + pkin(5);
t234 = -t144 * t228 - t198;
t161 = pkin(4) * t144 + t198;
t182 = qJD(5) * t144;
t233 = qJD(4) * t161 - t182;
t232 = 2 * m(5);
t231 = 2 * m(7);
t230 = 0.2e1 * pkin(8);
t229 = m(5) / 0.2e1;
t227 = pkin(9) - pkin(10);
t108 = -t139 * qJ(5) - t143 * t228;
t78 = t143 * qJD(5) + qJD(6) * t108;
t222 = t78 * mrSges(7,2);
t109 = t143 * qJ(5) - t139 * t228;
t79 = -t139 * qJD(5) - qJD(6) * t109;
t221 = t79 * mrSges(7,1);
t138 = sin(pkin(6));
t146 = cos(qJ(2));
t189 = qJD(2) * t146;
t178 = t138 * t189;
t142 = sin(qJ(2));
t196 = t138 * t142;
t199 = cos(pkin(6));
t84 = t141 * t199 + t145 * t196;
t55 = qJD(3) * t84 + t141 * t178;
t83 = t141 * t196 - t145 * t199;
t220 = t83 * t55;
t219 = qJD(3) / 0.2e1;
t158 = t139 * t140 + t143 * t144;
t26 = t158 * t187 + t235 * t80;
t50 = t235 * t158;
t27 = t141 * t50 - t159 * t187;
t217 = Ifges(7,5) * t26 + Ifges(7,6) * t27;
t174 = t144 * t187;
t152 = t174 - t176;
t66 = mrSges(5,1) * t188 - mrSges(5,3) * t152;
t67 = mrSges(6,2) * t174 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t186) * t141;
t216 = -t66 + t67;
t68 = -mrSges(5,2) * t188 - mrSges(5,3) * t151;
t69 = -mrSges(6,2) * t151 + mrSges(6,3) * t188;
t215 = t68 + t69;
t75 = -Ifges(6,4) * t145 + t141 * t164;
t76 = -Ifges(5,5) * t145 + t141 * t165;
t214 = t75 + t76;
t209 = Ifges(5,4) * t144;
t207 = Ifges(6,5) * t144;
t206 = Ifges(5,6) * t144;
t179 = qJD(2) * t196;
t195 = t138 * t146;
t181 = t140 * t195;
t56 = -qJD(3) * t83 + t145 * t178;
t13 = -qJD(4) * t181 + t140 * t56 - t144 * t179 + t184 * t84;
t205 = t13 * t140;
t57 = t140 * t84 + t144 * t195;
t14 = -qJD(4) * t57 + t140 * t179 + t144 * t56;
t204 = t14 * t144;
t203 = t145 * Ifges(5,6);
t202 = t55 * t141;
t201 = t56 * t145;
t115 = -t144 * mrSges(5,1) + t140 * mrSges(5,2);
t200 = t115 - mrSges(4,1);
t197 = qJ(5) * t144;
t101 = mrSges(5,2) * t145 - mrSges(5,3) * t194;
t104 = -mrSges(6,2) * t194 - mrSges(6,3) * t145;
t191 = t101 + t104;
t102 = -mrSges(5,1) * t145 - mrSges(5,3) * t193;
t103 = mrSges(6,1) * t145 + mrSges(6,2) * t193;
t190 = -t102 + t103;
t183 = qJD(5) * t140;
t121 = t227 * t144;
t116 = -Ifges(6,3) * t144 + t208;
t117 = Ifges(5,2) * t144 + t210;
t173 = t116 / 0.2e1 - t117 / 0.2e1;
t118 = Ifges(6,1) * t140 - t207;
t119 = Ifges(5,1) * t140 + t209;
t172 = t118 / 0.2e1 + t119 / 0.2e1;
t71 = t113 * t144 - t128;
t162 = Ifges(6,3) * t140 + t207;
t73 = -Ifges(6,6) * t145 + t141 * t162;
t163 = -Ifges(5,2) * t140 + t209;
t74 = t141 * t163 - t203;
t170 = t73 - t74 + t203;
t58 = t144 * t84 - t181;
t19 = t139 * t57 + t143 * t58;
t3 = -qJD(6) * t19 + t13 * t143 - t139 * t14;
t18 = -t139 * t58 + t143 * t57;
t4 = qJD(6) * t18 + t13 * t139 + t14 * t143;
t169 = t3 * mrSges(7,1) - t4 * mrSges(7,2);
t167 = mrSges(5,1) * t140 + mrSges(5,2) * t144;
t114 = -t144 * mrSges(6,1) - t140 * mrSges(6,3);
t166 = mrSges(6,1) * t140 - mrSges(6,3) * t144;
t160 = pkin(4) * t140 - t197;
t120 = t227 * t140;
t60 = t120 * t143 - t121 * t139;
t61 = t120 * t139 + t121 * t143;
t157 = pkin(8) + t160;
t156 = t187 * t83 + t202;
t155 = -t140 * t228 + t197;
t105 = t227 * t186;
t106 = qJD(4) * t121;
t28 = qJD(6) * t60 - t105 * t143 + t106 * t139;
t29 = -qJD(6) * t61 + t105 * t139 + t106 * t143;
t51 = t235 * t159;
t48 = Ifges(7,6) * t51;
t49 = Ifges(7,5) * t50;
t154 = t29 * mrSges(7,1) - t28 * mrSges(7,2) - t48 + t49;
t153 = -pkin(8) + t155;
t150 = -t151 * Ifges(6,6) - t174 * t218 + t239 * t188 + t217;
t33 = (-t144 * t188 - t145 * t186) * pkin(8) + t213;
t136 = Ifges(6,4) * t184;
t135 = Ifges(5,5) * t184;
t133 = Ifges(6,6) * t186;
t110 = -pkin(3) - t161;
t98 = t163 * qJD(4);
t97 = t162 * qJD(4);
t96 = (mrSges(4,1) * t141 + mrSges(4,2) * t145) * qJD(3);
t95 = t167 * qJD(4);
t94 = t166 * qJD(4);
t90 = pkin(3) - t234;
t87 = t167 * t141;
t86 = t166 * t141;
t82 = qJD(4) * t160 - t183;
t81 = t158 * t141;
t77 = t157 * t141;
t70 = qJD(4) * t155 + t183;
t65 = mrSges(7,1) * t145 - t81 * mrSges(7,3);
t64 = -mrSges(7,2) * t145 - t80 * mrSges(7,3);
t63 = t137 - t71;
t59 = t153 * t141;
t54 = -Ifges(7,1) * t159 - Ifges(7,4) * t158;
t53 = -Ifges(7,4) * t159 - Ifges(7,2) * t158;
t52 = mrSges(7,1) * t158 - mrSges(7,2) * t159;
t46 = mrSges(5,1) * t151 + mrSges(5,2) * t152;
t45 = mrSges(6,1) * t151 - mrSges(6,3) * t152;
t43 = mrSges(7,1) * t80 + mrSges(7,2) * t81;
t42 = -t119 * t185 + (Ifges(5,5) * t141 + t145 * t165) * qJD(3);
t41 = -t118 * t185 + (Ifges(6,4) * t141 + t145 * t164) * qJD(3);
t40 = -t117 * t185 + (Ifges(5,6) * t141 + t145 * t163) * qJD(3);
t39 = -t116 * t185 + (Ifges(6,6) * t141 + t145 * t162) * qJD(3);
t38 = Ifges(7,1) * t81 - Ifges(7,4) * t80 + Ifges(7,5) * t145;
t37 = Ifges(7,4) * t81 - Ifges(7,2) * t80 + Ifges(7,6) * t145;
t34 = t188 * t223 - t168;
t32 = t141 * t233 + t157 * t187;
t31 = t180 * t188 + t168;
t30 = -qJD(5) * t145 + t131 + t33;
t23 = (qJD(4) * t234 + t182) * t141 + t153 * t187;
t22 = Ifges(7,1) * t50 - Ifges(7,4) * t51;
t21 = Ifges(7,4) * t50 - Ifges(7,2) * t51;
t20 = mrSges(7,1) * t51 + mrSges(7,2) * t50;
t17 = mrSges(7,2) * t188 + t27 * mrSges(7,3);
t16 = -mrSges(7,1) * t188 - t26 * mrSges(7,3);
t9 = pkin(9) * t204;
t7 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t6 = Ifges(7,1) * t26 + Ifges(7,4) * t27 - Ifges(7,5) * t188;
t5 = Ifges(7,4) * t26 + Ifges(7,2) * t27 - Ifges(7,6) * t188;
t8 = [0.2e1 * m(7) * (t18 * t3 + t19 * t4 + t220) + 0.2e1 * m(4) * (-t138 ^ 2 * t142 * t189 + t84 * t56 + t220) + 0.2e1 * (m(6) + m(5)) * (t13 * t57 + t58 * t14 + t220); t18 * t16 + t19 * t17 + t3 * t65 + t4 * t64 + t215 * t58 + t216 * t57 + t191 * t14 + t190 * t13 + (t45 + t46 - t7) * t83 + (-t43 + t86 + t87) * t55 + (-t146 * t96 + (-t146 * mrSges(3,2) + (-mrSges(4,1) * t145 + mrSges(4,2) * t141 - mrSges(3,1)) * t142) * qJD(2)) * t138 + (t202 + t201 + (-t141 * t84 + t145 * t83) * qJD(3)) * mrSges(4,3) + m(6) * (t13 * t63 + t14 * t62 + t30 * t58 + t31 * t57 + t32 * t83 + t55 * t77) + m(7) * (t1 * t19 + t10 * t3 + t11 * t4 + t18 * t2 - t23 * t83 - t55 * t59) + m(5) * (-t71 * t13 + t72 * t14 + t33 * t58 - t34 * t57) - m(4) * pkin(2) * t179 + (t156 * t229 + m(4) * (-t84 * t188 + t156 + t201) / 0.2e1) * t230; (t46 * t230 + (t41 + t42) * t144 + (t39 - t40) * t140 + (t170 * t144 + (t145 * t218 - t214) * t140) * qJD(4) + (-Ifges(7,5) * t81 + Ifges(7,6) * t80 + (-(2 * Ifges(4,4)) + t218 * t144 + (-Ifges(5,6) + Ifges(6,6)) * t140) * t141 + (pkin(8) ^ 2 * t232 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(7,3)) + t239) * t145) * qJD(3)) * t141 + ((0.2e1 * Ifges(4,4) * t145 + t140 * t170 + t144 * t214 + t230 * t87) * qJD(3) + t150) * t145 + (t1 * t11 + t10 * t2 + t23 * t59) * t231 + (t72 * t33 + t71 * t34) * t232 + 0.2e1 * m(6) * (t30 * t62 + t31 * t63 + t32 * t77) + 0.2e1 * t10 * t16 + 0.2e1 * t11 * t17 + t27 * t37 + t26 * t38 + 0.2e1 * t23 * t43 + 0.2e1 * t59 * t7 + 0.2e1 * t1 * t64 + 0.2e1 * t2 * t65 + 0.2e1 * t63 * t67 + 0.2e1 * t62 * t69 + 0.2e1 * t71 * t66 + 0.2e1 * t72 * t68 + 0.2e1 * t77 * t45 - t80 * t5 + t81 * t6 + 0.2e1 * t32 * t86 - 0.2e1 * pkin(2) * t96 + 0.2e1 * t33 * t101 + 0.2e1 * t34 * t102 + 0.2e1 * t31 * t103 + 0.2e1 * t30 * t104; -t56 * mrSges(4,2) + (-t20 + t94 + t95) * t83 + (t114 - t52 + t200) * t55 + m(6) * (t110 * t55 + t82 * t83 + t9) + m(7) * (t29 * t18 + t28 * t19 + t60 * t3 + t61 * t4 - t55 * t90 - t70 * t83) + m(5) * (-pkin(3) * t55 + t9) + 0.2e1 * (m(6) / 0.2e1 + t229) * (t184 * t57 - t186 * t58 + t205) * pkin(9) + (-t158 * t4 + t159 * t3 - t18 * t50 - t19 * t51) * mrSges(7,3) + (mrSges(5,3) + mrSges(6,2)) * (t205 + t204 + (-t140 * t58 + t144 * t57) * qJD(4)); m(7) * (t61 * t1 + t29 * t10 + t28 * t11 + t60 * t2 + t23 * t90 + t70 * t59) + (t30 * mrSges(6,2) + t33 * mrSges(5,3) - t39 / 0.2e1 + t40 / 0.2e1 + t172 * t187 + (t63 * mrSges(6,2) - t71 * mrSges(5,3) + t75 / 0.2e1 + t76 / 0.2e1) * qJD(4) + (t190 * qJD(4) + m(6) * (qJD(4) * t63 + t30) + m(5) * (-qJD(4) * t71 + t33) + t215) * pkin(9)) * t144 + (t31 * mrSges(6,2) - t34 * mrSges(5,3) + t41 / 0.2e1 + t42 / 0.2e1 + t173 * t187 + (-t62 * mrSges(6,2) - t72 * mrSges(5,3) + t73 / 0.2e1 - t74 / 0.2e1 + t203 / 0.2e1) * qJD(4) + (-t191 * qJD(4) + m(6) * (-qJD(4) * t62 + t31) + m(5) * (-qJD(4) * t72 - t34) + t216) * pkin(9)) * t140 + (-t1 * t158 - t10 * t50 - t11 * t51 + t159 * t2) * mrSges(7,3) - t158 * t5 / 0.2e1 + (-qJD(3) * (-Ifges(7,5) * t159 - Ifges(7,6) * t158) / 0.2e1 - Ifges(4,6) * qJD(3) + t206 * t219 + (mrSges(4,2) * qJD(3) + t95) * pkin(8) + (t97 / 0.2e1 - t98 / 0.2e1 - t172 * qJD(4) + t218 * t219) * t140 + (t238 / 0.2e1 - Ifges(6,6) * t219 + t173 * qJD(4)) * t144) * t141 - t159 * t6 / 0.2e1 + m(6) * (t110 * t32 + t82 * t77) + (-t135 / 0.2e1 + t49 / 0.2e1 - t48 / 0.2e1 - t136 / 0.2e1 - t133 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t200) * pkin(8)) * qJD(3)) * t145 - pkin(3) * t46 + t50 * t38 / 0.2e1 - t51 * t37 / 0.2e1 + t23 * t52 + t27 * t53 / 0.2e1 + t26 * t54 / 0.2e1 + t59 * t20 + t60 * t16 + t61 * t17 + t28 * t64 + t29 * t65 + t70 * t43 - t80 * t21 / 0.2e1 + t81 * t22 / 0.2e1 + t82 * t86 + t90 * t7 + t77 * t94 + t110 * t45 + t32 * t114; (t61 * t28 + t60 * t29 + t70 * t90) * t231 + 0.2e1 * t70 * t52 + 0.2e1 * t90 * t20 + t50 * t54 - t159 * t22 - t51 * t53 - t158 * t21 + 0.2e1 * t82 * t114 - 0.2e1 * pkin(3) * t95 + (t98 - t97) * t144 + t238 * t140 + ((t118 + t119) * t144 + (t116 - t117) * t140) * qJD(4) + 0.2e1 * (m(6) * t82 + t94) * t110 + 0.2e1 * (-t158 * t28 + t159 * t29 - t50 * t60 - t51 * t61) * mrSges(7,3); (-mrSges(5,2) + mrSges(6,3)) * t14 + (-mrSges(5,1) - mrSges(6,1)) * t13 + m(6) * (-pkin(4) * t13 + qJ(5) * t14 + qJD(5) * t58) + m(7) * (t108 * t3 + t109 * t4 + t79 * t18 + t78 * t19) - t169; -t150 + (Ifges(7,3) * qJD(3) + (-t140 * t218 - t206) * qJD(4)) * t141 + m(6) * (-pkin(4) * t31 + qJ(5) * t30 + qJD(5) * t62) + m(7) * (t1 * t109 + t79 * t10 + t108 * t2 + t78 * t11) - Ifges(5,6) * t177 + t30 * mrSges(6,3) - t31 * mrSges(6,1) - t33 * mrSges(5,2) + t34 * mrSges(5,1) - pkin(4) * t67 + qJ(5) * t69 + t78 * t64 + t79 * t65 + qJD(5) * t104 + t108 * t16 + t109 * t17 + t236; m(7) * (t108 * t29 + t109 * t28 + t79 * t60 + t78 * t61) + t136 + t133 + t135 - Ifges(5,6) * t186 - t233 * mrSges(6,2) + (-t108 * t50 - t109 * t51 - t158 * t78 + t159 * t79) * mrSges(7,3) + (m(6) * t182 + (-m(6) * t161 + t114 + t115) * qJD(4)) * pkin(9) - t154; (t108 * t79 + t109 * t78) * t231 + 0.2e1 * t222 - 0.2e1 * t221 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t13 + m(7) * (t139 * t4 + t143 * t3 + (-t139 * t18 + t143 * t19) * qJD(6)); t139 * t17 + t143 * t16 + (-t139 * t65 + t143 * t64) * qJD(6) + m(7) * (t1 * t139 + t143 * t2 + (-t10 * t139 + t11 * t143) * qJD(6)) + m(6) * t31 + t67; m(7) * (t139 * t28 + t143 * t29 + (-t139 * t60 + t143 * t61) * qJD(6)) + (m(6) * pkin(9) + mrSges(6,2)) * t184 + (-t139 * t51 - t143 * t50 + (-t139 * t159 - t143 * t158) * qJD(6)) * mrSges(7,3); m(7) * (t139 * t78 + t143 * t79) + (m(7) * (-t108 * t139 + t109 * t143) - t237) * qJD(6); 0; t169; -Ifges(7,3) * t188 + t217 - t236; t154; t221 - t222; t237 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
