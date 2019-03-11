% Calculate time derivative of joint inertia matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:47
% EndTime: 2019-03-09 03:57:52
% DurationCPUTime: 2.99s
% Computational Cost: add. (4682->363), mult. (9548->549), div. (0->0), fcn. (9049->8), ass. (0->159)
t116 = sin(qJ(5));
t159 = qJD(5) * t116;
t114 = sin(pkin(10));
t117 = sin(qJ(3));
t163 = cos(pkin(10));
t187 = cos(qJ(3));
t136 = t163 * t187;
t91 = t114 * t117 - t136;
t148 = t91 * t159;
t119 = cos(qJ(5));
t92 = t114 * t187 + t117 * t163;
t86 = t92 * qJD(3);
t166 = t119 * t86;
t125 = t148 - t166;
t73 = t91 * t86;
t115 = sin(qJ(6));
t118 = cos(qJ(6));
t94 = t115 * t119 + t116 * t118;
t54 = t94 * t91;
t204 = -Ifges(4,1) + Ifges(4,2);
t128 = t115 * t116 - t118 * t119;
t197 = qJD(5) + qJD(6);
t69 = t197 * t128;
t70 = t197 * t94;
t36 = t70 * mrSges(7,1) - t69 * mrSges(7,2);
t134 = mrSges(6,1) * t116 + mrSges(6,2) * t119;
t95 = t134 * qJD(5);
t203 = t36 + t95;
t55 = t128 * t92;
t174 = Ifges(6,4) * t119;
t102 = Ifges(6,1) * t116 + t174;
t132 = -Ifges(6,2) * t116 + t174;
t202 = (t102 + t132) * qJD(5);
t120 = -pkin(1) - pkin(7);
t162 = qJ(4) - t120;
t158 = qJD(5) * t119;
t169 = t116 * t86;
t126 = t91 * t158 + t169;
t143 = qJD(3) * t187;
t161 = qJD(3) * t117;
t201 = -mrSges(4,1) * t161 - mrSges(4,2) * t143;
t85 = -qJD(3) * t136 + t114 * t161;
t39 = -mrSges(6,1) * t85 - mrSges(6,3) * t125;
t40 = mrSges(6,2) * t85 + mrSges(6,3) * t126;
t200 = -t116 * t39 + t119 * t40;
t199 = t114 * t85 + t163 * t86;
t123 = -qJD(4) * t187 + t161 * t162;
t124 = t162 * t187;
t79 = -qJD(3) * t124 - t117 * qJD(4);
t50 = t114 * t123 + t163 * t79;
t103 = pkin(3) * t143 + qJD(2);
t51 = -pkin(4) * t85 + pkin(8) * t86 + t103;
t108 = t117 * pkin(3) + qJ(2);
t62 = pkin(4) * t92 + pkin(8) * t91 + t108;
t98 = t162 * t117;
t72 = -t114 * t124 - t163 * t98;
t11 = t116 * t51 + t119 * t50 + t62 * t158 - t159 * t72;
t140 = -t116 * t50 + t119 * t51;
t65 = t119 * t72;
t35 = t116 * t62 + t65;
t12 = -qJD(5) * t35 + t140;
t198 = t11 * t119 - t12 * t116;
t138 = (-t116 ^ 2 - t119 ^ 2) * t85;
t100 = -mrSges(6,1) * t119 + mrSges(6,2) * t116;
t107 = -pkin(3) * t163 - pkin(4);
t196 = m(6) * t107 - mrSges(5,1) + t100;
t195 = 2 * m(7);
t194 = -2 * mrSges(5,3);
t193 = 0.2e1 * qJD(2);
t192 = m(5) * pkin(3);
t191 = m(7) * pkin(5);
t49 = t114 * t79 - t163 * t123;
t71 = -t114 * t98 + t163 * t124;
t186 = t49 * t71;
t185 = t50 * t92;
t184 = t69 * mrSges(7,3);
t183 = t70 * mrSges(7,3);
t182 = t71 * t86;
t180 = t128 * mrSges(7,3);
t179 = t94 * mrSges(7,3);
t106 = pkin(3) * t114 + pkin(8);
t178 = pkin(9) + t106;
t177 = -Ifges(7,5) * t69 - Ifges(7,6) * t70;
t175 = Ifges(6,4) * t116;
t173 = Ifges(6,6) * t116;
t168 = t116 * t91;
t165 = t119 * t91;
t157 = qJD(6) * t115;
t156 = qJD(6) * t118;
t155 = 0.2e1 * t85 * mrSges(5,3);
t154 = t91 * t194;
t19 = t128 * t86 + t197 * t54;
t21 = -t69 * t91 + t86 * t94;
t153 = Ifges(7,5) * t19 + Ifges(7,6) * t21 - Ifges(7,3) * t85;
t151 = pkin(5) * t159;
t145 = -t85 * mrSges(5,1) - t86 * mrSges(5,2);
t18 = t128 * t85 - t70 * t92;
t20 = t197 * t55 + t94 * t85;
t144 = t20 * mrSges(7,1) - t18 * mrSges(7,2);
t141 = t158 / 0.2e1;
t34 = -t116 * t72 + t119 * t62;
t139 = qJD(5) * t178;
t135 = t49 * t91 + t182;
t133 = Ifges(6,1) * t119 - t175;
t22 = pkin(5) * t92 + pkin(9) * t165 + t34;
t27 = pkin(9) * t168 + t35;
t9 = -t115 * t27 + t118 * t22;
t10 = t115 * t22 + t118 * t27;
t88 = t178 * t116;
t89 = t178 * t119;
t60 = -t115 * t89 - t118 * t88;
t61 = -t115 * t88 + t118 * t89;
t131 = t116 * t34 - t119 * t35;
t63 = -mrSges(6,2) * t92 + mrSges(6,3) * t168;
t64 = mrSges(6,1) * t92 + mrSges(6,3) * t165;
t130 = t116 * t64 - t119 * t63;
t83 = t116 * t139;
t84 = t119 * t139;
t32 = qJD(6) * t60 - t115 * t84 - t118 * t83;
t33 = -qJD(6) * t61 + t115 * t83 - t118 * t84;
t129 = t33 * mrSges(7,1) - t32 * mrSges(7,2) + t177;
t7 = pkin(9) * t166 - pkin(5) * t85 + (-t65 + (-pkin(9) * t91 - t62) * t116) * qJD(5) + t140;
t8 = pkin(9) * t126 + t11;
t2 = qJD(6) * t9 + t115 * t7 + t118 * t8;
t3 = -qJD(6) * t10 - t115 * t8 + t118 * t7;
t127 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t153;
t122 = (-t116 * t35 - t119 * t34) * qJD(5) + t198;
t121 = t125 * Ifges(6,5) + Ifges(6,6) * t126 - Ifges(6,3) * t85;
t109 = Ifges(6,5) * t158;
t101 = Ifges(6,2) * t119 + t175;
t99 = -t119 * pkin(5) + t107;
t97 = t133 * qJD(5);
t90 = (-mrSges(7,1) * t115 - mrSges(7,2) * t118) * qJD(6) * pkin(5);
t76 = Ifges(7,1) * t94 - Ifges(7,4) * t128;
t75 = Ifges(7,4) * t94 - Ifges(7,2) * t128;
t74 = mrSges(7,1) * t128 + mrSges(7,2) * t94;
t59 = t134 * t91;
t56 = t128 * t91;
t53 = t94 * t92;
t46 = -pkin(5) * t168 + t71;
t44 = t92 * Ifges(6,5) - t133 * t91;
t43 = t92 * Ifges(6,6) - t132 * t91;
t42 = mrSges(7,1) * t92 - mrSges(7,3) * t56;
t41 = -mrSges(7,2) * t92 + mrSges(7,3) * t54;
t38 = -Ifges(7,1) * t69 - Ifges(7,4) * t70;
t37 = -Ifges(7,4) * t69 - Ifges(7,2) * t70;
t31 = -mrSges(6,1) * t126 + mrSges(6,2) * t125;
t29 = -mrSges(7,1) * t54 + mrSges(7,2) * t56;
t28 = -pkin(5) * t126 + t49;
t26 = Ifges(7,1) * t56 + Ifges(7,4) * t54 + Ifges(7,5) * t92;
t25 = Ifges(7,4) * t56 + Ifges(7,2) * t54 + Ifges(7,6) * t92;
t24 = Ifges(6,1) * t125 + Ifges(6,4) * t126 - t85 * Ifges(6,5);
t23 = Ifges(6,4) * t125 + Ifges(6,2) * t126 - t85 * Ifges(6,6);
t14 = mrSges(7,2) * t85 + mrSges(7,3) * t21;
t13 = -mrSges(7,1) * t85 - mrSges(7,3) * t19;
t6 = -mrSges(7,1) * t21 + mrSges(7,2) * t19;
t5 = Ifges(7,1) * t19 + Ifges(7,4) * t21 - t85 * Ifges(7,5);
t4 = Ifges(7,4) * t19 + Ifges(7,2) * t21 - t85 * Ifges(7,6);
t1 = [((m(4) + m(3)) * t193 + 0.2e1 * (mrSges(4,1) * t187 - mrSges(4,2) * t117) * qJD(3)) * qJ(2) + 0.2e1 * Ifges(5,1) * t73 + (-0.2e1 * Ifges(4,4) * t187 + t204 * t117) * t143 + (0.2e1 * Ifges(4,4) * t117 + t204 * t187) * t161 + t125 * t44 + t126 * t43 + (t10 * t2 + t28 * t46 + t3 * t9) * t195 + (t117 * mrSges(4,1) + mrSges(4,2) * t187 + mrSges(3,3)) * t193 + 0.2e1 * (-t85 * t91 + t86 * t92) * Ifges(5,4) + 0.2e1 * m(6) * (t11 * t35 + t12 * t34 + t186) + 0.2e1 * m(5) * (t103 * t108 + t50 * t72 + t186) - t24 * t165 + t92 * t153 + 0.2e1 * t108 * t145 + t23 * t168 + t92 * t121 + 0.2e1 * t103 * (mrSges(5,1) * t92 - mrSges(5,2) * t91) + 0.2e1 * t71 * t31 + t56 * t5 + 0.2e1 * t11 * t63 + 0.2e1 * t12 * t64 + t54 * t4 + 0.2e1 * t34 * t39 + 0.2e1 * t35 * t40 + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t42 + 0.2e1 * t46 * t6 + t21 * t25 + t19 * t26 + 0.2e1 * t28 * t29 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + (-(-Ifges(6,5) * t119 + t173) * t91 - Ifges(7,5) * t56 - Ifges(7,6) * t54 + (-(2 * Ifges(5,2)) - Ifges(6,3) - Ifges(7,3)) * t92) * t85 + t72 * t155 + (t185 + t182) * t194 + (-0.2e1 * t59 + t154) * t49; -t53 * t13 - t55 * t14 + t18 * t41 + t20 * t42 + (t31 + t6) * t91 + t130 * t85 + (t29 - t59 + t154) * t86 + (t155 + (-t116 * t63 - t119 * t64) * qJD(5) + t200) * t92 + m(7) * (t10 * t18 - t2 * t55 + t20 * t9 + t28 * t91 - t3 * t53 + t46 * t86) + m(6) * (t122 * t92 + t131 * t85 + t135) + m(5) * (-t72 * t85 + t135 + t185); 0.2e1 * m(7) * (-t18 * t55 - t20 * t53 + t73) + 0.2e1 * m(6) * (t138 * t92 + t73) + 0.2e1 * m(5) * (-t85 * t92 + t73); -t128 * t4 / 0.2e1 - (Ifges(6,5) * t116 + Ifges(7,5) * t94 + Ifges(6,6) * t119 - Ifges(7,6) * t128) * t85 / 0.2e1 + t199 * mrSges(5,3) * pkin(3) + (-t163 * t192 + t196) * t49 + (m(6) * t122 - t158 * t64 - t159 * t63 + t200) * t106 + t201 * t120 + t202 * t168 / 0.2e1 - t2 * t180 - t10 * t183 + (-t158 * t34 - t159 * t35 + t198) * mrSges(6,3) + t9 * t184 - t3 * t179 - t97 * t165 / 0.2e1 - t102 * t166 / 0.2e1 - Ifges(4,5) * t161 - t43 * t159 / 0.2e1 + m(7) * (t10 * t32 + t151 * t46 + t2 * t61 + t28 * t99 + t3 * t60 + t33 * t9) - Ifges(4,6) * t143 + t44 * t141 + (t114 * t192 - mrSges(5,2)) * t50 + t119 * t23 / 0.2e1 + t116 * t24 / 0.2e1 + t94 * t5 / 0.2e1 + t71 * t95 + t99 * t6 + t107 * t31 - Ifges(5,5) * t86 + Ifges(5,6) * t85 - t69 * t26 / 0.2e1 - t70 * t25 / 0.2e1 + t28 * t74 + t21 * t75 / 0.2e1 + t19 * t76 / 0.2e1 + t56 * t38 / 0.2e1 + t60 * t13 + t61 * t14 + t54 * t37 / 0.2e1 + t32 * t41 + t33 * t42 + t46 * t36 + (t169 / 0.2e1 + t91 * t141) * t101 + t29 * t151 + (-Ifges(6,6) * t159 + t109 + t177) * t92 / 0.2e1; m(7) * (pkin(5) * t148 + t18 * t61 + t20 * t60 - t32 * t55 - t33 * t53) - t18 * t180 + t55 * t183 - t20 * t179 - t53 * t184 - t199 * t192 + t85 * mrSges(5,2) + t203 * t91 + (m(7) * t99 + t196 + t74) * t86 + t201 + (m(6) * t106 + mrSges(6,3)) * t138; (t151 * t99 + t32 * t61 + t33 * t60) * t195 - t70 * t75 - t128 * t37 + 0.2e1 * t74 * t151 + 0.2e1 * t99 * t36 - t69 * t76 + t94 * t38 + t116 * t97 - t101 * t159 + 0.2e1 * t107 * t95 + t202 * t119 + 0.2e1 * (-t128 * t32 - t33 * t94 + t60 * t69 - t61 * t70) * mrSges(7,3); t116 * t40 + t119 * t39 - t128 * t13 + t94 * t14 - t69 * t41 - t70 * t42 - t130 * qJD(5) + m(7) * (-t10 * t69 - t128 * t3 + t2 * t94 - t70 * t9) + m(6) * (-qJD(5) * t131 + t11 * t116 + t119 * t12) + m(5) * t103 + t145; m(7) * (-t128 * t20 + t18 * t94 + t53 * t70 + t55 * t69); m(7) * (-t128 * t33 + t32 * t94 - t60 * t70 - t61 * t69); (t128 * t70 - t69 * t94) * t195; t12 * mrSges(6,1) - t11 * mrSges(6,2) + (t41 * t156 + t115 * t14 - t42 * t157 + t118 * t13 + m(7) * (t10 * t156 + t115 * t2 + t118 * t3 - t157 * t9)) * pkin(5) + t121 + t127; (t119 * t85 + t159 * t92) * mrSges(6,2) + (t116 * t85 - t158 * t92) * mrSges(6,1) + (t115 * t18 + t118 * t20 + (t115 * t53 - t118 * t55) * qJD(6)) * t191 + t144; t109 + (t100 * t106 - t173) * qJD(5) + (m(7) * (t115 * t32 + t118 * t33 + (-t115 * t60 + t118 * t61) * qJD(6)) + (-t115 * t70 + t118 * t69 + (t115 * t94 - t118 * t128) * qJD(6)) * mrSges(7,3)) * pkin(5) + t129; (-t115 * t69 - t118 * t70 + (t115 * t128 + t118 * t94) * qJD(6)) * t191 - t203; 0.2e1 * t90; t127; t144; t129; -t36; t90; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
