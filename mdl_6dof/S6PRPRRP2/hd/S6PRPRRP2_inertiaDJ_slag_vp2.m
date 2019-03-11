% Calculate time derivative of joint inertia matrix for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:00:01
% EndTime: 2019-03-08 20:00:05
% DurationCPUTime: 2.09s
% Computational Cost: add. (1669->321), mult. (4607->466), div. (0->0), fcn. (4005->10), ass. (0->149)
t172 = m(6) / 0.2e1;
t180 = t172 + m(7) / 0.2e1;
t188 = 0.2e1 * t180;
t155 = Ifges(7,4) + Ifges(6,5);
t187 = -Ifges(7,2) - Ifges(6,3);
t92 = sin(qJ(5));
t95 = cos(qJ(5));
t143 = t92 ^ 2 + t95 ^ 2;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t94 = sin(qJ(2));
t97 = cos(qJ(2));
t112 = t89 * t94 - t91 * t97;
t90 = sin(pkin(6));
t181 = t112 * t90;
t40 = (t89 * t97 + t91 * t94) * t90;
t38 = qJD(2) * t40;
t186 = t181 * t38;
t96 = cos(qJ(4));
t140 = qJD(4) * t96;
t127 = t92 * t140;
t137 = qJD(5) * t95;
t93 = sin(qJ(4));
t128 = t93 * t137;
t103 = t127 + t128;
t80 = pkin(2) * t89 + pkin(8);
t159 = t80 * t96;
t81 = -pkin(2) * t91 - pkin(3);
t51 = -pkin(4) * t96 - pkin(9) * t93 + t81;
t182 = t159 * t95 + t51 * t92;
t185 = qJD(5) * t182;
t184 = m(6) + m(7);
t183 = mrSges(7,2) + mrSges(6,3);
t179 = m(7) * qJ(6) + mrSges(7,3);
t64 = (pkin(4) * t93 - pkin(9) * t96) * qJD(4);
t178 = -t64 * t95 + t185;
t114 = pkin(5) * t95 + qJ(6) * t92;
t136 = qJD(6) * t95;
t177 = qJD(5) * t114 - t136;
t139 = qJD(5) * t92;
t141 = qJD(4) * t93;
t151 = t137 * t51 + t64 * t92;
t10 = (-t139 * t96 - t141 * t95) * t80 + t151;
t157 = t93 * t95;
t60 = -mrSges(6,1) * t96 - mrSges(6,3) * t157;
t61 = mrSges(7,1) * t96 + mrSges(7,2) * t157;
t146 = -t60 + t61;
t35 = -mrSges(6,2) * t141 - mrSges(6,3) * t103;
t36 = -mrSges(7,2) * t103 + mrSges(7,3) * t141;
t152 = t35 + t36;
t123 = t80 * t92 + pkin(5);
t161 = t51 * t95;
t22 = t123 * t96 - t161;
t26 = -t159 * t92 + t161;
t7 = (-t139 * t80 - qJD(6)) * t96 + (-t80 * t95 + qJ(6)) * t141 + t151;
t176 = t146 * qJD(5) + m(7) * (qJD(5) * t22 + t7) + m(6) * (-qJD(5) * t26 + t10) + t152;
t129 = t80 * t141;
t11 = t129 * t92 - t178;
t158 = t92 * t93;
t59 = mrSges(6,2) * t96 - mrSges(6,3) * t158;
t62 = -mrSges(7,2) * t158 - mrSges(7,3) * t96;
t147 = t59 + t62;
t126 = t95 * t140;
t138 = qJD(5) * t93;
t104 = -t138 * t92 + t126;
t33 = mrSges(6,1) * t141 - mrSges(6,3) * t104;
t34 = mrSges(7,2) * t126 + (-mrSges(7,1) * qJD(4) - mrSges(7,2) * t139) * t93;
t153 = -t33 + t34;
t21 = -qJ(6) * t96 + t182;
t8 = -t123 * t141 + t178;
t175 = -t147 * qJD(5) + m(7) * (-qJD(5) * t21 + t8) + m(6) * (-t11 - t185) + t153;
t174 = 0.2e1 * m(6);
t173 = 0.2e1 * t80;
t142 = cos(pkin(6));
t105 = t142 * t96 - t40 * t93;
t39 = qJD(2) * t181;
t13 = qJD(4) * t105 - t39 * t96;
t29 = t142 * t93 + t40 * t96;
t15 = t181 * t92 + t29 * t95;
t4 = qJD(5) * t15 + t13 * t92 - t38 * t95;
t170 = t4 * t92;
t14 = -t181 * t95 + t29 * t92;
t5 = -qJD(5) * t14 + t13 * t95 + t38 * t92;
t169 = t5 * t95;
t168 = Ifges(6,4) * t92;
t167 = Ifges(6,4) * t95;
t166 = Ifges(7,5) * t92;
t165 = Ifges(7,5) * t95;
t12 = qJD(4) * t29 - t39 * t93;
t164 = t12 * t93;
t163 = t12 * t96;
t162 = t13 * t96;
t6 = t105 * t12;
t156 = t96 * Ifges(6,6);
t23 = mrSges(7,1) * t103 - mrSges(7,3) * t104;
t24 = mrSges(6,1) * t103 + mrSges(6,2) * t104;
t154 = t23 + t24;
t117 = Ifges(7,1) * t95 + t166;
t44 = -Ifges(7,4) * t96 + t117 * t93;
t118 = Ifges(6,1) * t95 - t168;
t45 = -Ifges(6,5) * t96 + t118 * t93;
t150 = t44 + t45;
t119 = t92 * mrSges(7,1) - t95 * mrSges(7,3);
t49 = t119 * t93;
t120 = t92 * mrSges(6,1) + t95 * mrSges(6,2);
t50 = t120 * t93;
t149 = t49 + t50;
t52 = t119 * qJD(5);
t53 = t120 * qJD(5);
t148 = t52 + t53;
t67 = -t95 * mrSges(6,1) + t92 * mrSges(6,2);
t145 = t67 - mrSges(5,1);
t144 = t143 * pkin(9) * t140;
t66 = -t95 * mrSges(7,1) - t92 * mrSges(7,3);
t131 = t66 + t145;
t68 = -Ifges(7,3) * t95 + t166;
t69 = Ifges(6,2) * t95 + t168;
t125 = t68 / 0.2e1 - t69 / 0.2e1;
t70 = Ifges(7,1) * t92 - t165;
t71 = Ifges(6,1) * t92 + t167;
t124 = t70 / 0.2e1 + t71 / 0.2e1;
t115 = Ifges(7,3) * t92 + t165;
t42 = -Ifges(7,6) * t96 + t115 * t93;
t116 = -Ifges(6,2) * t92 + t167;
t43 = t116 * t93 - t156;
t121 = t42 - t43 + t156;
t113 = pkin(5) * t92 - qJ(6) * t95;
t111 = -Ifges(7,6) * t103 - t126 * t155 + t141 * t187;
t109 = t113 + t80;
t108 = -t105 * t140 + t164;
t98 = m(7) * t136 + (-m(7) * t114 + t66 + t67) * qJD(5);
t86 = Ifges(7,4) * t137;
t85 = Ifges(6,5) * t137;
t83 = Ifges(7,6) * t139;
t65 = -pkin(4) - t114;
t58 = t118 * qJD(5);
t57 = t117 * qJD(5);
t56 = t116 * qJD(5);
t55 = t115 * qJD(5);
t54 = (mrSges(5,1) * t93 + mrSges(5,2) * t96) * qJD(4);
t47 = qJD(5) * t113 - qJD(6) * t92;
t32 = t109 * t93;
t20 = -t71 * t138 + (Ifges(6,5) * t93 + t118 * t96) * qJD(4);
t19 = -t70 * t138 + (Ifges(7,4) * t93 + t117 * t96) * qJD(4);
t18 = -t69 * t138 + (Ifges(6,6) * t93 + t116 * t96) * qJD(4);
t17 = -t68 * t138 + (Ifges(7,6) * t93 + t115 * t96) * qJD(4);
t16 = t109 * t140 + t177 * t93;
t3 = pkin(9) * t169;
t1 = [0.2e1 * m(5) * (t29 * t13 + t186 - t6) + 0.2e1 * m(4) * (-t40 * t39 + t186) + 0.2e1 * t184 * (t14 * t4 + t15 * t5 - t6); t39 * mrSges(4,2) + t147 * t5 + t146 * t4 + (-mrSges(5,1) * t96 + mrSges(5,2) * t93 - mrSges(4,1)) * t38 - t154 * t105 + t152 * t15 + t153 * t14 + t149 * t12 + (t112 * t54 + (-mrSges(3,1) * t94 - mrSges(3,2) * t97) * qJD(2)) * t90 + (t164 + t162 + (-t105 * t96 - t29 * t93) * qJD(4)) * mrSges(5,3) + m(6) * (t10 * t15 - t11 * t14 + t182 * t5 - t26 * t4) + m(5) * t38 * t81 + m(7) * (-t105 * t16 + t12 * t32 + t14 * t8 + t15 * t7 + t21 * t5 + t22 * t4) + (t108 * t172 + m(5) * (-t29 * t141 + t108 + t162) / 0.2e1) * t173 + m(4) * (-t38 * t91 - t39 * t89) * pkin(2); 0.2e1 * t10 * t59 + 0.2e1 * t11 * t60 + 0.2e1 * t16 * t49 + 0.2e1 * t21 * t36 + 0.2e1 * t22 * t34 + 0.2e1 * t32 * t23 + 0.2e1 * t26 * t33 + 0.2e1 * t182 * t35 + 0.2e1 * t81 * t54 + 0.2e1 * t8 * t61 + 0.2e1 * t7 * t62 + (t10 * t182 + t11 * t26) * t174 + 0.2e1 * m(7) * (t16 * t32 + t21 * t7 + t22 * t8) + ((0.2e1 * Ifges(5,4) * t96 + t121 * t92 + t150 * t95 + t173 * t50) * qJD(4) + t111) * t96 + (t24 * t173 + (t19 + t20) * t95 + (t17 - t18) * t92 + (t121 * t95 + (t155 * t96 - t150) * t92) * qJD(5) + ((-0.2e1 * Ifges(5,4) + t155 * t95 + (-Ifges(6,6) + Ifges(7,6)) * t92) * t93 + (t174 * t80 ^ 2 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) + t187) * t96) * qJD(4)) * t93; m(5) * (-t163 + t13 * t93 + (-t105 * t93 + t29 * t96) * qJD(4)) + t184 * (-t105 * t141 + t126 * t15 + t128 * t14 + t157 * t5 - t163) + t92 * (-t138 * t15 + t14 * t140 + t4 * t93) * t188; (-m(7) * t16 + (t147 * t95 + t146 * t92 + m(7) * (t21 * t95 + t22 * t92) + (t182 * t95 - t26 * t92 - t159) * m(6)) * qJD(4) - t154) * t96 + (m(6) * t129 + (m(7) * t32 + t149) * qJD(4) + t176 * t95 + t175 * t92) * t93; 0.4e1 * t180 * (-0.1e1 + t143) * t93 * t140; -t13 * mrSges(5,2) - t148 * t105 + t131 * t12 + m(6) * (-pkin(4) * t12 + t3) + m(7) * (-t105 * t47 + t12 * t65 + t3) + t183 * (t170 + t169 + (t14 * t95 - t15 * t92) * qJD(5)) + pkin(9) * (t137 * t14 - t139 * t15 + t170) * t188; t47 * t49 + t32 * t52 + t65 * t23 + t16 * t66 - pkin(4) * t24 + m(7) * (t16 * t65 + t32 * t47) + (-t85 / 0.2e1 - t86 / 0.2e1 - t83 / 0.2e1 + (Ifges(5,5) + (-m(6) * pkin(4) + t145) * t80) * qJD(4)) * t96 + (t8 * mrSges(7,2) - t11 * mrSges(6,3) + t19 / 0.2e1 + t20 / 0.2e1 + t125 * t140 + (t156 / 0.2e1 + t42 / 0.2e1 - t43 / 0.2e1 - t21 * mrSges(7,2) - t182 * mrSges(6,3)) * qJD(5) + t175 * pkin(9)) * t92 + (-t17 / 0.2e1 + t18 / 0.2e1 + t7 * mrSges(7,2) + t10 * mrSges(6,3) + t124 * t140 + (t44 / 0.2e1 + t45 / 0.2e1 + t22 * mrSges(7,2) - t26 * mrSges(6,3)) * qJD(5) + t176 * pkin(9)) * t95 + (t80 * t53 + (t57 / 0.2e1 + t58 / 0.2e1) * t95 + (t55 / 0.2e1 - t56 / 0.2e1) * t92 + (t80 * mrSges(5,2) - Ifges(5,6) + (Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t95 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t92) * qJD(4) + (-t124 * t92 + t125 * t95) * qJD(5)) * t93; t131 * t141 + m(6) * (-pkin(4) * t141 + t144) + m(7) * (t141 * t65 + t144) + (-m(7) * t47 + (t143 * t183 - mrSges(5,2)) * qJD(4) - t148) * t96; -0.2e1 * pkin(4) * t53 + 0.2e1 * t52 * t65 + (-t55 + t56) * t95 + (t57 + t58) * t92 + 0.2e1 * (m(7) * t65 + t66) * t47 + ((t70 + t71) * t95 + (t68 - t69) * t92) * qJD(5); m(7) * qJD(6) * t15 + (-mrSges(6,2) + t179) * t5 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t4; -Ifges(6,6) * t127 + m(7) * (-pkin(5) * t8 + qJ(6) * t7 + qJD(6) * t21) + qJD(6) * t62 + qJ(6) * t36 - t10 * mrSges(6,2) + t11 * mrSges(6,1) - t8 * mrSges(7,1) - pkin(5) * t34 + t7 * mrSges(7,3) + (-Ifges(6,6) * t95 - t155 * t92) * t138 - t111; (-m(7) * t113 - t119 - t120) * t140 + t98 * t93; -mrSges(7,2) * t177 - Ifges(6,6) * t139 + pkin(9) * t98 + t83 + t85 + t86; 0.2e1 * t179 * qJD(6); m(7) * t4; m(7) * t8 + t34; t103 * m(7); (m(7) * pkin(9) + mrSges(7,2)) * t137; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
