% Calculate time derivative of joint inertia matrix for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:18
% EndTime: 2019-03-08 19:55:23
% DurationCPUTime: 2.19s
% Computational Cost: add. (1657->342), mult. (4584->495), div. (0->0), fcn. (3975->10), ass. (0->160)
t152 = Ifges(6,5) + Ifges(7,5);
t184 = -Ifges(6,3) - Ifges(7,3);
t95 = cos(qJ(5));
t130 = qJD(5) * t95;
t93 = sin(qJ(4));
t119 = t93 * t130;
t96 = cos(qJ(4));
t133 = qJD(4) * t96;
t92 = sin(qJ(5));
t99 = t133 * t92 + t119;
t136 = cos(pkin(6));
t89 = sin(pkin(11));
t90 = sin(pkin(6));
t91 = cos(pkin(11));
t94 = sin(qJ(2));
t97 = cos(qJ(2));
t34 = (t89 * t97 + t91 * t94) * t90;
t101 = t136 * t96 - t34 * t93;
t183 = t101 * qJD(4);
t173 = m(6) / 0.2e1;
t127 = t173 + m(7) / 0.2e1;
t182 = 0.2e1 * t127;
t80 = -pkin(5) * t95 - pkin(4);
t181 = m(7) * t80;
t132 = qJD(5) * t92;
t51 = mrSges(7,1) * t132 + mrSges(7,2) * t130;
t164 = mrSges(6,2) * t95;
t109 = mrSges(6,1) * t92 + t164;
t52 = t109 * qJD(5);
t180 = t51 + t52;
t79 = -pkin(2) * t91 - pkin(3);
t49 = -pkin(4) * t96 - pkin(9) * t93 + t79;
t154 = t95 * t96;
t78 = pkin(2) * t89 + pkin(8);
t62 = t78 * t154;
t22 = t49 * t92 + t62;
t156 = t92 * t96;
t41 = t95 * t49;
t21 = -t156 * t78 + t41;
t118 = t96 * t132;
t134 = qJD(4) * t93;
t167 = pkin(9) * t96;
t168 = pkin(4) * t93;
t63 = (-t167 + t168) * qJD(4);
t146 = t130 * t49 + t63 * t92;
t7 = (-t134 * t95 - t118) * t78 + t146;
t179 = -qJD(5) * t21 + t7;
t139 = t92 ^ 2 + t95 ^ 2;
t178 = -m(6) * pkin(9) - mrSges(6,3);
t177 = 0.2e1 * m(6);
t176 = 0.2e1 * m(7);
t175 = -2 * mrSges(7,3);
t174 = 0.2e1 * t78;
t172 = m(6) * pkin(4);
t170 = m(7) * pkin(5);
t26 = pkin(5) * t99 + t133 * t78;
t169 = m(7) * t26;
t24 = t136 * t93 + t34 * t96;
t104 = t89 * t94 - t91 * t97;
t135 = qJD(2) * t90;
t32 = t104 * t135;
t9 = qJD(4) * t24 - t32 * t93;
t166 = t101 * t9;
t165 = t9 * t93;
t163 = Ifges(6,4) * t92;
t162 = Ifges(6,4) * t95;
t161 = Ifges(7,4) * t92;
t160 = Ifges(7,4) * t95;
t10 = -t32 * t96 + t183;
t159 = t10 * t96;
t31 = qJD(2) * t34;
t33 = t104 * t90;
t158 = t31 * t33;
t157 = t92 * t93;
t155 = t93 * t95;
t153 = mrSges(6,2) + mrSges(7,2);
t151 = Ifges(6,6) + Ifges(7,6);
t150 = -qJ(6) - pkin(9);
t131 = qJD(5) * t93;
t120 = t92 * t131;
t121 = t95 * t133;
t125 = -mrSges(7,1) * t99 - mrSges(7,2) * t121;
t19 = -mrSges(7,2) * t120 - t125;
t100 = -t120 + t121;
t20 = mrSges(6,1) * t99 + mrSges(6,2) * t100;
t149 = t19 + t20;
t27 = mrSges(7,1) * t134 - mrSges(7,3) * t100;
t28 = mrSges(6,1) * t134 - mrSges(6,3) * t100;
t148 = t27 + t28;
t29 = -mrSges(7,2) * t134 - mrSges(7,3) * t99;
t30 = -mrSges(6,2) * t134 - mrSges(6,3) * t99;
t147 = t29 + t30;
t107 = Ifges(7,1) * t95 - t161;
t38 = -Ifges(7,5) * t96 + t107 * t93;
t108 = Ifges(6,1) * t95 - t163;
t39 = -Ifges(6,5) * t96 + t108 * t93;
t145 = t38 + t39;
t123 = t78 * t134;
t144 = t123 * t92 + t63 * t95;
t47 = (mrSges(7,1) * t92 + mrSges(7,2) * t95) * t93;
t48 = t109 * t93;
t143 = t47 + t48;
t58 = mrSges(7,2) * t96 - mrSges(7,3) * t157;
t59 = mrSges(6,2) * t96 - mrSges(6,3) * t157;
t142 = t58 + t59;
t60 = -mrSges(7,1) * t96 - mrSges(7,3) * t155;
t61 = -mrSges(6,1) * t96 - mrSges(6,3) * t155;
t141 = t60 + t61;
t140 = -mrSges(6,1) * t95 + mrSges(6,2) * t92 - mrSges(5,1);
t138 = qJ(6) * t93;
t137 = qJ(6) * t95;
t129 = qJD(6) * t95;
t65 = -mrSges(7,1) * t95 + mrSges(7,2) * t92;
t126 = t65 + t140;
t124 = pkin(5) * t132;
t117 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t68 = Ifges(7,2) * t95 + t161;
t69 = Ifges(6,2) * t95 + t163;
t116 = -t68 / 0.2e1 - t69 / 0.2e1;
t70 = Ifges(7,1) * t92 + t160;
t71 = Ifges(6,1) * t92 + t162;
t115 = t70 / 0.2e1 + t71 / 0.2e1;
t114 = -mrSges(6,1) - t170;
t113 = mrSges(7,1) + t170;
t112 = t151 * t92;
t111 = qJD(5) * t150;
t110 = -t121 * t152 + t134 * t184;
t106 = -Ifges(6,2) * t92 + t162;
t105 = -Ifges(7,2) * t92 + t160;
t12 = t24 * t95 + t33 * t92;
t11 = -t24 * t92 + t33 * t95;
t36 = -t96 * Ifges(7,6) + t105 * t93;
t37 = -t96 * Ifges(6,6) + t106 * t93;
t103 = t151 * t96 - t36 - t37;
t102 = -t101 * t133 + t165;
t3 = -qJD(5) * t12 - t10 * t92 + t31 * t95;
t4 = qJD(5) * t11 + t10 * t95 + t31 * t92;
t98 = -t3 * t92 + t4 * t95 + (-t11 * t95 - t12 * t92) * qJD(5);
t86 = Ifges(6,5) * t130;
t85 = Ifges(7,5) * t130;
t67 = t150 * t95;
t64 = t150 * t92;
t57 = t108 * qJD(5);
t56 = t107 * qJD(5);
t55 = t106 * qJD(5);
t54 = t105 * qJD(5);
t53 = (t93 * mrSges(5,1) + t96 * mrSges(5,2)) * qJD(4);
t44 = (pkin(5) * t92 + t78) * t93;
t43 = -qJD(6) * t92 + t111 * t95;
t42 = t111 * t92 + t129;
t18 = -t138 * t92 + t22;
t17 = -t71 * t131 + (Ifges(6,5) * t93 + t108 * t96) * qJD(4);
t16 = -t70 * t131 + (Ifges(7,5) * t93 + t107 * t96) * qJD(4);
t15 = -t69 * t131 + (Ifges(6,6) * t93 + t106 * t96) * qJD(4);
t14 = -t68 * t131 + (Ifges(7,6) * t93 + t105 * t96) * qJD(4);
t13 = -t93 * t137 + t41 + (-t78 * t92 - pkin(5)) * t96;
t8 = -qJD(5) * t22 + t144;
t6 = (-qJ(6) * qJD(5) - qJD(4) * t78) * t155 + (-qJD(6) * t93 + (-qJ(6) * qJD(4) - qJD(5) * t78) * t96) * t92 + t146;
t5 = -t93 * t129 + (pkin(5) * t93 - t137 * t96) * qJD(4) + (-t62 + (-t49 + t138) * t92) * qJD(5) + t144;
t1 = [0.2e1 * m(5) * (t10 * t24 + t158 - t166) + 0.2e1 * m(4) * (-t32 * t34 + t158) + 0.2e1 * (t11 * t3 + t12 * t4 - t166) * t182; t32 * mrSges(4,2) + t33 * t53 + t143 * t9 + t142 * t4 + (-mrSges(5,1) * t96 + mrSges(5,2) * t93 - mrSges(4,1)) * t31 + t141 * t3 - t149 * t101 + t147 * t12 + t148 * t11 + (-mrSges(3,1) * t94 - mrSges(3,2) * t97) * t135 + (t159 + t165 + (-t101 * t96 - t24 * t93) * qJD(4)) * mrSges(5,3) + m(6) * (t11 * t8 + t12 * t7 + t21 * t3 + t22 * t4) + m(5) * t31 * t79 + m(7) * (-t101 * t26 + t11 * t5 + t12 * t6 + t13 * t3 + t18 * t4 + t44 * t9) + (t102 * t173 + m(5) * (-t24 * t134 + t102 + t159) / 0.2e1) * t174 + m(4) * (-t31 * t91 - t32 * t89) * pkin(2); 0.2e1 * t13 * t27 + 0.2e1 * t18 * t29 + 0.2e1 * t44 * t19 + 0.2e1 * t21 * t28 + 0.2e1 * t22 * t30 + 0.2e1 * t26 * t47 + 0.2e1 * t5 * t60 + 0.2e1 * t79 * t53 + 0.2e1 * t6 * t58 + 0.2e1 * t7 * t59 + 0.2e1 * t8 * t61 + (t21 * t8 + t22 * t7) * t177 + (t13 * t5 + t18 * t6 + t26 * t44) * t176 + ((0.2e1 * Ifges(5,4) * t96 + t103 * t92 + t145 * t95 + t174 * t48) * qJD(4) + t110) * t96 + (t20 * t174 + (t16 + t17) * t95 + (-t14 - t15) * t92 + (t103 * t95 + (t152 * t96 - t145) * t92) * qJD(5) + ((t152 * t95 - 0.2e1 * Ifges(5,4) - t112) * t93 + (t177 * t78 ^ 2 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) + t184) * t96) * qJD(4)) * t93; m(5) * (t93 * t10 - t96 * t9 + (-t101 * t93 + t24 * t96) * qJD(4)) + ((-t9 + (-t11 * t92 + t12 * t95) * qJD(4)) * t96 + (t98 - t183) * t93) * t182; (-t169 + (t142 * t95 - t141 * t92 + m(7) * (-t13 * t92 + t18 * t95) + (-t21 * t92 + t22 * t95 - t78 * t96) * m(6)) * qJD(4) - t149) * t96 + (t147 * t95 - t148 * t92 + t143 * qJD(4) + (-t141 * t95 - t142 * t92) * qJD(5) + m(7) * (qJD(4) * t44 - t13 * t130 - t132 * t18 - t5 * t92 + t6 * t95) + (-t22 * t132 + t179 * t95 - t8 * t92 + t123) * m(6)) * t93; 0.4e1 * t127 * (-0.1e1 + t139) * t93 * t133; -t10 * mrSges(5,2) - t180 * t101 + m(7) * (-t101 * t124 + t11 * t43 + t12 * t42 + t3 * t64 - t4 * t67) + (t126 - t172 + t181) * t9 + (mrSges(7,3) - t178) * t98; m(7) * (t13 * t43 + t18 * t42 + t26 * t80 + t5 * t64 - t6 * t67) + t80 * t19 + t42 * t58 + t43 * t60 + t64 * t27 + t26 * t65 - t67 * t29 + t44 * t51 - pkin(4) * t20 + (-t85 / 0.2e1 - t86 / 0.2e1 + (Ifges(5,5) + (t140 - t172) * t78) * qJD(4)) * t96 + (-t5 * mrSges(7,3) - t8 * mrSges(6,3) + t16 / 0.2e1 + t17 / 0.2e1 + t116 * t133 + (-m(6) * t8 - t28) * pkin(9) + (pkin(5) * t47 - pkin(9) * t59 - t18 * mrSges(7,3) - t36 / 0.2e1 - t37 / 0.2e1 + t117 * t96 + t44 * t170 + t178 * t22) * qJD(5)) * t92 + (t6 * mrSges(7,3) + t7 * mrSges(6,3) + t14 / 0.2e1 + t15 / 0.2e1 + t115 * t133 + (-t13 * mrSges(7,3) - t21 * mrSges(6,3) + t39 / 0.2e1 + t38 / 0.2e1) * qJD(5) + (m(6) * t179 - qJD(5) * t61 + t30) * pkin(9)) * t95 + (t78 * t52 + (t56 / 0.2e1 + t57 / 0.2e1) * t95 + (-t54 / 0.2e1 - t55 / 0.2e1) * t92 + (t78 * mrSges(5,2) - Ifges(5,6) + t117 * t95 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t92) * qJD(4) + (-t115 * t92 + t116 * t95) * qJD(5)) * t93; m(7) * (-pkin(5) * t118 - t119 * t64 + t120 * t67 + t155 * t42 - t157 * t43) + (t126 * t93 + m(6) * (t139 * t167 - t168) + m(7) * (-t154 * t67 - t156 * t64 + t80 * t93)) * qJD(4) + ((-mrSges(5,2) + (mrSges(6,3) + mrSges(7,3)) * t139) * qJD(4) - t180) * t96; 0.2e1 * t80 * t51 - 0.2e1 * pkin(4) * t52 + (-t42 * t67 + t43 * t64) * t176 + (t43 * t175 + t56 + t57 + (-t67 * t175 - t68 - t69 + 0.2e1 * (t65 + t181) * pkin(5)) * qJD(5)) * t92 + (0.2e1 * t42 * mrSges(7,3) + t54 + t55 + (t175 * t64 + t70 + t71) * qJD(5)) * t95; -t153 * t4 + (mrSges(6,1) + t113) * t3; mrSges(6,1) * t8 + mrSges(7,1) * t5 - mrSges(6,2) * t7 - mrSges(7,2) * t6 - t112 * t133 + (m(7) * t5 + t27) * pkin(5) + (-t151 * t95 - t152 * t92) * t131 - t110; (t114 * t92 - t164) * t133 + (t114 * t95 + t153 * t92) * t131 + t125; -mrSges(7,2) * t42 + t85 + t86 + t113 * t43 + ((-mrSges(6,1) * pkin(9) - mrSges(7,3) * pkin(5)) * t95 + (mrSges(6,2) * pkin(9) - t151) * t92) * qJD(5); 0; m(7) * t9; t19 + t169; m(7) * t134; m(7) * t124 + t51; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
