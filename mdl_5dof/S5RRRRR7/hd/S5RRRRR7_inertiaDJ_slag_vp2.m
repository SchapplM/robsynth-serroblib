% Calculate time derivative of joint inertia matrix for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:24
% EndTime: 2019-12-31 22:21:30
% DurationCPUTime: 2.32s
% Computational Cost: add. (5227->282), mult. (11563->437), div. (0->0), fcn. (10985->8), ass. (0->141)
t106 = sin(qJ(5));
t104 = t106 ^ 2;
t110 = cos(qJ(5));
t105 = t110 ^ 2;
t188 = t104 + t105;
t139 = qJD(5) * t110;
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t108 = sin(qJ(3));
t109 = sin(qJ(2));
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t146 = t112 * t113;
t81 = -t108 * t109 + t146;
t82 = t108 * t113 + t112 * t109;
t58 = t107 * t82 - t111 * t81;
t185 = qJD(2) + qJD(3);
t63 = t185 * t81;
t64 = t185 * t82;
t29 = -qJD(4) * t58 - t107 * t64 + t111 * t63;
t159 = t106 * t29;
t59 = t107 * t81 + t111 * t82;
t122 = t59 * t139 + t159;
t101 = -pkin(2) * t113 - pkin(1);
t67 = -pkin(3) * t81 + t101;
t34 = pkin(4) * t58 - pkin(9) * t59 + t67;
t176 = -pkin(7) - pkin(6);
t94 = t176 * t109;
t95 = t176 * t113;
t65 = t108 * t95 + t112 * t94;
t123 = -pkin(8) * t82 + t65;
t83 = t108 * t94;
t66 = -t112 * t95 + t83;
t52 = pkin(8) * t81 + t66;
t36 = t107 * t123 + t111 * t52;
t17 = t106 * t34 + t110 * t36;
t145 = t17 * qJD(5);
t118 = (t146 * t176 - t83) * qJD(2);
t143 = qJD(3) * t112;
t144 = qJD(3) * t108;
t174 = t63 * pkin(8);
t41 = qJD(2) * t176 * t82 + t94 * t143 + t144 * t95;
t33 = -pkin(8) * t64 + t41;
t35 = t107 * t52 - t111 * t123;
t11 = -t35 * qJD(4) + t111 * t33 + (t143 * t95 - t144 * t94 + t118 - t174) * t107;
t30 = qJD(4) * t59 + t107 * t63 + t111 * t64;
t53 = qJD(2) * t109 * pkin(2) + pkin(3) * t64;
t15 = pkin(4) * t30 - pkin(9) * t29 + t53;
t3 = -t106 * t11 + t110 * t15 - t145;
t187 = -t3 - t145;
t186 = t188 * t111;
t127 = mrSges(6,1) * t106 + mrSges(6,2) * t110;
t86 = t127 * qJD(5);
t99 = -pkin(3) * t111 - pkin(4);
t68 = t99 * t86;
t142 = qJD(4) * t107;
t166 = mrSges(6,1) * t110;
t89 = mrSges(6,2) * t106 - t166;
t78 = pkin(3) * t89 * t142;
t141 = qJD(4) * t111;
t137 = pkin(3) * t141;
t129 = mrSges(6,3) * t137;
t92 = t104 * t129;
t93 = t105 * t129;
t184 = t68 + t78 + t92 + t93;
t42 = -qJD(3) * t66 + t118;
t183 = t42 * mrSges(4,1) - t41 * mrSges(4,2) + Ifges(4,5) * t63 - Ifges(4,6) * t64;
t182 = 2 * m(5);
t181 = 2 * m(6);
t12 = t107 * t33 - t111 * (t42 - t174) + t36 * qJD(4);
t180 = 0.2e1 * t12;
t179 = 0.2e1 * t53;
t178 = 0.2e1 * t101;
t175 = pkin(4) * t86;
t16 = -t106 * t36 + t110 * t34;
t154 = qJD(5) * t16;
t2 = t106 * t15 + t11 * t110 + t154;
t173 = t110 * t2;
t172 = t12 * t35;
t171 = t3 * t106;
t100 = pkin(2) * t112 + pkin(3);
t149 = t108 * t111;
t56 = t100 * t142 + (t108 * t141 + (t107 * t112 + t149) * qJD(3)) * pkin(2);
t170 = t35 * t56;
t150 = t107 * t108;
t55 = t100 * t141 + (-t108 * t142 + (t111 * t112 - t150) * qJD(3)) * pkin(2);
t168 = t55 * mrSges(5,2);
t156 = t110 * t29;
t167 = Ifges(6,5) * t156 + Ifges(6,3) * t30;
t72 = pkin(2) * t149 + t100 * t107;
t165 = Ifges(6,4) * t106;
t164 = Ifges(6,4) * t110;
t163 = Ifges(6,6) * t106;
t162 = pkin(3) * qJD(4);
t161 = t104 * t55;
t160 = t105 * t55;
t158 = t106 * t59;
t155 = t110 * t59;
t70 = pkin(9) + t72;
t153 = qJD(5) * t70;
t140 = qJD(5) * t106;
t138 = 0.2e1 * t113;
t135 = t59 * t140;
t121 = t135 - t156;
t8 = mrSges(6,1) * t122 - mrSges(6,2) * t121;
t136 = m(6) * t12 + t8;
t131 = -(2 * Ifges(5,4)) - t163;
t130 = t188 * t55;
t128 = -t107 * mrSges(5,1) - t111 * mrSges(5,2);
t126 = Ifges(6,1) * t110 - t165;
t125 = -Ifges(6,2) * t106 + t164;
t124 = Ifges(6,5) * t106 + Ifges(6,6) * t110;
t71 = -pkin(2) * t150 + t100 * t111;
t87 = t125 * qJD(5);
t88 = t126 * qJD(5);
t90 = Ifges(6,2) * t110 + t165;
t91 = Ifges(6,1) * t106 + t164;
t120 = t106 * t88 + t110 * t87 + t139 * t91 - t140 * t90;
t119 = (-mrSges(4,1) * t108 - mrSges(4,2) * t112) * qJD(3) * pkin(2);
t48 = t56 * t89;
t50 = mrSges(6,3) * t161;
t51 = mrSges(6,3) * t160;
t54 = t56 * mrSges(5,1);
t69 = -pkin(4) - t71;
t62 = t69 * t86;
t117 = t120 + t48 + t50 + t51 - t54 + t62 - t168;
t102 = Ifges(6,5) * t139;
t23 = Ifges(6,6) * t58 + t125 * t59;
t24 = Ifges(6,5) * t58 + t126 * t59;
t6 = -Ifges(6,4) * t121 - Ifges(6,2) * t122 + Ifges(6,6) * t30;
t7 = -Ifges(6,1) * t121 - Ifges(6,4) * t122 + Ifges(6,5) * t30;
t116 = -t11 * mrSges(5,2) + mrSges(6,3) * t173 + t106 * t7 / 0.2e1 + t91 * t156 / 0.2e1 + t24 * t139 / 0.2e1 + Ifges(5,5) * t29 + t35 * t86 - t87 * t158 / 0.2e1 + t88 * t155 / 0.2e1 + t58 * (-Ifges(6,6) * t140 + t102) / 0.2e1 + t110 * t6 / 0.2e1 + (t124 / 0.2e1 - Ifges(5,6)) * t30 - t122 * t90 / 0.2e1 - (t59 * t91 + t23) * t140 / 0.2e1 + (-mrSges(5,1) + t89) * t12;
t13 = mrSges(6,1) * t30 + mrSges(6,3) * t121;
t14 = -mrSges(6,2) * t30 - mrSges(6,3) * t122;
t38 = -mrSges(6,2) * t58 - mrSges(6,3) * t158;
t39 = mrSges(6,1) * t58 - mrSges(6,3) * t155;
t115 = -t39 * t139 - t38 * t140 - t106 * t13 + t110 * t14 + m(6) * (-t139 * t16 - t140 * t17 - t171 + t173);
t114 = t116 + (-t171 + (-t106 * t17 - t110 * t16) * qJD(5)) * mrSges(6,3);
t98 = pkin(3) * t107 + pkin(9);
t37 = t127 * t59;
t1 = [(mrSges(4,1) * t64 + mrSges(4,2) * t63) * t178 - 0.2e1 * t81 * Ifges(4,2) * t64 + 0.2e1 * t63 * t82 * Ifges(4,1) + 0.2e1 * t67 * (mrSges(5,1) * t30 + mrSges(5,2) * t29) + 0.2e1 * t35 * t8 + t37 * t180 + 0.2e1 * t2 * t38 + 0.2e1 * t3 * t39 + 0.2e1 * t16 * t13 + 0.2e1 * t17 * t14 - t23 * t159 + t24 * t156 + (t16 * t3 + t17 * t2 + t172) * t181 + (t11 * t36 + t53 * t67 + t172) * t182 + 0.2e1 * m(4) * (t41 * t66 + t42 * t65) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t113) * t138 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t81 + mrSges(4,2) * t82) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t178 - 0.2e1 * Ifges(3,4) * t109 + (Ifges(3,1) - Ifges(3,2)) * t138) * t109) * qJD(2) + (mrSges(5,1) * t179 - 0.2e1 * t11 * mrSges(5,3) + ((2 * Ifges(5,2)) + Ifges(6,3)) * t30 + t131 * t29 + t167) * t58 + (mrSges(5,2) * t179 + mrSges(5,3) * t180 + 0.2e1 * Ifges(5,1) * t29 - t106 * t6 + t110 * t7 + (Ifges(6,5) * t110 + t131) * t30 + (-t106 * t24 - t110 * t23 - t124 * t58) * qJD(5)) * t59 + 0.2e1 * (t29 * t35 - t30 * t36) * mrSges(5,3) + 0.2e1 * (t63 * t81 - t64 * t82) * Ifges(4,4) + 0.2e1 * (t41 * t81 - t42 * t82 - t63 * t65 - t64 * t66) * mrSges(4,3); (m(4) * (t108 * t41 + t112 * t42 + (-t108 * t65 + t112 * t66) * qJD(3)) + (-t108 * t64 - t112 * t63 + (t108 * t82 + t112 * t81) * qJD(3)) * mrSges(4,3)) * pkin(2) + (-t38 * t153 + t187 * mrSges(6,3) + (-m(6) * t16 - t39) * t55 + (m(6) * t187 - t13) * t70) * t106 + (-mrSges(6,3) * t154 - t39 * t153 + t55 * t38 + t70 * t14 + m(6) * (-t153 * t16 + t17 * t55 + t2 * t70)) * t110 + (-t29 * t71 - t30 * t72 - t55 * t58 + t56 * t59) * mrSges(5,3) + t116 + (Ifges(3,5) * t113 - Ifges(3,6) * t109 + (-mrSges(3,1) * t113 + mrSges(3,2) * t109) * pkin(6)) * qJD(2) + m(6) * (t12 * t69 + t170) + t69 * t8 + t56 * t37 + m(5) * (t11 * t72 - t12 * t71 + t36 * t55 + t170) + t183; -0.2e1 * t168 + 0.2e1 * t48 + 0.2e1 * t50 + 0.2e1 * t51 - 0.2e1 * t54 + 0.2e1 * t62 + 0.2e1 * t119 + (t130 * t70 + t56 * t69) * t181 + (t55 * t72 - t56 * t71) * t182 + t120; t115 * t98 + t136 * t99 + (m(5) * (t107 * t11 - t111 * t12) + (-t107 * t30 - t111 * t29) * mrSges(5,3) + ((-t58 * mrSges(5,3) - t106 * t39 + t110 * t38 + m(5) * t36 + m(6) * (-t106 * t16 + t110 * t17)) * t111 + (t59 * mrSges(5,3) + t37 + (m(5) + m(6)) * t35) * t107) * qJD(4)) * pkin(3) + t114 + t183; m(6) * (t56 * t99 + (t160 + t161) * t98) + (m(5) * (t107 * t55 - t111 * t56) + (m(6) * (t107 * t69 + t186 * t70) + m(5) * (-t107 * t71 + t111 * t72) + t128) * qJD(4)) * pkin(3) + t117 + t119 + t184; 0.2e1 * t68 + 0.2e1 * t78 + 0.2e1 * t92 + 0.2e1 * t93 + 0.2e1 * (m(6) * (t107 * t99 + t186 * t98) + t128) * t162 + t120; -pkin(4) * t136 + pkin(9) * t115 + t114; m(6) * (-pkin(4) * t56 + pkin(9) * t130) - t175 + t117; -t175 + (m(6) * (-pkin(4) * t107 + pkin(9) * t186) + t128) * t162 + t120 + t184; t120 - 0.2e1 * t175; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t135 - Ifges(6,6) * t122 + t167; t102 - t127 * t55 + (-t70 * t166 + (mrSges(6,2) * t70 - Ifges(6,6)) * t106) * qJD(5); t102 - t127 * t137 + (-t98 * t166 + (mrSges(6,2) * t98 - Ifges(6,6)) * t106) * qJD(5); t102 + (pkin(9) * t89 - t163) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
