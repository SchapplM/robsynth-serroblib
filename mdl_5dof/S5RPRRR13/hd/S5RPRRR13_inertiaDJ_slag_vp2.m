% Calculate time derivative of joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:38
% DurationCPUTime: 1.89s
% Computational Cost: add. (1849->299), mult. (4225->466), div. (0->0), fcn. (3349->6), ass. (0->135)
t145 = 2 * qJD(2);
t96 = cos(qJ(3));
t139 = pkin(7) * t96;
t93 = sin(qJ(3));
t140 = pkin(3) * t93;
t76 = qJ(2) - t139 + t140;
t97 = -pkin(1) - pkin(6);
t131 = t93 * t97;
t95 = cos(qJ(4));
t84 = t95 * t131;
t92 = sin(qJ(4));
t51 = t92 * t76 + t84;
t151 = qJD(4) * t51;
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t103 = t91 * t92 - t94 * t95;
t148 = qJD(4) + qJD(5);
t150 = t148 * t103;
t127 = qJD(3) * t93;
t114 = t92 * t127;
t126 = qJD(3) * t96;
t149 = Ifges(5,6) * t114 + Ifges(5,3) * t126;
t128 = t92 ^ 2 + t95 ^ 2;
t147 = 2 * m(6);
t146 = -0.2e1 * t97;
t144 = -t103 / 0.2e1;
t68 = t91 * t95 + t92 * t94;
t143 = t68 / 0.2e1;
t142 = -t92 / 0.2e1;
t141 = -pkin(8) - pkin(7);
t138 = pkin(8) * t96;
t137 = Ifges(5,4) * t92;
t136 = Ifges(5,4) * t95;
t135 = Ifges(5,5) * t92;
t134 = Ifges(5,6) * t92;
t133 = Ifges(5,6) * t95;
t132 = t93 * Ifges(5,6);
t130 = t96 * mrSges(5,3);
t78 = -mrSges(5,1) * t95 + mrSges(5,2) * t92;
t129 = t78 - mrSges(4,1);
t125 = qJD(3) * t97;
t124 = qJD(4) * t92;
t123 = qJD(4) * t95;
t122 = qJD(4) * t96;
t121 = qJD(5) * t91;
t120 = qJD(5) * t94;
t38 = t148 * t68;
t20 = t103 * t127 - t38 * t96;
t22 = t68 * t127 + t150 * t96;
t119 = Ifges(6,5) * t20 + Ifges(6,6) * t22 + Ifges(6,3) * t126;
t118 = pkin(4) * t124;
t117 = t96 * t125;
t116 = t92 * t122;
t115 = t95 * t122;
t113 = t93 * t124;
t55 = t68 * t96;
t57 = t103 * t96;
t32 = mrSges(6,1) * t55 - mrSges(6,2) * t57;
t66 = (pkin(4) * t92 - t97) * t96;
t112 = m(6) * t66 + t32;
t111 = -t92 * t97 + pkin(4);
t110 = qJD(4) * t141;
t19 = -qJD(3) * t57 - t38 * t93;
t21 = -qJD(3) * t55 + t150 * t93;
t109 = t21 * mrSges(6,1) - t19 * mrSges(6,2);
t108 = -Ifges(5,5) * t95 + (2 * Ifges(4,4));
t65 = qJD(2) + (pkin(3) * t96 + pkin(7) * t93) * qJD(3);
t26 = -t97 * t113 + t95 * t117 + t76 * t123 + t92 * t65;
t64 = t95 * t76;
t50 = -t92 * t131 + t64;
t107 = -qJD(4) * t50 + t26;
t106 = mrSges(5,1) * t92 + mrSges(5,2) * t95;
t105 = Ifges(5,1) * t95 - t137;
t80 = Ifges(5,1) * t92 + t136;
t104 = -Ifges(5,2) * t92 + t136;
t79 = Ifges(5,2) * t95 + t137;
t34 = t111 * t93 - t95 * t138 + t64;
t39 = -t92 * t138 + t51;
t8 = t34 * t94 - t39 * t91;
t9 = t34 * t91 + t39 * t94;
t81 = t141 * t92;
t82 = t141 * t95;
t43 = t81 * t94 + t82 * t91;
t44 = t81 * t91 - t82 * t94;
t74 = t92 * t110;
t75 = t95 * t110;
t24 = t43 * qJD(5) + t74 * t94 + t75 * t91;
t25 = -t44 * qJD(5) - t74 * t91 + t75 * t94;
t35 = Ifges(6,6) * t38;
t36 = Ifges(6,5) * t150;
t102 = t25 * mrSges(6,1) - t24 * mrSges(6,2) - t35 - t36;
t99 = t114 - t115;
t12 = t99 * pkin(8) + t26;
t59 = t95 * t65;
t7 = t59 + (-t84 + (-t76 + t138) * t92) * qJD(4) + (pkin(8) * t93 * t95 + t111 * t96) * qJD(3);
t2 = t8 * qJD(5) + t12 * t94 + t7 * t91;
t3 = -t9 * qJD(5) - t12 * t91 + t7 * t94;
t101 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t119;
t100 = t95 * t127 + t116;
t88 = Ifges(5,5) * t123;
t85 = -pkin(4) * t95 - pkin(3);
t73 = mrSges(5,1) * t93 - t95 * t130;
t72 = -mrSges(5,2) * t93 - t92 * t130;
t71 = t105 * qJD(4);
t70 = t104 * qJD(4);
t69 = t106 * qJD(4);
t62 = (-mrSges(6,1) * t91 - mrSges(6,2) * t94) * qJD(5) * pkin(4);
t61 = t106 * t96;
t56 = t103 * t93;
t54 = t68 * t93;
t53 = Ifges(5,5) * t93 + t105 * t96;
t52 = t104 * t96 + t132;
t49 = -mrSges(5,2) * t126 + t99 * mrSges(5,3);
t48 = mrSges(5,1) * t126 + t100 * mrSges(5,3);
t47 = -t99 * pkin(4) + t93 * t125;
t46 = mrSges(6,1) * t93 + mrSges(6,3) * t57;
t45 = -mrSges(6,2) * t93 - mrSges(6,3) * t55;
t42 = Ifges(6,1) * t68 - Ifges(6,4) * t103;
t41 = Ifges(6,4) * t68 - Ifges(6,2) * t103;
t40 = mrSges(6,1) * t103 + mrSges(6,2) * t68;
t33 = -t99 * mrSges(5,1) - t100 * mrSges(5,2);
t31 = -t80 * t122 + (Ifges(5,5) * t96 - t105 * t93) * qJD(3);
t30 = -t79 * t122 + (Ifges(5,6) * t96 - t104 * t93) * qJD(3);
t29 = -Ifges(6,1) * t57 - Ifges(6,4) * t55 + Ifges(6,5) * t93;
t28 = -Ifges(6,4) * t57 - Ifges(6,2) * t55 + Ifges(6,6) * t93;
t27 = -t92 * t117 - t151 + t59;
t15 = -Ifges(6,1) * t150 - Ifges(6,4) * t38;
t14 = -Ifges(6,4) * t150 - Ifges(6,2) * t38;
t13 = mrSges(6,1) * t38 - mrSges(6,2) * t150;
t11 = -mrSges(6,2) * t126 + mrSges(6,3) * t22;
t10 = mrSges(6,1) * t126 - mrSges(6,3) * t20;
t6 = -mrSges(6,1) * t22 + mrSges(6,2) * t20;
t5 = Ifges(6,1) * t20 + Ifges(6,4) * t22 + Ifges(6,5) * t126;
t4 = Ifges(6,4) * t20 + Ifges(6,2) * t22 + Ifges(6,6) * t126;
t1 = [0.2e1 * t8 * t10 + 0.2e1 * t9 * t11 + 0.2e1 * t2 * t45 + t20 * t29 + t22 * t28 + 0.2e1 * t26 * t72 + 0.2e1 * t27 * t73 + 0.2e1 * t3 * t46 + 0.2e1 * t47 * t32 - t55 * t4 + 0.2e1 * t50 * t48 + 0.2e1 * t51 * t49 - t57 * t5 + 0.2e1 * t66 * t6 + (t2 * t9 + t3 * t8 + t47 * t66) * t147 + 0.2e1 * m(5) * (t51 * t26 + t50 * t27) + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t145 + (mrSges(4,1) * t145 + (-0.2e1 * qJ(2) * mrSges(4,2) + t108 * t93 + t92 * t52 - t95 * t53 + 0.2e1 * t97 * t61) * qJD(3) + t119 + t149) * t93 + (mrSges(4,2) * t145 - t92 * t30 + t95 * t31 + t33 * t146 + (-t95 * t52 - t92 * t53 + t93 * (-t133 - t135)) * qJD(4) + (-Ifges(6,5) * t57 - Ifges(6,6) * t55 + 0.2e1 * qJ(2) * mrSges(4,1) + (-t108 - t134) * t96 + (-0.2e1 * m(5) * t97 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t93) * qJD(3)) * t96; m(6) * (t19 * t9 - t2 * t56 + t21 * t8 - t3 * t54 - t47 * t96) + t19 * t45 - t56 * t11 + t21 * t46 - t54 * t10 - t96 * t6 - t96 * t33 + (t95 * t72 - t92 * t73 + m(5) * (-t50 * t92 + t51 * t95)) * t126 + (-t72 * t124 + t95 * t49 - t73 * t123 - t92 * t48 + m(5) * (-t50 * t123 - t51 * t124 + t26 * t95 - t27 * t92) + (m(5) * t96 * t146 + t112 + t61) * qJD(3)) * t93; (-t19 * t56 - t21 * t54) * t147 + 0.4e1 * (m(5) * (-0.1e1 + t128) / 0.2e1 - m(6) / 0.2e1) * t93 * t126; m(6) * (t2 * t44 + t24 * t9 + t25 * t8 + t3 * t43 + t47 * t85) + t85 * t6 + t66 * t13 + t4 * t144 + t5 * t143 - t55 * t14 / 0.2e1 - t57 * t15 / 0.2e1 - t150 * t29 / 0.2e1 - t38 * t28 / 0.2e1 + t22 * t41 / 0.2e1 + t20 * t42 / 0.2e1 + t43 * t10 + t44 * t11 + t24 * t45 + t25 * t46 + t47 * t40 - pkin(3) * t33 + (t88 / 0.2e1 - t36 / 0.2e1 - t35 / 0.2e1 + (-Ifges(4,5) + (-m(5) * pkin(3) + t129) * t97) * qJD(3)) * t93 + (t79 * t127 / 0.2e1 - t27 * mrSges(5,3) + t31 / 0.2e1 + (-t52 / 0.2e1 - t132 / 0.2e1 - t51 * mrSges(5,3) + t112 * pkin(4)) * qJD(4) + (-qJD(4) * t72 - t48 + m(5) * (-t27 - t151)) * pkin(7)) * t92 + (-t103 * t2 + t150 * t8 - t3 * t68 - t38 * t9) * mrSges(6,3) + (-t80 * t127 / 0.2e1 + qJD(4) * t53 / 0.2e1 + t30 / 0.2e1 + t107 * mrSges(5,3) + (m(5) * t107 - qJD(4) * t73 + t49) * pkin(7)) * t95 + (t95 * t71 / 0.2e1 + t70 * t142 - t97 * t69 + (-t95 * t79 / 0.2e1 + t80 * t142) * qJD(4) + (-t97 * mrSges(4,2) + t135 / 0.2e1 + t133 / 0.2e1 + Ifges(6,5) * t143 + Ifges(6,6) * t144 - Ifges(4,6)) * qJD(3)) * t96; m(6) * (-pkin(4) * t116 + t19 * t44 + t21 * t43 - t24 * t56 - t25 * t54) + (-t103 * t19 - t150 * t54 - t21 * t68 + t38 * t56) * mrSges(6,3) + (m(5) * (t128 * t139 - t140) + (m(6) * t85 + t129 + t40) * t93) * qJD(3) + (-t69 - t13 + (t128 * mrSges(5,3) - mrSges(4,2)) * qJD(3)) * t96; 0.2e1 * t40 * t118 + 0.2e1 * t85 * t13 - t150 * t42 + t68 * t15 - t38 * t41 - t103 * t14 + (t85 * t118 + t24 * t44 + t25 * t43) * t147 - 0.2e1 * pkin(3) * t69 + t92 * t71 - t79 * t124 + (qJD(4) * t80 + t70) * t95 + 0.2e1 * (-t103 * t24 + t150 * t43 - t25 * t68 - t38 * t44) * mrSges(6,3); -Ifges(5,6) * t115 + t27 * mrSges(5,1) - t26 * mrSges(5,2) - t100 * Ifges(5,5) + (m(6) * (t9 * t120 - t8 * t121 + t2 * t91 + t3 * t94) + t45 * t120 + t91 * t11 - t46 * t121 + t94 * t10) * pkin(4) + t101 + t149; (-t95 * t126 + t113) * mrSges(5,2) + (-t93 * t123 - t92 * t126) * mrSges(5,1) + m(6) * (t19 * t91 + t21 * t94 + (t54 * t91 - t56 * t94) * qJD(5)) * pkin(4) + t109; t88 + (t78 * pkin(7) - t134) * qJD(4) + (m(6) * (t24 * t91 + t25 * t94 + (-t43 * t91 + t44 * t94) * qJD(5)) + (t94 * t150 - t91 * t38 + (-t103 * t94 + t68 * t91) * qJD(5)) * mrSges(6,3)) * pkin(4) + t102; 0.2e1 * t62; t101; t109; t102; t62; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
