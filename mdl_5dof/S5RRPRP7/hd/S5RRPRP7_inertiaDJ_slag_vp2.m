% Calculate time derivative of joint inertia matrix for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:56
% EndTime: 2019-12-31 20:00:01
% DurationCPUTime: 1.74s
% Computational Cost: add. (1624->254), mult. (3719->363), div. (0->0), fcn. (3192->6), ass. (0->111)
t119 = Ifges(6,4) + Ifges(5,5);
t140 = Ifges(6,2) + Ifges(5,3);
t83 = cos(qJ(4));
t109 = qJD(4) * t83;
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t55 = t79 * t82 - t80 * t84;
t53 = t55 * qJD(2);
t81 = sin(qJ(4));
t122 = t81 * t53;
t56 = t79 * t84 + t80 * t82;
t87 = t56 * t109 - t122;
t110 = qJD(4) * t81;
t120 = t83 * t53;
t86 = t56 * t110 + t120;
t74 = -pkin(2) * t84 - pkin(1);
t33 = t55 * pkin(3) - t56 * pkin(7) + t74;
t118 = -qJ(3) - pkin(6);
t65 = t118 * t82;
t66 = t118 * t84;
t40 = t65 * t79 - t66 * t80;
t138 = t81 * t33 + t83 * t40;
t139 = qJD(4) * t138;
t137 = Ifges(6,6) * t110 + t119 * t109;
t99 = qJD(2) * t118;
t51 = qJD(3) * t84 + t82 * t99;
t85 = -t82 * qJD(3) + t84 * t99;
t28 = t80 * t51 + t79 * t85;
t106 = pkin(2) * qJD(2) * t82;
t52 = t56 * qJD(2);
t29 = t52 * pkin(3) + t53 * pkin(7) + t106;
t4 = -t28 * t81 + t29 * t83 - t139;
t108 = qJD(5) * t83;
t92 = pkin(4) * t83 + qJ(5) * t81;
t136 = t92 * qJD(4) - t108;
t135 = -2 * mrSges(4,3);
t134 = -2 * Ifges(4,4);
t39 = -t80 * t65 - t66 * t79;
t132 = 0.2e1 * t39;
t130 = Ifges(5,4) * t81;
t129 = Ifges(5,4) * t83;
t128 = Ifges(6,5) * t81;
t127 = Ifges(6,5) * t83;
t126 = Ifges(5,6) * t55;
t27 = t51 * t79 - t80 * t85;
t125 = t39 * t27;
t124 = t56 * t81;
t123 = t56 * t83;
t17 = t52 * mrSges(5,1) + t86 * mrSges(5,3);
t18 = -t52 * mrSges(6,1) - t86 * mrSges(6,2);
t117 = t17 - t18;
t19 = -t52 * mrSges(5,2) - t87 * mrSges(5,3);
t20 = -t87 * mrSges(6,2) + t52 * mrSges(6,3);
t116 = t19 + t20;
t93 = Ifges(6,3) * t81 + t127;
t21 = Ifges(6,6) * t55 + t93 * t56;
t94 = -Ifges(5,2) * t81 + t129;
t22 = t94 * t56 + t126;
t115 = t21 - t22;
t95 = Ifges(6,1) * t83 + t128;
t23 = Ifges(6,4) * t55 + t95 * t56;
t96 = Ifges(5,1) * t83 - t130;
t24 = Ifges(5,5) * t55 + t96 * t56;
t114 = t23 + t24;
t34 = -t55 * mrSges(5,2) - mrSges(5,3) * t124;
t37 = -mrSges(6,2) * t124 + t55 * mrSges(6,3);
t113 = t34 + t37;
t35 = t55 * mrSges(5,1) - mrSges(5,3) * t123;
t36 = -t55 * mrSges(6,1) + mrSges(6,2) * t123;
t112 = -t35 + t36;
t67 = -Ifges(6,3) * t83 + t128;
t68 = Ifges(5,2) * t83 + t130;
t103 = t67 / 0.2e1 - t68 / 0.2e1;
t69 = Ifges(6,1) * t81 - t127;
t70 = Ifges(5,1) * t81 + t129;
t102 = t69 / 0.2e1 + t70 / 0.2e1;
t73 = -pkin(2) * t80 - pkin(3);
t101 = t52 * mrSges(4,1) - t53 * mrSges(4,2);
t100 = 0.2e1 * t106;
t98 = t81 * mrSges(5,1) + t83 * mrSges(5,2);
t64 = -t83 * mrSges(6,1) - t81 * mrSges(6,3);
t97 = t81 * mrSges(6,1) - t83 * mrSges(6,3);
t91 = pkin(4) * t81 - qJ(5) * t83;
t14 = t33 * t83 - t40 * t81;
t88 = t87 * Ifges(6,6) - t119 * t120 + t140 * t52;
t3 = t33 * t109 - t40 * t110 + t83 * t28 + t81 * t29;
t72 = pkin(2) * t79 + pkin(7);
t63 = t96 * qJD(4);
t62 = t95 * qJD(4);
t61 = t94 * qJD(4);
t60 = t93 * qJD(4);
t59 = t98 * qJD(4);
t58 = t97 * qJD(4);
t54 = t73 - t92;
t50 = -pkin(4) * t110 + qJ(5) * t109 + qJD(5) * t81;
t32 = t97 * t56;
t16 = t91 * t56 + t39;
t13 = t87 * mrSges(5,1) - t86 * mrSges(5,2);
t12 = t87 * mrSges(6,1) + t86 * mrSges(6,3);
t11 = -t55 * pkin(4) - t14;
t10 = qJ(5) * t55 + t138;
t9 = -t86 * Ifges(5,1) - t87 * Ifges(5,4) + Ifges(5,5) * t52;
t8 = -t86 * Ifges(6,1) + Ifges(6,4) * t52 + t87 * Ifges(6,5);
t7 = -t86 * Ifges(5,4) - t87 * Ifges(5,2) + Ifges(5,6) * t52;
t6 = -t86 * Ifges(6,5) + Ifges(6,6) * t52 + t87 * Ifges(6,3);
t5 = t136 * t56 - t91 * t53 + t27;
t2 = -t52 * pkin(4) - t4;
t1 = t52 * qJ(5) + t55 * qJD(5) + t3;
t15 = [0.2e1 * t74 * t101 + t13 * t132 + 0.2e1 * t5 * t32 + 0.2e1 * t3 * t34 + 0.2e1 * t4 * t35 + 0.2e1 * t2 * t36 + 0.2e1 * t1 * t37 + 0.2e1 * t16 * t12 + 0.2e1 * t14 * t17 + 0.2e1 * t11 * t18 + 0.2e1 * t138 * t19 + 0.2e1 * t10 * t20 + t40 * t52 * t135 + 0.2e1 * m(4) * (t74 * t106 + t40 * t28 + t125) + 0.2e1 * m(5) * (t138 * t3 + t14 * t4 + t125) + 0.2e1 * m(6) * (t1 * t10 + t11 * t2 + t16 * t5) - (mrSges(4,3) * t132 + t114 * t83 + t115 * t81) * t53 + (mrSges(4,1) * t100 + t28 * t135 - (-Ifges(5,6) * t81 + t134) * t53 + ((2 * Ifges(4,2)) + t140) * t52 + t88) * t55 + (mrSges(4,2) * t100 - 0.2e1 * Ifges(4,1) * t53 + t134 * t52 + (t8 + t9 + t119 * t52 + (t115 - t126) * qJD(4)) * t83 + (t6 - t7 + (-Ifges(5,6) + Ifges(6,6)) * t52 + (-t119 * t55 - t114) * qJD(4)) * t81 + 0.2e1 * (t98 + mrSges(4,3)) * t27) * t56 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t82 + mrSges(3,2) * t84) + (-Ifges(3,2) + Ifges(3,1)) * t82 * t84 + (-t82 ^ 2 + t84 ^ 2) * Ifges(3,4)) * qJD(2); t5 * t64 + t16 * t58 + t39 * t59 - t50 * t32 - Ifges(4,6) * t52 - Ifges(4,5) * t53 + t54 * t12 - t27 * mrSges(4,1) - t28 * mrSges(4,2) + m(6) * (-t50 * t16 + t54 * t5) + (Ifges(3,5) * t84 - Ifges(3,6) * t82 + (-mrSges(3,1) * t84 + mrSges(3,2) * t82) * pkin(6)) * qJD(2) + (m(4) * (-t27 * t80 + t28 * t79) + (-t52 * t79 + t53 * t80) * mrSges(4,3)) * pkin(2) + (t27 * mrSges(5,2) + t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t4 * mrSges(5,3) + (t60 / 0.2e1 - t61 / 0.2e1) * t56 - t103 * t53 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t52 + (-t126 / 0.2e1 + t21 / 0.2e1 - t22 / 0.2e1 - t10 * mrSges(6,2) - t138 * mrSges(5,3) - t102 * t56) * qJD(4) + (-t113 * qJD(4) + m(6) * (-qJD(4) * t10 + t2) + m(5) * (-t4 - t139) - t117) * t72) * t81 + (-t27 * mrSges(5,1) - t6 / 0.2e1 + t7 / 0.2e1 + t1 * mrSges(6,2) + t3 * mrSges(5,3) + (t62 / 0.2e1 + t63 / 0.2e1) * t56 - t102 * t53 + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t52 + (t23 / 0.2e1 + t24 / 0.2e1 + t11 * mrSges(6,2) - t14 * mrSges(5,3) + t103 * t56) * qJD(4) + (t112 * qJD(4) + m(6) * (qJD(4) * t11 + t1) + m(5) * (-qJD(4) * t14 + t3) + t116) * t72) * t83 + (m(5) * t27 + t13) * t73 + t137 * t55 / 0.2e1; 0.2e1 * t54 * t58 + 0.2e1 * t59 * t73 + (-t60 + t61) * t83 + (t62 + t63) * t81 + 0.2e1 * (-m(6) * t54 - t64) * t50 + ((t69 + t70) * t83 + (t67 - t68) * t81) * qJD(4); m(4) * t106 + t117 * t83 + t116 * t81 + (t112 * t81 + t113 * t83) * qJD(4) + m(6) * (t1 * t81 - t2 * t83 + (t10 * t83 + t11 * t81) * qJD(4)) + m(5) * (t3 * t81 + t4 * t83 + (t138 * t83 - t14 * t81) * qJD(4)) + t101; 0; 0; Ifges(5,6) * t122 - t3 * mrSges(5,2) - pkin(4) * t18 + qJD(5) * t37 + qJ(5) * t20 + m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t10) + t1 * mrSges(6,3) - t2 * mrSges(6,1) + t4 * mrSges(5,1) + (-Ifges(5,6) * t83 - t119 * t81) * t56 * qJD(4) + t88; -Ifges(5,6) * t110 - t136 * mrSges(6,2) + (m(6) * t108 + (-m(6) * t92 - t83 * mrSges(5,1) + t81 * mrSges(5,2) + t64) * qJD(4)) * t72 + t137; m(6) * t50 + ((-mrSges(5,2) + mrSges(6,3)) * t83 + (-mrSges(5,1) - mrSges(6,1)) * t81) * qJD(4); 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t2 + t18; (m(6) * t72 + mrSges(6,2)) * t109; m(6) * t110; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
