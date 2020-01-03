% Calculate time derivative of joint inertia matrix for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:13
% EndTime: 2019-12-31 19:03:18
% DurationCPUTime: 1.81s
% Computational Cost: add. (1915->291), mult. (4506->447), div. (0->0), fcn. (3636->8), ass. (0->130)
t106 = -cos(pkin(9)) * pkin(1) - pkin(2);
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t61 = -pkin(3) * t94 - t91 * pkin(7) + t106;
t93 = cos(qJ(4));
t122 = t93 * t94;
t81 = sin(pkin(9)) * pkin(1) + pkin(6);
t71 = t81 * t122;
t90 = sin(qJ(4));
t37 = t90 * t61 + t71;
t146 = qJD(4) * t37;
t145 = 0.2e1 * t106;
t114 = qJD(4) * t93;
t117 = qJD(3) * t94;
t95 = t91 * t114 + t90 * t117;
t115 = qJD(4) * t91;
t108 = t90 * t115;
t109 = t93 * t117;
t96 = -t108 + t109;
t33 = t95 * mrSges(5,1) + t96 * mrSges(5,2);
t142 = qJD(4) + qJD(5);
t89 = sin(qJ(5));
t92 = cos(qJ(5));
t65 = t89 * t93 + t90 * t92;
t39 = t142 * t65;
t99 = t89 * t90 - t92 * t93;
t21 = -t99 * t117 - t39 * t91;
t54 = t99 * t91;
t22 = -t65 * t117 + t142 * t54;
t6 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t144 = -t33 - t6;
t118 = qJD(3) * t91;
t143 = -Ifges(5,5) * t109 - Ifges(5,3) * t118;
t119 = t90 ^ 2 + t93 ^ 2;
t141 = 2 * m(5);
t140 = 2 * m(6);
t139 = 0.2e1 * t81;
t138 = -t99 / 0.2e1;
t137 = t65 / 0.2e1;
t131 = Ifges(5,4) * t90;
t76 = Ifges(5,2) * t93 + t131;
t136 = -t76 / 0.2e1;
t135 = -t90 / 0.2e1;
t134 = -pkin(8) - pkin(7);
t133 = pkin(3) * t91;
t132 = pkin(7) * t94;
t130 = Ifges(5,4) * t93;
t129 = Ifges(5,5) * t90;
t128 = Ifges(5,6) * t90;
t127 = Ifges(5,6) * t93;
t126 = Ifges(5,6) * t94;
t125 = t81 * t90;
t124 = t90 * t91;
t123 = t91 * t93;
t74 = (-t132 + t133) * qJD(3);
t121 = t118 * t125 + t93 * t74;
t75 = -mrSges(5,1) * t93 + mrSges(5,2) * t90;
t120 = t75 - mrSges(4,1);
t116 = qJD(4) * t90;
t113 = qJD(5) * t89;
t112 = qJD(5) * t92;
t111 = -Ifges(6,5) * t21 - Ifges(6,6) * t22 - Ifges(6,3) * t118;
t110 = pkin(4) * t116;
t107 = t94 * t116;
t105 = qJD(4) * t134;
t104 = (2 * Ifges(4,4)) + t128;
t19 = t61 * t114 + t90 * t74 + (-t93 * t118 - t107) * t81;
t56 = t93 * t61;
t36 = -t94 * t125 + t56;
t103 = -qJD(4) * t36 + t19;
t102 = mrSges(5,1) * t90 + mrSges(5,2) * t93;
t101 = Ifges(5,1) * t93 - t131;
t77 = Ifges(5,1) * t90 + t130;
t100 = -Ifges(5,2) * t90 + t130;
t28 = -pkin(8) * t123 + t56 + (-pkin(4) - t125) * t94;
t32 = -pkin(8) * t124 + t37;
t7 = t28 * t92 - t32 * t89;
t8 = t28 * t89 + t32 * t92;
t78 = t134 * t90;
t79 = t134 * t93;
t44 = t78 * t92 + t79 * t89;
t45 = t78 * t89 - t79 * t92;
t72 = t90 * t105;
t73 = t93 * t105;
t24 = t44 * qJD(5) + t72 * t92 + t73 * t89;
t25 = -t45 * qJD(5) - t72 * t89 + t73 * t92;
t34 = Ifges(6,6) * t39;
t38 = t142 * t99;
t35 = Ifges(6,5) * t38;
t98 = t25 * mrSges(6,1) - t24 * mrSges(6,2) - t34 - t35;
t10 = -pkin(8) * t95 + t19;
t9 = (pkin(4) * t91 - pkin(8) * t122) * qJD(3) + (-t71 + (pkin(8) * t91 - t61) * t90) * qJD(4) + t121;
t2 = t7 * qJD(5) + t10 * t92 + t89 * t9;
t3 = -t8 * qJD(5) - t10 * t89 + t9 * t92;
t97 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t111;
t85 = Ifges(5,5) * t114;
t82 = -pkin(4) * t93 - pkin(3);
t70 = -mrSges(5,1) * t94 - mrSges(5,3) * t123;
t69 = mrSges(5,2) * t94 - mrSges(5,3) * t124;
t68 = t101 * qJD(4);
t67 = t100 * qJD(4);
t66 = t102 * qJD(4);
t62 = (-mrSges(6,1) * t89 - mrSges(6,2) * t92) * qJD(5) * pkin(4);
t60 = t102 * t91;
t57 = (pkin(4) * t90 + t81) * t91;
t53 = t65 * t91;
t52 = -t94 * Ifges(5,5) + t101 * t91;
t51 = t100 * t91 - t126;
t49 = -mrSges(5,2) * t118 - mrSges(5,3) * t95;
t48 = mrSges(5,1) * t118 - mrSges(5,3) * t96;
t47 = -mrSges(6,1) * t94 + t54 * mrSges(6,3);
t46 = mrSges(6,2) * t94 - t53 * mrSges(6,3);
t43 = pkin(4) * t95 + t81 * t117;
t42 = Ifges(6,1) * t65 - Ifges(6,4) * t99;
t41 = Ifges(6,4) * t65 - Ifges(6,2) * t99;
t40 = mrSges(6,1) * t99 + mrSges(6,2) * t65;
t31 = mrSges(6,1) * t53 - mrSges(6,2) * t54;
t30 = -t77 * t115 + (Ifges(5,5) * t91 + t101 * t94) * qJD(3);
t29 = -t76 * t115 + (Ifges(5,6) * t91 + t100 * t94) * qJD(3);
t27 = -Ifges(6,1) * t54 - Ifges(6,4) * t53 - Ifges(6,5) * t94;
t26 = -Ifges(6,4) * t54 - Ifges(6,2) * t53 - Ifges(6,6) * t94;
t20 = t121 - t146;
t15 = -Ifges(6,1) * t38 - Ifges(6,4) * t39;
t14 = -Ifges(6,4) * t38 - Ifges(6,2) * t39;
t13 = mrSges(6,1) * t39 - mrSges(6,2) * t38;
t12 = -mrSges(6,2) * t118 + mrSges(6,3) * t22;
t11 = mrSges(6,1) * t118 - mrSges(6,3) * t21;
t5 = Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t118;
t4 = Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t118;
t1 = [0.2e1 * t7 * t11 + 0.2e1 * t8 * t12 + 0.2e1 * t19 * t69 + 0.2e1 * t2 * t46 + 0.2e1 * t20 * t70 + t21 * t27 + t22 * t26 + 0.2e1 * t3 * t47 + 0.2e1 * t43 * t31 + 0.2e1 * t36 * t48 + 0.2e1 * t37 * t49 - t53 * t4 - t54 * t5 + 0.2e1 * t57 * t6 + (t37 * t19 + t36 * t20) * t141 + (t2 * t8 + t3 * t7 + t43 * t57) * t140 + ((mrSges(4,2) * t145 + t104 * t94 + t60 * t139 - t90 * t51 + t93 * t52) * qJD(3) + t111 + t143) * t94 + (-t90 * t29 + t93 * t30 + t33 * t139 + (-t94 * (-t127 - t129) - t93 * t51 - t90 * t52) * qJD(4) + (mrSges(4,1) * t145 - Ifges(6,5) * t54 - Ifges(6,6) * t53 + (Ifges(5,5) * t93 - t104) * t91 + (t81 ^ 2 * t141 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) - Ifges(6,3)) * t94) * qJD(3)) * t91; -t53 * t11 - t54 * t12 + t21 * t46 + t22 * t47 + t144 * t94 + (t31 + t60) * t118 + m(6) * (t57 * t118 - t2 * t54 + t8 * t21 + t7 * t22 - t3 * t53 - t43 * t94) + m(5) * (t91 ^ 2 - t94 ^ 2) * qJD(3) * t81 + (t69 * t117 + t91 * t49 - t70 * t115 + m(5) * (-t36 * t115 + t37 * t117 + t19 * t91)) * t93 + (-t69 * t115 - t70 * t117 - t91 * t48 + m(5) * (-t37 * t115 - t36 * t117 - t20 * t91)) * t90; (-t54 * t21 - t53 * t22) * t140 + 0.4e1 * (-m(6) / 0.2e1 + m(5) * (-0.1e1 + t119) / 0.2e1) * t91 * t117; t82 * t6 + t57 * t13 + t4 * t138 + t5 * t137 + m(6) * (t2 * t45 + t24 * t8 + t25 * t7 + t3 * t44 + t43 * t82) + t44 * t11 + t45 * t12 + t24 * t46 + t25 * t47 - t53 * t14 / 0.2e1 - t54 * t15 / 0.2e1 - pkin(3) * t33 - t38 * t27 / 0.2e1 - t39 * t26 / 0.2e1 + t22 * t41 / 0.2e1 + t21 * t42 / 0.2e1 + t43 * t40 + (-t85 / 0.2e1 + t35 / 0.2e1 + t34 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t120) * t81) * qJD(3)) * t94 + (t30 / 0.2e1 - t20 * mrSges(5,3) + t117 * t136 + (t126 / 0.2e1 - t51 / 0.2e1 - t37 * mrSges(5,3) + (m(6) * t57 + t31) * pkin(4)) * qJD(4) + (-qJD(4) * t69 - t48 + m(5) * (-t20 - t146)) * pkin(7)) * t90 + (-t2 * t99 - t3 * t65 + t38 * t7 - t39 * t8) * mrSges(6,3) + (t29 / 0.2e1 + qJD(4) * t52 / 0.2e1 + t77 * t117 / 0.2e1 + t103 * mrSges(5,3) + (m(5) * t103 - qJD(4) * t70 + t49) * pkin(7)) * t93 + (t93 * t68 / 0.2e1 + t67 * t135 + t81 * t66 + (t77 * t135 + t93 * t136) * qJD(4) + (t129 / 0.2e1 + t127 / 0.2e1 + Ifges(6,5) * t137 + Ifges(6,6) * t138 - Ifges(4,6) + t81 * mrSges(4,2)) * qJD(3)) * t91; m(6) * (-pkin(4) * t107 + t45 * t21 + t44 * t22 - t24 * t54 - t25 * t53) + (-t21 * t99 - t22 * t65 - t38 * t53 + t39 * t54) * mrSges(6,3) + (m(5) * (t119 * t132 - t133) + (m(6) * t82 + t120 + t40) * t91) * qJD(3) + (-t13 - t66 + (t119 * mrSges(5,3) - mrSges(4,2)) * qJD(3)) * t94; -t38 * t42 + t65 * t15 - t39 * t41 - t99 * t14 + (t82 * t110 + t24 * t45 + t25 * t44) * t140 + 0.2e1 * t40 * t110 + 0.2e1 * t82 * t13 - 0.2e1 * pkin(3) * t66 + t90 * t68 - t76 * t116 + (qJD(4) * t77 + t67) * t93 + 0.2e1 * (-t24 * t99 - t25 * t65 + t38 * t44 - t39 * t45) * mrSges(6,3); -Ifges(5,5) * t108 + t20 * mrSges(5,1) - t19 * mrSges(5,2) - t95 * Ifges(5,6) + (m(6) * (t8 * t112 - t7 * t113 + t2 * t89 + t3 * t92) + t46 * t112 + t89 * t12 - t47 * t113 + t92 * t11) * pkin(4) + t97 - t143; m(6) * (t21 * t89 + t22 * t92 + (t53 * t89 - t54 * t92) * qJD(5)) * pkin(4) + t144; t85 + (t75 * pkin(7) - t128) * qJD(4) + (m(6) * (t24 * t89 + t25 * t92 + (-t44 * t89 + t45 * t92) * qJD(5)) + (t92 * t38 - t89 * t39 + (t65 * t89 - t92 * t99) * qJD(5)) * mrSges(6,3)) * pkin(4) + t98; 0.2e1 * t62; t97; -t6; t98; t62; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
