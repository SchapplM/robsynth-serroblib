% Calculate time derivative of joint inertia matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:03
% EndTime: 2019-12-05 16:35:10
% DurationCPUTime: 1.40s
% Computational Cost: add. (1227->281), mult. (3382->464), div. (0->0), fcn. (2918->10), ass. (0->134)
t144 = 2 * qJD(4);
t82 = sin(pkin(5));
t90 = cos(qJ(2));
t120 = t82 * t90;
t87 = sin(qJ(2));
t121 = t82 * t87;
t84 = cos(pkin(5));
t86 = sin(qJ(3));
t89 = cos(qJ(3));
t56 = t89 * t121 + t84 * t86;
t81 = sin(pkin(10));
t83 = cos(pkin(10));
t31 = -t81 * t120 + t56 * t83;
t55 = t86 * t121 - t84 * t89;
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t10 = t31 * t88 + t55 * t85;
t9 = -t31 * t85 + t55 * t88;
t143 = (t10 * t88 - t85 * t9) * qJD(5);
t142 = 2 * m(5);
t141 = 2 * m(6);
t140 = 2 * pkin(7);
t80 = t83 ^ 2;
t139 = t81 / 0.2e1;
t138 = t83 / 0.2e1;
t137 = -t85 / 0.2e1;
t118 = t83 * t89;
t62 = (mrSges(5,1) * t86 - mrSges(5,3) * t118) * qJD(3);
t114 = t88 * t89;
t117 = t85 * t86;
t57 = -t83 * t117 - t114;
t21 = t57 * qJD(5) + (t83 * t114 + t117) * qJD(3);
t115 = t86 * t88;
t116 = t85 * t89;
t58 = t83 * t115 - t116;
t22 = -t58 * qJD(5) + (-t83 * t116 + t115) * qJD(3);
t8 = -mrSges(6,1) * t22 + mrSges(6,2) * t21;
t136 = t8 - t62;
t135 = mrSges(5,2) * t83;
t134 = mrSges(6,3) * t81;
t133 = Ifges(5,4) * t81;
t132 = Ifges(5,4) * t83;
t131 = Ifges(6,4) * t85;
t130 = Ifges(6,4) * t88;
t101 = qJD(2) * t121;
t110 = qJD(2) * t90;
t100 = t82 * t110;
t33 = -t55 * qJD(3) + t89 * t100;
t15 = -t83 * t101 + t33 * t81;
t30 = t83 * t120 + t56 * t81;
t129 = t15 * t30;
t32 = t56 * qJD(3) + t86 * t100;
t128 = t32 * t86;
t127 = t33 * t89;
t54 = -qJD(4) * t86 + (pkin(3) * t86 - qJ(4) * t89) * qJD(3);
t126 = t54 * t83;
t125 = t55 * t32;
t71 = -pkin(3) * t89 - qJ(4) * t86 - pkin(2);
t124 = t71 * t83;
t123 = t81 * t86;
t122 = t81 * t89;
t119 = t83 * mrSges(5,3);
t113 = -mrSges(5,1) * t83 + mrSges(5,2) * t81 - mrSges(4,1);
t25 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t67 = -mrSges(5,1) * t89 - t86 * t119;
t112 = t25 - t67;
t108 = qJD(3) * t89;
t99 = t81 * t108;
t49 = mrSges(5,1) * t99 + t108 * t135;
t44 = pkin(7) * t118 + t81 * t71;
t111 = qJ(4) * t83;
t109 = qJD(3) * t86;
t107 = qJD(4) * t83;
t106 = qJD(5) * t81;
t105 = qJD(5) * t85;
t104 = qJD(5) * t88;
t103 = qJ(4) * qJD(4);
t5 = Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t99;
t102 = pkin(7) * t109;
t98 = pkin(7) * t81 + pkin(4);
t96 = mrSges(6,1) * t85 + mrSges(6,2) * t88;
t95 = -Ifges(5,5) * t83 + Ifges(5,6) * t81;
t35 = -pkin(8) * t89 + t44;
t92 = pkin(4) * t81 - pkin(8) * t83 + pkin(7);
t48 = t92 * t86;
t11 = -t35 * t85 + t48 * t88;
t12 = t35 * t88 + t48 * t85;
t94 = t11 * t85 - t12 * t88;
t69 = -pkin(4) * t83 - pkin(8) * t81 - pkin(3);
t41 = -t85 * t111 + t69 * t88;
t42 = t88 * t111 + t69 * t85;
t93 = t41 * t85 - t42 * t88;
t91 = t55 * t108 + t128;
t51 = (-Ifges(6,5) * t85 - Ifges(6,6) * t88) * t106;
t79 = t81 ^ 2;
t77 = t79 * t103;
t68 = (mrSges(4,1) * t86 + mrSges(4,2) * t89) * qJD(3);
t66 = mrSges(5,2) * t89 - mrSges(5,3) * t123;
t65 = -mrSges(6,1) * t83 - t88 * t134;
t64 = mrSges(6,2) * t83 - t85 * t134;
t61 = (-mrSges(5,2) * t86 - mrSges(5,3) * t122) * qJD(3);
t60 = t96 * t81;
t59 = (mrSges(5,1) * t81 + t135) * t86;
t53 = (-Ifges(6,1) * t85 - t130) * t106;
t52 = (-Ifges(6,2) * t88 - t131) * t106;
t50 = (mrSges(6,1) * t88 - mrSges(6,2) * t85) * t106;
t47 = -Ifges(6,5) * t83 + (Ifges(6,1) * t88 - t131) * t81;
t46 = -Ifges(6,6) * t83 + (-Ifges(6,2) * t85 + t130) * t81;
t45 = t81 * t54;
t43 = -pkin(7) * t122 + t124;
t40 = t92 * t108;
t39 = (Ifges(5,5) * t86 + (Ifges(5,1) * t83 - t133) * t89) * qJD(3);
t38 = (Ifges(5,6) * t86 + (-Ifges(5,2) * t81 + t132) * t89) * qJD(3);
t37 = mrSges(6,1) * t123 - mrSges(6,3) * t58;
t36 = -mrSges(6,2) * t123 + mrSges(6,3) * t57;
t34 = t98 * t89 - t124;
t29 = -t83 * t102 + t45;
t28 = t81 * t102 + t126;
t27 = -t42 * qJD(5) - t85 * t107;
t26 = t41 * qJD(5) + t88 * t107;
t24 = t45 + (-pkin(7) * t83 + pkin(8)) * t109;
t23 = -t98 * t109 - t126;
t18 = Ifges(6,1) * t58 + Ifges(6,4) * t57 + Ifges(6,5) * t123;
t17 = Ifges(6,4) * t58 + Ifges(6,2) * t57 + Ifges(6,6) * t123;
t16 = t81 * t101 + t33 * t83;
t14 = -mrSges(6,2) * t99 + mrSges(6,3) * t22;
t13 = mrSges(6,1) * t99 - mrSges(6,3) * t21;
t7 = Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t99;
t6 = Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t99;
t4 = -t12 * qJD(5) - t24 * t85 + t40 * t88;
t3 = t11 * qJD(5) + t24 * t88 + t40 * t85;
t2 = t9 * qJD(5) + t16 * t88 + t32 * t85;
t1 = -t10 * qJD(5) - t16 * t85 + t32 * t88;
t19 = [0.2e1 * m(5) * (t16 * t31 + t125 + t129) + 0.2e1 * m(6) * (t1 * t9 + t10 * t2 + t129) + 0.2e1 * m(4) * (-t82 ^ 2 * t87 * t110 + t56 * t33 + t125); t1 * t37 + t10 * t14 + t9 * t13 + t16 * t66 + t2 * t36 + t31 * t61 + t32 * t59 + t55 * t49 + t136 * t30 + t112 * t15 + (-t90 * t68 + (-t90 * mrSges(3,2) + (-t89 * mrSges(4,1) + t86 * mrSges(4,2) - mrSges(3,1)) * t87) * qJD(2)) * t82 + m(5) * (-t15 * t43 + t16 * t44 - t28 * t30 + t29 * t31) + m(6) * (t1 * t11 + t10 * t3 + t12 * t2 + t15 * t34 + t23 * t30 + t4 * t9) - m(4) * pkin(2) * t101 + (m(5) * t91 / 0.2e1 + m(4) * (-t56 * t109 + t127 + t91) / 0.2e1) * t140 + (t128 + t127 + (t55 * t89 - t56 * t86) * qJD(3)) * mrSges(4,3); -0.2e1 * pkin(2) * t68 + 0.2e1 * t11 * t13 + 0.2e1 * t12 * t14 + t22 * t17 + t21 * t18 + 0.2e1 * t23 * t25 + 0.2e1 * t28 * t67 + 0.2e1 * t29 * t66 + 0.2e1 * t3 * t36 + 0.2e1 * t34 * t8 + 0.2e1 * t4 * t37 + 0.2e1 * t43 * t62 + 0.2e1 * t44 * t61 + t57 * t6 + t58 * t7 + (t28 * t43 + t29 * t44) * t142 + (t11 * t4 + t12 * t3 + t23 * t34) * t141 + (t49 * t140 + t83 * t39 + (-t38 + t5) * t81 + (-(2 * Ifges(4,4)) - t95) * t109) * t86 + (t81 * (Ifges(6,5) * t58 + Ifges(6,6) * t57) + t59 * t140 + 0.2e1 * (Ifges(4,4) + t95) * t89 + (-(2 * Ifges(5,3)) - (2 * Ifges(4,2)) + (2 * Ifges(4,1)) + (pkin(7) ^ 2 * t142) + Ifges(5,1) * t80 + (-0.2e1 * t132 + (Ifges(6,3) + Ifges(5,2)) * t81) * t81) * t86) * t108; t16 * t119 - t33 * mrSges(4,2) + t1 * t65 + t15 * t60 + t2 * t64 + t30 * t50 + t113 * t32 + (t15 * mrSges(5,3) - mrSges(6,3) * t143) * t81 + m(5) * (-pkin(3) * t32 + (t30 * t81 + t31 * t83) * qJD(4) + (t15 * t81 + t16 * t83) * qJ(4)) + m(6) * (t1 * t41 + t10 * t26 + t2 * t42 + t27 * t9 + (qJ(4) * t15 + qJD(4) * t30) * t81); t57 * t52 / 0.2e1 + t58 * t53 / 0.2e1 + t23 * t60 + t3 * t64 + t4 * t65 + t41 * t13 + t42 * t14 + t22 * t46 / 0.2e1 + t21 * t47 / 0.2e1 - pkin(3) * t49 + t34 * t50 + t26 * t36 + t27 * t37 + m(6) * (t11 * t27 + t12 * t26 + t3 * t42 + t4 * t41) + (t38 / 0.2e1 - t5 / 0.2e1 + qJ(4) * t61 + qJD(4) * t66 + t29 * mrSges(5,3) + m(5) * (qJ(4) * t29 + qJD(4) * t44)) * t83 + (t39 / 0.2e1 + t88 * t7 / 0.2e1 + t6 * t137 + t86 * t51 / 0.2e1 - t28 * mrSges(5,3) + t112 * qJD(4) + t136 * qJ(4) + m(6) * (qJ(4) * t23 + qJD(4) * t34) + m(5) * (-qJ(4) * t28 - qJD(4) * t43) + (-t88 * t17 / 0.2e1 + t18 * t137 + t94 * mrSges(6,3)) * qJD(5)) * t81 + ((pkin(7) * mrSges(4,2) + Ifges(5,5) * t139 + Ifges(5,6) * t138 - Ifges(4,6)) * t86 + (Ifges(4,5) + (-Ifges(6,3) * t83 + (Ifges(6,5) * t88 - Ifges(6,6) * t85) * t81) * t139 + (Ifges(5,1) * t81 + t132) * t138 - t81 * (Ifges(5,2) * t83 + t133) / 0.2e1 + (-m(5) * pkin(3) + t113) * pkin(7)) * t89) * qJD(3); 0.2e1 * t26 * t64 + 0.2e1 * t27 * t65 - t83 * t51 + (t79 + t80) * mrSges(5,3) * t144 + (t26 * t42 + t27 * t41 + t77) * t141 + (t80 * t103 + t77) * t142 + (0.2e1 * qJ(4) * t50 + t60 * t144 - t85 * t52 + t88 * t53 + (0.2e1 * mrSges(6,3) * t93 - t46 * t88 - t47 * t85) * qJD(5)) * t81; m(5) * t32 + m(6) * (t1 * t88 + t2 * t85 + t143); m(6) * (-t94 * qJD(5) + t3 * t85 + t4 * t88) + t36 * t104 + t85 * t14 - t37 * t105 + t88 * t13 + m(5) * pkin(7) * t108 + t49; m(6) * (-t93 * qJD(5) + t26 * t85 + t27 * t88) + t64 * t104 - t65 * t105; 0; mrSges(6,1) * t1 - mrSges(6,2) * t2; mrSges(6,1) * t4 - mrSges(6,2) * t3 + t5; mrSges(6,1) * t27 - mrSges(6,2) * t26 + t51; -t96 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
