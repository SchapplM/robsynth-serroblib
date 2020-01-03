% Calculate time derivative of joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:56
% EndTime: 2019-12-31 18:51:00
% DurationCPUTime: 1.41s
% Computational Cost: add. (1536->244), mult. (3592->347), div. (0->0), fcn. (3193->6), ass. (0->112)
t118 = Ifges(5,5) + Ifges(6,5);
t139 = Ifges(5,3) + Ifges(6,3);
t85 = cos(qJ(4));
t110 = qJD(4) * t85;
t81 = sin(pkin(8));
t82 = cos(pkin(8));
t84 = sin(qJ(3));
t86 = cos(qJ(3));
t57 = t81 * t86 + t84 * t82;
t107 = t57 * t110;
t56 = t81 * t84 - t86 * t82;
t53 = t56 * qJD(3);
t83 = sin(qJ(4));
t89 = -t83 * t53 + t107;
t111 = qJD(4) * t83;
t119 = t85 * t53;
t88 = t57 * t111 + t119;
t103 = -pkin(2) * t82 - pkin(1);
t35 = t56 * pkin(3) - t57 * pkin(7) + t103;
t116 = pkin(6) + qJ(2);
t64 = t116 * t81;
t65 = t116 * t82;
t42 = -t84 * t64 + t65 * t86;
t40 = t85 * t42;
t15 = t83 * t35 + t40;
t137 = -t86 * t64 - t65 * t84;
t117 = -Ifges(5,6) - Ifges(6,6);
t136 = (t117 * t85 - t118 * t83) * qJD(4);
t135 = 2 * m(6);
t134 = -2 * mrSges(4,3);
t133 = -2 * mrSges(6,3);
t131 = -0.2e1 * t137;
t130 = m(6) * pkin(4);
t129 = mrSges(5,2) * t85;
t128 = Ifges(5,4) * t83;
t127 = Ifges(5,4) * t85;
t126 = Ifges(6,4) * t83;
t125 = Ifges(6,4) * t85;
t26 = t57 * qJD(2) + t42 * qJD(3);
t124 = t26 * t137;
t123 = t57 * t83;
t122 = t57 * t85;
t115 = -qJ(5) - pkin(7);
t91 = -Ifges(6,2) * t83 + t125;
t20 = Ifges(6,6) * t56 + t91 * t57;
t92 = -Ifges(5,2) * t83 + t127;
t21 = Ifges(5,6) * t56 + t92 * t57;
t114 = -t20 - t21;
t93 = Ifges(6,1) * t85 - t126;
t22 = Ifges(6,5) * t56 + t93 * t57;
t94 = Ifges(5,1) * t85 - t128;
t23 = Ifges(5,5) * t56 + t94 * t57;
t113 = t22 + t23;
t58 = mrSges(6,1) * t111 + mrSges(6,2) * t110;
t112 = qJ(5) * t57;
t25 = -t56 * qJD(2) + t137 * qJD(3);
t54 = t57 * qJD(3);
t34 = pkin(3) * t54 + pkin(7) * t53;
t109 = t35 * t110 + t85 * t25 + t83 * t34;
t106 = -Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t69 = Ifges(6,2) * t85 + t126;
t70 = Ifges(5,2) * t85 + t128;
t105 = -t69 / 0.2e1 - t70 / 0.2e1;
t71 = Ifges(6,1) * t83 + t125;
t72 = Ifges(5,1) * t83 + t127;
t104 = t71 / 0.2e1 + t72 / 0.2e1;
t102 = t117 * t83;
t101 = t54 * mrSges(4,1) - t53 * mrSges(4,2);
t100 = -t25 * t83 + t85 * t34;
t14 = t85 * t35 - t42 * t83;
t99 = qJD(4) * t115;
t97 = -t118 * t119 + t139 * t54;
t96 = -(2 * Ifges(4,4)) + t102;
t95 = mrSges(5,1) * t83 + t129;
t90 = qJ(5) * t53 - qJD(5) * t57;
t12 = t89 * mrSges(6,1) - t88 * mrSges(6,2);
t78 = Ifges(5,5) * t110;
t77 = Ifges(6,5) * t110;
t74 = -pkin(4) * t85 - pkin(3);
t68 = t115 * t85;
t67 = -mrSges(6,1) * t85 + mrSges(6,2) * t83;
t66 = t115 * t83;
t63 = t94 * qJD(4);
t62 = t93 * qJD(4);
t61 = t92 * qJD(4);
t60 = t91 * qJD(4);
t59 = t95 * qJD(4);
t52 = -qJD(5) * t83 + t85 * t99;
t51 = qJD(5) * t85 + t83 * t99;
t39 = t56 * mrSges(5,1) - mrSges(5,3) * t122;
t38 = t56 * mrSges(6,1) - mrSges(6,3) * t122;
t37 = -t56 * mrSges(5,2) - mrSges(5,3) * t123;
t36 = -t56 * mrSges(6,2) - mrSges(6,3) * t123;
t33 = (mrSges(6,1) * t83 + mrSges(6,2) * t85) * t57;
t27 = pkin(4) * t123 - t137;
t19 = -t54 * mrSges(5,2) - t89 * mrSges(5,3);
t18 = -t54 * mrSges(6,2) - t89 * mrSges(6,3);
t17 = t54 * mrSges(5,1) + t88 * mrSges(5,3);
t16 = t54 * mrSges(6,1) + t88 * mrSges(6,3);
t13 = t89 * mrSges(5,1) - t88 * mrSges(5,2);
t11 = t89 * pkin(4) + t26;
t10 = -t83 * t112 + t15;
t9 = -t88 * Ifges(5,1) - t89 * Ifges(5,4) + Ifges(5,5) * t54;
t8 = -t88 * Ifges(6,1) - t89 * Ifges(6,4) + Ifges(6,5) * t54;
t7 = -t88 * Ifges(5,4) - t89 * Ifges(5,2) + Ifges(5,6) * t54;
t6 = -t88 * Ifges(6,4) - t89 * Ifges(6,2) + Ifges(6,6) * t54;
t5 = t56 * pkin(4) - t85 * t112 + t14;
t4 = -t15 * qJD(4) + t100;
t3 = -t42 * t111 + t109;
t2 = -qJ(5) * t107 + (-qJD(4) * t42 + t90) * t83 + t109;
t1 = t54 * pkin(4) + t90 * t85 + (-t40 + (-t35 + t112) * t83) * qJD(4) + t100;
t24 = [0.2e1 * t103 * t101 + t13 * t131 + 0.2e1 * t27 * t12 + 0.2e1 * t11 * t33 + 0.2e1 * t2 * t36 + 0.2e1 * t3 * t37 + 0.2e1 * t1 * t38 + 0.2e1 * t4 * t39 + 0.2e1 * t5 * t16 + 0.2e1 * t14 * t17 + 0.2e1 * t10 * t18 + 0.2e1 * t15 * t19 + t42 * t54 * t134 - (mrSges(4,3) * t131 + t113 * t85 + t114 * t83) * t53 + (t1 * t5 + t10 * t2 + t11 * t27) * t135 + 0.2e1 * m(5) * (t14 * t4 + t15 * t3 - t124) + 0.2e1 * m(4) * (t25 * t42 - t124) + (t25 * t134 + ((2 * Ifges(4,2)) + t139) * t54 - t96 * t53 + t97) * t56 + (t56 * t136 - 0.2e1 * Ifges(4,1) * t53 + t96 * t54 + (t114 * qJD(4) + t118 * t54 + t8 + t9) * t85 + (-t113 * qJD(4) - t6 - t7) * t83 + 0.2e1 * (mrSges(4,3) + t95) * t26) * t57 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t81 ^ 2 + t82 ^ 2) * qJD(2); (t16 + t17) * t85 + (t18 + t19) * t83 + ((t36 + t37) * t85 + (-t38 - t39) * t83) * qJD(4) + m(5) * (t3 * t83 + t4 * t85 + (-t14 * t83 + t15 * t85) * qJD(4)) + m(6) * (t1 * t85 + t2 * t83 + (t10 * t85 - t5 * t83) * qJD(4)) + t101; 0; -t26 * mrSges(4,1) - t25 * mrSges(4,2) - Ifges(4,5) * t53 - Ifges(4,6) * t54 + t11 * t67 + t74 * t12 + t66 * t16 - t68 * t18 + t27 * t58 + t51 * t36 + t52 * t38 - t137 * t59 + (t77 / 0.2e1 + t78 / 0.2e1) * t56 + m(6) * (t1 * t66 + t51 * t10 + t11 * t74 - t2 * t68 + t52 * t5) + (-t1 * mrSges(6,3) - t4 * mrSges(5,3) + t26 * mrSges(5,2) + t8 / 0.2e1 + t9 / 0.2e1 + (-t60 / 0.2e1 - t61 / 0.2e1) * t57 + (Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t54 - t105 * t53 + (-m(5) * t4 - t17) * pkin(7) + (-t20 / 0.2e1 - t21 / 0.2e1 - t10 * mrSges(6,3) - t15 * mrSges(5,3) + pkin(4) * t33 - t104 * t57 + t106 * t56 + t27 * t130 + (-m(5) * t15 - t37) * pkin(7)) * qJD(4)) * t83 + (t2 * mrSges(6,3) + t3 * mrSges(5,3) - t26 * mrSges(5,1) + t6 / 0.2e1 + t7 / 0.2e1 + (t62 / 0.2e1 + t63 / 0.2e1) * t57 - t106 * t54 - t104 * t53 + (m(5) * t3 + t19) * pkin(7) + (t22 / 0.2e1 + t23 / 0.2e1 - t14 * mrSges(5,3) - t5 * mrSges(6,3) + t105 * t57 + (-m(5) * t14 - t39) * pkin(7)) * qJD(4)) * t85 + (-m(5) * t26 - t13) * pkin(3); m(6) * (t51 * t83 + t52 * t85 + (-t66 * t83 - t68 * t85) * qJD(4)); (-t51 * t68 + t52 * t66) * t135 + 0.2e1 * t74 * t58 - 0.2e1 * pkin(3) * t59 + (t52 * t133 + t62 + t63 + (-t68 * t133 - t69 - t70 + 0.2e1 * (m(6) * t74 + t67) * pkin(4)) * qJD(4)) * t83 + (0.2e1 * t51 * mrSges(6,3) + t60 + t61 + (t66 * t133 + t71 + t72) * qJD(4)) * t85; t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) - t53 * t102 + (m(6) * t1 + t16) * pkin(4) + t57 * t136 + t97; (-t129 + (-mrSges(5,1) - t130) * t83) * qJD(4) - t58; -t51 * mrSges(6,2) + t77 + t78 + (mrSges(6,1) + t130) * t52 + ((-mrSges(5,1) * pkin(7) - (mrSges(6,3) * pkin(4))) * t85 + (mrSges(5,2) * pkin(7) + t117) * t83) * qJD(4); 0; m(6) * t11 + t12; 0; t111 * t130 + t58; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t24(1), t24(2), t24(4), t24(7), t24(11); t24(2), t24(3), t24(5), t24(8), t24(12); t24(4), t24(5), t24(6), t24(9), t24(13); t24(7), t24(8), t24(9), t24(10), t24(14); t24(11), t24(12), t24(13), t24(14), t24(15);];
Mq = res;
