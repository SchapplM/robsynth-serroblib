% Calculate time derivative of joint inertia matrix for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:28
% EndTime: 2019-12-31 18:53:31
% DurationCPUTime: 1.42s
% Computational Cost: add. (1533->234), mult. (3554->333), div. (0->0), fcn. (3162->6), ass. (0->104)
t113 = Ifges(6,4) + Ifges(5,5);
t136 = Ifges(6,2) + Ifges(5,3);
t80 = cos(qJ(4));
t103 = qJD(4) * t80;
t76 = sin(pkin(8));
t77 = cos(pkin(8));
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t53 = t76 * t79 - t81 * t77;
t50 = t53 * qJD(3);
t78 = sin(qJ(4));
t115 = t78 * t50;
t54 = t76 * t81 + t79 * t77;
t83 = t54 * t103 - t115;
t104 = qJD(4) * t78;
t114 = t80 * t50;
t82 = t54 * t104 + t114;
t97 = -pkin(2) * t77 - pkin(1);
t33 = t53 * pkin(3) - t54 * pkin(7) + t97;
t112 = pkin(6) + qJ(2);
t62 = t112 * t76;
t63 = t112 * t77;
t40 = -t79 * t62 + t63 * t81;
t133 = t78 * t33 + t80 * t40;
t135 = qJD(4) * t133;
t132 = -t81 * t62 - t63 * t79;
t131 = Ifges(6,6) * t104 + t103 * t113;
t26 = -t53 * qJD(2) + t132 * qJD(3);
t51 = t54 * qJD(3);
t32 = pkin(3) * t51 + pkin(7) * t50;
t4 = -t26 * t78 + t32 * t80 - t135;
t102 = qJD(5) * t80;
t88 = pkin(4) * t80 + qJ(5) * t78;
t130 = t88 * qJD(4) - t102;
t129 = -2 * mrSges(4,3);
t128 = -2 * Ifges(4,4);
t126 = -0.2e1 * t132;
t124 = Ifges(5,4) * t78;
t123 = Ifges(5,4) * t80;
t122 = Ifges(6,5) * t78;
t121 = Ifges(6,5) * t80;
t120 = Ifges(5,6) * t53;
t27 = t54 * qJD(2) + t40 * qJD(3);
t119 = t27 * t132;
t118 = t54 * t78;
t117 = t54 * t80;
t17 = t51 * mrSges(5,1) + t82 * mrSges(5,3);
t18 = -t51 * mrSges(6,1) - t82 * mrSges(6,2);
t111 = t17 - t18;
t19 = -t51 * mrSges(5,2) - t83 * mrSges(5,3);
t20 = -t83 * mrSges(6,2) + t51 * mrSges(6,3);
t110 = t19 + t20;
t89 = Ifges(6,3) * t78 + t121;
t21 = Ifges(6,6) * t53 + t89 * t54;
t90 = -Ifges(5,2) * t78 + t123;
t22 = t90 * t54 + t120;
t109 = t21 - t22;
t91 = Ifges(6,1) * t80 + t122;
t23 = Ifges(6,4) * t53 + t91 * t54;
t92 = Ifges(5,1) * t80 - t124;
t24 = Ifges(5,5) * t53 + t92 * t54;
t108 = t23 + t24;
t34 = -t53 * mrSges(5,2) - mrSges(5,3) * t118;
t37 = -mrSges(6,2) * t118 + t53 * mrSges(6,3);
t107 = t34 + t37;
t35 = t53 * mrSges(5,1) - mrSges(5,3) * t117;
t36 = -t53 * mrSges(6,1) + mrSges(6,2) * t117;
t106 = -t35 + t36;
t65 = -Ifges(6,3) * t80 + t122;
t66 = Ifges(5,2) * t80 + t124;
t99 = t65 / 0.2e1 - t66 / 0.2e1;
t67 = Ifges(6,1) * t78 - t121;
t68 = Ifges(5,1) * t78 + t123;
t98 = t67 / 0.2e1 + t68 / 0.2e1;
t96 = t51 * mrSges(4,1) - t50 * mrSges(4,2);
t94 = t78 * mrSges(5,1) + t80 * mrSges(5,2);
t64 = -t80 * mrSges(6,1) - t78 * mrSges(6,3);
t93 = t78 * mrSges(6,1) - t80 * mrSges(6,3);
t87 = pkin(4) * t78 - qJ(5) * t80;
t14 = t33 * t80 - t40 * t78;
t84 = t83 * Ifges(6,6) - t113 * t114 + t136 * t51;
t3 = t33 * t103 - t40 * t104 + t80 * t26 + t78 * t32;
t61 = -pkin(3) - t88;
t60 = t92 * qJD(4);
t59 = t91 * qJD(4);
t58 = t90 * qJD(4);
t57 = t89 * qJD(4);
t56 = t94 * qJD(4);
t55 = t93 * qJD(4);
t49 = -pkin(4) * t104 + qJ(5) * t103 + qJD(5) * t78;
t31 = t93 * t54;
t16 = t87 * t54 - t132;
t13 = t83 * mrSges(5,1) - t82 * mrSges(5,2);
t12 = t83 * mrSges(6,1) + t82 * mrSges(6,3);
t11 = -t53 * pkin(4) - t14;
t10 = qJ(5) * t53 + t133;
t9 = -t82 * Ifges(5,1) - t83 * Ifges(5,4) + Ifges(5,5) * t51;
t8 = -t82 * Ifges(6,1) + Ifges(6,4) * t51 + t83 * Ifges(6,5);
t7 = -t82 * Ifges(5,4) - t83 * Ifges(5,2) + Ifges(5,6) * t51;
t6 = -t82 * Ifges(6,5) + Ifges(6,6) * t51 + t83 * Ifges(6,3);
t5 = t130 * t54 - t87 * t50 + t27;
t2 = -t51 * pkin(4) - t4;
t1 = t51 * qJ(5) + t53 * qJD(5) + t3;
t15 = [t40 * t51 * t129 + 0.2e1 * t97 * t96 + 0.2e1 * t3 * t34 + 0.2e1 * t4 * t35 + 0.2e1 * t2 * t36 + 0.2e1 * t1 * t37 + t13 * t126 + 0.2e1 * t11 * t18 + 0.2e1 * t133 * t19 + 0.2e1 * t10 * t20 + 0.2e1 * t5 * t31 + 0.2e1 * t16 * t12 + 0.2e1 * t14 * t17 - (mrSges(4,3) * t126 + t108 * t80 + t109 * t78) * t50 + 0.2e1 * m(5) * (t133 * t3 + t14 * t4 - t119) + 0.2e1 * m(4) * (t26 * t40 - t119) + 0.2e1 * m(6) * (t1 * t10 + t11 * t2 + t16 * t5) + (t26 * t129 - (-Ifges(5,6) * t78 + t128) * t50 + ((2 * Ifges(4,2)) + t136) * t51 + t84) * t53 + (-0.2e1 * Ifges(4,1) * t50 + t128 * t51 + (t8 + t9 + t113 * t51 + (t109 - t120) * qJD(4)) * t80 + (t6 - t7 + (-Ifges(5,6) + Ifges(6,6)) * t51 + (-t113 * t53 - t108) * qJD(4)) * t78 + 0.2e1 * (mrSges(4,3) + t94) * t27) * t54 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t76 ^ 2 + t77 ^ 2) * qJD(2); t111 * t80 + t110 * t78 + (t106 * t78 + t107 * t80) * qJD(4) + m(6) * (t1 * t78 - t2 * t80 + (t10 * t80 + t11 * t78) * qJD(4)) + m(5) * (t3 * t78 + t4 * t80 + (t133 * t80 - t14 * t78) * qJD(4)) + t96; 0; t61 * t12 + t5 * t64 + t16 * t55 - t132 * t56 - t49 * t31 - Ifges(4,5) * t50 - Ifges(4,6) * t51 - t26 * mrSges(4,2) - t27 * mrSges(4,1) + m(6) * (-t49 * t16 + t61 * t5) + (t2 * mrSges(6,2) - t4 * mrSges(5,3) + t27 * mrSges(5,2) + t8 / 0.2e1 + t9 / 0.2e1 + (t57 / 0.2e1 - t58 / 0.2e1) * t54 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t51 - t99 * t50 + (t21 / 0.2e1 - t22 / 0.2e1 - t120 / 0.2e1 - t10 * mrSges(6,2) - t133 * mrSges(5,3) - t98 * t54) * qJD(4) + (-t107 * qJD(4) + m(6) * (-t10 * qJD(4) + t2) + m(5) * (-t4 - t135) - t111) * pkin(7)) * t78 + (t3 * mrSges(5,3) + t1 * mrSges(6,2) - t27 * mrSges(5,1) - t6 / 0.2e1 + t7 / 0.2e1 + (t59 / 0.2e1 + t60 / 0.2e1) * t54 + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t51 - t98 * t50 + (t23 / 0.2e1 + t24 / 0.2e1 + t11 * mrSges(6,2) - t14 * mrSges(5,3) + t99 * t54) * qJD(4) + (t106 * qJD(4) + m(6) * (t11 * qJD(4) + t1) + m(5) * (-t14 * qJD(4) + t3) + t110) * pkin(7)) * t80 + t131 * t53 / 0.2e1 + (-m(5) * t27 - t13) * pkin(3); 0; -0.2e1 * pkin(3) * t56 + 0.2e1 * t61 * t55 + (-t57 + t58) * t80 + (t59 + t60) * t78 + 0.2e1 * (-m(6) * t61 - t64) * t49 + ((t67 + t68) * t80 + (t65 - t66) * t78) * qJD(4); Ifges(5,6) * t115 - pkin(4) * t18 + m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t10) + qJD(5) * t37 + qJ(5) * t20 + t1 * mrSges(6,3) - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t2 * mrSges(6,1) + (-Ifges(5,6) * t80 - t113 * t78) * t54 * qJD(4) + t84; m(6) * t49 + ((-mrSges(5,2) + mrSges(6,3)) * t80 + (-mrSges(5,1) - mrSges(6,1)) * t78) * qJD(4); -Ifges(5,6) * t104 - t130 * mrSges(6,2) + (m(6) * t102 + (-m(6) * t88 - t80 * mrSges(5,1) + t78 * mrSges(5,2) + t64) * qJD(4)) * pkin(7) + t131; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t2 + t18; m(6) * t104; (m(6) * pkin(7) + mrSges(6,2)) * t103; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
