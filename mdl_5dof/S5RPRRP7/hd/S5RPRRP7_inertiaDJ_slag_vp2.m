% Calculate time derivative of joint inertia matrix for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:34
% EndTime: 2019-12-31 18:44:38
% DurationCPUTime: 1.26s
% Computational Cost: add. (867->237), mult. (2089->339), div. (0->0), fcn. (1393->6), ass. (0->108)
t115 = Ifges(6,4) + Ifges(5,5);
t135 = -Ifges(6,2) - Ifges(5,3);
t68 = sin(qJ(4));
t70 = cos(qJ(4));
t106 = t68 ^ 2 + t70 ^ 2;
t101 = qJD(4) * t70;
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t104 = qJD(3) * t71;
t92 = t68 * t104;
t73 = t69 * t101 + t92;
t59 = sin(pkin(8)) * pkin(1) + pkin(6);
t121 = t59 * t71;
t88 = -cos(pkin(8)) * pkin(1) - pkin(2);
t32 = -pkin(3) * t71 - t69 * pkin(7) + t88;
t132 = t70 * t121 + t68 * t32;
t134 = qJD(4) * t132;
t133 = 0.2e1 * t88;
t100 = qJD(5) * t70;
t79 = pkin(4) * t70 + qJ(5) * t68;
t131 = t79 * qJD(4) - t100;
t44 = (pkin(3) * t69 - pkin(7) * t71) * qJD(3);
t130 = -t70 * t44 + t134;
t103 = qJD(4) * t68;
t105 = qJD(3) * t69;
t112 = t32 * t101 + t68 * t44;
t1 = (-t59 * t103 - qJD(5)) * t71 + (-t59 * t70 + qJ(5)) * t105 + t112;
t119 = t69 * t70;
t40 = -t71 * mrSges(5,1) - mrSges(5,3) * t119;
t41 = t71 * mrSges(6,1) + mrSges(6,2) * t119;
t109 = -t40 + t41;
t118 = t70 * t32;
t87 = t59 * t68 + pkin(4);
t11 = t87 * t71 - t118;
t14 = -t68 * t121 + t118;
t20 = -mrSges(5,2) * t105 - t73 * mrSges(5,3);
t21 = -t73 * mrSges(6,2) + mrSges(6,3) * t105;
t3 = (-t71 * t103 - t70 * t105) * t59 + t112;
t129 = t109 * qJD(4) + m(5) * (-t14 * qJD(4) + t3) + m(6) * (t11 * qJD(4) + t1) + t20 + t21;
t10 = -qJ(5) * t71 + t132;
t120 = t68 * t69;
t39 = mrSges(5,2) * t71 - mrSges(5,3) * t120;
t42 = -mrSges(6,2) * t120 - mrSges(6,3) * t71;
t110 = t39 + t42;
t102 = qJD(4) * t69;
t91 = t70 * t104;
t74 = -t68 * t102 + t91;
t18 = mrSges(5,1) * t105 - t74 * mrSges(5,3);
t19 = mrSges(6,2) * t91 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t103) * t69;
t2 = -t87 * t105 + t130;
t94 = t59 * t105;
t4 = t68 * t94 - t130;
t128 = -t110 * qJD(4) + m(5) * (-t4 - t134) + m(6) * (-t10 * qJD(4) + t2) - t18 + t19;
t127 = 0.2e1 * m(5);
t126 = 0.2e1 * t59;
t125 = Ifges(5,4) * t68;
t124 = Ifges(5,4) * t70;
t123 = Ifges(6,5) * t68;
t122 = Ifges(6,5) * t70;
t116 = t71 * Ifges(5,6);
t82 = Ifges(6,1) * t70 + t123;
t25 = -Ifges(6,4) * t71 + t82 * t69;
t83 = Ifges(5,1) * t70 - t125;
t26 = -Ifges(5,5) * t71 + t83 * t69;
t111 = t25 + t26;
t47 = -t70 * mrSges(5,1) + t68 * mrSges(5,2);
t108 = t47 - mrSges(4,1);
t107 = t106 * pkin(7) * t104;
t48 = -Ifges(6,3) * t70 + t123;
t49 = Ifges(5,2) * t70 + t125;
t90 = t48 / 0.2e1 - t49 / 0.2e1;
t50 = Ifges(6,1) * t68 - t122;
t51 = Ifges(5,1) * t68 + t124;
t89 = t50 / 0.2e1 + t51 / 0.2e1;
t80 = Ifges(6,3) * t68 + t122;
t23 = -Ifges(6,6) * t71 + t80 * t69;
t81 = -Ifges(5,2) * t68 + t124;
t24 = t81 * t69 - t116;
t86 = t23 - t24 + t116;
t85 = t68 * mrSges(5,1) + t70 * mrSges(5,2);
t46 = -t70 * mrSges(6,1) - t68 * mrSges(6,3);
t84 = t68 * mrSges(6,1) - t70 * mrSges(6,3);
t78 = pkin(4) * t68 - qJ(5) * t70;
t77 = -t73 * Ifges(6,6) + t135 * t105 - t115 * t91;
t75 = t59 + t78;
t72 = m(6) * t100 + (-m(6) * t79 + t46 + t47) * qJD(4);
t64 = Ifges(6,4) * t101;
t63 = Ifges(5,5) * t101;
t61 = Ifges(6,6) * t103;
t45 = -pkin(3) - t79;
t38 = t83 * qJD(4);
t37 = t82 * qJD(4);
t36 = t81 * qJD(4);
t35 = t80 * qJD(4);
t34 = t85 * qJD(4);
t33 = t84 * qJD(4);
t31 = t85 * t69;
t30 = t84 * t69;
t28 = t78 * qJD(4) - qJD(5) * t68;
t17 = t75 * t69;
t13 = t73 * mrSges(5,1) + t74 * mrSges(5,2);
t12 = t73 * mrSges(6,1) - t74 * mrSges(6,3);
t9 = -t51 * t102 + (Ifges(5,5) * t69 + t83 * t71) * qJD(3);
t8 = -t50 * t102 + (Ifges(6,4) * t69 + t82 * t71) * qJD(3);
t7 = -t49 * t102 + (Ifges(5,6) * t69 + t81 * t71) * qJD(3);
t6 = -t48 * t102 + (Ifges(6,6) * t69 + t80 * t71) * qJD(3);
t5 = t75 * t104 + t131 * t69;
t15 = [0.2e1 * t1 * t42 + 0.2e1 * t10 * t21 + 0.2e1 * t11 * t19 + 0.2e1 * t17 * t12 + 0.2e1 * t14 * t18 + 0.2e1 * t132 * t20 + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t39 + 0.2e1 * t5 * t30 + 0.2e1 * t4 * t40 + (t132 * t3 + t14 * t4) * t127 + 0.2e1 * m(6) * (t1 * t10 + t11 * t2 + t17 * t5) + ((mrSges(4,2) * t133 + 0.2e1 * Ifges(4,4) * t71 + t111 * t70 + t31 * t126 + t86 * t68) * qJD(3) + t77) * t71 + (t13 * t126 + (t8 + t9) * t70 + (t6 - t7) * t68 + (t86 * t70 + (t115 * t71 - t111) * t68) * qJD(4) + (mrSges(4,1) * t133 + (-0.2e1 * Ifges(4,4) + t115 * t70 + (-Ifges(5,6) + Ifges(6,6)) * t68) * t69 + (t59 ^ 2 * t127 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t135) * t71) * qJD(3)) * t69; (-t13 - t12 - m(6) * t5 + (t110 * t70 + t109 * t68 + m(6) * (t10 * t70 + t11 * t68) + (t132 * t70 - t14 * t68 - t121) * m(5)) * qJD(3)) * t71 + (m(5) * t94 + (m(6) * t17 + t30 + t31) * qJD(3) + t129 * t70 + t128 * t68) * t69; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t106) * t69 * t104; t45 * t12 + t5 * t46 + t28 * t30 + t17 * t33 - pkin(3) * t13 + m(6) * (t17 * t28 + t45 * t5) + (-t63 / 0.2e1 - t64 / 0.2e1 - t61 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t108) * t59) * qJD(3)) * t71 + (t2 * mrSges(6,2) - t4 * mrSges(5,3) + t8 / 0.2e1 + t9 / 0.2e1 + t90 * t104 + (t23 / 0.2e1 - t24 / 0.2e1 + t116 / 0.2e1 - t10 * mrSges(6,2) - t132 * mrSges(5,3)) * qJD(4) + t128 * pkin(7)) * t68 + (t1 * mrSges(6,2) + t3 * mrSges(5,3) - t6 / 0.2e1 + t7 / 0.2e1 + t89 * t104 + (t25 / 0.2e1 + t26 / 0.2e1 + t11 * mrSges(6,2) - t14 * mrSges(5,3)) * qJD(4) + t129 * pkin(7)) * t70 + (t59 * t34 + (t37 / 0.2e1 + t38 / 0.2e1) * t70 + (t35 / 0.2e1 - t36 / 0.2e1) * t68 + (t59 * mrSges(4,2) - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t70 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t68) * qJD(3) + (-t89 * t68 + t90 * t70) * qJD(4)) * t69; (t46 + t108) * t105 + m(5) * (-pkin(3) * t105 + t107) + m(6) * (t45 * t105 + t107) + (-t34 - t33 - m(6) * t28 + (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t106) * qJD(3)) * t71; -0.2e1 * pkin(3) * t34 + 0.2e1 * t33 * t45 + (-t35 + t36) * t70 + (t37 + t38) * t68 + 0.2e1 * (m(6) * t45 + t46) * t28 + ((t50 + t51) * t70 + (t48 - t49) * t68) * qJD(4); -Ifges(5,6) * t92 + m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t10) - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t19 + t1 * mrSges(6,3) + qJD(5) * t42 + qJ(5) * t21 + (-Ifges(5,6) * t70 - t115 * t68) * t102 - t77; (-m(6) * t78 - t84 - t85) * t104 + t72 * t69; -t131 * mrSges(6,2) - Ifges(5,6) * t103 + t72 * pkin(7) + t61 + t63 + t64; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t2 + t19; t73 * m(6); (m(6) * pkin(7) + mrSges(6,2)) * t101; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
