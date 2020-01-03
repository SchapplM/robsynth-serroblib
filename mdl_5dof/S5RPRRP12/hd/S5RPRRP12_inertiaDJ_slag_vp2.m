% Calculate time derivative of joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:17
% EndTime: 2019-12-31 18:56:21
% DurationCPUTime: 1.73s
% Computational Cost: add. (853->257), mult. (1997->378), div. (0->0), fcn. (1295->4), ass. (0->120)
t130 = 2 * qJD(2);
t140 = Ifges(5,6) + Ifges(6,6);
t139 = Ifges(5,3) + Ifges(6,3);
t70 = sin(qJ(3));
t108 = qJD(3) * t70;
t69 = sin(qJ(4));
t94 = t69 * t108;
t72 = cos(qJ(3));
t101 = qJD(4) * t72;
t71 = cos(qJ(4));
t95 = t71 * t101;
t75 = t94 - t95;
t138 = Ifges(5,5) + Ifges(6,5);
t73 = -pkin(1) - pkin(6);
t121 = t70 * t73;
t126 = pkin(7) * t72;
t127 = pkin(3) * t70;
t46 = qJ(2) - t126 + t127;
t20 = t71 * t121 + t69 * t46;
t110 = t69 ^ 2 + t71 ^ 2;
t105 = qJD(3) * t73;
t14 = -pkin(4) * t75 + t105 * t70;
t107 = qJD(3) * t71;
t92 = t70 * t107;
t96 = t69 * t101;
t76 = t92 + t96;
t9 = -t75 * mrSges(6,1) - mrSges(6,2) * t76;
t137 = m(6) * t14 + t9;
t106 = qJD(3) * t72;
t16 = mrSges(5,1) * t106 + mrSges(5,3) * t76;
t34 = qJD(2) + (pkin(3) * t72 + pkin(7) * t70) * qJD(3);
t74 = -t20 * qJD(4) + t71 * t34;
t97 = t72 * t105;
t4 = -t69 * t97 + t74;
t136 = -m(5) * t4 - t16;
t100 = qJD(5) * t71;
t104 = qJD(4) * t69;
t87 = -t69 * t73 + pkin(4);
t1 = qJ(5) * t92 + (qJ(5) * t104 + qJD(3) * t87 - t100) * t72 + t74;
t15 = mrSges(6,1) * t106 + mrSges(6,3) * t76;
t135 = m(6) * t1 + t15;
t49 = -mrSges(6,1) * t71 + mrSges(6,2) * t69;
t60 = -pkin(4) * t71 - pkin(3);
t134 = m(6) * t60 + t49;
t18 = -mrSges(5,2) * t106 + mrSges(5,3) * t75;
t33 = t71 * t46;
t19 = -t121 * t69 + t33;
t103 = qJD(4) * t70;
t93 = t69 * t103;
t102 = qJD(4) * t71;
t98 = t46 * t102 + t69 * t34 + t71 * t97;
t3 = -t73 * t93 + t98;
t133 = m(5) * (-t19 * qJD(4) + t3) + t18;
t132 = 0.2e1 * m(6);
t131 = -0.2e1 * mrSges(6,3);
t129 = m(6) * pkin(4);
t125 = Ifges(5,4) * t69;
t124 = Ifges(5,4) * t71;
t123 = Ifges(6,4) * t69;
t122 = Ifges(6,4) * t71;
t120 = t72 * mrSges(5,3);
t119 = t72 * mrSges(6,3);
t118 = mrSges(5,2) + mrSges(6,2);
t115 = -qJ(5) - pkin(7);
t78 = -Ifges(6,2) * t69 + t122;
t21 = Ifges(6,6) * t70 + t72 * t78;
t79 = -Ifges(5,2) * t69 + t124;
t22 = Ifges(5,6) * t70 + t72 * t79;
t114 = t21 + t22;
t42 = -mrSges(6,2) * t70 - t119 * t69;
t43 = -mrSges(5,2) * t70 - t120 * t69;
t113 = t42 + t43;
t44 = mrSges(6,1) * t70 - t119 * t71;
t45 = mrSges(5,1) * t70 - t120 * t71;
t112 = -t44 - t45;
t111 = -mrSges(5,1) * t71 + mrSges(5,2) * t69 - mrSges(4,1);
t36 = mrSges(6,1) * t104 + mrSges(6,2) * t102;
t109 = qJ(5) * t72;
t91 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t52 = Ifges(6,2) * t71 + t123;
t53 = Ifges(5,2) * t71 + t125;
t90 = t53 / 0.2e1 + t52 / 0.2e1;
t54 = Ifges(6,1) * t69 + t122;
t55 = Ifges(5,1) * t69 + t124;
t89 = -t54 / 0.2e1 - t55 / 0.2e1;
t88 = mrSges(6,1) + t129;
t86 = t138 * t70;
t85 = qJD(4) * t115;
t84 = t139 * t106 + t140 * t94;
t83 = -mrSges(5,1) - t88;
t82 = mrSges(5,1) * t69 + mrSges(5,2) * t71;
t81 = Ifges(5,1) * t71 - t125;
t80 = Ifges(6,1) * t71 - t123;
t23 = Ifges(6,5) * t70 + t72 * t80;
t24 = Ifges(5,5) * t70 + t72 * t81;
t77 = -t23 - t24 - t86;
t66 = Ifges(5,5) * t102;
t65 = Ifges(6,5) * t102;
t51 = t115 * t71;
t48 = t115 * t69;
t41 = t81 * qJD(4);
t40 = t80 * qJD(4);
t39 = t79 * qJD(4);
t38 = t78 * qJD(4);
t37 = t82 * qJD(4);
t35 = (pkin(4) * t69 - t73) * t72;
t31 = t82 * t72;
t30 = (mrSges(6,1) * t69 + mrSges(6,2) * t71) * t72;
t26 = -qJD(5) * t69 + t71 * t85;
t25 = t69 * t85 + t100;
t17 = -mrSges(6,2) * t106 + mrSges(6,3) * t75;
t12 = -t109 * t69 + t20;
t11 = -t109 * t71 + t70 * t87 + t33;
t10 = -mrSges(5,1) * t75 - mrSges(5,2) * t76;
t8 = -t55 * t101 + (Ifges(5,5) * t72 - t70 * t81) * qJD(3);
t7 = -t54 * t101 + (Ifges(6,5) * t72 - t70 * t80) * qJD(3);
t6 = -t53 * t101 + (Ifges(5,6) * t72 - t70 * t79) * qJD(3);
t5 = -t52 * t101 + (Ifges(6,6) * t72 - t70 * t78) * qJD(3);
t2 = -qJ(5) * t95 + (-qJD(5) * t72 + (qJ(5) * qJD(3) - qJD(4) * t73) * t70) * t69 + t98;
t13 = [0.2e1 * t1 * t44 + 0.2e1 * t11 * t15 + 0.2e1 * t12 * t17 + 0.2e1 * t14 * t30 + 0.2e1 * t19 * t16 + 0.2e1 * t20 * t18 + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43 + 0.2e1 * t35 * t9 + 0.2e1 * t4 * t45 + 0.2e1 * m(5) * (t19 * t4 + t20 * t3) + (t1 * t11 + t12 * t2 + t14 * t35) * t132 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t130 + (mrSges(4,1) * t130 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * Ifges(4,4) * t70 + t114 * t69 + 0.2e1 * t73 * t31 + t71 * t77) * qJD(3) + t84) * t70 + (mrSges(4,2) * t130 - 0.2e1 * t73 * t10 + (t7 + t8) * t71 + (-t5 - t6) * t69 + ((-t140 * t70 - t114) * t71 + t77 * t69) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + (t138 * t71 - t140 * t69 - 0.2e1 * Ifges(4,4)) * t72 + (-0.2e1 * m(5) * t73 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + t139) * t70) * qJD(3)) * t72; (-t10 + (t113 * t71 + t112 * t69 + m(6) * (-t11 * t69 + t12 * t71) + m(5) * (-t19 * t69 + t20 * t71)) * qJD(3) - t137) * t72 + ((t30 + t31) * qJD(3) + m(6) * (qJD(3) * t35 - t102 * t11 - t104 * t12) + m(5) * (-t104 * t20 - 0.2e1 * t97) + (-t113 * qJD(4) - t135 + t136) * t69 + (m(6) * t2 + t112 * qJD(4) + t133 + t17) * t71) * t70; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t110) * t70 * t106; m(6) * (t1 * t48 + t11 * t26 + t12 * t25 + t14 * t60 - t2 * t51) + t60 * t9 + t25 * t42 + t26 * t44 + t48 * t15 + t14 * t49 - t51 * t17 + t35 * t36 - pkin(3) * t10 + (t65 / 0.2e1 + t66 / 0.2e1 + (-Ifges(4,5) + (-m(5) * pkin(3) + t111) * t73) * qJD(3)) * t70 + (-t4 * mrSges(5,3) - t1 * mrSges(6,3) + t7 / 0.2e1 + t8 / 0.2e1 + t90 * t108 + (-t21 / 0.2e1 - t22 / 0.2e1 - t20 * mrSges(5,3) + pkin(4) * t30 - t12 * mrSges(6,3) - t91 * t70 + t35 * t129) * qJD(4) + ((-m(5) * t20 - t43) * qJD(4) + t136) * pkin(7)) * t69 + (t3 * mrSges(5,3) + t2 * mrSges(6,3) + t5 / 0.2e1 + t6 / 0.2e1 + t89 * t108 + (t23 / 0.2e1 + t24 / 0.2e1 - t11 * mrSges(6,3) - t19 * mrSges(5,3)) * qJD(4) + (-qJD(4) * t45 + t133) * pkin(7)) * t71 + (-t73 * t37 + (t40 / 0.2e1 + t41 / 0.2e1) * t71 + (-t38 / 0.2e1 - t39 / 0.2e1) * t69 + (-t73 * mrSges(4,2) - Ifges(4,6) + t91 * t71 + (Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1) * t69) * qJD(3) + (t69 * t89 - t71 * t90) * qJD(4)) * t72; m(6) * (-pkin(4) * t96 + t51 * t93) + m(5) * (t110 * t126 - t127) * qJD(3) + (m(6) * (-t102 * t48 + t25 * t71 - t26 * t69) + (t111 + t134) * qJD(3)) * t70 + (-t37 - t36 + (-mrSges(4,2) + m(6) * (-t48 * t69 - t51 * t71) + (mrSges(5,3) + mrSges(6,3)) * t110) * qJD(3)) * t72; (-t25 * t51 + t26 * t48) * t132 + 0.2e1 * t60 * t36 - 0.2e1 * pkin(3) * t37 + (t26 * t131 + t40 + t41 + (0.2e1 * t134 * pkin(4) - t51 * t131 - t52 - t53) * qJD(4)) * t69 + (0.2e1 * t25 * mrSges(6,3) + t38 + t39 + (t131 * t48 + t54 + t55) * qJD(4)) * t71; mrSges(5,1) * t4 + mrSges(6,1) * t1 - mrSges(5,2) * t3 - mrSges(6,2) * t2 - t86 * t107 + t135 * pkin(4) + (-t138 * t69 - t140 * t71) * t101 + t84; (t118 * t69 + t71 * t83) * t103 + (-t118 * t71 + t69 * t83) * t106; -mrSges(6,2) * t25 + t65 + t66 + t88 * t26 + ((-mrSges(5,1) * pkin(7) - mrSges(6,3) * pkin(4)) * t71 + (mrSges(5,2) * pkin(7) - t140) * t69) * qJD(4); 0; t137; m(6) * t108; t104 * t129 + t36; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
