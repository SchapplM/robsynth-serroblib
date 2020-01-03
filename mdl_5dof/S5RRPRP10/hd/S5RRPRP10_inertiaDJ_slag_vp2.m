% Calculate time derivative of joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:26
% EndTime: 2019-12-31 20:09:31
% DurationCPUTime: 1.64s
% Computational Cost: add. (911->248), mult. (2062->347), div. (0->0), fcn. (1336->4), ass. (0->110)
t136 = Ifges(5,5) + Ifges(6,5);
t78 = sin(qJ(2));
t138 = t136 * t78;
t137 = Ifges(5,6) + Ifges(6,6);
t123 = pkin(3) + pkin(6);
t135 = Ifges(5,3) + Ifges(6,3);
t80 = cos(qJ(2));
t79 = cos(qJ(4));
t119 = Ifges(6,4) * t79;
t77 = sin(qJ(4));
t90 = Ifges(6,1) * t77 + t119;
t121 = Ifges(5,4) * t79;
t91 = Ifges(5,1) * t77 + t121;
t114 = (-t90 - t91) * t80 + t138;
t134 = t114 * t77;
t120 = Ifges(6,4) * t77;
t88 = Ifges(6,2) * t79 + t120;
t122 = Ifges(5,4) * t77;
t89 = Ifges(5,2) * t79 + t122;
t115 = (-t88 - t89) * t80 + t137 * t78;
t133 = t115 * t79;
t107 = qJD(4) * t80;
t105 = t77 * t107;
t111 = qJD(2) * t78;
t83 = t79 * t111 + t105;
t103 = t77 * t111;
t104 = t79 * t107;
t82 = t103 - t104;
t129 = -2 * mrSges(6,3);
t81 = -pkin(2) - pkin(7);
t96 = -qJ(3) * t78 - pkin(1);
t35 = t81 * t80 + t96;
t60 = t123 * t78;
t14 = t79 * t35 + t77 * t60;
t132 = -t136 * t77 - t137 * t79;
t131 = 2 * m(6);
t130 = -0.2e1 * pkin(1);
t110 = qJD(2) * t80;
t94 = pkin(2) * t111 - qJD(3) * t78;
t29 = -qJ(3) * t110 + t94;
t128 = 0.2e1 * t29;
t53 = -pkin(2) * t80 + t96;
t127 = -0.2e1 * t53;
t118 = t77 * t80;
t117 = t79 * t80;
t61 = t123 * t80;
t113 = qJ(5) * t80;
t112 = qJ(5) - t81;
t109 = qJD(4) * t77;
t108 = qJD(4) * t79;
t106 = t79 * qJD(5);
t56 = -Ifges(6,2) * t77 + t119;
t57 = -Ifges(5,2) * t77 + t121;
t101 = t56 / 0.2e1 + t57 / 0.2e1;
t58 = Ifges(6,1) * t79 - t120;
t59 = Ifges(5,1) * t79 - t122;
t100 = t58 / 0.2e1 + t59 / 0.2e1;
t99 = m(4) * pkin(6) + mrSges(4,1);
t98 = m(6) * pkin(4) + mrSges(6,1);
t97 = m(5) * t81 - mrSges(5,3);
t52 = t112 * t79;
t95 = -t35 + t113;
t39 = mrSges(6,1) * t108 - mrSges(6,2) * t109;
t10 = -t79 * t113 + t14;
t38 = t79 * t60;
t5 = pkin(4) * t78 + t95 * t77 + t38;
t93 = t10 * t79 - t5 * t77;
t92 = mrSges(5,1) * t79 - mrSges(5,2) * t77;
t13 = -t35 * t77 + t38;
t85 = t13 * t77 - t14 * t79;
t21 = (pkin(7) * t78 - qJ(3) * t80) * qJD(2) + t94;
t50 = t123 * t110;
t3 = t60 * t108 - t35 * t109 + t79 * t21 + t77 * t50;
t84 = t136 * t103 + t135 * t110 + t137 * t83;
t11 = -t83 * mrSges(6,1) + t82 * mrSges(6,2);
t70 = pkin(4) * t77 + qJ(3);
t69 = pkin(4) * t108 + qJD(3);
t55 = mrSges(5,1) * t77 + mrSges(5,2) * t79;
t54 = mrSges(6,1) * t77 + mrSges(6,2) * t79;
t51 = t112 * t77;
t49 = t123 * t111;
t48 = -mrSges(5,2) * t78 - mrSges(5,3) * t117;
t47 = -mrSges(6,2) * t78 - mrSges(6,3) * t117;
t46 = mrSges(5,1) * t78 + mrSges(5,3) * t118;
t45 = mrSges(6,1) * t78 + mrSges(6,3) * t118;
t44 = t91 * qJD(4);
t43 = t90 * qJD(4);
t42 = t89 * qJD(4);
t41 = t88 * qJD(4);
t40 = t92 * qJD(4);
t34 = t92 * t80;
t33 = (mrSges(6,1) * t79 - mrSges(6,2) * t77) * t80;
t32 = t79 * t50;
t30 = pkin(4) * t117 + t61;
t27 = -qJD(4) * t52 - t77 * qJD(5);
t26 = t112 * t109 - t106;
t20 = mrSges(5,1) * t110 - t82 * mrSges(5,3);
t19 = mrSges(6,1) * t110 - t82 * mrSges(6,3);
t18 = -mrSges(5,2) * t110 + t83 * mrSges(5,3);
t17 = -mrSges(6,2) * t110 + t83 * mrSges(6,3);
t15 = -pkin(4) * t105 + (-pkin(4) * t79 - t123) * t111;
t12 = -t83 * mrSges(5,1) + t82 * mrSges(5,2);
t9 = -t59 * t107 + (t80 * Ifges(5,5) + t91 * t78) * qJD(2);
t8 = -t58 * t107 + (t80 * Ifges(6,5) + t90 * t78) * qJD(2);
t7 = -t57 * t107 + (t80 * Ifges(5,6) + t89 * t78) * qJD(2);
t6 = -t56 * t107 + (t80 * Ifges(6,6) + t88 * t78) * qJD(2);
t4 = -t14 * qJD(4) - t21 * t77 + t32;
t2 = t83 * qJ(5) - t80 * t106 + t3;
t1 = pkin(4) * t110 + t32 + t95 * t108 + (-qJ(5) * t111 - qJD(4) * t60 + qJD(5) * t80 - t21) * t77;
t16 = [m(4) * t53 * t128 + 0.2e1 * t1 * t45 + 0.2e1 * t10 * t17 + 0.2e1 * t30 * t11 + 0.2e1 * t61 * t12 + 0.2e1 * t13 * t20 + 0.2e1 * t14 * t18 + 0.2e1 * t15 * t33 + 0.2e1 * t5 * t19 + 0.2e1 * t2 * t47 + 0.2e1 * t3 * t48 - 0.2e1 * t49 * t34 + 0.2e1 * t4 * t46 + 0.2e1 * m(5) * (t13 * t4 + t14 * t3 - t49 * t61) + (t1 * t5 + t10 * t2 + t15 * t30) * t131 + (-0.2e1 * t29 * mrSges(4,3) + (mrSges(3,1) * t130 + mrSges(4,2) * t127 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6)) * t78 + t133 + t134) * qJD(2) + t84) * t78 + (mrSges(4,2) * t128 + (-t6 - t7) * t79 + (-t8 - t9) * t77 + (t115 * t77 + (-t114 - t138) * t79) * qJD(4) + (mrSges(3,2) * t130 + mrSges(4,3) * t127 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t132) * t80 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + t135) * t78) * qJD(2)) * t80; qJ(3) * t12 + qJD(3) * t34 + t70 * t11 + t15 * t54 - t51 * t17 - t52 * t19 + t26 * t45 + t27 * t47 + t30 * t39 + t69 * t33 + t61 * t40 - t49 * t55 + m(5) * (-qJ(3) * t49 + qJD(3) * t61) + m(6) * (-t52 * t1 + t27 * t10 + t15 * t70 - t51 * t2 + t26 * t5 + t30 * t69) + (t81 * t20 - t1 * mrSges(6,3) + t8 / 0.2e1 + t9 / 0.2e1 + t97 * t4) * t79 + (t81 * t18 - t2 * mrSges(6,3) - t6 / 0.2e1 - t7 / 0.2e1 + t97 * t3) * t77 + ((t41 / 0.2e1 + t42 / 0.2e1) * t79 + (t43 / 0.2e1 + t44 / 0.2e1) * t77 + t99 * qJD(3)) * t80 + ((-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) + (Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t79 + (-Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t77 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(6)) * t80 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t101 * t79 + t100 * t77 + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(6)) * t78) * qJD(2) + (-t93 * mrSges(6,3) + t85 * mrSges(5,3) + (-m(5) * t85 - t77 * t46 + t79 * t48) * t81 + (-t100 * t79 + t101 * t77) * t80 - t134 / 0.2e1 + t132 * t78 / 0.2e1 - t133 / 0.2e1) * qJD(4); 0.2e1 * qJ(3) * t40 + 0.2e1 * t69 * t54 + 0.2e1 * t70 * t39 + (-t52 * t26 - t51 * t27 + t69 * t70) * t131 + (t26 * t129 - t43 - t44) * t79 + (t27 * t129 + t41 + t42) * t77 + 0.2e1 * (mrSges(4,3) + t55 + (m(4) + m(5)) * qJ(3)) * qJD(3) + ((-t51 * t129 - t56 - t57) * t79 + (t52 * t129 - t58 - t59) * t77) * qJD(4); (t19 + t20) * t79 + (t17 + t18) * t77 + t99 * t110 + ((t47 + t48) * t79 + (-t45 - t46) * t77) * qJD(4) + m(6) * (t93 * qJD(4) + t1 * t79 + t2 * t77) + m(5) * (-t85 * qJD(4) + t3 * t77 + t4 * t79); m(6) * (t26 * t79 + t27 * t77 + (-t51 * t79 + t52 * t77) * qJD(4)); 0; t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) - t136 * t104 + (m(6) * t1 + t19) * pkin(4) + t84; -t27 * mrSges(6,2) + t98 * t26 + ((-mrSges(5,2) * t81 - t137) * t79 + (-mrSges(5,1) * t81 + (mrSges(6,3) * pkin(4)) - t136) * t77) * qJD(4); ((-mrSges(5,2) - mrSges(6,2)) * t79 + (-mrSges(5,1) - t98) * t77) * qJD(4); 0; m(6) * t15 + t11; m(6) * t69 + t39; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
