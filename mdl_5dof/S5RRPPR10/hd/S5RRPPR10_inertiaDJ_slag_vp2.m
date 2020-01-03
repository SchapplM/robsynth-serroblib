% Calculate time derivative of joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:43:08
% DurationCPUTime: 1.62s
% Computational Cost: add. (1155->278), mult. (2839->419), div. (0->0), fcn. (2200->6), ass. (0->118)
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t138 = qJD(3) * (t91 ^ 2 + t92 ^ 2);
t137 = Ifges(5,1) + Ifges(4,1);
t128 = pkin(3) + pkin(4);
t136 = t128 * t91;
t94 = sin(qJ(2));
t110 = qJD(2) * t94;
t96 = cos(qJ(2));
t135 = qJ(4) * t110 - qJD(4) * t96;
t134 = 2 * m(4);
t133 = 2 * m(5);
t132 = 2 * m(6);
t131 = -2 * pkin(1);
t130 = 2 * pkin(6);
t109 = qJD(2) * t96;
t103 = t92 * t109;
t117 = t92 * t94;
t112 = -qJ(4) * t103 - qJD(4) * t117;
t126 = pkin(3) * t91;
t29 = (pkin(6) + t126) * t109 + t112;
t127 = m(5) * t29;
t125 = mrSges(5,3) * t92;
t124 = Ifges(4,4) * t91;
t123 = Ifges(4,4) * t92;
t122 = Ifges(5,5) * t91;
t121 = Ifges(5,5) * t92;
t120 = t91 * t94;
t119 = t91 * t96;
t56 = -t94 * qJD(3) + (pkin(2) * t94 - qJ(3) * t96) * qJD(2);
t118 = t92 * t56;
t116 = t92 * t96;
t115 = -pkin(7) + qJ(3);
t93 = sin(qJ(5));
t95 = cos(qJ(5));
t98 = t91 * t93 + t92 * t95;
t57 = t98 * qJD(5);
t67 = t91 * t95 - t92 * t93;
t18 = t67 * t109 - t94 * t57;
t50 = t67 * t94;
t19 = qJD(5) * t50 + t98 * t109;
t114 = Ifges(6,5) * t19 + Ifges(6,6) * t18;
t58 = t67 * qJD(5);
t113 = -Ifges(6,5) * t57 - Ifges(6,6) * t58;
t104 = t91 * t109;
t53 = mrSges(4,1) * t104 + mrSges(4,2) * t103;
t73 = -pkin(2) * t96 - t94 * qJ(3) - pkin(1);
t47 = pkin(6) * t116 + t91 * t73;
t111 = qJ(3) * t138;
t108 = qJD(4) * t91;
t105 = pkin(6) * t110;
t102 = -pkin(6) * t91 - pkin(3);
t101 = qJ(4) * t91 + pkin(2);
t100 = qJ(4) * t92 - pkin(6);
t83 = pkin(6) * t119;
t46 = t92 * t73 - t83;
t38 = -qJ(4) * t96 + t47;
t63 = -mrSges(5,1) * t110 + mrSges(5,2) * t103;
t5 = -t18 * mrSges(6,1) + t19 * mrSges(6,2);
t24 = mrSges(6,1) * t58 - mrSges(6,2) * t57;
t88 = t96 * pkin(3);
t23 = pkin(4) * t96 + t83 + t88 + (-pkin(7) * t94 - t73) * t92;
t28 = pkin(7) * t120 + t38;
t6 = t23 * t95 - t28 * t93;
t7 = t23 * t93 + t28 * t95;
t74 = t115 * t91;
t76 = t115 * t92;
t35 = t74 * t95 - t76 * t93;
t36 = t74 * t93 + t76 * t95;
t48 = t91 * t56;
t34 = -t92 * t105 + t48;
t97 = (-Ifges(5,4) - Ifges(4,5)) * t92 + (Ifges(4,6) - Ifges(5,6)) * t91;
t78 = mrSges(5,1) * t104;
t75 = -mrSges(5,1) * t92 - mrSges(5,3) * t91;
t72 = -pkin(3) * t92 - t101;
t71 = -mrSges(5,2) * t120 - mrSges(5,3) * t96;
t70 = mrSges(5,1) * t96 + mrSges(5,2) * t117;
t69 = -mrSges(4,1) * t96 - mrSges(4,3) * t117;
t68 = mrSges(4,2) * t96 - mrSges(4,3) * t120;
t64 = (-mrSges(5,2) * t119 + mrSges(5,3) * t94) * qJD(2);
t62 = (mrSges(4,1) * t94 - mrSges(4,3) * t116) * qJD(2);
t61 = (-mrSges(4,2) * t94 - mrSges(4,3) * t119) * qJD(2);
t60 = (mrSges(5,1) * t91 - t125) * t94;
t59 = t128 * t92 + t101;
t52 = -mrSges(5,3) * t103 + t78;
t51 = t98 * t94;
t49 = (-t100 + t126) * t94;
t45 = (Ifges(4,5) * t94 + (Ifges(4,1) * t92 - t124) * t96) * qJD(2);
t44 = (Ifges(5,4) * t94 + (Ifges(5,1) * t92 + t122) * t96) * qJD(2);
t43 = (Ifges(4,6) * t94 + (-Ifges(4,2) * t91 + t123) * t96) * qJD(2);
t42 = (Ifges(5,6) * t94 + (Ifges(5,3) * t91 + t121) * t96) * qJD(2);
t41 = mrSges(6,1) * t96 - t51 * mrSges(6,3);
t40 = -mrSges(6,2) * t96 + t50 * mrSges(6,3);
t39 = -t46 + t88;
t37 = (t100 - t136) * t94;
t33 = t91 * t105 + t118;
t32 = Ifges(6,1) * t67 - Ifges(6,4) * t98;
t31 = Ifges(6,4) * t67 - Ifges(6,2) * t98;
t30 = mrSges(6,1) * t98 + mrSges(6,2) * t67;
t27 = t102 * t110 - t118;
t26 = -Ifges(6,1) * t57 - Ifges(6,4) * t58;
t25 = -Ifges(6,4) * t57 - Ifges(6,2) * t58;
t22 = (pkin(6) + t136) * t109 + t112;
t21 = t34 + t135;
t20 = -mrSges(6,1) * t50 + mrSges(6,2) * t51;
t15 = Ifges(6,1) * t51 + Ifges(6,4) * t50 + Ifges(6,5) * t96;
t14 = Ifges(6,4) * t51 + Ifges(6,2) * t50 + Ifges(6,6) * t96;
t13 = t67 * qJD(3) - t36 * qJD(5);
t12 = t98 * qJD(3) + t35 * qJD(5);
t11 = t48 + (-pkin(6) * t117 + pkin(7) * t119) * qJD(2) + t135;
t10 = -t118 + (-pkin(7) * t116 + (-pkin(4) + t102) * t94) * qJD(2);
t9 = -mrSges(6,1) * t110 - mrSges(6,3) * t19;
t8 = mrSges(6,2) * t110 + mrSges(6,3) * t18;
t4 = Ifges(6,1) * t19 + Ifges(6,4) * t18 - Ifges(6,5) * t110;
t3 = Ifges(6,4) * t19 + Ifges(6,2) * t18 - Ifges(6,6) * t110;
t2 = -t7 * qJD(5) + t10 * t95 - t11 * t93;
t1 = t6 * qJD(5) + t10 * t93 + t11 * t95;
t16 = [0.2e1 * t29 * t60 + 0.2e1 * t47 * t61 + 0.2e1 * t46 * t62 + 0.2e1 * t39 * t63 + 0.2e1 * t38 * t64 + 0.2e1 * t34 * t68 + 0.2e1 * t33 * t69 + 0.2e1 * t27 * t70 + 0.2e1 * t21 * t71 + t50 * t3 + t51 * t4 + 0.2e1 * t49 * t52 + 0.2e1 * t37 * t5 + 0.2e1 * t1 * t40 + 0.2e1 * t2 * t41 + t19 * t15 - 0.2e1 * t22 * t20 + 0.2e1 * t6 * t9 + t18 * t14 + 0.2e1 * t7 * t8 + (((mrSges(3,1) * t131) - Ifges(6,5) * t51 - Ifges(6,6) * t50 + (-(2 * Ifges(3,4)) - t97) * t94) * t94 + ((mrSges(3,2) * t131) + 0.2e1 * (Ifges(3,4) + t97) * t96 + (-(2 * Ifges(6,3)) + (pkin(6) ^ 2 * t134) + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(5,2)) - (2 * Ifges(4,3)) + (mrSges(4,2) * t130 + t137 * t92) * t92 + (mrSges(4,1) * t130 + (Ifges(5,3) + Ifges(4,2)) * t91 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t92) * t91) * t94) * t96) * qJD(2) + (t53 * t130 + (t44 + t45) * t92 + (t42 - t43) * t91) * t94 + (t1 * t7 + t2 * t6 - t22 * t37) * t132 + (t21 * t38 + t27 * t39 + t29 * t49) * t133 + (t46 * t33 + t47 * t34) * t134 + t96 * t114; (-t1 * t98 - t2 * t67 + t6 * t57 - t7 * t58) * mrSges(6,3) + m(6) * (t1 * t36 + t12 * t7 + t13 * t6 + t2 * t35 - t22 * t59) + (t34 * mrSges(4,3) + t21 * mrSges(5,2) - t42 / 0.2e1 + t43 / 0.2e1 + (t68 + t71) * qJD(3) + (t61 + t64) * qJ(3) + m(5) * (qJ(3) * t21 + qJD(3) * t38) + m(4) * (qJ(3) * t34 + qJD(3) * t47)) * t92 + t29 * t75 - t98 * t3 / 0.2e1 + t67 * t4 / 0.2e1 + t50 * t25 / 0.2e1 + t51 * t26 / 0.2e1 - pkin(2) * t53 - t57 * t15 / 0.2e1 - t58 * t14 / 0.2e1 + t59 * t5 + t18 * t31 / 0.2e1 + t19 * t32 / 0.2e1 + t35 * t9 + t36 * t8 + t37 * t24 + t12 * t40 + t13 * t41 - t22 * t30 + ((-Ifges(6,5) * t67 / 0.2e1 + Ifges(6,6) * t98 / 0.2e1 - Ifges(3,6) + (pkin(6) * mrSges(3,2)) + (-Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1) * t92 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t91) * t94 + (Ifges(3,5) - t91 * (Ifges(4,2) * t92 + t124) / 0.2e1 + t91 * (-Ifges(5,3) * t92 + t122) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t92 + mrSges(4,2) * t91 - mrSges(3,1)) * pkin(6) + (t137 * t91 - t121 + t123) * t92 / 0.2e1) * t96) * qJD(2) + (t27 * mrSges(5,2) - t33 * mrSges(4,3) + t44 / 0.2e1 + t45 / 0.2e1 + (-t69 + t70) * qJD(3) + (-t62 + t63) * qJ(3) + m(5) * (qJ(3) * t27 + qJD(3) * t39) + m(4) * (-qJ(3) * t33 - qJD(3) * t46) + (-m(5) * t49 + m(6) * t37 + t20 - t60) * qJD(4)) * t91 + t96 * t113 / 0.2e1 + (t52 + t127) * t72; 0.2e1 * t59 * t24 - t98 * t25 + t67 * t26 - t58 * t31 - t57 * t32 + (t59 * t108 + t12 * t36 + t13 * t35) * t132 + t111 * t134 + (-t72 * t108 + t111) * t133 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t138 + 0.2e1 * (t30 - t75) * t108 + 0.2e1 * (-t12 * t98 - t13 * t67 + t35 * t57 - t36 * t58) * mrSges(6,3); t78 + ((m(4) * pkin(6)) - t125) * t109 + t127 + m(6) * t22 - t5 + t53; (-m(5) - m(6)) * t108 - t24; 0; t93 * t8 + t95 * t9 + (t95 * t40 - t93 * t41) * qJD(5) + m(6) * (t1 * t93 + t2 * t95 + (-t6 * t93 + t7 * t95) * qJD(5)) + m(5) * t27 + t63; m(6) * (t12 * t93 + t13 * t95 + (-t35 * t93 + t36 * t95) * qJD(5)) + m(5) * t91 * qJD(3) + (t95 * t57 - t93 * t58 + (t67 * t93 - t95 * t98) * qJD(5)) * mrSges(6,3); 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,3) * t110 + t114; mrSges(6,1) * t13 - mrSges(6,2) * t12 + t113; 0; (-mrSges(6,1) * t93 - mrSges(6,2) * t95) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
