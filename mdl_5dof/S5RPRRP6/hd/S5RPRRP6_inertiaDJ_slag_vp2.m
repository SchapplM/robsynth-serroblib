% Calculate time derivative of joint inertia matrix for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:10
% EndTime: 2019-12-31 18:42:14
% DurationCPUTime: 1.50s
% Computational Cost: add. (875->264), mult. (2134->378), div. (0->0), fcn. (1438->6), ass. (0->118)
t117 = Ifges(5,5) + Ifges(6,5);
t138 = -Ifges(5,3) - Ifges(6,3);
t76 = cos(qJ(3));
t104 = qJD(3) * t76;
t73 = sin(qJ(4));
t94 = t73 * t104;
t75 = cos(qJ(4));
t101 = qJD(4) * t75;
t74 = sin(qJ(3));
t95 = t74 * t101;
t77 = t94 + t95;
t87 = -cos(pkin(8)) * pkin(1) - pkin(2);
t137 = 0.2e1 * t87;
t120 = t73 * t76;
t35 = -pkin(3) * t76 - t74 * pkin(7) + t87;
t27 = t75 * t35;
t62 = sin(pkin(8)) * pkin(1) + pkin(6);
t13 = -t120 * t62 + t27;
t105 = qJD(3) * t74;
t127 = pkin(7) * t76;
t128 = pkin(3) * t74;
t48 = (-t127 + t128) * qJD(3);
t114 = t35 * t101 + t73 * t48;
t103 = qJD(4) * t73;
t92 = t76 * t103;
t3 = (-t105 * t75 - t92) * t62 + t114;
t136 = -t13 * qJD(4) + t3;
t116 = Ifges(5,6) + Ifges(6,6);
t118 = t75 * t76;
t47 = t62 * t118;
t14 = t73 * t35 + t47;
t108 = t73 ^ 2 + t75 ^ 2;
t102 = qJD(4) * t74;
t96 = t73 * t102;
t93 = t75 * t104;
t98 = -t77 * mrSges(6,1) - mrSges(6,2) * t93;
t11 = -mrSges(6,2) * t96 - t98;
t16 = pkin(4) * t77 + t104 * t62;
t135 = m(6) * t16 + t11;
t134 = 2 * m(5);
t133 = 0.2e1 * m(6);
t132 = -2 * mrSges(6,3);
t131 = 0.2e1 * t62;
t130 = m(6) * pkin(4);
t126 = mrSges(5,2) * t75;
t125 = Ifges(5,4) * t73;
t124 = Ifges(5,4) * t75;
t123 = Ifges(6,4) * t73;
t122 = Ifges(6,4) * t75;
t121 = t73 * t74;
t119 = t74 * t75;
t115 = -qJ(5) - pkin(7);
t82 = Ifges(6,1) * t75 - t123;
t24 = -Ifges(6,5) * t76 + t74 * t82;
t83 = Ifges(5,1) * t75 - t125;
t25 = -Ifges(5,5) * t76 + t74 * t83;
t113 = t24 + t25;
t97 = t62 * t105;
t112 = t75 * t48 + t73 * t97;
t43 = mrSges(6,2) * t76 - mrSges(6,3) * t121;
t44 = mrSges(5,2) * t76 - mrSges(5,3) * t121;
t111 = t43 + t44;
t45 = -mrSges(6,1) * t76 - mrSges(6,3) * t119;
t46 = -mrSges(5,1) * t76 - mrSges(5,3) * t119;
t110 = -t45 - t46;
t109 = -mrSges(5,1) * t75 + mrSges(5,2) * t73 - mrSges(4,1);
t37 = mrSges(6,1) * t103 + mrSges(6,2) * t101;
t107 = qJ(5) * t74;
t106 = qJ(5) * t75;
t100 = qJD(5) * t75;
t91 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t53 = Ifges(6,2) * t75 + t123;
t54 = Ifges(5,2) * t75 + t125;
t90 = -t53 / 0.2e1 - t54 / 0.2e1;
t55 = Ifges(6,1) * t73 + t122;
t56 = Ifges(5,1) * t73 + t124;
t89 = t55 / 0.2e1 + t56 / 0.2e1;
t88 = -mrSges(5,1) - t130;
t86 = qJD(4) * t115;
t85 = t138 * t105 - t117 * t93;
t84 = mrSges(5,1) * t73 + t126;
t81 = -Ifges(5,2) * t73 + t124;
t80 = -Ifges(6,2) * t73 + t122;
t22 = -t76 * Ifges(6,6) + t74 * t80;
t23 = -t76 * Ifges(5,6) + t74 * t81;
t79 = t116 * t76 - t22 - t23;
t78 = t93 - t96;
t69 = Ifges(5,5) * t101;
t68 = Ifges(6,5) * t101;
t63 = -pkin(4) * t75 - pkin(3);
t52 = t115 * t75;
t50 = -mrSges(6,1) * t75 + mrSges(6,2) * t73;
t49 = t115 * t73;
t42 = t83 * qJD(4);
t41 = t82 * qJD(4);
t40 = t81 * qJD(4);
t39 = t80 * qJD(4);
t38 = t84 * qJD(4);
t34 = t84 * t74;
t33 = (mrSges(6,1) * t73 + mrSges(6,2) * t75) * t74;
t30 = (pkin(4) * t73 + t62) * t74;
t29 = -qJD(5) * t73 + t75 * t86;
t28 = t73 * t86 + t100;
t20 = -mrSges(5,2) * t105 - mrSges(5,3) * t77;
t19 = -mrSges(6,2) * t105 - mrSges(6,3) * t77;
t18 = mrSges(5,1) * t105 - mrSges(5,3) * t78;
t17 = mrSges(6,1) * t105 - mrSges(6,3) * t78;
t12 = mrSges(5,1) * t77 + mrSges(5,2) * t78;
t10 = -t107 * t73 + t14;
t9 = -t56 * t102 + (Ifges(5,5) * t74 + t76 * t83) * qJD(3);
t8 = -t55 * t102 + (Ifges(6,5) * t74 + t76 * t82) * qJD(3);
t7 = -t54 * t102 + (Ifges(5,6) * t74 + t76 * t81) * qJD(3);
t6 = -t53 * t102 + (Ifges(6,6) * t74 + t76 * t80) * qJD(3);
t5 = -t74 * t106 + t27 + (-t62 * t73 - pkin(4)) * t76;
t4 = -t14 * qJD(4) + t112;
t2 = (-qJ(5) * qJD(4) - qJD(3) * t62) * t119 + (-qJD(5) * t74 + (-qJ(5) * qJD(3) - qJD(4) * t62) * t76) * t73 + t114;
t1 = -t74 * t100 + (pkin(4) * t74 - t106 * t76) * qJD(3) + (-t47 + (-t35 + t107) * t73) * qJD(4) + t112;
t15 = [0.2e1 * t1 * t45 + 0.2e1 * t10 * t19 + 0.2e1 * t30 * t11 + 0.2e1 * t13 * t18 + 0.2e1 * t14 * t20 + 0.2e1 * t16 * t33 + 0.2e1 * t5 * t17 + 0.2e1 * t2 * t43 + 0.2e1 * t3 * t44 + 0.2e1 * t4 * t46 + (t13 * t4 + t14 * t3) * t134 + (t1 * t5 + t10 * t2 + t16 * t30) * t133 + ((mrSges(4,2) * t137 + 0.2e1 * Ifges(4,4) * t76 + t113 * t75 + t131 * t34 + t73 * t79) * qJD(3) + t85) * t76 + (t12 * t131 + (t8 + t9) * t75 + (-t6 - t7) * t73 + (t79 * t75 + (t117 * t76 - t113) * t73) * qJD(4) + (mrSges(4,1) * t137 + (-t116 * t73 + t117 * t75 - 0.2e1 * Ifges(4,4)) * t74 + (t134 * t62 ^ 2 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t138) * t76) * qJD(3)) * t74; (-t12 + (t111 * t75 + t110 * t73 + m(6) * (t10 * t75 - t5 * t73) + (-t13 * t73 + t14 * t75 - t62 * t76) * m(5)) * qJD(3) - t135) * t76 + ((t19 + t20) * t75 + (-t17 - t18) * t73 + (t33 + t34) * qJD(3) + (t110 * t75 - t111 * t73) * qJD(4) + m(6) * (qJD(3) * t30 - t1 * t73 - t10 * t103 - t101 * t5 + t2 * t75) + (-t103 * t14 + t136 * t75 - t4 * t73 + t97) * m(5)) * t74; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t108) * t74 * t104; m(6) * (t1 * t49 + t10 * t28 + t16 * t63 - t2 * t52 + t29 * t5) + t63 * t11 + t49 * t17 + t16 * t50 - t52 * t19 + t30 * t37 + t28 * t43 + t29 * t45 - pkin(3) * t12 + (-t68 / 0.2e1 - t69 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t109) * t62) * qJD(3)) * t76 + (-t1 * mrSges(6,3) - t4 * mrSges(5,3) + t8 / 0.2e1 + t9 / 0.2e1 + t90 * t104 + (-t22 / 0.2e1 - t23 / 0.2e1 + pkin(4) * t33 - t10 * mrSges(6,3) - t14 * mrSges(5,3) + t91 * t76 + t30 * t130) * qJD(4) + (-m(5) * t4 - t18 + (-m(5) * t14 - t44) * qJD(4)) * pkin(7)) * t73 + (t2 * mrSges(6,3) + t3 * mrSges(5,3) + t6 / 0.2e1 + t7 / 0.2e1 + t89 * t104 + (t24 / 0.2e1 + t25 / 0.2e1 - t5 * mrSges(6,3) - t13 * mrSges(5,3)) * qJD(4) + (m(5) * t136 - qJD(4) * t46 + t20) * pkin(7)) * t75 + (t62 * t38 + (t41 / 0.2e1 + t42 / 0.2e1) * t75 + (-t39 / 0.2e1 - t40 / 0.2e1) * t73 + (t62 * mrSges(4,2) - Ifges(4,6) + t91 * t75 + (Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1) * t73) * qJD(3) + (-t73 * t89 + t75 * t90) * qJD(4)) * t74; m(6) * (-pkin(4) * t92 + t28 * t119 - t29 * t121 - t49 * t95 + t52 * t96) + ((t50 + t109) * t74 + m(5) * (t108 * t127 - t128) + m(6) * (-t118 * t52 - t120 * t49 + t63 * t74)) * qJD(3) + (-t38 - t37 + (-mrSges(4,2) + (mrSges(5,3) + mrSges(6,3)) * t108) * qJD(3)) * t76; -0.2e1 * pkin(3) * t38 + 0.2e1 * t63 * t37 + (-t28 * t52 + t29 * t49) * t133 + (t29 * t132 + t41 + t42 + (-t52 * t132 - t53 - t54 + 0.2e1 * (m(6) * t63 + t50) * pkin(4)) * qJD(4)) * t73 + (0.2e1 * t28 * mrSges(6,3) + t39 + t40 + (t132 * t49 + t55 + t56) * qJD(4)) * t75; t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) - t116 * t94 + (m(6) * t1 + t17) * pkin(4) + (-t116 * t75 - t117 * t73) * t102 - t85; (t73 * t88 - t126) * t104 + (t88 * t75 + (mrSges(5,2) + mrSges(6,2)) * t73) * t102 + t98; -mrSges(6,2) * t28 + t68 + t69 + (mrSges(6,1) + t130) * t29 + ((-mrSges(5,1) * pkin(7) - mrSges(6,3) * pkin(4)) * t75 + (mrSges(5,2) * pkin(7) - t116) * t73) * qJD(4); 0; t135; m(6) * t105; t103 * t130 + t37; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
