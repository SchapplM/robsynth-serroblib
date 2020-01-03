% Calculate time derivative of joint inertia matrix for
% S5RRPRP6
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:47
% EndTime: 2019-12-31 19:56:51
% DurationCPUTime: 1.60s
% Computational Cost: add. (1634->261), mult. (3760->374), div. (0->0), fcn. (3232->6), ass. (0->120)
t125 = Ifges(5,5) + Ifges(6,5);
t145 = Ifges(5,3) + Ifges(6,3);
t88 = cos(qJ(4));
t117 = qJD(4) * t88;
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t87 = sin(qJ(2));
t89 = cos(qJ(2));
t60 = t84 * t89 + t85 * t87;
t112 = t60 * t117;
t59 = t84 * t87 - t85 * t89;
t56 = t59 * qJD(2);
t86 = sin(qJ(4));
t93 = -t86 * t56 + t112;
t118 = qJD(4) * t86;
t126 = t88 * t56;
t92 = t60 * t118 + t126;
t79 = -pkin(2) * t89 - pkin(1);
t35 = pkin(3) * t59 - pkin(7) * t60 + t79;
t123 = -qJ(3) - pkin(6);
t70 = t123 * t87;
t71 = t123 * t89;
t42 = t70 * t84 - t71 * t85;
t40 = t88 * t42;
t15 = t86 * t35 + t40;
t77 = pkin(2) * t84 + pkin(7);
t108 = m(5) * t77 + mrSges(5,3);
t124 = -Ifges(5,6) - Ifges(6,6);
t144 = (t124 * t88 - t125 * t86) * qJD(4);
t143 = 2 * m(6);
t142 = -2 * mrSges(4,3);
t141 = -2 * mrSges(6,3);
t41 = -t70 * t85 - t71 * t84;
t139 = 0.2e1 * t41;
t138 = m(6) * pkin(4);
t136 = mrSges(5,2) * t88;
t135 = Ifges(5,4) * t86;
t134 = Ifges(5,4) * t88;
t133 = Ifges(6,4) * t86;
t132 = Ifges(6,4) * t88;
t103 = qJD(2) * t123;
t54 = qJD(3) * t89 + t103 * t87;
t91 = -t87 * qJD(3) + t103 * t89;
t28 = t54 * t84 - t85 * t91;
t131 = t28 * t41;
t130 = t60 * t86;
t129 = t60 * t88;
t95 = -Ifges(6,2) * t86 + t132;
t20 = t59 * Ifges(6,6) + t60 * t95;
t96 = -Ifges(5,2) * t86 + t134;
t21 = t59 * Ifges(5,6) + t60 * t96;
t122 = -t20 - t21;
t97 = Ifges(6,1) * t88 - t133;
t22 = Ifges(6,5) * t59 + t60 * t97;
t98 = Ifges(5,1) * t88 - t135;
t23 = Ifges(5,5) * t59 + t60 * t98;
t121 = t22 + t23;
t62 = mrSges(6,1) * t118 + mrSges(6,2) * t117;
t120 = qJ(5) * t60;
t119 = qJ(5) + t77;
t29 = t85 * t54 + t84 * t91;
t114 = pkin(2) * qJD(2) * t87;
t55 = t60 * qJD(2);
t30 = pkin(3) * t55 + pkin(7) * t56 + t114;
t116 = t117 * t35 + t29 * t88 + t30 * t86;
t111 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t72 = Ifges(6,2) * t88 + t133;
t73 = Ifges(5,2) * t88 + t135;
t110 = -t72 / 0.2e1 - t73 / 0.2e1;
t74 = Ifges(6,1) * t86 + t132;
t75 = Ifges(5,1) * t86 + t134;
t109 = t74 / 0.2e1 + t75 / 0.2e1;
t78 = -pkin(2) * t85 - pkin(3);
t107 = t124 * t86;
t106 = t55 * mrSges(4,1) - mrSges(4,2) * t56;
t105 = -t29 * t86 + t30 * t88;
t14 = t35 * t88 - t42 * t86;
t104 = 0.2e1 * t114;
t102 = -t125 * t126 + t145 * t55;
t101 = qJD(4) * t119;
t100 = -(2 * Ifges(4,4)) + t107;
t99 = mrSges(5,1) * t86 + t136;
t94 = qJ(5) * t56 - qJD(5) * t60;
t12 = mrSges(6,1) * t93 - mrSges(6,2) * t92;
t83 = Ifges(5,5) * t117;
t82 = Ifges(6,5) * t117;
t69 = -mrSges(6,1) * t88 + mrSges(6,2) * t86;
t68 = -pkin(4) * t88 + t78;
t67 = t98 * qJD(4);
t66 = t97 * qJD(4);
t65 = t96 * qJD(4);
t64 = t95 * qJD(4);
t63 = t99 * qJD(4);
t58 = t119 * t88;
t57 = t119 * t86;
t44 = -qJD(5) * t86 - t101 * t88;
t43 = qJD(5) * t88 - t101 * t86;
t39 = mrSges(5,1) * t59 - mrSges(5,3) * t129;
t38 = mrSges(6,1) * t59 - mrSges(6,3) * t129;
t37 = -mrSges(5,2) * t59 - mrSges(5,3) * t130;
t36 = -mrSges(6,2) * t59 - mrSges(6,3) * t130;
t34 = (mrSges(6,1) * t86 + mrSges(6,2) * t88) * t60;
t27 = pkin(4) * t130 + t41;
t19 = -mrSges(5,2) * t55 - mrSges(5,3) * t93;
t18 = -mrSges(6,2) * t55 - mrSges(6,3) * t93;
t17 = mrSges(5,1) * t55 + mrSges(5,3) * t92;
t16 = mrSges(6,1) * t55 + mrSges(6,3) * t92;
t13 = mrSges(5,1) * t93 - mrSges(5,2) * t92;
t11 = pkin(4) * t93 + t28;
t10 = -t120 * t86 + t15;
t9 = -Ifges(5,1) * t92 - Ifges(5,4) * t93 + t55 * Ifges(5,5);
t8 = -Ifges(6,1) * t92 - Ifges(6,4) * t93 + t55 * Ifges(6,5);
t7 = -Ifges(5,4) * t92 - Ifges(5,2) * t93 + t55 * Ifges(5,6);
t6 = -Ifges(6,4) * t92 - Ifges(6,2) * t93 + t55 * Ifges(6,6);
t5 = pkin(4) * t59 - t120 * t88 + t14;
t4 = -qJD(4) * t15 + t105;
t3 = -t118 * t42 + t116;
t2 = -qJ(5) * t112 + (-qJD(4) * t42 + t94) * t86 + t116;
t1 = pkin(4) * t55 + t94 * t88 + (-t40 + (-t35 + t120) * t86) * qJD(4) + t105;
t24 = [0.2e1 * t79 * t106 + 0.2e1 * t11 * t34 + 0.2e1 * t2 * t36 + 0.2e1 * t3 * t37 + 0.2e1 * t1 * t38 + 0.2e1 * t4 * t39 + t13 * t139 + 0.2e1 * t5 * t16 + 0.2e1 * t14 * t17 + 0.2e1 * t10 * t18 + 0.2e1 * t15 * t19 + 0.2e1 * t27 * t12 + t42 * t55 * t142 + 0.2e1 * m(4) * (t114 * t79 + t29 * t42 + t131) + (t1 * t5 + t10 * t2 + t11 * t27) * t143 + 0.2e1 * m(5) * (t14 * t4 + t15 * t3 + t131) - (mrSges(4,3) * t139 + t121 * t88 + t122 * t86) * t56 + (mrSges(4,1) * t104 + t29 * t142 - t100 * t56 + ((2 * Ifges(4,2)) + t145) * t55 + t102) * t59 + (t59 * t144 + mrSges(4,2) * t104 - 0.2e1 * Ifges(4,1) * t56 + t100 * t55 + (qJD(4) * t122 + t125 * t55 + t8 + t9) * t88 + (-qJD(4) * t121 - t6 - t7) * t86 + 0.2e1 * (t99 + mrSges(4,3)) * t28) * t60 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t87 + mrSges(3,2) * t89) + (-Ifges(3,2) + Ifges(3,1)) * t87 * t89 + (-t87 ^ 2 + t89 ^ 2) * Ifges(3,4)) * qJD(2); -t28 * mrSges(4,1) - t29 * mrSges(4,2) - Ifges(4,5) * t56 - Ifges(4,6) * t55 + t11 * t69 + t68 * t12 - t57 * t16 + t58 * t18 + t27 * t62 + t43 * t36 + t44 * t38 + t41 * t63 + (t82 / 0.2e1 + t83 / 0.2e1) * t59 + m(6) * (-t1 * t57 + t10 * t43 + t11 * t68 + t2 * t58 + t44 * t5) + (Ifges(3,5) * t89 - Ifges(3,6) * t87 + (-mrSges(3,1) * t89 + mrSges(3,2) * t87) * pkin(6)) * qJD(2) + (m(4) * (-t28 * t85 + t29 * t84) + (-t55 * t84 + t56 * t85) * mrSges(4,3)) * pkin(2) + (t28 * mrSges(5,2) + t8 / 0.2e1 + t9 / 0.2e1 - t77 * t17 - t1 * mrSges(6,3) + (-t64 / 0.2e1 - t65 / 0.2e1) * t60 - t110 * t56 + (Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1) * t55 - t108 * t4 + (-t77 * t37 - t10 * mrSges(6,3) + pkin(4) * t34 - t20 / 0.2e1 - t21 / 0.2e1 - t109 * t60 - t111 * t59 + t27 * t138 - t108 * t15) * qJD(4)) * t86 + (-t28 * mrSges(5,1) + t6 / 0.2e1 + t7 / 0.2e1 + t77 * t19 + t2 * mrSges(6,3) + (t66 / 0.2e1 + t67 / 0.2e1) * t60 - t109 * t56 + t111 * t55 + t108 * t3 + (-t77 * t39 - t5 * mrSges(6,3) + t22 / 0.2e1 + t23 / 0.2e1 + t110 * t60 - t108 * t14) * qJD(4)) * t88 + (m(5) * t28 + t13) * t78; (t43 * t58 - t44 * t57) * t143 + 0.2e1 * t78 * t63 + 0.2e1 * t68 * t62 + (t44 * t141 + t66 + t67 + (t58 * t141 - t72 - t73 + 0.2e1 * (m(6) * t68 + t69) * pkin(4)) * qJD(4)) * t86 + (0.2e1 * t43 * mrSges(6,3) + t64 + t65 + (-t141 * t57 + t74 + t75) * qJD(4)) * t88; m(4) * t114 + (t16 + t17) * t88 + (t18 + t19) * t86 + ((t36 + t37) * t88 + (-t38 - t39) * t86) * qJD(4) + m(5) * (t3 * t86 + t4 * t88 + (-t14 * t86 + t15 * t88) * qJD(4)) + m(6) * (t1 * t88 + t2 * t86 + (t10 * t88 - t5 * t86) * qJD(4)) + t106; m(6) * (t43 * t86 + t44 * t88 + (t57 * t86 + t58 * t88) * qJD(4)); 0; mrSges(5,1) * t4 + mrSges(6,1) * t1 - mrSges(5,2) * t3 - mrSges(6,2) * t2 - t56 * t107 + (m(6) * t1 + t16) * pkin(4) + t60 * t144 + t102; -mrSges(6,2) * t43 + t82 + t83 + (mrSges(6,1) + t138) * t44 + ((-mrSges(5,1) * t77 - (mrSges(6,3) * pkin(4))) * t88 + (mrSges(5,2) * t77 + t124) * t86) * qJD(4); (-t136 + (-mrSges(5,1) - t138) * t86) * qJD(4) - t62; 0; m(6) * t11 + t12; t118 * t138 + t62; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t24(1), t24(2), t24(4), t24(7), t24(11); t24(2), t24(3), t24(5), t24(8), t24(12); t24(4), t24(5), t24(6), t24(9), t24(13); t24(7), t24(8), t24(9), t24(10), t24(14); t24(11), t24(12), t24(13), t24(14), t24(15);];
Mq = res;
