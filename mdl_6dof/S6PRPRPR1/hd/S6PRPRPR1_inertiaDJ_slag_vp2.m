% Calculate time derivative of joint inertia matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:46
% EndTime: 2019-03-08 19:24:49
% DurationCPUTime: 1.35s
% Computational Cost: add. (2107->267), mult. (5406->412), div. (0->0), fcn. (5293->12), ass. (0->126)
t116 = cos(pkin(12));
t71 = sin(pkin(12));
t77 = sin(qJ(4));
t123 = t71 * t77;
t80 = cos(qJ(4));
t85 = t116 * t80 - t123;
t49 = t85 * qJD(4);
t76 = sin(qJ(6));
t79 = cos(qJ(6));
t101 = (t76 ^ 2 + t79 ^ 2) * t49;
t135 = m(6) * pkin(4);
t139 = t71 * t135 - mrSges(6,2);
t58 = -mrSges(7,1) * t79 + mrSges(7,2) * t76;
t100 = t116 * pkin(4);
t64 = -t100 - pkin(5);
t138 = m(7) * t64 - t116 * t135 - mrSges(6,1) + t58;
t72 = sin(pkin(11));
t63 = pkin(2) * t72 + pkin(8);
t117 = qJ(5) + t63;
t96 = qJD(4) * t117;
t39 = qJD(5) * t80 - t77 * t96;
t84 = -qJD(5) * t77 - t80 * t96;
t21 = -t116 * t84 + t39 * t71;
t137 = 0.2e1 * t21;
t50 = t117 * t80;
t97 = t116 * t77;
t28 = t117 * t97 + t50 * t71;
t136 = 0.2e1 * t28;
t134 = t79 / 0.2e1;
t133 = pkin(4) * t71;
t73 = sin(pkin(6));
t74 = cos(pkin(11));
t78 = sin(qJ(2));
t81 = cos(qJ(2));
t44 = (t72 * t81 + t74 * t78) * t73;
t75 = cos(pkin(6));
t36 = -t44 * t77 + t75 * t80;
t37 = t44 * t80 + t75 * t77;
t14 = -t116 * t36 + t37 * t71;
t115 = qJD(2) * t73;
t88 = t72 * t78 - t74 * t81;
t41 = t88 * t115;
t19 = -t37 * qJD(4) + t41 * t77;
t20 = t36 * qJD(4) - t41 * t80;
t5 = -t116 * t19 + t20 * t71;
t132 = t14 * t5;
t131 = Ifges(7,4) * t76;
t130 = Ifges(7,4) * t79;
t129 = Ifges(7,6) * t76;
t128 = t21 * t28;
t40 = qJD(2) * t44;
t43 = t88 * t73;
t25 = t43 * t40;
t52 = t71 * t80 + t97;
t48 = t52 * qJD(4);
t127 = t48 * mrSges(6,3);
t126 = t48 * t85;
t124 = t49 * t79;
t122 = t76 * mrSges(7,3);
t59 = Ifges(7,2) * t79 + t131;
t121 = t76 * t59;
t120 = t79 * mrSges(7,3);
t60 = Ifges(7,1) * t76 + t130;
t119 = t79 * t60;
t118 = Ifges(7,5) * t124 + Ifges(7,3) * t48;
t114 = qJD(4) * t77;
t113 = qJD(4) * t80;
t112 = qJD(6) * t52;
t111 = qJD(6) * t76;
t110 = qJD(6) * t79;
t65 = -pkin(2) * t74 - pkin(3);
t57 = -pkin(4) * t80 + t65;
t27 = -pkin(5) * t85 - pkin(9) * t52 + t57;
t29 = t116 * t50 - t117 * t123;
t11 = t27 * t79 - t29 * t76;
t109 = t11 * qJD(6);
t12 = t27 * t76 + t29 * t79;
t108 = t12 * qJD(6);
t107 = 0.2e1 * t77;
t106 = pkin(4) * t114;
t105 = t52 * t111;
t104 = t52 * t110;
t103 = mrSges(7,3) * t111;
t102 = mrSges(7,3) * t110;
t99 = -t111 / 0.2e1;
t30 = t48 * mrSges(6,1) + t49 * mrSges(6,2);
t98 = -(2 * Ifges(6,4)) - t129;
t15 = t116 * t37 + t71 * t36;
t10 = t15 * t79 + t43 * t76;
t9 = -t15 * t76 + t43 * t79;
t95 = t10 * t79 - t76 * t9;
t94 = t14 * t21 + t28 * t5;
t93 = t14 * t48 - t5 * t85;
t92 = Ifges(7,1) * t79 - t131;
t91 = -Ifges(7,2) * t76 + t130;
t90 = Ifges(7,5) * t76 + Ifges(7,6) * t79;
t89 = -t21 * t85 + t28 * t48;
t87 = t49 * t76 + t104;
t86 = t105 - t124;
t6 = t116 * t20 + t71 * t19;
t1 = t9 * qJD(6) + t40 * t76 + t6 * t79;
t2 = -t10 * qJD(6) + t40 * t79 - t6 * t76;
t83 = t1 * t79 - t2 * t76 + (-t10 * t76 - t79 * t9) * qJD(6);
t82 = -t19 * t77 + t20 * t80 + (-t36 * t80 - t37 * t77) * qJD(4);
t68 = Ifges(7,5) * t110;
t62 = pkin(9) + t133;
t56 = t92 * qJD(6);
t55 = t91 * qJD(6);
t54 = (mrSges(5,1) * t77 + mrSges(5,2) * t80) * qJD(4);
t53 = -mrSges(7,1) * t111 - mrSges(7,2) * t110;
t35 = -mrSges(6,1) * t85 + mrSges(6,2) * t52;
t33 = -mrSges(7,1) * t85 - t52 * t120;
t32 = mrSges(7,2) * t85 - t52 * t122;
t31 = (mrSges(7,1) * t76 + mrSges(7,2) * t79) * t52;
t26 = pkin(5) * t48 - pkin(9) * t49 + t106;
t24 = -Ifges(7,5) * t85 + t92 * t52;
t23 = -Ifges(7,6) * t85 + t91 * t52;
t22 = t116 * t39 + t71 * t84;
t18 = -mrSges(7,2) * t48 - t87 * mrSges(7,3);
t17 = mrSges(7,1) * t48 + t86 * mrSges(7,3);
t13 = t87 * mrSges(7,1) - t86 * mrSges(7,2);
t8 = -t86 * Ifges(7,1) - t87 * Ifges(7,4) + Ifges(7,5) * t48;
t7 = -t86 * Ifges(7,4) - t87 * Ifges(7,2) + Ifges(7,6) * t48;
t4 = -t22 * t76 + t26 * t79 - t108;
t3 = t22 * t79 + t26 * t76 + t109;
t16 = [0.2e1 * m(7) * (t1 * t10 + t2 * t9 + t132) + 0.2e1 * m(5) * (t19 * t36 + t20 * t37 + t25) + 0.2e1 * m(6) * (t15 * t6 + t132 + t25) + 0.2e1 * m(4) * (-t41 * t44 + t25); t41 * mrSges(4,2) + t1 * t32 + t10 * t18 + t14 * t13 + t9 * t17 + t2 * t33 + t5 * t31 + (t30 + t54) * t43 + (-mrSges(3,1) * t78 - mrSges(3,2) * t81) * t115 + (-t80 * mrSges(5,1) + t77 * mrSges(5,2) - mrSges(4,1) + t35) * t40 + (t14 * t49 - t15 * t48 + t5 * t52 + t6 * t85) * mrSges(6,3) + t82 * mrSges(5,3) + m(6) * (t43 * t106 + t15 * t22 + t29 * t6 + t40 * t57 + t94) + m(7) * (t1 * t12 + t10 * t3 + t11 * t2 + t4 * t9 + t94) + m(4) * (-t40 * t74 - t41 * t72) * pkin(2) + (t40 * t65 + t82 * t63) * m(5); -0.2e1 * t29 * t127 + 0.2e1 * t11 * t17 + 0.2e1 * t12 * t18 + t13 * t136 + t31 * t137 + 0.2e1 * t3 * t32 + 0.2e1 * t57 * t30 + 0.2e1 * t4 * t33 + 0.2e1 * t65 * t54 + (-Ifges(5,4) * t77 + pkin(4) * t35) * qJD(4) * t107 + 0.2e1 * m(6) * (t57 * t106 + t22 * t29 + t128) + 0.2e1 * m(7) * (t11 * t4 + t12 * t3 + t128) + (mrSges(6,3) * t136 - t76 * t23 + t79 * t24) * t49 + (0.2e1 * Ifges(5,4) * t80 + (Ifges(5,1) - Ifges(5,2)) * t107) * t113 - (-0.2e1 * t22 * mrSges(6,3) + t98 * t49 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t48 + t118) * t85 + (mrSges(6,3) * t137 + 0.2e1 * Ifges(6,1) * t49 - t76 * t7 + t79 * t8 + (Ifges(7,5) * t79 + t98) * t48 + (-t79 * t23 - t76 * t24 + t85 * t90) * qJD(6)) * t52; m(7) * (t95 * t49 + t83 * t52 + t93) + m(5) * (t19 * t80 + t20 * t77 + (-t36 * t77 + t37 * t80) * qJD(4)) + m(6) * (t15 * t49 + t52 * t6 + t93); -t85 * t13 + t48 * t31 + m(7) * t89 + m(6) * (t22 * t52 + t29 * t49 + t89) + (-t33 * t112 + (m(7) * t12 + t32) * t49 + (m(7) * (-t109 + t3) + t18) * t52) * t79 + (-t32 * t112 + (-m(7) * t11 - t33) * t49 + (m(7) * (-t108 - t4) - t17) * t52) * t76; 0.2e1 * m(6) * (t49 * t52 - t126) + 0.2e1 * m(7) * (t52 * t101 - t126); m(7) * t83 * t62 + t19 * mrSges(5,1) - t20 * mrSges(5,2) + t1 * t120 - t10 * t103 - t9 * t102 - t2 * t122 + t138 * t5 + t139 * t6 - t14 * t53; t3 * t120 + t7 * t134 - t28 * t53 + t64 * t13 + t76 * t8 / 0.2e1 + t23 * t99 - t11 * t102 - t12 * t103 - t59 * t104 / 0.2e1 + t24 * t110 / 0.2e1 - t85 * (-Ifges(7,6) * t111 + t68) / 0.2e1 - t4 * t122 - t127 * t133 + (t90 / 0.2e1 - Ifges(6,6)) * t48 + t139 * t22 + (t63 * mrSges(5,2) - Ifges(5,6)) * t114 + (-t63 * mrSges(5,1) + Ifges(5,5)) * t113 + (t56 * t134 + t60 * t99 - t76 * t55 / 0.2e1) * t52 + t138 * t21 + (-t33 * t110 - t32 * t111 + m(7) * (t3 * t79 - t4 * t76 + (-t11 * t79 - t12 * t76) * qJD(6)) + t79 * t18 - t76 * t17) * t62 + (-mrSges(6,3) * t100 + Ifges(6,5) + t119 / 0.2e1 - t121 / 0.2e1) * t49; -mrSges(5,2) * t113 - mrSges(5,1) * t114 + (-t116 * t48 + t49 * t71) * t135 + t48 * t58 + t85 * t53 + m(7) * (t62 * t101 + t64 * t48) - t30 + mrSges(7,3) * t101; -0.2e1 * t53 * t64 + t55 * t79 + t56 * t76 + (t119 - t121) * qJD(6); m(7) * (t95 * qJD(6) + t1 * t76 + t2 * t79) + m(6) * t40; m(7) * (t3 * t76 + t4 * t79 + (-t11 * t76 + t12 * t79) * qJD(6)) - t33 * t111 + t79 * t17 + t32 * t110 + t76 * t18 + m(6) * t106 + t30; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1; mrSges(7,1) * t4 - mrSges(7,2) * t3 - Ifges(7,5) * t105 - t87 * Ifges(7,6) + t118; -t13; t68 + (t58 * t62 - t129) * qJD(6); t53; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
