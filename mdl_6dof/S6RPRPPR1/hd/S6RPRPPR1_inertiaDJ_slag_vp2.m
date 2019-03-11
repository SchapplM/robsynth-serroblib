% Calculate time derivative of joint inertia matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:06
% EndTime: 2019-03-09 02:38:09
% DurationCPUTime: 1.51s
% Computational Cost: add. (3006->287), mult. (6593->447), div. (0->0), fcn. (6207->10), ass. (0->117)
t91 = sin(pkin(11));
t93 = cos(pkin(11));
t96 = sin(qJ(6));
t98 = cos(qJ(6));
t101 = t91 * t96 - t93 * t98;
t73 = t101 * qJD(6);
t136 = 2 * m(6);
t135 = 2 * m(7);
t90 = t93 ^ 2;
t87 = sin(pkin(9)) * pkin(1) + pkin(7);
t116 = qJ(4) + t87;
t106 = qJD(3) * t116;
t97 = sin(qJ(3));
t99 = cos(qJ(3));
t100 = -t97 * qJD(4) - t106 * t99;
t63 = qJD(4) * t99 - t106 * t97;
t92 = sin(pkin(10));
t94 = cos(pkin(10));
t37 = -t94 * t100 + t63 * t92;
t134 = 0.2e1 * t37;
t112 = -cos(pkin(9)) * pkin(1) - pkin(2);
t83 = -pkin(3) * t99 + t112;
t133 = 0.2e1 * t83;
t132 = m(5) * pkin(3);
t85 = pkin(3) * t92 + qJ(5);
t131 = pkin(8) + t85;
t78 = t92 * t97 - t94 * t99;
t72 = t78 * qJD(3);
t120 = t93 * t72;
t121 = t91 * t72;
t41 = -mrSges(6,1) * t121 - mrSges(6,2) * t120;
t81 = t91 * t98 + t93 * t96;
t74 = t81 * qJD(6);
t80 = t92 * t99 + t94 * t97;
t20 = t101 * t72 - t74 * t80;
t21 = t72 * t81 + t73 * t80;
t7 = -t21 * mrSges(7,1) + t20 * mrSges(7,2);
t130 = t41 + t7;
t129 = Ifges(6,4) * t91;
t128 = Ifges(6,4) * t93;
t107 = t116 * t97;
t77 = t116 * t99;
t50 = t94 * t107 + t77 * t92;
t127 = t37 * t50;
t71 = t80 * qJD(3);
t126 = t71 * Ifges(6,5);
t125 = t71 * Ifges(6,6);
t61 = t78 * t71;
t124 = t80 * t91;
t123 = t80 * t93;
t122 = t91 * Ifges(6,2);
t119 = -mrSges(6,1) * t93 + mrSges(6,2) * t91 - mrSges(5,1);
t113 = pkin(3) * qJD(3) * t97;
t30 = pkin(4) * t71 + qJ(5) * t72 - qJD(5) * t80 + t113;
t38 = t100 * t92 + t94 * t63;
t12 = t91 * t30 + t93 * t38;
t43 = t78 * pkin(4) - t80 * qJ(5) + t83;
t51 = -t107 * t92 + t94 * t77;
t23 = t91 * t43 + t93 * t51;
t118 = -Ifges(7,5) * t73 - Ifges(7,6) * t74;
t117 = t91 ^ 2 + t90;
t115 = 0.2e1 * t99;
t114 = Ifges(7,5) * t20 + Ifges(7,6) * t21 + Ifges(7,3) * t71;
t88 = -pkin(3) * t94 - pkin(4);
t111 = mrSges(6,3) * t117;
t110 = t117 * t80;
t109 = t117 * t85;
t67 = t72 * mrSges(5,2);
t108 = t71 * mrSges(5,1) - t67;
t11 = t93 * t30 - t38 * t91;
t22 = t93 * t43 - t51 * t91;
t105 = Ifges(6,5) * t93 - Ifges(6,6) * t91;
t10 = pkin(5) * t78 - pkin(8) * t123 + t22;
t15 = -pkin(8) * t124 + t23;
t3 = t10 * t98 - t15 * t96;
t4 = t10 * t96 + t15 * t98;
t104 = -t11 * t91 + t12 * t93;
t103 = -t22 * t91 + t23 * t93;
t102 = t37 * t78 + t50 * t71;
t75 = t131 * t91;
t76 = t131 * t93;
t48 = -t75 * t98 - t76 * t96;
t49 = -t75 * t96 + t76 * t98;
t82 = -pkin(5) * t93 + t88;
t60 = Ifges(7,1) * t81 - Ifges(7,4) * t101;
t59 = Ifges(7,4) * t81 - Ifges(7,2) * t101;
t58 = mrSges(7,1) * t101 + mrSges(7,2) * t81;
t57 = mrSges(6,1) * t78 - mrSges(6,3) * t123;
t56 = -mrSges(6,2) * t78 - mrSges(6,3) * t124;
t55 = -Ifges(7,1) * t73 - Ifges(7,4) * t74;
t54 = -Ifges(7,4) * t73 - Ifges(7,2) * t74;
t53 = mrSges(7,1) * t74 - mrSges(7,2) * t73;
t52 = (mrSges(6,1) * t91 + mrSges(6,2) * t93) * t80;
t47 = mrSges(6,1) * t71 + mrSges(6,3) * t120;
t46 = -mrSges(6,2) * t71 + mrSges(6,3) * t121;
t45 = t101 * t80;
t44 = t81 * t80;
t36 = pkin(5) * t124 + t50;
t34 = -qJD(5) * t81 - qJD(6) * t49;
t33 = -qJD(5) * t101 + qJD(6) * t48;
t32 = mrSges(7,1) * t78 + mrSges(7,3) * t45;
t31 = -mrSges(7,2) * t78 - mrSges(7,3) * t44;
t29 = t126 - (t93 * Ifges(6,1) - t129) * t72;
t28 = t125 - (-t122 + t128) * t72;
t25 = -pkin(5) * t121 + t37;
t24 = mrSges(7,1) * t44 - mrSges(7,2) * t45;
t17 = -Ifges(7,1) * t45 - Ifges(7,4) * t44 + Ifges(7,5) * t78;
t16 = -Ifges(7,4) * t45 - Ifges(7,2) * t44 + Ifges(7,6) * t78;
t14 = -mrSges(7,2) * t71 + mrSges(7,3) * t21;
t13 = mrSges(7,1) * t71 - mrSges(7,3) * t20;
t9 = pkin(8) * t121 + t12;
t8 = pkin(5) * t71 + pkin(8) * t120 + t11;
t6 = Ifges(7,1) * t20 + Ifges(7,4) * t21 + t71 * Ifges(7,5);
t5 = Ifges(7,4) * t20 + Ifges(7,2) * t21 + t71 * Ifges(7,6);
t2 = -qJD(6) * t4 + t8 * t98 - t9 * t96;
t1 = qJD(6) * t3 + t8 * t96 + t9 * t98;
t18 = [t108 * t133 + t52 * t134 + (t1 * t4 + t2 * t3 + t25 * t36) * t135 + (t11 * t22 + t12 * t23 + t127) * t136 + (-0.2e1 * t38 * mrSges(5,3) + ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t71 + t114) * t78 - 0.2e1 * (-Ifges(5,4) + t105) * t72 * t78 + t71 * (-Ifges(7,5) * t45 - Ifges(7,6) * t44) + 0.2e1 * (-t50 * t72 - t51 * t71) * mrSges(5,3) + (mrSges(5,3) * t134 - t91 * t28 + t93 * t29 + (-0.2e1 * Ifges(5,4) + t105) * t71 - (Ifges(6,1) * t90 + (2 * Ifges(5,1)) + (t122 - 0.2e1 * t128) * t91) * t72) * t80 + ((t112 * mrSges(4,2) + Ifges(4,4) * t99) * t115 + (t132 * t133 + 0.2e1 * pkin(3) * (mrSges(5,1) * t78 + mrSges(5,2) * t80) + 0.2e1 * t112 * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t97 + (Ifges(4,1) - Ifges(4,2)) * t115) * t97) * qJD(3) + 0.2e1 * m(5) * (t38 * t51 + t127) + 0.2e1 * t3 * t13 + 0.2e1 * t4 * t14 + t20 * t17 + t21 * t16 + 0.2e1 * t25 * t24 + 0.2e1 * t1 * t31 + 0.2e1 * t2 * t32 + 0.2e1 * t36 * t7 - t44 * t5 - t45 * t6 + 0.2e1 * t23 * t46 + 0.2e1 * t22 * t47 + 0.2e1 * t50 * t41 + 0.2e1 * t12 * t56 + 0.2e1 * t11 * t57; -t44 * t13 - t45 * t14 + t20 * t31 + t21 * t32 + (t93 * t46 - t91 * t47) * t80 + t130 * t78 - (t93 * t56 - t91 * t57) * t72 + (t24 + t52) * t71 + m(7) * (-t1 * t45 - t2 * t44 + t20 * t4 + t21 * t3 + t25 * t78 + t36 * t71) + m(5) * (t38 * t80 - t51 * t72 + t102) + m(6) * (-t103 * t72 + t104 * t80 + t102); 0.2e1 * m(7) * (-t20 * t45 - t21 * t44 + t61) + 0.2e1 * m(5) * (-t72 * t80 + t61) + 0.2e1 * m(6) * (-t110 * t72 + t61); (Ifges(4,5) * t99 - Ifges(4,6) * t97 + (-mrSges(4,1) * t99 + mrSges(4,2) * t97) * t87) * qJD(3) + (m(5) * (-t37 * t94 + t38 * t92) + (-t92 * t71 + t94 * t72) * mrSges(5,3)) * pkin(3) - (t93 * (Ifges(6,1) * t91 + t128) / 0.2e1 - t91 * (Ifges(6,2) * t93 + t129) / 0.2e1 + Ifges(5,5)) * t72 + (qJD(5) * t56 + t12 * mrSges(6,3) + t85 * t46 + t28 / 0.2e1 + t125 / 0.2e1) * t93 + (-qJD(5) * t57 - t11 * mrSges(6,3) - t85 * t47 + t29 / 0.2e1 + t126 / 0.2e1) * t91 + t78 * t118 / 0.2e1 + t119 * t37 + m(6) * (t103 * qJD(5) + t104 * t85 + t37 * t88) + m(7) * (t1 * t49 + t2 * t48 + t25 * t82 + t3 * t34 + t33 * t4) + (-t1 * t101 - t2 * t81 + t3 * t73 - t4 * t74) * mrSges(7,3) + t71 * (Ifges(7,5) * t81 - Ifges(7,6) * t101) / 0.2e1 - t101 * t5 / 0.2e1 + t33 * t31 + t34 * t32 - t38 * mrSges(5,2) + t48 * t13 + t49 * t14 + t36 * t53 - t44 * t54 / 0.2e1 - t45 * t55 / 0.2e1 + t25 * t58 + t21 * t59 / 0.2e1 + t20 * t60 / 0.2e1 - Ifges(5,6) * t71 - t73 * t17 / 0.2e1 - t74 * t16 / 0.2e1 + t81 * t6 / 0.2e1 + t82 * t7 + t88 * t41; t78 * t53 + t67 + (-mrSges(4,1) * t97 - mrSges(4,2) * t99) * qJD(3) - t72 * t111 + (t58 + t119) * t71 + m(7) * (t20 * t49 + t21 * t48 - t33 * t45 - t34 * t44 + t71 * t82) + m(6) * (qJD(5) * t110 - t109 * t72 + t71 * t88) + (-t71 * t94 - t72 * t92) * t132 + (-t101 * t20 - t21 * t81 - t44 * t73 + t45 * t74) * mrSges(7,3); 0.2e1 * t82 * t53 - t73 * t60 + t81 * t55 - t74 * t59 - t101 * t54 + (t33 * t49 + t34 * t48) * t135 + 0.2e1 * (-t101 * t33 - t34 * t81 + t48 * t73 - t49 * t74) * mrSges(7,3) + (t109 * t136 + 0.2e1 * t111) * qJD(5); m(5) * t113 - t101 * t13 + t81 * t14 - t73 * t31 - t74 * t32 + t91 * t46 + t93 * t47 + m(7) * (t1 * t81 - t101 * t2 - t3 * t74 - t4 * t73) + m(6) * (t11 * t93 + t12 * t91) + t108; m(7) * (-t101 * t21 + t20 * t81 + t44 * t74 + t45 * t73); m(7) * (-t101 * t34 + t33 * t81 - t48 * t74 - t49 * t73); (t101 * t74 - t73 * t81) * t135; m(6) * t37 + m(7) * t25 + t130; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t71; t53; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t114; -t7; mrSges(7,1) * t34 - mrSges(7,2) * t33 + t118; -t53; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
