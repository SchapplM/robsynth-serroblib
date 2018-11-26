% Calculate joint inertia matrix for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:28:53
% EndTime: 2018-11-23 15:28:54
% DurationCPUTime: 1.01s
% Computational Cost: add. (1295->291), mult. (2828->398), div. (0->0), fcn. (2878->10), ass. (0->116)
t145 = 2 * pkin(8);
t101 = sin(qJ(5));
t105 = cos(qJ(5));
t137 = -mrSges(6,2) - mrSges(7,2);
t144 = (t105 * mrSges(6,1) + t137 * t101) * pkin(4);
t100 = cos(pkin(6));
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t104 = sin(qJ(2));
t99 = sin(pkin(6));
t128 = t104 * t99;
t56 = -t100 * t107 + t103 * t128;
t55 = t56 ^ 2;
t143 = 2 * mrSges(7,1);
t142 = m(7) * pkin(5);
t141 = -pkin(10) - pkin(9);
t140 = m(7) * t101;
t139 = pkin(4) * t105;
t138 = pkin(8) * t107;
t93 = t103 * pkin(8);
t136 = Ifges(6,3) + Ifges(7,3);
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t123 = t103 * t106;
t75 = -pkin(3) * t107 - pkin(9) * t103 - pkin(2);
t66 = t106 * t75;
t24 = -pkin(10) * t123 + t66 + (-pkin(8) * t102 - pkin(4)) * t107;
t124 = t102 * t103;
t44 = t102 * t75 + t106 * t138;
t35 = -pkin(10) * t124 + t44;
t8 = t101 * t24 + t105 * t35;
t68 = t101 * t106 + t102 * t105;
t53 = t68 * t103;
t39 = mrSges(7,2) * t107 - mrSges(7,3) * t53;
t40 = mrSges(6,2) * t107 - mrSges(6,3) * t53;
t135 = t39 + t40;
t80 = t141 * t102;
t81 = t141 * t106;
t38 = t101 * t80 - t105 * t81;
t76 = -mrSges(5,1) * t106 + mrSges(5,2) * t102;
t134 = t76 - mrSges(4,1);
t74 = pkin(4) * t124 + t93;
t133 = t102 ^ 2 + t106 ^ 2;
t132 = Ifges(5,4) * t102;
t131 = Ifges(5,4) * t106;
t67 = -t101 * t102 + t105 * t106;
t130 = t101 * t67;
t129 = t103 * t56;
t58 = t100 * t103 + t107 * t128;
t126 = t107 * t58;
t108 = cos(qJ(2));
t125 = t108 * t99;
t122 = Ifges(5,3) + t136;
t54 = t67 * t103;
t7 = -t101 * t35 + t105 * t24;
t3 = -pkin(5) * t107 - qJ(6) * t54 + t7;
t41 = -mrSges(7,1) * t107 - mrSges(7,3) * t54;
t121 = m(7) * t3 + t41;
t88 = -pkin(4) * t106 - pkin(3);
t20 = t53 * mrSges(7,1) + t54 * mrSges(7,2);
t27 = -t67 * mrSges(7,1) + t68 * mrSges(7,2);
t37 = t101 * t81 + t105 * t80;
t119 = (-Ifges(6,5) - Ifges(7,5)) * t54 + (Ifges(6,6) + Ifges(7,6)) * t53;
t18 = -qJ(6) * t68 + t37;
t118 = m(7) * t18 - t68 * mrSges(7,3);
t116 = mrSges(5,1) * t102 + mrSges(5,2) * t106;
t33 = -t102 * t58 - t106 * t125;
t34 = -t102 * t125 + t106 * t58;
t115 = -t102 * t33 + t106 * t34;
t11 = -t101 * t34 + t105 * t33;
t12 = t101 * t33 + t105 * t34;
t114 = t12 * t137 + (mrSges(6,1) + mrSges(7,1)) * t11;
t4 = -qJ(6) * t53 + t8;
t113 = t7 * mrSges(6,1) + t3 * mrSges(7,1) - t8 * mrSges(6,2) - t4 * mrSges(7,2) - t119;
t19 = qJ(6) * t67 + t38;
t61 = Ifges(7,6) * t67;
t62 = Ifges(6,6) * t67;
t63 = Ifges(7,5) * t68;
t64 = Ifges(6,5) * t68;
t112 = t37 * mrSges(6,1) + t18 * mrSges(7,1) - t38 * mrSges(6,2) - t19 * mrSges(7,2) + t61 + t62 + t63 + t64;
t111 = pkin(4) ^ 2;
t110 = pkin(8) ^ 2;
t98 = t107 ^ 2;
t96 = t103 ^ 2;
t94 = t99 ^ 2;
t92 = t96 * t110;
t91 = t101 ^ 2 * t111;
t90 = Ifges(5,5) * t102;
t89 = Ifges(5,6) * t106;
t87 = pkin(5) + t139;
t85 = t94 * t108 ^ 2;
t82 = Ifges(5,5) * t123;
t79 = Ifges(5,1) * t102 + t131;
t78 = Ifges(5,2) * t106 + t132;
t77 = -mrSges(4,1) * t107 + mrSges(4,2) * t103;
t73 = -mrSges(5,1) * t107 - mrSges(5,3) * t123;
t72 = mrSges(5,2) * t107 - mrSges(5,3) * t124;
t59 = t116 * t103;
t52 = -Ifges(5,5) * t107 + (Ifges(5,1) * t106 - t132) * t103;
t51 = -Ifges(5,6) * t107 + (-Ifges(5,2) * t102 + t131) * t103;
t45 = -pkin(5) * t67 + t88;
t43 = -t102 * t138 + t66;
t42 = -mrSges(6,1) * t107 - mrSges(6,3) * t54;
t32 = Ifges(6,1) * t68 + Ifges(6,4) * t67;
t31 = Ifges(7,1) * t68 + Ifges(7,4) * t67;
t30 = Ifges(6,4) * t68 + Ifges(6,2) * t67;
t29 = Ifges(7,4) * t68 + Ifges(7,2) * t67;
t28 = -mrSges(6,1) * t67 + mrSges(6,2) * t68;
t25 = pkin(5) * t53 + t74;
t21 = mrSges(6,1) * t53 + mrSges(6,2) * t54;
t16 = Ifges(6,1) * t54 - Ifges(6,4) * t53 - Ifges(6,5) * t107;
t15 = Ifges(7,1) * t54 - Ifges(7,4) * t53 - Ifges(7,5) * t107;
t14 = Ifges(6,4) * t54 - Ifges(6,2) * t53 - Ifges(6,6) * t107;
t13 = Ifges(7,4) * t54 - Ifges(7,2) * t53 - Ifges(7,6) * t107;
t5 = t101 * pkin(4) * t12;
t1 = [m(2) + m(5) * (t33 ^ 2 + t34 ^ 2 + t55) + m(4) * (t58 ^ 2 + t55 + t85) + m(3) * (t104 ^ 2 * t94 + t100 ^ 2 + t85) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t11 ^ 2 + t12 ^ 2 + t55); mrSges(4,3) * t126 + t33 * t73 + t34 * t72 + t135 * t12 + (t41 + t42) * t11 + (-t104 * mrSges(3,2) + (mrSges(3,1) - t77) * t108) * t99 + (t103 * mrSges(4,3) + t20 + t21 + t59) * t56 + m(6) * (t11 * t7 + t12 * t8 + t56 * t74) + m(7) * (t11 * t3 + t12 * t4 + t25 * t56) + m(5) * (pkin(8) * t129 + t33 * t43 + t34 * t44) + m(4) * (pkin(2) * t125 + (t126 + t129) * pkin(8)); -0.2e1 * pkin(2) * t77 + 0.2e1 * t25 * t20 + 0.2e1 * t74 * t21 + 0.2e1 * t3 * t41 + 0.2e1 * t4 * t39 + 0.2e1 * t8 * t40 + 0.2e1 * t7 * t42 + 0.2e1 * t43 * t73 + 0.2e1 * t44 * t72 + Ifges(3,3) + (t15 + t16) * t54 - (t13 + t14) * t53 + (t96 + t98) * mrSges(4,3) * t145 + (Ifges(4,1) * t103 - t102 * t51 + t106 * t52 + t59 * t145) * t103 + m(4) * (pkin(2) ^ 2 + t110 * t98 + t92) + m(5) * (t43 ^ 2 + t44 ^ 2 + t92) + m(6) * (t7 ^ 2 + t74 ^ 2 + t8 ^ 2) + m(7) * (t25 ^ 2 + t3 ^ 2 + t4 ^ 2) + (-t82 + (Ifges(4,2) + t122) * t107 + (Ifges(5,6) * t102 + (2 * Ifges(4,4))) * t103 + t119) * t107; -t58 * mrSges(4,2) + t115 * mrSges(5,3) + (t27 + t28 + t134) * t56 + m(6) * (t11 * t37 + t12 * t38 + t56 * t88) + m(7) * (t11 * t18 + t12 * t19 + t45 * t56) + m(5) * (-pkin(3) * t56 + pkin(9) * t115) + (mrSges(7,3) + mrSges(6,3)) * (-t11 * t68 + t12 * t67); -pkin(3) * t59 + t18 * t41 + t19 * t39 + t45 * t20 + t88 * t21 + t25 * t27 + t74 * t28 + t37 * t42 + t38 * t40 + (t31 / 0.2e1 + t32 / 0.2e1) * t54 - (t29 / 0.2e1 + t30 / 0.2e1) * t53 + (-t90 / 0.2e1 - t89 / 0.2e1 - t63 / 0.2e1 - t61 / 0.2e1 - t64 / 0.2e1 - t62 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t107 + (t51 / 0.2e1 + pkin(9) * t72 + t44 * mrSges(5,3)) * t106 + (t52 / 0.2e1 - pkin(9) * t73 - t43 * mrSges(5,3)) * t102 + (Ifges(4,5) + t106 * t79 / 0.2e1 - t102 * t78 / 0.2e1 + t134 * pkin(8)) * t103 + m(5) * (-pkin(3) * t93 + (-t43 * t102 + t44 * t106) * pkin(9)) + m(6) * (t37 * t7 + t38 * t8 + t74 * t88) + m(7) * (t18 * t3 + t19 * t4 + t25 * t45) + (t15 / 0.2e1 + t16 / 0.2e1 - t7 * mrSges(6,3) - t3 * mrSges(7,3)) * t68 + (t13 / 0.2e1 + t14 / 0.2e1 + t8 * mrSges(6,3) + t4 * mrSges(7,3)) * t67; -0.2e1 * pkin(3) * t76 + t102 * t79 + t106 * t78 + 0.2e1 * t45 * t27 + 0.2e1 * t88 * t28 + Ifges(4,3) + 0.2e1 * t133 * pkin(9) * mrSges(5,3) + m(7) * (t18 ^ 2 + t19 ^ 2 + t45 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2 + t88 ^ 2) + m(5) * (pkin(9) ^ 2 * t133 + pkin(3) ^ 2) + (-0.2e1 * mrSges(6,3) * t37 - 0.2e1 * mrSges(7,3) * t18 + t31 + t32) * t68 + (0.2e1 * mrSges(6,3) * t38 + 0.2e1 * mrSges(7,3) * t19 + t29 + t30) * t67; t33 * mrSges(5,1) - t34 * mrSges(5,2) + m(6) * (t11 * t139 + t5) + m(7) * (t11 * t87 + t5) + t114; -Ifges(5,6) * t124 + t43 * mrSges(5,1) - t44 * mrSges(5,2) + t82 + t121 * t87 + (t105 * t42 + t135 * t101 + t4 * t140 + m(6) * (t101 * t8 + t105 * t7)) * pkin(4) - t122 * t107 + t113; t89 + t90 + t118 * t87 - t116 * pkin(9) + (mrSges(7,3) * t130 + (-t105 * t68 + t130) * mrSges(6,3) + t19 * t140 + m(6) * (t101 * t38 + t105 * t37)) * pkin(4) + t112; t87 * t143 + m(7) * (t87 ^ 2 + t91) + m(6) * (t105 ^ 2 * t111 + t91) + 0.2e1 * t144 + t122; t11 * t142 + t114; pkin(5) * t121 - t107 * t136 + t113; pkin(5) * t118 + t112; t87 * t142 + (pkin(5) + t87) * mrSges(7,1) + t144 + t136; (t143 + t142) * pkin(5) + t136; m(7) * t56; m(7) * t25 + t20; m(7) * t45 + t27; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
