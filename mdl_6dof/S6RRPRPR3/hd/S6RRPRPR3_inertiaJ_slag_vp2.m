% Calculate joint inertia matrix for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:53
% EndTime: 2019-03-09 10:15:56
% DurationCPUTime: 1.20s
% Computational Cost: add. (2818->285), mult. (5371->423), div. (0->0), fcn. (6127->10), ass. (0->103)
t147 = Ifges(5,3) + Ifges(6,3);
t146 = m(6) * pkin(4);
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t106 = sin(pkin(11));
t108 = cos(pkin(11));
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t84 = -t106 * t111 + t108 * t114;
t107 = sin(pkin(10));
t109 = cos(pkin(10));
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t87 = t107 * t115 + t109 * t112;
t46 = t84 * t87;
t129 = t114 * t87;
t85 = t107 * t112 - t109 * t115;
t99 = -pkin(2) * t115 - pkin(1);
t55 = pkin(3) * t85 - pkin(8) * t87 + t99;
t135 = -qJ(3) - pkin(7);
t123 = t135 * t112;
t91 = t135 * t115;
t66 = t107 * t123 - t109 * t91;
t33 = -t111 * t66 + t114 * t55;
t18 = pkin(4) * t85 - qJ(5) * t129 + t33;
t130 = t111 * t87;
t34 = t111 * t55 + t114 * t66;
t25 = -qJ(5) * t130 + t34;
t6 = -t106 * t25 + t108 * t18;
t4 = pkin(5) * t85 - pkin(9) * t46 + t6;
t86 = t106 * t114 + t108 * t111;
t45 = t86 * t87;
t7 = t106 * t18 + t108 * t25;
t5 = -pkin(9) * t45 + t7;
t2 = -t110 * t5 + t113 * t4;
t3 = t110 * t4 + t113 * t5;
t145 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t144 = Ifges(5,5) * t111 + Ifges(6,5) * t86 + Ifges(5,6) * t114 + Ifges(6,6) * t84;
t64 = -t107 * t91 - t109 * t123;
t143 = t64 ^ 2;
t142 = 0.2e1 * t64;
t141 = 0.2e1 * t99;
t138 = pkin(4) * t106;
t97 = pkin(4) * t108 + pkin(5);
t73 = -t110 * t138 + t113 * t97;
t137 = t73 * mrSges(7,1);
t74 = t110 * t97 + t113 * t138;
t136 = t74 * mrSges(7,2);
t56 = -t110 * t86 + t113 * t84;
t57 = t110 * t84 + t113 * t86;
t134 = Ifges(7,5) * t57 + Ifges(7,6) * t56;
t96 = pkin(2) * t107 + pkin(8);
t128 = qJ(5) + t96;
t82 = t128 * t111;
t83 = t128 * t114;
t50 = -t106 * t82 + t108 * t83;
t132 = Ifges(5,4) * t111;
t131 = Ifges(5,4) * t114;
t126 = t111 ^ 2 + t114 ^ 2;
t125 = t112 ^ 2 + t115 ^ 2;
t26 = -t110 * t46 - t113 * t45;
t27 = -t110 * t45 + t113 * t46;
t124 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t85;
t98 = -pkin(2) * t109 - pkin(3);
t29 = t45 * mrSges(6,1) + t46 * mrSges(6,2);
t60 = -t84 * mrSges(6,1) + t86 * mrSges(6,2);
t10 = -t26 * mrSges(7,1) + t27 * mrSges(7,2);
t30 = -t56 * mrSges(7,1) + t57 * mrSges(7,2);
t49 = -t106 * t83 - t108 * t82;
t41 = pkin(4) * t130 + t64;
t90 = -t114 * mrSges(5,1) + t111 * mrSges(5,2);
t122 = mrSges(5,1) * t111 + mrSges(5,2) * t114;
t37 = -pkin(9) * t86 + t49;
t38 = pkin(9) * t84 + t50;
t12 = -t110 * t38 + t113 * t37;
t13 = t110 * t37 + t113 * t38;
t121 = t12 * mrSges(7,1) - t13 * mrSges(7,2) + t134;
t89 = -pkin(4) * t114 + t98;
t120 = -t30 - t60;
t119 = Ifges(5,5) * t129 + Ifges(6,5) * t46 - Ifges(6,6) * t45 + t147 * t85 + t124;
t93 = Ifges(5,1) * t111 + t131;
t92 = Ifges(5,2) * t114 + t132;
t76 = t87 * mrSges(4,2);
t67 = -pkin(5) * t84 + t89;
t62 = Ifges(6,1) * t86 + Ifges(6,4) * t84;
t61 = Ifges(6,4) * t86 + Ifges(6,2) * t84;
t59 = mrSges(5,1) * t85 - mrSges(5,3) * t129;
t58 = -mrSges(5,2) * t85 - mrSges(5,3) * t130;
t51 = t122 * t87;
t40 = Ifges(5,5) * t85 + (Ifges(5,1) * t114 - t132) * t87;
t39 = Ifges(5,6) * t85 + (-Ifges(5,2) * t111 + t131) * t87;
t36 = mrSges(6,1) * t85 - mrSges(6,3) * t46;
t35 = -mrSges(6,2) * t85 - mrSges(6,3) * t45;
t32 = Ifges(7,1) * t57 + Ifges(7,4) * t56;
t31 = Ifges(7,4) * t57 + Ifges(7,2) * t56;
t28 = pkin(5) * t45 + t41;
t20 = Ifges(6,1) * t46 - Ifges(6,4) * t45 + Ifges(6,5) * t85;
t19 = Ifges(6,4) * t46 - Ifges(6,2) * t45 + Ifges(6,6) * t85;
t17 = mrSges(7,1) * t85 - mrSges(7,3) * t27;
t16 = -mrSges(7,2) * t85 + mrSges(7,3) * t26;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t85;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t85;
t1 = [m(3) * (pkin(7) ^ 2 * t125 + pkin(1) ^ 2) + m(6) * (t41 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(7) * (t2 ^ 2 + t28 ^ 2 + t3 ^ 2) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t115 + mrSges(3,2) * t112) + t112 * (Ifges(3,1) * t112 + Ifges(3,4) * t115) + t115 * (Ifges(3,4) * t112 + Ifges(3,2) * t115) + 0.2e1 * t34 * t58 + 0.2e1 * t33 * t59 + 0.2e1 * t41 * t29 - t45 * t19 + t46 * t20 + t27 * t9 + 0.2e1 * t28 * t10 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t36 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + t26 * t8 + t76 * t141 + t51 * t142 + 0.2e1 * t125 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t141 - 0.2e1 * t66 * mrSges(4,3) + Ifges(4,2) * t85 + (-Ifges(5,6) * t111 - (2 * Ifges(4,4))) * t87 + t119) * t85 + (mrSges(4,3) * t142 + Ifges(4,1) * t87 - t111 * t39 + t114 * t40) * t87 + m(4) * (t66 ^ 2 + t99 ^ 2 + t143) + m(5) * (t33 ^ 2 + t34 ^ 2 + t143) + Ifges(2,3); (m(4) * (t107 * t66 - t109 * t64) + (-t107 * t85 - t109 * t87) * mrSges(4,3)) * pkin(2) + (t134 + t144) * t85 / 0.2e1 + Ifges(3,6) * t115 + Ifges(3,5) * t112 + t89 * t29 + t98 * t51 + t84 * t19 / 0.2e1 - Ifges(4,6) * t85 + t86 * t20 / 0.2e1 + Ifges(4,5) * t87 + t41 * t60 - t45 * t61 / 0.2e1 + t46 * t62 / 0.2e1 - t66 * mrSges(4,2) + t67 * t10 + t49 * t36 + t50 * t35 + t56 * t8 / 0.2e1 + t57 * t9 / 0.2e1 + t28 * t30 + t26 * t31 / 0.2e1 + t27 * t32 / 0.2e1 + t13 * t16 + t12 * t17 + (t90 - mrSges(4,1)) * t64 + (t87 * t93 / 0.2e1 + t96 * t58 + t34 * mrSges(5,3) + t39 / 0.2e1) * t114 + (-t87 * t92 / 0.2e1 - t96 * t59 - t33 * mrSges(5,3) + t40 / 0.2e1) * t111 + m(6) * (t41 * t89 + t49 * t6 + t50 * t7) + m(7) * (t12 * t2 + t13 * t3 + t28 * t67) + m(5) * (t64 * t98 + (-t33 * t111 + t34 * t114) * t96) + (-t112 * mrSges(3,1) - t115 * mrSges(3,2)) * pkin(7) + (-t2 * t57 + t3 * t56) * mrSges(7,3) + (-t6 * t86 + t7 * t84) * mrSges(6,3); t111 * t93 + t114 * t92 + 0.2e1 * t67 * t30 + t56 * t31 + t57 * t32 + 0.2e1 * t89 * t60 + t84 * t61 + t86 * t62 + 0.2e1 * t98 * t90 + Ifges(3,3) + Ifges(4,3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t67 ^ 2) + m(6) * (t49 ^ 2 + t50 ^ 2 + t89 ^ 2) + m(5) * (t126 * t96 ^ 2 + t98 ^ 2) + m(4) * (t107 ^ 2 + t109 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t109 * mrSges(4,1) - t107 * mrSges(4,2)) * pkin(2) + 0.2e1 * (-t12 * t57 + t13 * t56) * mrSges(7,3) + 0.2e1 * (-t49 * t86 + t50 * t84) * mrSges(6,3) + 0.2e1 * t126 * t96 * mrSges(5,3); t85 * mrSges(4,1) + t111 * t58 + t114 * t59 + t57 * t16 + t56 * t17 + t86 * t35 + t84 * t36 + t76 + m(7) * (t2 * t56 + t3 * t57) + m(6) * (t6 * t84 + t7 * t86) + m(5) * (t111 * t34 + t114 * t33) + m(4) * t99; m(7) * (t12 * t56 + t13 * t57) + m(6) * (t49 * t84 + t50 * t86); m(4) + m(5) * t126 + m(6) * (t84 ^ 2 + t86 ^ 2) + m(7) * (t56 ^ 2 + t57 ^ 2); m(7) * (t2 * t73 + t3 * t74) - Ifges(5,6) * t130 + t74 * t16 + t73 * t17 + t33 * mrSges(5,1) - t34 * mrSges(5,2) + t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(6) * (t106 * t7 + t108 * t6) + t106 * t35 + t108 * t36) * pkin(4) + t119 + t145; m(7) * (t12 * t73 + t13 * t74) - t50 * mrSges(6,2) + t49 * mrSges(6,1) - t122 * t96 + (t56 * t74 - t57 * t73) * mrSges(7,3) + (m(6) * (t106 * t50 + t108 * t49) + (t106 * t84 - t108 * t86) * mrSges(6,3)) * pkin(4) + t121 + t144; m(7) * (t56 * t73 + t57 * t74) + (t106 * t86 + t108 * t84) * t146 + t120 - t90; 0.2e1 * t137 - 0.2e1 * t136 + Ifges(7,3) + m(7) * (t73 ^ 2 + t74 ^ 2) + (0.2e1 * mrSges(6,1) * t108 - 0.2e1 * mrSges(6,2) * t106 + (t106 ^ 2 + t108 ^ 2) * t146) * pkin(4) + t147; m(6) * t41 + m(7) * t28 + t10 + t29; m(6) * t89 + m(7) * t67 - t120; 0; 0; m(6) + m(7); t124 + t145; t121; -t30; Ifges(7,3) - t136 + t137; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
