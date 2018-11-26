% Calculate joint inertia matrix for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:07:57
% EndTime: 2018-11-23 18:07:58
% DurationCPUTime: 1.11s
% Computational Cost: add. (1522->331), mult. (2960->428), div. (0->0), fcn. (2772->6), ass. (0->114)
t147 = 2 * pkin(7);
t146 = (mrSges(7,2) + mrSges(6,3));
t106 = sin(qJ(4));
t90 = pkin(3) * t106 + qJ(5);
t145 = t90 ^ 2;
t144 = -2 * mrSges(7,3);
t143 = -pkin(9) - pkin(8);
t111 = cos(qJ(2));
t142 = pkin(7) * t111;
t108 = sin(qJ(2));
t99 = t108 * pkin(7);
t107 = sin(qJ(3));
t109 = cos(qJ(4));
t110 = cos(qJ(3));
t129 = t109 * t110;
t74 = t106 * t107 - t129;
t141 = t74 * t90;
t139 = mrSges(6,2) - mrSges(7,3);
t138 = pkin(4) + qJ(6);
t130 = t108 * t110;
t80 = -pkin(2) * t111 - pkin(8) * t108 - pkin(1);
t73 = t110 * t80;
t27 = -pkin(9) * t130 + t73 + (-pkin(7) * t107 - pkin(3)) * t111;
t131 = t107 * t108;
t53 = t107 * t80 + t110 * t142;
t39 = -pkin(9) * t131 + t53;
t10 = t106 * t27 + t109 * t39;
t75 = t106 * t110 + t107 * t109;
t63 = t75 * t108;
t47 = mrSges(6,1) * t63 + mrSges(6,3) * t111;
t48 = -t63 * mrSges(7,1) - t111 * mrSges(7,2);
t137 = -t47 + t48;
t64 = -t106 * t131 + t108 * t129;
t46 = t64 * mrSges(7,1) + t111 * mrSges(7,3);
t79 = pkin(3) * t131 + t99;
t136 = Ifges(4,4) * t107;
t135 = Ifges(4,4) * t110;
t134 = qJ(5) * t74;
t133 = qJ(5) * t90;
t128 = t107 ^ 2 + t110 ^ 2;
t127 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t125 = t143 * t107;
t84 = t143 * t110;
t42 = -t106 * t84 - t109 * t125;
t44 = t106 * t125 - t109 * t84;
t126 = t42 ^ 2 + t44 ^ 2;
t94 = -pkin(3) * t110 - pkin(2);
t93 = -pkin(3) * t109 - pkin(4);
t9 = -t106 * t39 + t109 * t27;
t124 = Ifges(4,3) + t127;
t123 = 2 * t146;
t122 = -qJ(5) * t64 + t79;
t5 = qJ(5) * t111 - t10;
t6 = t111 * pkin(4) - t9;
t49 = t64 * mrSges(6,1) - t111 * mrSges(6,2);
t121 = mrSges(4,1) * t107 + mrSges(4,2) * t110;
t120 = (Ifges(6,4) - Ifges(5,5) - Ifges(7,5)) * t64 + (-Ifges(7,4) - Ifges(6,5) + Ifges(5,6)) * t63;
t119 = -qJ(5) * t75 + t94;
t118 = (mrSges(5,1) * t109 - mrSges(5,2) * t106) * pkin(3);
t23 = pkin(5) * t75 + t42;
t24 = -t74 * pkin(5) + t44;
t66 = Ifges(5,6) * t74;
t67 = Ifges(7,5) * t75;
t68 = Ifges(6,5) * t74;
t69 = Ifges(5,5) * t75;
t70 = Ifges(7,4) * t74;
t71 = Ifges(6,4) * t75;
t117 = t24 * mrSges(7,2) + t67 + t68 + t69 + t70 + (-mrSges(5,2) + mrSges(6,3)) * t44 - t66 - t71 - t23 * mrSges(7,3) + (mrSges(6,2) - mrSges(5,1)) * t42;
t1 = pkin(5) * t64 + qJ(6) * t111 + t6;
t3 = -pkin(5) * t63 - t5;
t116 = t9 * mrSges(5,1) - t10 * mrSges(5,2) + t6 * mrSges(6,2) + t3 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3) - t120;
t114 = pkin(7) ^ 2;
t112 = qJ(5) ^ 2;
t104 = t111 ^ 2;
t102 = t108 ^ 2;
t98 = t102 * t114;
t97 = Ifges(4,5) * t107;
t96 = Ifges(4,6) * t110;
t87 = Ifges(4,5) * t130;
t86 = -qJ(6) + t93;
t83 = Ifges(4,1) * t107 + t135;
t82 = Ifges(4,2) * t110 + t136;
t81 = -mrSges(4,1) * t110 + mrSges(4,2) * t107;
t78 = -mrSges(4,1) * t111 - mrSges(4,3) * t130;
t77 = mrSges(4,2) * t111 - mrSges(4,3) * t131;
t65 = t121 * t108;
t62 = -Ifges(4,5) * t111 + (Ifges(4,1) * t110 - t136) * t108;
t61 = -Ifges(4,6) * t111 + (-Ifges(4,2) * t107 + t135) * t108;
t52 = -t107 * t142 + t73;
t51 = -mrSges(5,1) * t111 - mrSges(5,3) * t64;
t50 = mrSges(5,2) * t111 - mrSges(5,3) * t63;
t38 = Ifges(5,1) * t75 - Ifges(5,4) * t74;
t37 = Ifges(5,4) * t75 - Ifges(5,2) * t74;
t36 = -Ifges(6,2) * t75 + Ifges(6,6) * t74;
t35 = Ifges(7,2) * t74 + Ifges(7,6) * t75;
t34 = -Ifges(6,6) * t75 + Ifges(6,3) * t74;
t33 = Ifges(7,6) * t74 + Ifges(7,3) * t75;
t32 = -mrSges(6,2) * t74 - mrSges(6,3) * t75;
t31 = mrSges(5,1) * t74 + mrSges(5,2) * t75;
t30 = -mrSges(7,2) * t75 + mrSges(7,3) * t74;
t26 = pkin(4) * t74 + t119;
t22 = -mrSges(6,2) * t63 - mrSges(6,3) * t64;
t21 = mrSges(5,1) * t63 + mrSges(5,2) * t64;
t20 = -mrSges(7,2) * t64 + mrSges(7,3) * t63;
t18 = Ifges(5,1) * t64 - Ifges(5,4) * t63 - Ifges(5,5) * t111;
t17 = Ifges(5,4) * t64 - Ifges(5,2) * t63 - Ifges(5,6) * t111;
t16 = -Ifges(6,4) * t111 - Ifges(6,2) * t64 + Ifges(6,6) * t63;
t15 = -Ifges(7,4) * t111 + Ifges(7,2) * t63 + Ifges(7,6) * t64;
t14 = -Ifges(6,5) * t111 - Ifges(6,6) * t64 + Ifges(6,3) * t63;
t13 = -Ifges(7,5) * t111 + Ifges(7,6) * t63 + Ifges(7,3) * t64;
t12 = t138 * t74 + t119;
t11 = pkin(4) * t63 + t122;
t7 = t138 * t63 + t122;
t2 = [0.2e1 * t1 * t46 + 0.2e1 * t10 * t50 + 0.2e1 * t11 * t22 + 0.2e1 * t7 * t20 + 0.2e1 * t79 * t21 + 0.2e1 * t3 * t48 + 0.2e1 * t5 * t47 + 0.2e1 * t6 * t49 + 0.2e1 * t9 * t51 + 0.2e1 * t52 * t78 + 0.2e1 * t53 * t77 + Ifges(2,3) + (t102 + t104) * mrSges(3,3) * t147 + (-t16 + t18 + t13) * t64 + (t15 + t14 - t17) * t63 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t108 - t107 * t61 + t110 * t62 + t65 * t147) * t108 + m(3) * (pkin(1) ^ 2 + t104 * t114 + t98) + m(4) * (t52 ^ 2 + t53 ^ 2 + t98) + m(5) * (t10 ^ 2 + t79 ^ 2 + t9 ^ 2) + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t3 ^ 2 + t7 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t87 + (Ifges(3,2) + t124) * t111 + (Ifges(4,6) * t107 + (2 * Ifges(3,4))) * t108 + t120) * t111; -pkin(2) * t65 + t11 * t32 + t12 * t20 + t94 * t21 + t26 * t22 + t23 * t46 + t24 * t48 + t7 * t30 + t79 * t31 + (-t47 + t50) * t44 + (t49 - t51) * t42 + (-pkin(7) * mrSges(3,2) - t97 / 0.2e1 - t96 / 0.2e1 - t69 / 0.2e1 + t66 / 0.2e1 - t70 / 0.2e1 - t67 / 0.2e1 + t71 / 0.2e1 - t68 / 0.2e1 + Ifges(3,6)) * t111 + (-t36 / 0.2e1 + t33 / 0.2e1 + t38 / 0.2e1) * t64 + (t34 / 0.2e1 + t35 / 0.2e1 - t37 / 0.2e1) * t63 + (pkin(8) * t77 + t53 * mrSges(4,3) + t61 / 0.2e1) * t110 + (-pkin(8) * t78 - t52 * mrSges(4,3) + t62 / 0.2e1) * t107 + (t110 * t83 / 0.2e1 - t107 * t82 / 0.2e1 + Ifges(3,5) + (t81 - mrSges(3,1)) * pkin(7)) * t108 + m(4) * (-pkin(2) * t99 + (-t107 * t52 + t110 * t53) * pkin(8)) + m(5) * (t10 * t44 - t42 * t9 + t79 * t94) + m(6) * (t11 * t26 + t42 * t6 - t44 * t5) + m(7) * (t1 * t23 + t12 * t7 + t24 * t3) + (-t9 * mrSges(5,3) + t1 * mrSges(7,1) + t6 * mrSges(6,1) + t18 / 0.2e1 + t13 / 0.2e1 - t16 / 0.2e1) * t75 + (-t10 * mrSges(5,3) + t15 / 0.2e1 + t14 / 0.2e1 - t17 / 0.2e1 - t3 * mrSges(7,1) + t5 * mrSges(6,1)) * t74; -0.2e1 * pkin(2) * t81 + t107 * t83 + t110 * t82 + 0.2e1 * t12 * t30 + 0.2e1 * t26 * t32 + 0.2e1 * t94 * t31 + Ifges(3,3) + m(4) * (t128 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t94 ^ 2 + t126) + m(6) * (t26 ^ 2 + t126) + m(7) * (t12 ^ 2 + t23 ^ 2 + t24 ^ 2) + (0.2e1 * t23 * mrSges(7,1) + t33 - t36 + t38) * t75 + (-0.2e1 * t24 * mrSges(7,1) + t34 + t35 - t37) * t74 + 0.2e1 * t128 * pkin(8) * mrSges(4,3) + 0.2e1 * (t42 * t75 - t44 * t74) * (mrSges(6,1) + mrSges(5,3)); (m(5) * (t10 * t106 + t109 * t9) + t109 * t51 + t106 * t50) * pkin(3) - Ifges(4,6) * t131 + t87 + t93 * t49 + t86 * t46 + t52 * mrSges(4,1) - t53 * mrSges(4,2) + t116 - t124 * t111 + m(7) * (t1 * t86 + t3 * t90) + m(6) * (-t5 * t90 + t6 * t93) + t137 * t90; m(7) * (t23 * t86 + t24 * t90) + m(6) * (t42 * t93 + t44 * t90) - t121 * pkin(8) + (t75 * t86 - t141) * mrSges(7,1) + (t75 * t93 - t141) * mrSges(6,1) + (m(5) * (t106 * t44 - t109 * t42) + (-t106 * t74 - t109 * t75) * mrSges(5,3)) * pkin(3) + t96 + t97 + t117; 0.2e1 * t93 * mrSges(6,2) + t86 * t144 + t90 * t123 + 0.2e1 * t118 + m(6) * (t93 ^ 2 + t145) + m(7) * (t86 ^ 2 + t145) + m(5) * (t106 ^ 2 + t109 ^ 2) * pkin(3) ^ 2 + t124; t137 * qJ(5) + m(7) * (qJ(5) * t3 - t1 * t138) + m(6) * (-pkin(4) * t6 - qJ(5) * t5) - t127 * t111 - t138 * t46 - pkin(4) * t49 + t116; m(6) * (-pkin(4) * t42 + qJ(5) * t44) + m(7) * (qJ(5) * t24 - t138 * t23) + (-t138 * t75 - t134) * mrSges(7,1) + (-pkin(4) * t75 - t134) * mrSges(6,1) + t117; t118 + (t138 - t86) * mrSges(7,3) + (t93 - pkin(4)) * mrSges(6,2) + m(6) * (-pkin(4) * t93 + t133) + m(7) * (-t138 * t86 + t133) + t127 + t146 * (qJ(5) + t90); -0.2e1 * pkin(4) * mrSges(6,2) - t138 * t144 + qJ(5) * t123 + m(6) * (pkin(4) ^ 2 + t112) + m(7) * (t138 ^ 2 + t112) + t127; m(6) * t6 + m(7) * t1 + t46 + t49; (mrSges(6,1) + mrSges(7,1)) * t75 + m(6) * t42 + m(7) * t23; m(6) * t93 + m(7) * t86 + t139; -m(6) * pkin(4) - m(7) * t138 + t139; m(6) + m(7); m(7) * t3 + t48; m(7) * t24 - t74 * mrSges(7,1); m(7) * t90 + mrSges(7,2); m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
