% Calculate joint inertia matrix for
% S6RRRRPP5
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
% Datum: 2018-11-23 18:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:07:09
% EndTime: 2018-11-23 18:07:10
% DurationCPUTime: 1.09s
% Computational Cost: add. (1518->328), mult. (2960->422), div. (0->0), fcn. (2780->6), ass. (0->112)
t145 = 2 * pkin(7);
t144 = (mrSges(6,3) + mrSges(7,2));
t143 = -2 * mrSges(7,1);
t112 = -pkin(4) - pkin(5);
t142 = -pkin(9) - pkin(8);
t111 = cos(qJ(2));
t141 = pkin(7) * t111;
t108 = sin(qJ(2));
t100 = t108 * pkin(7);
t140 = -mrSges(6,1) - mrSges(7,1);
t138 = mrSges(6,2) - mrSges(7,3);
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t130 = t108 * t110;
t80 = -pkin(2) * t111 - pkin(8) * t108 - pkin(1);
t71 = t110 * t80;
t27 = -pkin(9) * t130 + t71 + (-pkin(7) * t107 - pkin(3)) * t111;
t131 = t107 * t108;
t53 = t107 * t80 + t110 * t141;
t39 = -pkin(9) * t131 + t53;
t10 = t106 * t27 + t109 * t39;
t73 = t106 * t110 + t107 * t109;
t62 = t73 * t108;
t46 = -mrSges(6,2) * t62 - mrSges(6,3) * t111;
t47 = -mrSges(7,2) * t111 + mrSges(7,3) * t62;
t137 = t46 + t47;
t72 = t106 * t107 - t109 * t110;
t63 = t72 * t108;
t51 = t111 * mrSges(6,1) - t63 * mrSges(6,2);
t126 = t142 * t107;
t84 = t142 * t110;
t45 = t106 * t126 - t109 * t84;
t136 = Ifges(4,4) * t107;
t135 = Ifges(4,4) * t110;
t134 = t111 * Ifges(7,5);
t133 = t111 * Ifges(7,6);
t79 = pkin(3) * t131 + t100;
t129 = t107 ^ 2 + t110 ^ 2;
t43 = -t106 * t84 - t109 * t126;
t128 = t43 ^ 2 + t45 ^ 2;
t127 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t94 = -pkin(3) * t110 - pkin(2);
t93 = -pkin(3) * t109 - pkin(4);
t23 = -t62 * mrSges(7,1) - t63 * mrSges(7,2);
t31 = -t72 * mrSges(7,1) + t73 * mrSges(7,2);
t49 = t111 * mrSges(7,1) + t63 * mrSges(7,3);
t9 = -t106 * t39 + t109 * t27;
t125 = Ifges(4,3) + t127;
t124 = 2 * t144;
t123 = -qJ(5) * t63 - t79;
t5 = -qJ(5) * t111 + t10;
t6 = t111 * pkin(4) - t9;
t122 = mrSges(4,1) * t107 + mrSges(4,2) * t110;
t121 = qJ(5) * t73 - t94;
t120 = (mrSges(5,1) * t109 - mrSges(5,2) * t106) * pkin(3);
t119 = (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t63 + (Ifges(5,6) - Ifges(6,6) + Ifges(7,6)) * t62;
t20 = -qJ(6) * t73 + t43;
t21 = qJ(6) * t72 + t45;
t66 = Ifges(6,6) * t72;
t67 = Ifges(5,6) * t72;
t68 = Ifges(5,5) * t73;
t69 = Ifges(6,4) * t73;
t118 = -t20 * mrSges(7,1) + t21 * mrSges(7,2) + t66 - t67 + t68 + t69 + (-mrSges(5,2) + mrSges(6,3)) * t45 + (-mrSges(5,1) - mrSges(6,1)) * t43;
t1 = pkin(5) * t111 + qJ(6) * t63 + t6;
t3 = qJ(6) * t62 + t5;
t117 = t9 * mrSges(5,1) - t6 * mrSges(6,1) - t1 * mrSges(7,1) - t10 * mrSges(5,2) + t3 * mrSges(7,2) + t5 * mrSges(6,3) - t119;
t115 = pkin(7) ^ 2;
t113 = qJ(5) ^ 2;
t105 = t111 ^ 2;
t103 = t108 ^ 2;
t99 = t103 * t115;
t98 = Ifges(4,5) * t107;
t97 = Ifges(4,6) * t110;
t91 = pkin(3) * t106 + qJ(5);
t90 = -pkin(5) + t93;
t89 = t91 ^ 2;
t86 = Ifges(4,5) * t130;
t85 = qJ(5) * t91;
t83 = Ifges(4,1) * t107 + t135;
t82 = Ifges(4,2) * t110 + t136;
t81 = -mrSges(4,1) * t110 + mrSges(4,2) * t107;
t78 = -mrSges(4,1) * t111 - mrSges(4,3) * t130;
t77 = mrSges(4,2) * t111 - mrSges(4,3) * t131;
t64 = t122 * t108;
t61 = -Ifges(4,5) * t111 + (Ifges(4,1) * t110 - t136) * t108;
t60 = -Ifges(4,6) * t111 + (-Ifges(4,2) * t107 + t135) * t108;
t52 = -t107 * t141 + t71;
t50 = -mrSges(5,1) * t111 + mrSges(5,3) * t63;
t48 = mrSges(5,2) * t111 - mrSges(5,3) * t62;
t38 = Ifges(5,1) * t73 - Ifges(5,4) * t72;
t37 = Ifges(6,1) * t73 + Ifges(6,5) * t72;
t36 = Ifges(7,1) * t73 + Ifges(7,4) * t72;
t35 = Ifges(5,4) * t73 - Ifges(5,2) * t72;
t34 = Ifges(7,4) * t73 + Ifges(7,2) * t72;
t33 = Ifges(6,5) * t73 + Ifges(6,3) * t72;
t32 = mrSges(5,1) * t72 + mrSges(5,2) * t73;
t30 = mrSges(6,1) * t72 - mrSges(6,3) * t73;
t26 = pkin(4) * t72 - t121;
t24 = mrSges(5,1) * t62 - mrSges(5,2) * t63;
t22 = mrSges(6,1) * t62 + mrSges(6,3) * t63;
t18 = -Ifges(5,1) * t63 - Ifges(5,4) * t62 - Ifges(5,5) * t111;
t17 = -Ifges(6,1) * t63 - Ifges(6,4) * t111 + Ifges(6,5) * t62;
t16 = -Ifges(7,1) * t63 + Ifges(7,4) * t62 + t134;
t15 = -Ifges(5,4) * t63 - Ifges(5,2) * t62 - Ifges(5,6) * t111;
t14 = -Ifges(7,4) * t63 + Ifges(7,2) * t62 + t133;
t13 = -Ifges(6,5) * t63 - Ifges(6,6) * t111 + Ifges(6,3) * t62;
t12 = t112 * t72 + t121;
t11 = pkin(4) * t62 - t123;
t7 = t112 * t62 + t123;
t2 = [0.2e1 * t1 * t49 + 0.2e1 * t10 * t48 + 0.2e1 * t11 * t22 + 0.2e1 * t7 * t23 + 0.2e1 * t79 * t24 + 0.2e1 * t3 * t47 + 0.2e1 * t5 * t46 + 0.2e1 * t9 * t50 + 0.2e1 * t6 * t51 + 0.2e1 * t52 * t78 + 0.2e1 * t53 * t77 + Ifges(2,3) + (t103 + t105) * mrSges(3,3) * t145 - (t16 + t17 + t18) * t63 + (t13 + t14 - t15) * t62 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t108 - t107 * t60 + t110 * t61 + t64 * t145) * t108 + m(3) * (pkin(1) ^ 2 + t105 * t115 + t99) + m(4) * (t52 ^ 2 + t53 ^ 2 + t99) + m(5) * (t10 ^ 2 + t79 ^ 2 + t9 ^ 2) + m(7) * (t1 ^ 2 + t3 ^ 2 + t7 ^ 2) + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t86 + (Ifges(3,2) + t125) * t111 + (Ifges(4,6) * t107 + (2 * Ifges(3,4))) * t108 + t119) * t111; -pkin(2) * t64 + t11 * t30 + t12 * t23 + t20 * t49 + t21 * t47 + t26 * t22 + t94 * t24 + t7 * t31 + t79 * t32 + (t46 + t48) * t45 + (t51 - t50) * t43 + (-pkin(7) * mrSges(3,2) - t98 / 0.2e1 - t97 / 0.2e1 - t68 / 0.2e1 + t67 / 0.2e1 - t69 / 0.2e1 - t66 / 0.2e1 + Ifges(3,6)) * t111 - (t36 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1) * t63 + (t34 / 0.2e1 + t33 / 0.2e1 - t35 / 0.2e1) * t62 + (pkin(8) * t77 + t53 * mrSges(4,3) + t60 / 0.2e1) * t110 + (-pkin(8) * t78 - t52 * mrSges(4,3) + t61 / 0.2e1) * t107 + (t134 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1 + t18 / 0.2e1 - t9 * mrSges(5,3) - t1 * mrSges(7,3) + t6 * mrSges(6,2)) * t73 + (t133 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 - t15 / 0.2e1 + t3 * mrSges(7,3) - t5 * mrSges(6,2) - t10 * mrSges(5,3)) * t72 + (-t107 * t82 / 0.2e1 + t110 * t83 / 0.2e1 + Ifges(3,5) + (t81 - mrSges(3,1)) * pkin(7)) * t108 + m(4) * (-pkin(2) * t100 + (-t52 * t107 + t53 * t110) * pkin(8)) + m(5) * (t10 * t45 - t43 * t9 + t79 * t94) + m(6) * (t11 * t26 + t43 * t6 + t45 * t5) + m(7) * (t1 * t20 + t12 * t7 + t21 * t3); -0.2e1 * pkin(2) * t81 + t107 * t83 + t110 * t82 + 0.2e1 * t12 * t31 + 0.2e1 * t26 * t30 + 0.2e1 * t94 * t32 + Ifges(3,3) + m(4) * (t129 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t94 ^ 2 + t128) + m(7) * (t12 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t26 ^ 2 + t128) + (-0.2e1 * t20 * mrSges(7,3) + t36 + t37 + t38) * t73 + (0.2e1 * t21 * mrSges(7,3) + t33 + t34 - t35) * t72 + 0.2e1 * t129 * pkin(8) * mrSges(4,3) + 0.2e1 * (t43 * t73 - t45 * t72) * (mrSges(6,2) + mrSges(5,3)); t137 * t91 + m(7) * (t1 * t90 + t3 * t91) + m(6) * (t5 * t91 + t6 * t93) - t125 * t111 + (m(5) * (t10 * t106 + t109 * t9) + t109 * t50 + t106 * t48) * pkin(3) - Ifges(4,6) * t131 + t86 + t90 * t49 + t93 * t51 + t52 * mrSges(4,1) - t53 * mrSges(4,2) + t117; (-t138 * t91 - Ifges(7,6)) * t72 + (mrSges(6,2) * t93 - mrSges(7,3) * t90 - Ifges(7,5)) * t73 + m(7) * (t20 * t90 + t21 * t91) + m(6) * (t43 * t93 + t45 * t91) - t122 * pkin(8) + (m(5) * (t106 * t45 - t109 * t43) + (-t106 * t72 - t109 * t73) * mrSges(5,3)) * pkin(3) + t97 + t98 + t118; -0.2e1 * t93 * mrSges(6,1) + t90 * t143 + t91 * t124 + 0.2e1 * t120 + m(6) * (t93 ^ 2 + t89) + m(7) * (t90 ^ 2 + t89) + m(5) * (t106 ^ 2 + t109 ^ 2) * pkin(3) ^ 2 + t125; t137 * qJ(5) + m(7) * (qJ(5) * t3 + t1 * t112) + m(6) * (-pkin(4) * t6 + qJ(5) * t5) - t127 * t111 + t112 * t49 - pkin(4) * t51 + t117; (-t138 * qJ(5) - Ifges(7,6)) * t72 + m(6) * (-pkin(4) * t43 + qJ(5) * t45) + m(7) * (qJ(5) * t21 + t112 * t20) + (-mrSges(6,2) * pkin(4) - mrSges(7,3) * t112 - Ifges(7,5)) * t73 + t118; t120 + (-t112 - t90) * mrSges(7,1) + (-t93 + pkin(4)) * mrSges(6,1) + m(6) * (-pkin(4) * t93 + t85) + m(7) * (t112 * t90 + t85) + t127 + t144 * (t91 + qJ(5)); 0.2e1 * pkin(4) * mrSges(6,1) + t112 * t143 + qJ(5) * t124 + m(6) * (pkin(4) ^ 2 + t113) + m(7) * (t112 ^ 2 + t113) + t127; m(6) * t6 + m(7) * t1 + t49 + t51; m(6) * t43 + m(7) * t20 + t138 * t73; m(6) * t93 + m(7) * t90 + t140; -m(6) * pkin(4) + m(7) * t112 + t140; m(6) + m(7); m(7) * t7 + t23; m(7) * t12 + t31; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
