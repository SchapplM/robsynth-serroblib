% Calculate joint inertia matrix for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2018-11-23 16:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:10
% EndTime: 2018-11-23 16:45:11
% DurationCPUTime: 0.94s
% Computational Cost: add. (1743->249), mult. (3335->341), div. (0->0), fcn. (3630->8), ass. (0->94)
t134 = Ifges(7,2) + Ifges(6,3);
t133 = m(7) + m(6);
t132 = -m(7) * pkin(5) - mrSges(7,1);
t116 = -qJ(3) - pkin(7);
t98 = sin(qJ(2));
t108 = t116 * t98;
t99 = cos(qJ(2));
t81 = t116 * t99;
t94 = sin(pkin(9));
t96 = cos(pkin(9));
t53 = -t96 * t108 - t81 * t94;
t131 = t53 ^ 2;
t130 = 0.2e1 * t53;
t87 = -pkin(2) * t99 - pkin(1);
t129 = 0.2e1 * t87;
t93 = sin(pkin(10));
t128 = t93 / 0.2e1;
t95 = cos(pkin(10));
t127 = t95 / 0.2e1;
t126 = pkin(2) * t94;
t125 = pkin(2) * t96;
t83 = qJ(4) + t126;
t124 = pkin(8) + t83;
t123 = cos(qJ(5));
t74 = t94 * t99 + t96 * t98;
t120 = t74 * t93;
t71 = t94 * t98 - t96 * t99;
t43 = pkin(3) * t71 - qJ(4) * t74 + t87;
t55 = t94 * t108 - t96 * t81;
t17 = t93 * t43 + t95 * t55;
t13 = -pkin(8) * t120 + t17;
t119 = t74 * t95;
t16 = t95 * t43 - t55 * t93;
t7 = pkin(4) * t71 - pkin(8) * t119 + t16;
t97 = sin(qJ(5));
t4 = t123 * t13 + t97 * t7;
t122 = Ifges(5,4) * t93;
t121 = Ifges(5,4) * t95;
t118 = t97 * t93;
t109 = t123 * t93;
t75 = t97 * t95 + t109;
t34 = t75 * t74;
t18 = -mrSges(6,2) * t71 - mrSges(6,3) * t34;
t21 = -mrSges(7,2) * t34 + mrSges(7,3) * t71;
t115 = t18 + t21;
t103 = t123 * t95 - t118;
t35 = t103 * t74;
t19 = mrSges(6,1) * t71 - mrSges(6,3) * t35;
t20 = -t71 * mrSges(7,1) + t35 * mrSges(7,2);
t114 = -t19 + t20;
t42 = mrSges(5,1) * t120 + mrSges(5,2) * t119;
t113 = t93 ^ 2 + t95 ^ 2;
t112 = t98 ^ 2 + t99 ^ 2;
t60 = t124 * t95;
t39 = t124 * t109 + t60 * t97;
t41 = -t124 * t118 + t123 * t60;
t111 = t39 ^ 2 + t41 ^ 2;
t86 = -pkin(3) - t125;
t78 = -t95 * mrSges(5,1) + t93 * mrSges(5,2);
t15 = t34 * mrSges(6,1) + t35 * mrSges(6,2);
t47 = -mrSges(6,1) * t103 + t75 * mrSges(6,2);
t14 = t34 * mrSges(7,1) - t35 * mrSges(7,3);
t46 = -mrSges(7,1) * t103 - t75 * mrSges(7,3);
t25 = pkin(4) * t120 + t53;
t106 = pkin(5) * t103 + qJ(6) * t75;
t105 = -t16 * t93 + t17 * t95;
t77 = -pkin(4) * t95 + t86;
t104 = t134 * t71 + (Ifges(7,4) + Ifges(6,5)) * t35 + (-Ifges(6,6) + Ifges(7,6)) * t34;
t3 = t123 * t7 - t97 * t13;
t102 = -t46 - t47;
t80 = Ifges(5,1) * t93 + t121;
t79 = Ifges(5,2) * t95 + t122;
t69 = Ifges(7,4) * t75;
t68 = Ifges(6,5) * t75;
t66 = Ifges(6,6) * t103;
t65 = Ifges(7,6) * t103;
t62 = t74 * mrSges(4,2);
t51 = Ifges(6,1) * t75 + Ifges(6,4) * t103;
t50 = Ifges(7,1) * t75 - Ifges(7,5) * t103;
t49 = Ifges(6,4) * t75 + Ifges(6,2) * t103;
t48 = Ifges(7,5) * t75 - Ifges(7,3) * t103;
t45 = mrSges(5,1) * t71 - mrSges(5,3) * t119;
t44 = -mrSges(5,2) * t71 - mrSges(5,3) * t120;
t33 = -t106 + t77;
t24 = t71 * Ifges(5,5) + (Ifges(5,1) * t95 - t122) * t74;
t23 = t71 * Ifges(5,6) + (-Ifges(5,2) * t93 + t121) * t74;
t11 = Ifges(6,1) * t35 - Ifges(6,4) * t34 + Ifges(6,5) * t71;
t10 = Ifges(7,1) * t35 + Ifges(7,4) * t71 + Ifges(7,5) * t34;
t9 = Ifges(6,4) * t35 - Ifges(6,2) * t34 + Ifges(6,6) * t71;
t8 = Ifges(7,5) * t35 + Ifges(7,6) * t71 + Ifges(7,3) * t34;
t5 = pkin(5) * t34 - qJ(6) * t35 + t25;
t2 = -t71 * pkin(5) - t3;
t1 = qJ(6) * t71 + t4;
t6 = [-0.2e1 * pkin(1) * (-t99 * mrSges(3,1) + t98 * mrSges(3,2)) + t99 * (Ifges(3,4) * t98 + Ifges(3,2) * t99) + t98 * (Ifges(3,1) * t98 + Ifges(3,4) * t99) + t62 * t129 + t42 * t130 + 0.2e1 * t17 * t44 + 0.2e1 * t16 * t45 + 0.2e1 * t25 * t15 + 0.2e1 * t5 * t14 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t1 * t21 + Ifges(2,3) + (t10 + t11) * t35 + (-t9 + t8) * t34 + 0.2e1 * t112 * pkin(7) * mrSges(3,3) + (mrSges(4,3) * t130 + Ifges(4,1) * t74 - t93 * t23 + t95 * t24) * t74 + (mrSges(4,1) * t129 - 0.2e1 * t55 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t71 + (Ifges(5,5) * t95 - Ifges(5,6) * t93 - (2 * Ifges(4,4))) * t74 + t104) * t71 + m(3) * (t112 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t55 ^ 2 + t87 ^ 2 + t131) + m(5) * (t16 ^ 2 + t17 ^ 2 + t131) + m(6) * (t25 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); (-t98 * mrSges(3,1) - t99 * mrSges(3,2)) * pkin(7) - (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t8 / 0.2e1 - t9 / 0.2e1) * t103 + (-mrSges(4,3) * t126 + Ifges(5,5) * t128 + Ifges(5,6) * t127 - Ifges(4,6) + t68 / 0.2e1 + t66 / 0.2e1 + t69 / 0.2e1 - t65 / 0.2e1) * t71 + (t95 * t44 - t93 * t45) * t83 + (t78 - mrSges(4,1)) * t53 + Ifges(3,5) * t98 + Ifges(3,6) * t99 + t77 * t15 + t86 * t42 + t5 * t46 + t25 * t47 - t55 * mrSges(4,2) + (-t93 * t79 / 0.2e1 + t80 * t127 - mrSges(4,3) * t125 + Ifges(4,5)) * t74 + t114 * t39 + t115 * t41 + t105 * mrSges(5,3) + m(5) * (t105 * t83 + t53 * t86) + t33 * t14 + t23 * t127 + t24 * t128 + (t48 / 0.2e1 - t49 / 0.2e1) * t34 + (t50 / 0.2e1 + t51 / 0.2e1) * t35 + m(4) * (-t53 * t96 + t55 * t94) * pkin(2) + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1) * t75 + m(6) * (t25 * t77 - t3 * t39 + t4 * t41) + m(7) * (t1 * t41 + t2 * t39 + t33 * t5); 0.2e1 * t33 * t46 + 0.2e1 * t77 * t47 + 0.2e1 * t86 * t78 + t95 * t79 + t93 * t80 + Ifges(3,3) + Ifges(4,3) + (t50 + t51) * t75 - (t48 - t49) * t103 + m(7) * (t33 ^ 2 + t111) + m(6) * (t77 ^ 2 + t111) + m(5) * (t113 * t83 ^ 2 + t86 ^ 2) + m(4) * (t94 ^ 2 + t96 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t96 - mrSges(4,2) * t94) * pkin(2) + 0.2e1 * t113 * t83 * mrSges(5,3) + 0.2e1 * (t103 * t41 + t39 * t75) * (mrSges(7,2) + mrSges(6,3)); t71 * mrSges(4,1) + t93 * t44 + t95 * t45 + t62 + t115 * t75 - t114 * t103 + m(7) * (t1 * t75 - t103 * t2) + m(6) * (t103 * t3 + t4 * t75) + m(5) * (t16 * t95 + t17 * t93) + m(4) * t87; t133 * (-t103 * t39 + t75 * t41); m(5) * t113 + m(4) + t133 * (t103 ^ 2 + t75 ^ 2); m(5) * t53 + m(6) * t25 + m(7) * t5 + t14 + t15 + t42; m(5) * t86 + m(6) * t77 + m(7) * t33 - t102 + t78; 0; m(5) + t133; qJ(6) * t21 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t20 + t104; t66 + t69 - t65 + t68 + (-pkin(5) * t75 + qJ(6) * t103) * mrSges(7,2) + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t41 + (-mrSges(6,1) + t132) * t39; m(7) * t106 + t102; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t134; m(7) * t2 + t20; m(7) * t39 + t75 * mrSges(7,2); -m(7) * t103; 0; t132; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
