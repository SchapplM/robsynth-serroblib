% Calculate Gravitation load on the joints for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:22
% EndTime: 2019-03-10 03:33:23
% DurationCPUTime: 0.87s
% Computational Cost: add. (700->157), mult. (579->206), div. (0->0), fcn. (521->12), ass. (0->91)
t125 = rSges(6,3) + pkin(10);
t50 = qJ(2) + qJ(3);
t46 = qJ(4) + t50;
t38 = sin(t46);
t39 = cos(t46);
t54 = cos(qJ(5));
t40 = pkin(5) * t54 + pkin(4);
t124 = t38 * rSges(7,3) + t39 * t40;
t55 = cos(qJ(2));
t47 = t55 * pkin(2);
t43 = sin(t50);
t45 = cos(t50);
t81 = t45 * rSges(4,1) - rSges(4,2) * t43;
t123 = t47 + t81;
t122 = t39 * rSges(5,1) - rSges(5,2) * t38;
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t121 = g(1) * t56 + g(2) * t53;
t69 = t39 * pkin(4) + t125 * t38;
t58 = -pkin(8) - pkin(7);
t49 = qJ(5) + qJ(6);
t42 = sin(t49);
t101 = t42 * t53;
t44 = cos(t49);
t98 = t44 * t56;
t5 = t39 * t101 + t98;
t100 = t42 * t56;
t99 = t44 * t53;
t6 = -t39 * t99 + t100;
t120 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t39 * t100 + t99;
t8 = t39 * t98 + t101;
t119 = t7 * rSges(7,1) - t8 * rSges(7,2);
t52 = sin(qJ(2));
t118 = pkin(2) * t52;
t117 = pkin(3) * t43;
t51 = sin(qJ(5));
t116 = pkin(5) * t51;
t113 = g(3) * t38;
t112 = rSges(3,3) + pkin(7);
t111 = rSges(6,1) * t54;
t110 = rSges(7,1) * t44;
t107 = rSges(6,2) * t51;
t106 = rSges(7,2) * t42;
t57 = -pkin(11) - pkin(10);
t105 = t38 * t57;
t104 = t39 * t53;
t103 = t39 * t56;
t102 = t39 * t57;
t97 = t51 * t53;
t96 = t51 * t56;
t95 = t53 * t54;
t94 = t54 * t56;
t93 = rSges(4,3) - t58;
t48 = -pkin(9) + t58;
t92 = rSges(5,3) - t48;
t87 = t38 * t106;
t91 = rSges(7,3) * t104 + t53 * t87;
t90 = rSges(7,3) * t103 + t56 * t87;
t37 = pkin(3) * t45;
t89 = t37 + t47;
t88 = t38 * t107;
t86 = t125 * t104 + t53 * t88;
t85 = t125 * t103 + t56 * t88;
t84 = -pkin(4) - t111;
t83 = -t48 + t116;
t82 = -t40 - t110;
t79 = t37 + t122;
t78 = rSges(3,1) * t55 - rSges(3,2) * t52;
t76 = -rSges(4,1) * t43 - rSges(4,2) * t45;
t74 = -rSges(5,1) * t38 - rSges(5,2) * t39;
t73 = -rSges(7,1) * t42 - rSges(7,2) * t44;
t72 = pkin(1) + t78;
t11 = -t39 * t96 + t95;
t9 = t39 * t97 + t94;
t71 = pkin(1) + t123;
t70 = t124 + (-t106 + t110) * t39;
t68 = t69 + (-t107 + t111) * t39;
t65 = -t105 + t124;
t64 = g(1) * t90 + g(2) * t91;
t63 = t37 + t68;
t62 = t70 - t105;
t60 = t121 * t84 * t38;
t23 = -t117 - t118;
t22 = pkin(1) + t89;
t17 = t56 * t23;
t16 = t53 * t23;
t13 = t56 * t22;
t12 = t39 * t94 + t97;
t10 = -t39 * t95 + t96;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t53 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 - rSges(2,2) * t53)) - m(3) * ((g(1) * t112 + g(2) * t72) * t56 + (-g(1) * t72 + g(2) * t112) * t53) - m(4) * ((g(1) * t93 + g(2) * t71) * t56 + (-g(1) * t71 + g(2) * t93) * t53) - m(5) * (g(2) * t13 + (g(1) * t92 + g(2) * t122) * t56 + (g(1) * (-t22 - t122) + g(2) * t92) * t53) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t11 + t13) + (-g(1) * t48 + g(2) * t69) * t56 + (g(1) * (-t22 - t69) - g(2) * t48) * t53) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t13) + (g(1) * t83 + g(2) * t65) * t56 + (g(1) * (-t22 - t65) + g(2) * t83) * t53) -m(3) * (g(3) * t78 + t121 * (-rSges(3,1) * t52 - rSges(3,2) * t55)) - m(4) * (g(3) * t123 + t121 * (t76 - t118)) - m(5) * (g(1) * (t56 * t74 + t17) + g(2) * (t53 * t74 + t16) + g(3) * (t47 + t79)) - m(6) * (g(1) * (t17 + t85) + g(2) * (t16 + t86) + g(3) * (t47 + t63) + t60) - m(7) * (g(1) * (-t102 * t56 + t17 + t90) + g(2) * (-t102 * t53 + t16 + t91) + g(3) * (t70 + t89) + (-g(3) * t57 + t121 * t82) * t38) -m(4) * (g(3) * t81 + t121 * t76) - m(5) * (g(3) * t79 + t121 * (t74 - t117)) - m(6) * (g(1) * (-t117 * t56 + t85) + g(2) * (-t117 * t53 + t86) + g(3) * t63 + t60) - m(7) * (g(3) * (t37 + t62) + t64 + t121 * (t38 * t82 - t102 - t117)) -m(5) * g(3) * t122 - m(6) * (g(1) * t85 + g(2) * t86 + g(3) * t68) - m(7) * (g(3) * t62 + t64) + t121 * ((m(5) * rSges(5,2) + m(7) * t57) * t39 + (m(5) * rSges(5,1) - m(6) * t84 - m(7) * t82) * t38) -m(6) * (g(1) * (rSges(6,1) * t11 - rSges(6,2) * t12) + g(2) * (-rSges(6,1) * t9 + rSges(6,2) * t10)) - m(7) * (g(1) * (pkin(5) * t11 + t119) + g(2) * (-pkin(5) * t9 + t120)) + (-m(6) * (-rSges(6,1) * t51 - rSges(6,2) * t54) - m(7) * (t73 - t116)) * t113, -m(7) * (g(1) * t119 + g(2) * t120 + t73 * t113)];
taug  = t1(:);
