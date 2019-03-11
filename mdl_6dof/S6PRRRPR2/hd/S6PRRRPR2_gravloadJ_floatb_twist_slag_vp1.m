% Calculate Gravitation load on the joints for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:34
% EndTime: 2019-03-08 23:05:36
% DurationCPUTime: 0.81s
% Computational Cost: add. (706->137), mult. (1098->195), div. (0->0), fcn. (1278->14), ass. (0->69)
t127 = rSges(7,3) + pkin(10) + qJ(5);
t60 = sin(pkin(6));
t64 = sin(qJ(2));
t104 = t60 * t64;
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t94 = cos(pkin(6));
t126 = -t104 * t63 + t65 * t94;
t125 = rSges(6,3) + qJ(5);
t66 = cos(qJ(2));
t92 = sin(pkin(11));
t82 = t94 * t92;
t93 = cos(pkin(11));
t42 = -t64 * t82 + t66 * t93;
t89 = t60 * t92;
t124 = -t42 * t63 + t65 * t89;
t59 = sin(pkin(12));
t61 = cos(pkin(12));
t123 = rSges(6,1) * t61 - rSges(6,2) * t59 + pkin(4);
t57 = pkin(12) + qJ(6);
t53 = sin(t57);
t54 = cos(t57);
t122 = -rSges(7,1) * t54 + rSges(7,2) * t53 - pkin(5) * t61 - pkin(4);
t121 = rSges(6,1) * t59 + rSges(6,2) * t61;
t83 = t94 * t93;
t39 = t64 * t92 - t66 * t83;
t115 = g(2) * t39;
t41 = t64 * t93 + t66 * t82;
t120 = g(1) * t41 + t115;
t58 = qJ(3) + qJ(4);
t55 = sin(t58);
t56 = cos(t58);
t119 = t123 * t56 + t125 * t55;
t118 = -t122 * t56 + t127 * t55;
t117 = -m(6) - m(7);
t40 = t64 * t83 + t66 * t92;
t114 = g(2) * t40;
t103 = t60 * t66;
t52 = pkin(3) * t65 + pkin(2);
t113 = g(3) * t52 * t103;
t112 = g(3) * t60;
t111 = rSges(4,3) + pkin(8);
t90 = t60 * t93;
t20 = t40 * t55 + t56 * t90;
t21 = t40 * t56 - t55 * t90;
t100 = -rSges(5,1) * t20 - rSges(5,2) * t21;
t22 = t42 * t55 - t56 * t89;
t23 = t42 * t56 + t55 * t89;
t99 = -rSges(5,1) * t22 - rSges(5,2) * t23;
t67 = -pkin(9) - pkin(8);
t98 = -t39 * t52 - t40 * t67;
t97 = -t41 * t52 - t42 * t67;
t35 = t104 * t55 - t56 * t94;
t36 = t104 * t56 + t55 * t94;
t96 = -rSges(5,1) * t35 - rSges(5,2) * t36;
t87 = t124 * pkin(3);
t85 = rSges(5,1) * t56 - rSges(5,2) * t55;
t84 = t126 * pkin(3);
t81 = rSges(4,1) * t65 - rSges(4,2) * t63 + pkin(2);
t78 = t122 * t20 + t127 * t21;
t77 = t122 * t22 + t127 * t23;
t76 = t122 * t35 + t127 * t36;
t75 = rSges(7,1) * t53 + rSges(7,2) * t54 + pkin(5) * t59;
t74 = -t40 * t63 - t65 * t90;
t73 = -t123 * t20 + t125 * t21;
t72 = -t123 * t22 + t125 * t23;
t71 = -t123 * t35 + t125 * t36;
t70 = t74 * pkin(3);
t1 = [(-m(2) - m(3) - m(4) - m(5) + t117) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t41 - rSges(3,2) * t42) + g(2) * (-rSges(3,1) * t39 - rSges(3,2) * t40) + (rSges(3,1) * t66 - rSges(3,2) * t64) * t112) - m(4) * (g(1) * (t111 * t42 - t41 * t81) + t111 * t114 - t81 * t115 + (t111 * t64 + t66 * t81) * t112) - m(5) * (g(1) * (rSges(5,3) * t42 - t41 * t85 + t97) + g(2) * (rSges(5,3) * t40 - t39 * t85 + t98) + t113 + (t85 * t66 + (rSges(5,3) - t67) * t64) * t112) - m(6) * (g(1) * (t121 * t42 + t97) + g(2) * (t121 * t40 + t98) + t113 + ((-t67 + t121) * t64 + t119 * t66) * t112 - t120 * t119) - m(7) * (g(2) * t98 + t113 + t75 * t114 + ((-t67 + t75) * t64 + t118 * t66) * t112 - t120 * t118 + (t42 * t75 + t97) * g(1)) -m(4) * (g(1) * (t124 * rSges(4,1) + (-t42 * t65 - t63 * t89) * rSges(4,2)) + g(2) * (t74 * rSges(4,1) + (-t40 * t65 + t63 * t90) * rSges(4,2)) + g(3) * (t126 * rSges(4,1) + (-t104 * t65 - t63 * t94) * rSges(4,2))) - m(5) * (g(1) * (t87 + t99) + g(2) * (t70 + t100) + g(3) * (t84 + t96)) - m(6) * (g(1) * (t72 + t87) + g(2) * (t70 + t73) + g(3) * (t71 + t84)) - m(7) * (g(1) * (t77 + t87) + g(2) * (t70 + t78) + g(3) * (t76 + t84)) -m(5) * (g(1) * t99 + g(2) * t100 + g(3) * t96) - m(6) * (g(1) * t72 + g(2) * t73 + g(3) * t71) - m(7) * (g(1) * t77 + g(2) * t78 + g(3) * t76) t117 * (g(1) * t22 + g(2) * t20 + g(3) * t35) -m(7) * (g(1) * ((-t23 * t53 + t41 * t54) * rSges(7,1) + (-t23 * t54 - t41 * t53) * rSges(7,2)) + g(2) * ((-t21 * t53 + t39 * t54) * rSges(7,1) + (-t21 * t54 - t39 * t53) * rSges(7,2)) + g(3) * ((-t103 * t54 - t36 * t53) * rSges(7,1) + (t103 * t53 - t36 * t54) * rSges(7,2)))];
taug  = t1(:);
