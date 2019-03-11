% Calculate Gravitation load on the joints for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:18
% EndTime: 2019-03-09 18:06:20
% DurationCPUTime: 0.82s
% Computational Cost: add. (586->150), mult. (501->194), div. (0->0), fcn. (452->12), ass. (0->86)
t118 = rSges(6,3) + pkin(9);
t51 = qJ(2) + qJ(3);
t43 = pkin(11) + t51;
t38 = sin(t43);
t39 = cos(t43);
t55 = cos(qJ(5));
t41 = pkin(5) * t55 + pkin(4);
t117 = t38 * rSges(7,3) + t39 * t41;
t72 = t39 * rSges(5,1) - rSges(5,2) * t38;
t45 = sin(t51);
t47 = cos(t51);
t77 = t47 * rSges(4,1) - rSges(4,2) * t45;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t116 = g(1) * t57 + g(2) * t54;
t67 = t39 * pkin(4) + t118 * t38;
t59 = -pkin(8) - pkin(7);
t50 = qJ(5) + qJ(6);
t46 = cos(t50);
t93 = t46 * t57;
t44 = sin(t50);
t96 = t44 * t54;
t5 = t39 * t96 + t93;
t94 = t46 * t54;
t95 = t44 * t57;
t6 = -t39 * t94 + t95;
t115 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t39 * t95 + t94;
t8 = t39 * t93 + t96;
t114 = t7 * rSges(7,1) - t8 * rSges(7,2);
t53 = sin(qJ(2));
t113 = pkin(2) * t53;
t112 = pkin(3) * t45;
t52 = sin(qJ(5));
t111 = pkin(5) * t52;
t108 = g(3) * t38;
t107 = rSges(3,3) + pkin(7);
t56 = cos(qJ(2));
t48 = t56 * pkin(2);
t42 = t48 + pkin(1);
t106 = rSges(6,1) * t55;
t105 = rSges(7,1) * t46;
t102 = rSges(6,2) * t52;
t101 = rSges(7,2) * t44;
t58 = -pkin(10) - pkin(9);
t100 = t38 * t58;
t99 = t39 * t54;
t98 = t39 * t57;
t97 = t39 * t58;
t92 = t52 * t54;
t91 = t52 * t57;
t90 = t54 * t55;
t89 = t55 * t57;
t88 = rSges(4,3) - t59;
t49 = -qJ(4) + t59;
t87 = rSges(5,3) - t49;
t83 = t38 * t101;
t86 = rSges(7,3) * t99 + t54 * t83;
t85 = rSges(7,3) * t98 + t57 * t83;
t84 = t38 * t102;
t82 = t118 * t99 + t54 * t84;
t81 = t118 * t98 + t57 * t84;
t79 = -t49 + t111;
t78 = -t41 - t105;
t40 = pkin(3) * t47;
t76 = t40 + t72;
t75 = rSges(3,1) * t56 - rSges(3,2) * t53;
t73 = -rSges(4,1) * t45 - rSges(4,2) * t47;
t71 = -rSges(5,1) * t38 - rSges(5,2) * t39;
t70 = -rSges(7,1) * t44 - rSges(7,2) * t46;
t69 = pkin(1) + t75;
t11 = -t39 * t91 + t90;
t9 = t39 * t92 + t89;
t68 = t42 + t77;
t66 = t117 + t40 + (-t101 + t105) * t39;
t63 = -t100 + t117;
t62 = t40 + t67 + (-t102 + t106) * t39;
t60 = t116 * (-pkin(4) - t106) * t38;
t23 = -t112 - t113;
t22 = t40 + t42;
t17 = t57 * t23;
t16 = t54 * t23;
t13 = t57 * t22;
t12 = t39 * t89 + t92;
t10 = -t39 * t90 + t91;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t54 - rSges(2,2) * t57) + g(2) * (rSges(2,1) * t57 - rSges(2,2) * t54)) - m(3) * ((g(1) * t107 + g(2) * t69) * t57 + (-g(1) * t69 + g(2) * t107) * t54) - m(4) * ((g(1) * t88 + g(2) * t68) * t57 + (-g(1) * t68 + g(2) * t88) * t54) - m(5) * (g(2) * t13 + (g(1) * t87 + g(2) * t72) * t57 + (g(1) * (-t22 - t72) + g(2) * t87) * t54) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t11 + t13) + (-g(1) * t49 + g(2) * t67) * t57 + (g(1) * (-t22 - t67) - g(2) * t49) * t54) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t13) + (g(1) * t79 + g(2) * t63) * t57 + (g(1) * (-t22 - t63) + g(2) * t79) * t54) -m(3) * (g(3) * t75 + t116 * (-rSges(3,1) * t53 - rSges(3,2) * t56)) - m(4) * (g(3) * (t48 + t77) + t116 * (t73 - t113)) - m(5) * (g(1) * (t71 * t57 + t17) + g(2) * (t71 * t54 + t16) + g(3) * (t48 + t76)) - m(6) * (g(1) * (t17 + t81) + g(2) * (t16 + t82) + g(3) * (t48 + t62) + t60) - m(7) * (g(1) * (-t57 * t97 + t17 + t85) + g(2) * (-t54 * t97 + t16 + t86) + g(3) * (t48 + t66) + (-g(3) * t58 + t116 * t78) * t38) -m(4) * (g(3) * t77 + t116 * t73) - m(5) * (g(3) * t76 + t116 * (t71 - t112)) - m(6) * (g(1) * (-t57 * t112 + t81) + g(2) * (-t54 * t112 + t82) + g(3) * t62 + t60) - m(7) * (g(1) * t85 + g(2) * t86 + g(3) * (t66 - t100) + t116 * (t78 * t38 - t112 - t97)) (-m(5) - m(6) - m(7)) * (g(1) * t54 - g(2) * t57) -m(6) * (g(1) * (rSges(6,1) * t11 - rSges(6,2) * t12) + g(2) * (-rSges(6,1) * t9 + rSges(6,2) * t10)) - m(7) * (g(1) * (t11 * pkin(5) + t114) + g(2) * (-t9 * pkin(5) + t115)) + (-m(6) * (-rSges(6,1) * t52 - rSges(6,2) * t55) - m(7) * (t70 - t111)) * t108, -m(7) * (g(1) * t114 + g(2) * t115 + t70 * t108)];
taug  = t1(:);
