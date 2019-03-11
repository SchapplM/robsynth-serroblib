% Calculate Gravitation load on the joints for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:06
% EndTime: 2019-03-10 01:12:09
% DurationCPUTime: 1.04s
% Computational Cost: add. (678->169), mult. (693->223), div. (0->0), fcn. (664->10), ass. (0->91)
t144 = rSges(5,3) + pkin(9);
t143 = rSges(7,1) + pkin(5);
t142 = rSges(7,3) + qJ(6);
t65 = sin(qJ(1));
t66 = cos(qJ(4));
t115 = t65 * t66;
t63 = sin(qJ(4));
t68 = cos(qJ(1));
t117 = t63 * t68;
t62 = qJ(2) + qJ(3);
t59 = cos(t62);
t17 = -t59 * t117 + t115;
t57 = sin(t62);
t141 = t59 * rSges(4,1) - rSges(4,2) * t57;
t140 = g(1) * t68 + g(2) * t65;
t78 = t59 * pkin(3) + t144 * t57;
t120 = t59 * t65;
t61 = qJ(4) + qJ(5);
t56 = sin(t61);
t58 = cos(t61);
t11 = t56 * t120 + t58 * t68;
t112 = t68 * t56;
t116 = t65 * t58;
t12 = t59 * t116 - t112;
t139 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t59 * t112 - t116;
t119 = t59 * t68;
t14 = t58 * t119 + t56 * t65;
t138 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t64 = sin(qJ(2));
t137 = pkin(2) * t64;
t136 = pkin(4) * t63;
t70 = -pkin(8) - pkin(7);
t133 = g(2) * t70;
t132 = g(3) * t57;
t69 = -pkin(10) - pkin(9);
t131 = g(3) * t69;
t129 = rSges(3,3) + pkin(7);
t128 = rSges(5,1) * t66;
t126 = rSges(5,2) * t63;
t125 = t56 * t57;
t124 = t56 * t59;
t50 = t57 * rSges(7,2);
t48 = t57 * rSges(6,3);
t122 = t57 * t68;
t121 = t58 * t59;
t54 = pkin(4) * t66 + pkin(3);
t26 = t59 * t54;
t118 = t63 * t65;
t114 = t65 * t69;
t113 = t66 * t68;
t111 = t68 * t69;
t110 = rSges(4,3) - t70;
t103 = rSges(6,2) * t125;
t109 = rSges(6,3) * t120 + t65 * t103;
t108 = rSges(6,3) * t119 + t68 * t103;
t107 = t142 * t57 * t58;
t106 = t65 * t137;
t105 = t68 * t137;
t104 = t57 * t126;
t101 = t59 * t114;
t100 = t59 * t111;
t99 = t65 * t104 + t144 * t120;
t98 = t68 * t104 + t144 * t119;
t96 = -rSges(6,1) * t58 - t54;
t67 = cos(qJ(2));
t60 = t67 * pkin(2);
t55 = t60 + pkin(1);
t94 = -t55 - t26;
t93 = rSges(7,2) * t120 - t101;
t92 = rSges(7,2) * t119 - t100;
t91 = pkin(4) * t117 + t57 * t114 - t68 * t70;
t90 = rSges(3,1) * t67 - rSges(3,2) * t64;
t87 = -rSges(4,1) * t57 - rSges(4,2) * t59;
t86 = -rSges(6,1) * t56 - rSges(6,2) * t58;
t85 = -t143 * t11 + t142 * t12;
t84 = t17 * pkin(4);
t83 = -t143 * t13 + t142 * t14;
t82 = t143 * t121 + t142 * t124 + t26 + t50;
t81 = pkin(1) + t90;
t15 = t59 * t118 + t113;
t80 = rSges(6,1) * t121 - rSges(6,2) * t124 + t26 + t48;
t33 = t68 * t55;
t79 = pkin(4) * t118 - t57 * t111 + t54 * t119 + t33;
t77 = t78 + (-t126 + t128) * t59;
t76 = t15 * pkin(4);
t72 = t140 * (-pkin(3) - t128) * t57;
t71 = (-t131 + t140 * (-t142 * t56 - t143 * t58 - t54)) * t57;
t18 = t59 * t113 + t118;
t16 = -t59 * t115 + t117;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t65 - rSges(2,2) * t68) + g(2) * (rSges(2,1) * t68 - rSges(2,2) * t65)) - m(3) * ((g(1) * t129 + g(2) * t81) * t68 + (-g(1) * t81 + g(2) * t129) * t65) - m(4) * (g(2) * t33 + (g(1) * t110 + g(2) * t141) * t68 + (g(1) * (-t55 - t141) + g(2) * t110) * t65) - m(5) * (g(1) * (t16 * rSges(5,1) + t15 * rSges(5,2)) + g(2) * (t18 * rSges(5,1) + t17 * rSges(5,2) + t33) + (-g(1) * t70 + g(2) * t78) * t68 + (g(1) * (-t55 - t78) - t133) * t65) - m(6) * (g(1) * (-t12 * rSges(6,1) + t11 * rSges(6,2) + t91) + g(2) * (t14 * rSges(6,1) - t13 * rSges(6,2) + rSges(6,3) * t122 + t79) + (g(1) * (t94 - t48) - t133) * t65) - m(7) * (g(1) * (-t142 * t11 - t143 * t12 + t91) + g(2) * (rSges(7,2) * t122 + t142 * t13 + t143 * t14 + t79) + (g(1) * (t94 - t50) - t133) * t65) -m(3) * (g(3) * t90 + t140 * (-rSges(3,1) * t64 - rSges(3,2) * t67)) - m(4) * (g(3) * (t60 + t141) + t140 * (t87 - t137)) - m(5) * (g(1) * (t98 - t105) + g(2) * (t99 - t106) + g(3) * (t60 + t77) + t72) - m(6) * (g(1) * t108 + g(2) * t109 + g(3) * (-t57 * t69 + t60 + t80) + t140 * (t96 * t57 - t59 * t69 - t137)) - m(7) * (g(1) * (t92 - t105) + g(2) * (t93 - t106) + g(3) * (t60 + t82) + t71) -m(4) * (g(3) * t141 + t140 * t87) - m(5) * (g(1) * t98 + g(2) * t99 + g(3) * t77 + t72) - m(6) * (g(1) * (-t100 + t108) + g(2) * (-t101 + t109) + g(3) * t80 + (t140 * t96 - t131) * t57) - m(7) * (g(1) * t92 + g(2) * t93 + g(3) * t82 + t71) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (t84 + t138) + g(2) * (-t76 + t139)) - m(7) * (g(1) * (t83 + t84) + g(2) * (-t76 + t85) + g(3) * t107) + (-m(5) * (-rSges(5,1) * t63 - rSges(5,2) * t66) - m(6) * (t86 - t136) - m(7) * (-t143 * t56 - t136)) * t132, -m(6) * (g(1) * t138 + g(2) * t139 + t86 * t132) - m(7) * (g(1) * t83 + g(2) * t85 + g(3) * (-t125 * t143 + t107)) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t125)];
taug  = t1(:);
