% Calculate Gravitation load on the joints for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:52
% EndTime: 2019-03-10 01:05:55
% DurationCPUTime: 1.02s
% Computational Cost: add. (614->174), mult. (627->234), div. (0->0), fcn. (584->10), ass. (0->88)
t131 = rSges(5,3) + pkin(9);
t56 = qJ(4) + qJ(5);
t51 = cos(t56);
t61 = cos(qJ(4));
t53 = t61 * pkin(4);
t26 = pkin(5) * t51 + t53;
t22 = pkin(3) + t26;
t57 = qJ(2) + qJ(3);
t50 = sin(t57);
t52 = cos(t57);
t130 = t50 * rSges(7,3) + t52 * t22;
t47 = t53 + pkin(3);
t129 = t50 * rSges(6,3) + t52 * t47;
t128 = t52 * rSges(4,1) - rSges(4,2) * t50;
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t127 = g(1) * t63 + g(2) * t60;
t73 = t52 * pkin(3) + t131 * t50;
t126 = t127 * t50;
t64 = -pkin(10) - pkin(9);
t105 = t52 * t60;
t49 = sin(t56);
t109 = t49 * t63;
t10 = -t105 * t51 + t109;
t9 = t105 * t49 + t51 * t63;
t125 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t124 = -t9 * rSges(7,1) + t10 * rSges(7,2);
t104 = t52 * t63;
t11 = -t104 * t49 + t51 * t60;
t110 = t49 * t60;
t12 = t104 * t51 + t110;
t123 = t11 * rSges(7,1) - t12 * rSges(7,2);
t122 = t11 * rSges(6,1) - t12 * rSges(6,2);
t59 = sin(qJ(2));
t121 = pkin(2) * t59;
t58 = sin(qJ(4));
t120 = pkin(4) * t58;
t119 = pkin(5) * t49;
t116 = g(3) * t50;
t115 = rSges(3,3) + pkin(7);
t114 = rSges(5,1) * t61;
t112 = rSges(5,2) * t58;
t111 = t49 * t52;
t55 = -qJ(6) + t64;
t108 = t50 * t55;
t107 = t50 * t64;
t106 = t51 * t52;
t103 = t52 * t64;
t102 = t58 * t60;
t101 = t58 * t63;
t100 = t60 * t61;
t99 = t61 * t63;
t65 = -pkin(8) - pkin(7);
t98 = rSges(4,3) - t65;
t91 = t50 * t110;
t97 = rSges(7,2) * t91 + rSges(7,3) * t105;
t96 = rSges(6,2) * t91 + rSges(6,3) * t105;
t90 = t50 * t109;
t95 = rSges(7,2) * t90 + rSges(7,3) * t104;
t94 = rSges(6,2) * t90 + rSges(6,3) * t104;
t25 = t119 + t120;
t93 = t25 - t65;
t92 = t50 * t112;
t89 = t105 * t131 + t60 * t92;
t88 = t104 * t131 + t63 * t92;
t86 = -t65 + t120;
t85 = -rSges(6,1) * t51 - t47;
t84 = -rSges(7,1) * t51 - t22;
t62 = cos(qJ(2));
t82 = rSges(3,1) * t62 - rSges(3,2) * t59;
t79 = -rSges(4,1) * t50 - rSges(4,2) * t52;
t78 = -rSges(6,1) * t49 - rSges(6,2) * t51;
t77 = -rSges(7,1) * t49 - rSges(7,2) * t51;
t76 = pkin(1) + t82;
t15 = -t101 * t52 + t100;
t13 = t102 * t52 + t99;
t75 = rSges(6,1) * t106 - rSges(6,2) * t111 + t129;
t74 = rSges(7,1) * t106 - rSges(7,2) * t111 + t130;
t72 = t73 + (-t112 + t114) * t52;
t70 = -t107 + t129;
t69 = -t108 + t130;
t66 = (-pkin(3) - t114) * t126;
t54 = t62 * pkin(2);
t48 = t54 + pkin(1);
t30 = t63 * t48;
t16 = t52 * t99 + t102;
t14 = -t100 * t52 + t101;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t60 - rSges(2,2) * t63) + g(2) * (rSges(2,1) * t63 - rSges(2,2) * t60)) - m(3) * ((g(1) * t115 + g(2) * t76) * t63 + (-g(1) * t76 + g(2) * t115) * t60) - m(4) * (g(2) * t30 + (g(1) * t98 + g(2) * t128) * t63 + (g(1) * (-t48 - t128) + g(2) * t98) * t60) - m(5) * (g(1) * (t14 * rSges(5,1) + t13 * rSges(5,2)) + g(2) * (t16 * rSges(5,1) + t15 * rSges(5,2) + t30) + (-g(1) * t65 + g(2) * t73) * t63 + (g(1) * (-t48 - t73) - g(2) * t65) * t60) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2)) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t30) + (g(1) * t86 + g(2) * t70) * t63 + (g(1) * (-t48 - t70) + g(2) * t86) * t60) - m(7) * (g(1) * (t10 * rSges(7,1) + t9 * rSges(7,2)) + g(2) * (t12 * rSges(7,1) + t11 * rSges(7,2) + t30) + (g(1) * t93 + g(2) * t69) * t63 + (g(1) * (-t48 - t69) + g(2) * t93) * t60) -m(3) * (g(3) * t82 + t127 * (-rSges(3,1) * t59 - rSges(3,2) * t62)) - m(4) * (g(3) * (t54 + t128) + t127 * (t79 - t121)) - m(5) * (g(1) * (-t121 * t63 + t88) + g(2) * (-t121 * t60 + t89) + g(3) * (t54 + t72) + t66) - m(6) * (g(1) * t94 + g(2) * t96 + g(3) * (t54 + t75 - t107) + t127 * (t85 * t50 - t103 - t121)) - m(7) * (g(1) * t95 + g(2) * t97 + g(3) * (t54 + t74 - t108) + t127 * (t84 * t50 - t52 * t55 - t121)) -m(4) * (g(3) * t128 + t127 * t79) - m(5) * (g(1) * t88 + g(2) * t89 + g(3) * t72 + t66) - m(6) * (g(1) * (-t103 * t63 + t94) + g(2) * (-t103 * t60 + t96) + g(3) * t75 + (-g(3) * t64 + t127 * t85) * t50) - m(7) * (g(1) * (-t104 * t55 + t95) + g(2) * (-t105 * t55 + t97) + g(3) * t74 + (-g(3) * t55 + t127 * t84) * t50) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (-rSges(5,1) * t13 + rSges(5,2) * t14)) - m(6) * (g(1) * (pkin(4) * t15 + t122) + g(2) * (-pkin(4) * t13 + t125)) - m(7) * (g(1) * (-t104 * t25 + t26 * t60 + t123) + g(2) * (-t105 * t25 - t26 * t63 + t124)) + (-m(5) * (-rSges(5,1) * t58 - rSges(5,2) * t61) - m(6) * (t78 - t120) - m(7) * (-t25 + t77)) * t116, -m(6) * (g(1) * t122 + g(2) * t125) - m(7) * (g(1) * (pkin(5) * t11 + t123) + g(2) * (-pkin(5) * t9 + t124)) + (-m(6) * t78 - m(7) * (t77 - t119)) * t116, -m(7) * (-g(3) * t52 + t126)];
taug  = t1(:);
