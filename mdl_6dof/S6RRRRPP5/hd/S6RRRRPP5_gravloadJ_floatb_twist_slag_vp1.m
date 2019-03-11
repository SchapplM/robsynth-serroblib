% Calculate Gravitation load on the joints for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:17
% EndTime: 2019-03-09 21:05:20
% DurationCPUTime: 1.08s
% Computational Cost: add. (568->174), mult. (766->224), div. (0->0), fcn. (775->8), ass. (0->74)
t98 = rSges(7,1) + pkin(5);
t48 = sin(qJ(3));
t50 = sin(qJ(1));
t51 = cos(qJ(3));
t52 = cos(qJ(2));
t53 = cos(qJ(1));
t89 = t52 * t53;
t26 = -t48 * t89 + t50 * t51;
t47 = qJ(3) + qJ(4);
t42 = sin(t47);
t43 = cos(t47);
t90 = t50 * t52;
t20 = t42 * t90 + t43 * t53;
t88 = t53 * t42;
t21 = t43 * t90 - t88;
t113 = -t20 * rSges(6,1) + t21 * rSges(6,3);
t49 = sin(qJ(2));
t101 = g(3) * t49;
t41 = pkin(3) * t51 + pkin(2);
t31 = t52 * t41;
t75 = -pkin(1) - t31;
t22 = -t50 * t43 + t52 * t88;
t23 = t42 * t50 + t43 * t89;
t112 = -t22 * rSges(6,1) + t23 * rSges(6,3);
t111 = rSges(7,2) + qJ(5);
t102 = g(2) * t50;
t110 = g(1) * t53 + t102;
t109 = t21 * rSges(7,2) - t98 * t20;
t108 = t23 * rSges(7,2) - t98 * t22;
t97 = rSges(4,3) + pkin(8);
t107 = t52 * pkin(2) + t97 * t49;
t106 = -t20 * rSges(5,1) - t21 * rSges(5,2);
t105 = pkin(3) * t48;
t104 = g(1) * t50;
t99 = -rSges(6,1) - pkin(4);
t96 = rSges(3,2) * t49;
t95 = rSges(5,2) * t43;
t94 = t43 * t49;
t93 = t48 * t50;
t92 = t48 * t53;
t54 = -pkin(9) - pkin(8);
t87 = rSges(6,2) - t54;
t86 = rSges(5,3) - t54;
t85 = -t22 * rSges(5,1) - t23 * rSges(5,2);
t84 = t53 * pkin(1) + t50 * pkin(7);
t83 = rSges(6,3) + qJ(5);
t82 = rSges(7,3) + qJ(6);
t81 = -pkin(4) - t98;
t45 = t53 * pkin(7);
t79 = t50 * t49 * t54 + pkin(3) * t92 + t45;
t78 = -t54 - t82;
t29 = qJ(5) * t94;
t77 = g(3) * (rSges(6,3) * t94 + t29);
t76 = g(3) * (rSges(7,2) * t94 + t29);
t74 = t87 * t53;
t73 = t86 * t53;
t72 = -t20 * pkin(4) + qJ(5) * t21;
t71 = -t22 * pkin(4) + qJ(5) * t23;
t70 = pkin(3) * t93 + t41 * t89 + t84;
t69 = g(3) * (t31 + (pkin(4) * t43 + qJ(5) * t42) * t52);
t68 = t78 * t53;
t67 = t23 * pkin(4) + t70;
t66 = rSges(3,1) * t52 - t96;
t64 = rSges(5,1) * t43 - rSges(5,2) * t42;
t63 = t26 * pkin(3);
t62 = rSges(4,1) * t51 - rSges(4,2) * t48 + pkin(2);
t24 = t48 * t90 + t51 * t53;
t60 = -t21 * pkin(4) - qJ(5) * t20 + t79;
t59 = t24 * pkin(3);
t58 = t63 + t71;
t55 = -t59 + t72;
t27 = t51 * t89 + t93;
t25 = -t51 * t90 + t92;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t50 - rSges(2,2) * t53) + g(2) * (rSges(2,1) * t53 - rSges(2,2) * t50)) - m(3) * (g(1) * (rSges(3,3) * t53 + t45) + g(2) * (rSges(3,1) * t89 - t53 * t96 + t84) + (g(1) * (-pkin(1) - t66) + g(2) * rSges(3,3)) * t50) - m(4) * (g(1) * (rSges(4,1) * t25 + rSges(4,2) * t24 + t45) + (-pkin(1) - t107) * t104 + (rSges(4,1) * t27 + rSges(4,2) * t26 + t107 * t53 + t84) * g(2)) - m(5) * (g(1) * (-rSges(5,1) * t21 + rSges(5,2) * t20 + t79) + g(2) * (t23 * rSges(5,1) - t22 * rSges(5,2) + t49 * t73 + t70) + (-t49 * rSges(5,3) + t75) * t104) - m(6) * (g(1) * (-rSges(6,1) * t21 - rSges(6,3) * t20 + t60) + g(2) * (t23 * rSges(6,1) + t83 * t22 + t49 * t74 + t67) + (-t49 * rSges(6,2) + t75) * t104) - m(7) * (g(1) * (-rSges(7,2) * t20 - t98 * t21 + t75 * t50 + t60) + g(2) * (t111 * t22 + t98 * t23 + t67) + (g(2) * t68 + t82 * t104) * t49) -m(3) * (g(3) * t66 + t110 * (-rSges(3,1) * t49 - rSges(3,2) * t52)) - m(4) * ((g(3) * t62 + t110 * t97) * t52 + (g(3) * t97 - t110 * t62) * t49) - m(5) * (g(3) * t31 + (g(1) * t73 + g(3) * t64 + t86 * t102) * t52 + (g(3) * t86 + t110 * (-t41 - t64)) * t49) - m(6) * (t69 + (g(3) * (rSges(6,1) * t43 + rSges(6,3) * t42) + g(1) * t74 + t87 * t102) * t52 + (g(3) * t87 + t110 * (-t83 * t42 + t99 * t43 - t41)) * t49) - m(7) * (t69 + (g(3) * (rSges(7,2) * t42 + t98 * t43) + g(1) * t68 + t78 * t102) * t52 + (g(3) * t78 + t110 * (-t111 * t42 + t81 * t43 - t41)) * t49) -m(4) * (g(1) * (rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-rSges(4,1) * t24 + rSges(4,2) * t25) + (-rSges(4,1) * t48 - rSges(4,2) * t51) * t101) - m(5) * (g(1) * (t63 + t85) + g(2) * (-t59 + t106) + (-rSges(5,1) * t42 - t105 - t95) * t101) - m(6) * (g(1) * (t58 + t112) + g(2) * (t55 + t113) + t77 + (t99 * t42 - t105) * t101) - m(7) * (g(1) * (t58 + t108) + g(2) * (t55 + t109) + t76 + (t81 * t42 - t105) * t101) -m(5) * (g(1) * t85 + g(2) * t106) - m(6) * (g(1) * (t71 + t112) + g(2) * (t72 + t113) + t77) - m(7) * (g(1) * (t71 + t108) + g(2) * (t72 + t109) + t76) + (m(5) * t95 + (m(5) * rSges(5,1) - m(6) * t99 - m(7) * t81) * t42) * t101 (-m(6) - m(7)) * (g(1) * t22 + g(2) * t20 + t42 * t101) -m(7) * (g(3) * t52 - t110 * t49)];
taug  = t1(:);
