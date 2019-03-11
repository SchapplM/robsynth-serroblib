% Calculate Gravitation load on the joints for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:40
% EndTime: 2019-03-09 20:58:42
% DurationCPUTime: 0.99s
% Computational Cost: add. (639->168), mult. (692->225), div. (0->0), fcn. (681->10), ass. (0->79)
t102 = rSges(7,1) + pkin(5);
t57 = qJ(3) + qJ(4);
t50 = sin(t57);
t51 = cos(t57);
t60 = sin(qJ(1));
t62 = cos(qJ(2));
t63 = cos(qJ(1));
t95 = t62 * t63;
t21 = -t50 * t95 + t51 * t60;
t87 = rSges(7,3) + qJ(6);
t64 = -pkin(9) - pkin(8);
t92 = rSges(5,3) - t64;
t106 = g(2) * t60;
t115 = g(1) * t63 + t106;
t101 = rSges(4,3) + pkin(8);
t59 = sin(qJ(2));
t114 = t62 * pkin(2) + t101 * t59;
t49 = pkin(10) + t57;
t45 = sin(t49);
t46 = cos(t49);
t113 = t102 * t46 + t87 * t45;
t96 = t60 * t62;
t11 = t45 * t96 + t46 * t63;
t94 = t63 * t45;
t12 = t46 * t96 - t94;
t112 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = -t60 * t46 + t62 * t94;
t14 = t45 * t60 + t46 * t95;
t111 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t58 = sin(qJ(3));
t110 = pkin(3) * t58;
t109 = pkin(4) * t50;
t108 = g(1) * t60;
t61 = cos(qJ(3));
t53 = t61 * pkin(3);
t38 = pkin(4) * t51 + t53;
t36 = pkin(2) + t38;
t31 = t62 * t36;
t105 = g(3) * t31;
t104 = g(3) * t59;
t100 = rSges(3,2) * t59;
t99 = t45 * t59;
t37 = t109 + t110;
t35 = t63 * t37;
t56 = -qJ(5) + t64;
t93 = rSges(7,2) - t56;
t91 = rSges(6,3) - t56;
t90 = -t62 * t35 + t60 * t38;
t89 = t87 * t46 * t59;
t88 = t63 * pkin(1) + t60 * pkin(7);
t54 = t63 * pkin(7);
t85 = t60 * t59 * t56 + t35 + t54;
t84 = -pkin(1) - t31;
t83 = t93 * t63;
t82 = t91 * t63;
t81 = -t37 * t96 - t38 * t63;
t80 = t36 * t95 + t60 * t37 + t88;
t79 = rSges(3,1) * t62 - t100;
t77 = -rSges(5,1) * t50 - rSges(5,2) * t51;
t76 = rSges(6,1) * t46 - rSges(6,2) * t45;
t75 = -rSges(6,1) * t45 - rSges(6,2) * t46;
t74 = -t102 * t11 + t87 * t12;
t73 = t21 * pkin(4);
t72 = -t102 * t13 + t87 * t14;
t71 = rSges(4,1) * t61 - rSges(4,2) * t58 + pkin(2);
t19 = t50 * t96 + t51 * t63;
t29 = -t58 * t95 + t60 * t61;
t27 = t58 * t96 + t61 * t63;
t48 = t53 + pkin(2);
t70 = rSges(5,1) * t51 - rSges(5,2) * t50 + t48;
t68 = t19 * pkin(4);
t67 = t62 * t48 + t59 * t92;
t20 = t50 * t63 - t51 * t96;
t22 = t50 * t60 + t51 * t95;
t66 = g(1) * (t21 * rSges(5,1) - t22 * rSges(5,2)) + g(2) * (-t19 * rSges(5,1) + t20 * rSges(5,2));
t32 = t59 * t37;
t30 = t58 * t60 + t61 * t95;
t28 = t58 * t63 - t61 * t96;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t60 - rSges(2,2) * t63) + g(2) * (rSges(2,1) * t63 - rSges(2,2) * t60)) - m(3) * (g(1) * (rSges(3,3) * t63 + t54) + g(2) * (rSges(3,1) * t95 - t63 * t100 + t88) + (g(1) * (-pkin(1) - t79) + g(2) * rSges(3,3)) * t60) - m(4) * (g(1) * (rSges(4,1) * t28 + rSges(4,2) * t27 + t54) + (-pkin(1) - t114) * t108 + (rSges(4,1) * t30 + rSges(4,2) * t29 + t114 * t63 + t88) * g(2)) - m(5) * (g(1) * (t20 * rSges(5,1) + t19 * rSges(5,2) + t54) + g(2) * (t22 * rSges(5,1) + t21 * rSges(5,2) + t88) + (g(1) * t110 + g(2) * t67) * t63 + (g(1) * (-pkin(1) - t67) + g(2) * t110) * t60) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + t85) + g(2) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t59 * t82 + t80) + (-rSges(6,3) * t59 + t84) * t108) - m(7) * (g(1) * (-t102 * t12 - t87 * t11 + t85) + g(2) * (t102 * t14 + t87 * t13 + t59 * t83 + t80) + (-rSges(7,2) * t59 + t84) * t108) -m(3) * (g(3) * t79 + t115 * (-rSges(3,1) * t59 - rSges(3,2) * t62)) - m(4) * ((g(3) * t71 + t101 * t115) * t62 + (g(3) * t101 - t115 * t71) * t59) - m(5) * ((g(3) * t70 + t115 * t92) * t62 + (g(3) * t92 - t115 * t70) * t59) - m(6) * (t105 + (g(1) * t82 + g(3) * t76 + t91 * t106) * t62 + (g(3) * t91 + t115 * (-t36 - t76)) * t59) - m(7) * (t105 + (g(1) * t83 + g(3) * t113 + t93 * t106) * t62 + (g(3) * t93 + t115 * (-t113 - t36)) * t59) -m(4) * (g(1) * (rSges(4,1) * t29 - rSges(4,2) * t30) + g(2) * (-rSges(4,1) * t27 + rSges(4,2) * t28) + (-rSges(4,1) * t58 - rSges(4,2) * t61) * t104) - m(5) * ((g(1) * t29 - g(2) * t27) * pkin(3) + (t77 - t110) * t104 + t66) - m(6) * (g(1) * (t90 + t111) + g(2) * (t81 + t112) + g(3) * (t75 * t59 - t32)) - m(7) * (g(1) * (t72 + t90) + g(2) * (t74 + t81) + g(3) * (-t102 * t99 - t32 + t89)) -m(5) * t66 - m(6) * (g(1) * (t73 + t111) + g(2) * (-t68 + t112)) - m(7) * (g(1) * (t72 + t73) + g(2) * (-t68 + t74) + g(3) * t89) + (-m(5) * t77 - m(6) * (t75 - t109) - m(7) * (-t102 * t45 - t109)) * t104 (-m(6) - m(7)) * (-g(3) * t62 + t115 * t59) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t99)];
taug  = t1(:);
