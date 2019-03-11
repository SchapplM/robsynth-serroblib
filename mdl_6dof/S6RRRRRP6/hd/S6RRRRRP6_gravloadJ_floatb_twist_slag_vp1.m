% Calculate Gravitation load on the joints for
% S6RRRRRP6
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:35
% EndTime: 2019-03-10 01:26:38
% DurationCPUTime: 1.00s
% Computational Cost: add. (723->171), mult. (748->230), div. (0->0), fcn. (743->10), ass. (0->81)
t103 = rSges(7,1) + pkin(5);
t56 = qJ(3) + qJ(4);
t48 = sin(t56);
t49 = cos(t56);
t59 = sin(qJ(1));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t96 = t61 * t62;
t21 = -t48 * t96 + t49 * t59;
t88 = rSges(7,3) + qJ(6);
t63 = -pkin(9) - pkin(8);
t93 = rSges(5,3) - t63;
t107 = g(2) * t59;
t116 = g(1) * t62 + t107;
t102 = rSges(4,3) + pkin(8);
t58 = sin(qJ(2));
t115 = t61 * pkin(2) + t102 * t58;
t50 = qJ(5) + t56;
t45 = sin(t50);
t46 = cos(t50);
t114 = t103 * t46 + t88 * t45;
t97 = t59 * t61;
t11 = t45 * t97 + t46 * t62;
t95 = t62 * t45;
t12 = t46 * t97 - t95;
t113 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = -t59 * t46 + t61 * t95;
t14 = t45 * t59 + t46 * t96;
t112 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t57 = sin(qJ(3));
t111 = pkin(3) * t57;
t110 = pkin(4) * t48;
t109 = g(1) * t59;
t60 = cos(qJ(3));
t52 = t60 * pkin(3);
t37 = pkin(4) * t49 + t52;
t35 = pkin(2) + t37;
t30 = t61 * t35;
t106 = g(3) * t30;
t105 = g(3) * t58;
t101 = rSges(3,2) * t58;
t100 = t45 * t58;
t36 = t110 + t111;
t34 = t62 * t36;
t55 = -pkin(10) + t63;
t94 = rSges(7,2) - t55;
t92 = rSges(6,3) - t55;
t91 = -t61 * t34 + t59 * t37;
t90 = t88 * t46 * t58;
t89 = t62 * pkin(1) + t59 * pkin(7);
t53 = t62 * pkin(7);
t86 = t59 * t58 * t55 + t34 + t53;
t85 = -pkin(1) - t30;
t84 = t94 * t62;
t83 = t92 * t62;
t82 = -t36 * t97 - t37 * t62;
t81 = t35 * t96 + t59 * t36 + t89;
t80 = rSges(3,1) * t61 - t101;
t78 = -rSges(5,1) * t48 - rSges(5,2) * t49;
t77 = rSges(6,1) * t46 - rSges(6,2) * t45;
t76 = -rSges(6,1) * t45 - rSges(6,2) * t46;
t75 = -t103 * t11 + t88 * t12;
t74 = t21 * pkin(4);
t73 = -t103 * t13 + t88 * t14;
t72 = rSges(4,1) * t60 - rSges(4,2) * t57 + pkin(2);
t19 = t48 * t97 + t49 * t62;
t28 = -t57 * t96 + t59 * t60;
t26 = t57 * t97 + t60 * t62;
t47 = t52 + pkin(2);
t71 = rSges(5,1) * t49 - rSges(5,2) * t48 + t47;
t69 = t76 * t58;
t68 = -t103 * t100 + t90;
t67 = t19 * pkin(4);
t66 = t61 * t47 + t58 * t93;
t20 = t48 * t62 - t49 * t97;
t22 = t48 * t59 + t49 * t96;
t65 = g(1) * (t21 * rSges(5,1) - t22 * rSges(5,2)) + g(2) * (-t19 * rSges(5,1) + t20 * rSges(5,2));
t31 = t58 * t36;
t29 = t57 * t59 + t60 * t96;
t27 = t57 * t62 - t60 * t97;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t59 - rSges(2,2) * t62) + g(2) * (rSges(2,1) * t62 - rSges(2,2) * t59)) - m(3) * (g(1) * (rSges(3,3) * t62 + t53) + g(2) * (rSges(3,1) * t96 - t62 * t101 + t89) + (g(1) * (-pkin(1) - t80) + g(2) * rSges(3,3)) * t59) - m(4) * (g(1) * (rSges(4,1) * t27 + rSges(4,2) * t26 + t53) + (-pkin(1) - t115) * t109 + (rSges(4,1) * t29 + rSges(4,2) * t28 + t115 * t62 + t89) * g(2)) - m(5) * (g(1) * (t20 * rSges(5,1) + t19 * rSges(5,2) + t53) + g(2) * (t22 * rSges(5,1) + t21 * rSges(5,2) + t89) + (g(1) * t111 + g(2) * t66) * t62 + (g(1) * (-pkin(1) - t66) + g(2) * t111) * t59) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + t86) + g(2) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t58 * t83 + t81) + (-rSges(6,3) * t58 + t85) * t109) - m(7) * (g(1) * (-t103 * t12 - t88 * t11 + t86) + g(2) * (t103 * t14 + t88 * t13 + t58 * t84 + t81) + (-rSges(7,2) * t58 + t85) * t109) -m(3) * (g(3) * t80 + t116 * (-rSges(3,1) * t58 - rSges(3,2) * t61)) - m(4) * ((g(3) * t72 + t102 * t116) * t61 + (g(3) * t102 - t116 * t72) * t58) - m(5) * ((g(3) * t71 + t116 * t93) * t61 + (g(3) * t93 - t116 * t71) * t58) - m(6) * (t106 + (g(1) * t83 + g(3) * t77 + t92 * t107) * t61 + (g(3) * t92 + t116 * (-t35 - t77)) * t58) - m(7) * (t106 + (g(1) * t84 + g(3) * t114 + t94 * t107) * t61 + (g(3) * t94 + t116 * (-t114 - t35)) * t58) -m(4) * (g(1) * (rSges(4,1) * t28 - rSges(4,2) * t29) + g(2) * (-rSges(4,1) * t26 + rSges(4,2) * t27) + (-rSges(4,1) * t57 - rSges(4,2) * t60) * t105) - m(5) * ((g(1) * t28 - g(2) * t26) * pkin(3) + (t78 - t111) * t105 + t65) - m(6) * (g(1) * (t91 + t112) + g(2) * (t82 + t113) + g(3) * (-t31 + t69)) - m(7) * (g(1) * (t73 + t91) + g(2) * (t75 + t82) + g(3) * (-t31 + t68)) -m(5) * t65 - m(6) * (g(1) * (t74 + t112) + g(2) * (-t67 + t113)) - m(7) * (g(1) * (t73 + t74) + g(2) * (-t67 + t75) + g(3) * t90) + (-m(5) * t78 - m(6) * (t76 - t110) - m(7) * (-t103 * t45 - t110)) * t105, -m(6) * (g(1) * t112 + g(2) * t113 + g(3) * t69) - m(7) * (g(1) * t73 + g(2) * t75 + g(3) * t68) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t100)];
taug  = t1(:);
