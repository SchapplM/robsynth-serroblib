% Calculate Gravitation load on the joints for
% S6RRRRRP2
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:41
% EndTime: 2019-03-10 01:00:44
% DurationCPUTime: 0.92s
% Computational Cost: add. (700->156), mult. (623->196), div. (0->0), fcn. (574->10), ass. (0->84)
t108 = rSges(7,1) + pkin(5);
t46 = cos(qJ(5));
t109 = -t108 * t46 - pkin(4);
t45 = sin(qJ(1));
t42 = qJ(2) + qJ(3);
t39 = qJ(4) + t42;
t34 = sin(t39);
t43 = sin(qJ(5));
t94 = t34 * t43;
t73 = rSges(6,2) * t94;
t35 = cos(t39);
t91 = t35 * t45;
t107 = rSges(6,3) * t91 + t45 * t73;
t48 = cos(qJ(1));
t89 = t35 * t48;
t106 = rSges(6,3) * t89 + t48 * t73;
t27 = t35 * rSges(5,1);
t64 = -rSges(5,2) * t34 + t27;
t79 = t35 * pkin(4) + t34 * pkin(10);
t37 = sin(t42);
t38 = cos(t42);
t71 = t38 * rSges(4,1) - rSges(4,2) * t37;
t77 = rSges(7,3) + qJ(6);
t105 = g(1) * t48 + g(2) * t45;
t104 = t105 * t34;
t49 = -pkin(8) - pkin(7);
t44 = sin(qJ(2));
t103 = pkin(2) * t44;
t102 = pkin(3) * t37;
t41 = -pkin(9) + t49;
t100 = g(2) * t41;
t97 = rSges(3,3) + pkin(7);
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t36 = t40 + pkin(1);
t26 = t34 * rSges(7,2);
t25 = t34 * rSges(6,3);
t93 = t34 * t48;
t92 = t35 * t43;
t90 = t35 * t46;
t88 = t41 * t48;
t87 = t43 * t45;
t86 = t45 * t46;
t85 = t46 * t48;
t84 = t48 * t43;
t83 = rSges(4,3) - t49;
t82 = rSges(5,3) - t41;
t14 = rSges(7,2) * t91;
t20 = pkin(10) * t91;
t81 = t14 + t20;
t19 = rSges(7,2) * t89;
t23 = pkin(10) * t89;
t80 = t19 + t23;
t78 = qJ(6) * t43;
t76 = t20 + t107;
t75 = t23 + t106;
t31 = pkin(3) * t38;
t10 = t31 + t36;
t5 = t48 * t10;
t74 = pkin(4) * t89 + pkin(10) * t93 + t5;
t72 = -rSges(6,1) * t46 - pkin(4);
t70 = -t45 * t102 + t20;
t69 = -t48 * t102 + t23;
t68 = t31 + t64;
t67 = rSges(3,1) * t47 - rSges(3,2) * t44;
t65 = -rSges(4,1) * t37 - rSges(4,2) * t38;
t63 = -rSges(5,1) * t34 - rSges(5,2) * t35;
t62 = -t10 - t79;
t61 = pkin(1) + t67;
t60 = t36 + t71;
t59 = rSges(7,3) * t92 + t108 * t90 + t35 * t78 + t26 + t79;
t58 = rSges(6,1) * t90 - rSges(6,2) * t92 + t25 + t79;
t55 = t31 + t59;
t54 = t31 + t58;
t52 = t72 * t104;
t50 = (-t77 * t43 + t109) * t104;
t11 = -t102 - t103;
t7 = t48 * t11;
t6 = t45 * t11;
t4 = t35 * t85 + t87;
t3 = t35 * t84 - t86;
t2 = t35 * t86 - t84;
t1 = t35 * t87 + t85;
t8 = [-m(2) * (g(1) * (-rSges(2,1) * t45 - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - rSges(2,2) * t45)) - m(3) * ((g(1) * t97 + g(2) * t61) * t48 + (-g(1) * t61 + g(2) * t97) * t45) - m(4) * ((g(1) * t83 + g(2) * t60) * t48 + (-g(1) * t60 + g(2) * t83) * t45) - m(5) * (g(2) * t5 + (g(1) * t82 + g(2) * t64) * t48 + (g(1) * (-t10 - t64) + g(2) * t82) * t45) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 - t88) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t93 + t74) + (g(1) * (t62 - t25) - t100) * t45) - m(7) * (g(1) * (-t77 * t1 - t108 * t2 - t88) + g(2) * (rSges(7,2) * t93 + t108 * t4 + t77 * t3 + t74) + (g(1) * (t62 - t26) - t100) * t45) -m(3) * (g(3) * t67 + t105 * (-rSges(3,1) * t44 - rSges(3,2) * t47)) - m(4) * (g(3) * (t40 + t71) + t105 * (t65 - t103)) - m(5) * (g(1) * (t48 * t63 + t7) + g(2) * (t45 * t63 + t6) + g(3) * (t40 + t68)) - m(6) * (g(1) * (t7 + t75) + g(2) * (t6 + t76) + g(3) * (t40 + t54) + t52) - m(7) * (g(1) * (t7 + t80) + g(2) * (t6 + t81) + g(3) * (t40 + t55) + t50) -m(4) * (g(3) * t71 + t105 * t65) - m(5) * (g(3) * t68 + t105 * (t63 - t102)) - m(6) * (g(1) * (t69 + t106) + g(2) * (t70 + t107) + g(3) * t54 + t52) - m(7) * (g(1) * (t19 + t69) + g(2) * (t14 + t70) + g(3) * t55 + t50) -m(5) * (g(3) * t27 + (-g(1) * t89 - g(2) * t91) * rSges(5,2)) - m(6) * (g(1) * t75 + g(2) * t76 + g(3) * t58) - m(7) * (g(1) * t80 + g(2) * t81 + g(3) * t59) + (m(5) * g(3) * rSges(5,2) + t105 * (m(5) * rSges(5,1) - m(6) * t72 - m(7) * (-rSges(7,3) * t43 + t109 - t78))) * t34, -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t108 * t3 + t77 * t4) + g(2) * (-t1 * t108 + t77 * t2)) + (-m(6) * (-rSges(6,1) * t43 - rSges(6,2) * t46) - m(7) * (-t108 * t43 + t77 * t46)) * g(3) * t34, -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t94)];
taug  = t8(:);
