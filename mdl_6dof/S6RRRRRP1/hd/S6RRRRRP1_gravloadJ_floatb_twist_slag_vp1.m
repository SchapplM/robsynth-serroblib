% Calculate Gravitation load on the joints for
% S6RRRRRP1
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:35
% EndTime: 2019-03-10 00:55:37
% DurationCPUTime: 0.87s
% Computational Cost: add. (641->150), mult. (561->197), div. (0->0), fcn. (502->10), ass. (0->75)
t105 = rSges(6,3) + pkin(10);
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t39 = qJ(2) + qJ(3);
t34 = sin(t39);
t35 = cos(t39);
t69 = t35 * rSges(4,1) - rSges(4,2) * t34;
t104 = t37 + t69;
t103 = rSges(7,1) + pkin(5);
t36 = qJ(4) + t39;
t30 = sin(t36);
t31 = cos(t36);
t44 = cos(qJ(5));
t32 = pkin(5) * t44 + pkin(4);
t102 = t30 * rSges(7,3) + t31 * t32;
t101 = t31 * rSges(5,1) - rSges(5,2) * t30;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t100 = g(1) * t46 + g(2) * t43;
t58 = t31 * pkin(4) + t105 * t30;
t99 = t100 * t30;
t47 = -pkin(8) - pkin(7);
t42 = sin(qJ(2));
t98 = pkin(2) * t42;
t97 = pkin(3) * t34;
t41 = sin(qJ(5));
t96 = pkin(5) * t41;
t93 = rSges(3,3) + pkin(7);
t84 = t41 * t43;
t75 = t30 * t84;
t87 = t31 * t43;
t92 = rSges(7,2) * t75 + rSges(7,3) * t87;
t40 = -qJ(6) - pkin(10);
t89 = t30 * t40;
t88 = t31 * t41;
t86 = t31 * t44;
t85 = t31 * t46;
t83 = t41 * t46;
t82 = t43 * t44;
t81 = t44 * t46;
t80 = rSges(4,3) - t47;
t38 = -pkin(9) + t47;
t79 = rSges(5,3) - t38;
t74 = t30 * t83;
t78 = rSges(7,2) * t74 + rSges(7,3) * t85;
t29 = pkin(3) * t35;
t77 = t29 + t37;
t76 = rSges(6,2) * t75 + t105 * t87;
t73 = rSges(6,2) * t74 + t105 * t85;
t72 = -rSges(6,1) * t44 - pkin(4);
t71 = -t38 + t96;
t70 = -rSges(7,1) * t44 - t32;
t67 = t29 + t101;
t66 = rSges(3,1) * t45 - rSges(3,2) * t42;
t64 = -rSges(4,1) * t34 - rSges(4,2) * t35;
t62 = -rSges(5,1) * t30 - rSges(5,2) * t31;
t61 = pkin(1) + t66;
t3 = -t31 * t83 + t82;
t1 = t31 * t84 + t81;
t60 = pkin(1) + t104;
t59 = rSges(7,1) * t86 - rSges(7,2) * t88 + t102;
t57 = rSges(6,1) * t86 - rSges(6,2) * t88 + t58;
t54 = -t89 + t102;
t53 = g(1) * t78 + g(2) * t92;
t52 = t29 + t57;
t51 = t59 - t89;
t49 = t72 * t99;
t14 = -t97 - t98;
t13 = pkin(1) + t77;
t7 = t46 * t14;
t6 = t43 * t14;
t5 = t46 * t13;
t4 = t31 * t81 + t84;
t2 = -t31 * t82 + t83;
t8 = [-m(2) * (g(1) * (-rSges(2,1) * t43 - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - rSges(2,2) * t43)) - m(3) * ((g(1) * t93 + g(2) * t61) * t46 + (-g(1) * t61 + g(2) * t93) * t43) - m(4) * ((g(1) * t80 + g(2) * t60) * t46 + (-g(1) * t60 + g(2) * t80) * t43) - m(5) * (g(2) * t5 + (g(1) * t79 + g(2) * t101) * t46 + (g(1) * (-t13 - t101) + g(2) * t79) * t43) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t3 + t5) + (-g(1) * t38 + g(2) * t58) * t46 + (g(1) * (-t13 - t58) - g(2) * t38) * t43) - m(7) * (g(1) * (rSges(7,1) * t2 + rSges(7,2) * t1) + g(2) * (rSges(7,1) * t4 + rSges(7,2) * t3 + t5) + (g(1) * t71 + g(2) * t54) * t46 + (g(1) * (-t13 - t54) + g(2) * t71) * t43) -m(3) * (g(3) * t66 + t100 * (-rSges(3,1) * t42 - rSges(3,2) * t45)) - m(4) * (g(3) * t104 + t100 * (t64 - t98)) - m(5) * (g(1) * (t46 * t62 + t7) + g(2) * (t43 * t62 + t6) + g(3) * (t37 + t67)) - m(6) * (g(1) * (t7 + t73) + g(2) * (t6 + t76) + g(3) * (t37 + t52) + t49) - m(7) * (g(1) * (-t40 * t85 + t7 + t78) + g(2) * (-t40 * t87 + t6 + t92) + g(3) * (t59 + t77) + (-g(3) * t40 + t100 * t70) * t30) -m(4) * (g(3) * t69 + t100 * t64) - m(5) * (g(3) * t67 + t100 * (t62 - t97)) - m(6) * (g(1) * (-t46 * t97 + t73) + g(2) * (-t43 * t97 + t76) + g(3) * t52 + t49) - m(7) * (g(3) * (t29 + t51) + t53 + t100 * (t30 * t70 - t31 * t40 - t97)) -m(5) * g(3) * t101 - m(6) * (g(1) * t73 + g(2) * t76 + g(3) * t57) - m(7) * (g(3) * t51 + t53) + t100 * ((m(5) * rSges(5,2) + m(7) * t40) * t31 + (m(5) * rSges(5,1) - m(6) * t72 - m(7) * t70) * t30) -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2)) - m(7) * (g(1) * (-rSges(7,2) * t4 + t103 * t3) + g(2) * (rSges(7,2) * t2 - t103 * t1)) + (-m(6) * (-rSges(6,1) * t41 - rSges(6,2) * t44) - m(7) * (-rSges(7,1) * t41 - rSges(7,2) * t44 - t96)) * g(3) * t30, -m(7) * (-g(3) * t31 + t99)];
taug  = t8(:);
