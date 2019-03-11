% Calculate Gravitation load on the joints for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:20
% EndTime: 2019-03-09 10:45:22
% DurationCPUTime: 0.91s
% Computational Cost: add. (448->162), mult. (707->202), div. (0->0), fcn. (733->10), ass. (0->81)
t46 = sin(qJ(2));
t50 = cos(qJ(2));
t81 = qJ(4) + pkin(10);
t75 = sin(t81);
t76 = cos(t81);
t12 = t46 * t75 + t50 * t76;
t13 = t46 * t76 - t50 * t75;
t47 = sin(qJ(1));
t3 = t13 * t47;
t44 = sin(qJ(6));
t48 = cos(qJ(6));
t51 = cos(qJ(1));
t62 = t51 * t75;
t63 = t51 * t76;
t6 = -t46 * t63 + t50 * t62;
t111 = (g(1) * t6 - g(2) * t3 + g(3) * t12) * (t48 * rSges(7,1) - t44 * rSges(7,2) + pkin(5));
t110 = -rSges(5,3) - pkin(8);
t109 = -pkin(9) - rSges(7,3);
t49 = cos(qJ(4));
t37 = pkin(4) * t49 + pkin(3);
t27 = t50 * t37;
t38 = t46 * qJ(3);
t84 = t50 * pkin(2) + t38;
t108 = t27 + t84;
t107 = g(1) * t51 + g(2) * t47;
t105 = t107 * t46;
t102 = -m(6) - m(7);
t100 = pkin(3) * t50;
t45 = sin(qJ(4));
t99 = pkin(4) * t45;
t98 = g(1) * t47;
t93 = -pkin(2) - t37;
t92 = rSges(4,1) * t50;
t91 = t45 * t46;
t90 = t46 * t49;
t89 = t46 * t51;
t88 = t47 * t50;
t87 = t49 * t50;
t86 = t50 * t51;
t41 = t51 * pkin(7);
t43 = -qJ(5) - pkin(8);
t85 = t51 * t43 + t41;
t83 = t51 * pkin(1) + t47 * pkin(7);
t82 = qJ(3) * t50;
t80 = pkin(4) * t87;
t33 = pkin(4) * t91;
t79 = t45 * t86;
t78 = t49 * t89;
t77 = pkin(2) * t86 + t51 * t38 + t83;
t67 = t45 * t50 - t90;
t8 = t67 * t47;
t66 = t87 + t91;
t9 = t66 * t47;
t74 = -rSges(5,1) * t8 - rSges(5,2) * t9;
t4 = t12 * t47;
t73 = -t4 * t48 - t44 * t51;
t72 = t4 * t44 - t48 * t51;
t71 = rSges(3,1) * t50 - rSges(3,2) * t46;
t10 = -t78 + t79;
t11 = t66 * t51;
t69 = -t10 * rSges(5,1) - t11 * rSges(5,2);
t68 = -rSges(5,1) * t66 + rSges(5,2) * t67;
t22 = t88 * t99;
t61 = -t109 * t4 - t22;
t25 = pkin(4) * t79;
t5 = -t46 * t62 - t50 * t63;
t60 = t109 * t5 - t25;
t59 = rSges(6,1) * t3 - rSges(6,2) * t4 - t22;
t58 = -t6 * rSges(6,1) + t5 * rSges(6,2) - t25;
t57 = -t109 * t13 - t33;
t55 = -pkin(1) - t84;
t54 = -rSges(6,1) * t12 - rSges(6,2) * t13 - t33;
t53 = t51 * t33 + t37 * t86 + t47 * t43 + t77;
t52 = t93 * t105;
t32 = t51 * t82;
t30 = t47 * t82;
t24 = pkin(4) * t78;
t21 = t47 * pkin(4) * t90;
t2 = -t44 * t47 - t48 * t5;
t1 = t44 * t5 - t47 * t48;
t7 = [-m(2) * (g(1) * (-t47 * rSges(2,1) - rSges(2,2) * t51) + g(2) * (rSges(2,1) * t51 - t47 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t51 + t41) + g(2) * (rSges(3,1) * t86 - rSges(3,2) * t89 + t83) + (g(1) * (-pkin(1) - t71) + g(2) * rSges(3,3)) * t47) - m(4) * (g(1) * (rSges(4,2) * t51 + t41) + g(2) * (rSges(4,1) * t86 + rSges(4,3) * t89 + t77) + (g(1) * (-rSges(4,3) * t46 + t55 - t92) + g(2) * rSges(4,2)) * t47) - m(5) * (g(1) * (-t9 * rSges(5,1) + t8 * rSges(5,2) + t110 * t51 + t41) + g(2) * (t11 * rSges(5,1) - t10 * rSges(5,2) + pkin(3) * t86 + t77) + (g(1) * (t55 - t100) + g(2) * t110) * t47) - m(6) * (g(1) * (-t4 * rSges(6,1) - t3 * rSges(6,2) - rSges(6,3) * t51 + t85) + g(2) * (-rSges(6,1) * t5 - rSges(6,2) * t6 + t53) + (g(1) * (-t33 + t55 - t27) - g(2) * rSges(6,3)) * t47) - m(7) * (g(1) * (rSges(7,1) * t73 + rSges(7,2) * t72 - t4 * pkin(5) - t109 * t3 + t85) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 - pkin(5) * t5 - t109 * t6 + t53) + (-pkin(1) + t93 * t50 + (-qJ(3) - t99) * t46) * t98) -m(3) * (g(3) * t71 + t107 * (-rSges(3,1) * t46 - rSges(3,2) * t50)) - m(4) * (g(1) * (rSges(4,3) * t86 + t32) + g(2) * (rSges(4,3) * t88 + t30) + g(3) * (t84 + t92) + (g(3) * rSges(4,3) + t107 * (-rSges(4,1) - pkin(2))) * t46) - m(5) * (g(1) * (t32 - t69) + g(2) * (t30 - t74) + g(3) * (-t68 + t84 + t100) + (-pkin(2) - pkin(3)) * t105) - m(6) * (g(1) * (t32 - t58) + g(2) * (t30 - t59) + g(3) * (-t54 + t108) + t52) - m(7) * (g(1) * (t32 - t60) + g(2) * (t30 - t61) + g(3) * (-t57 + t108) + t52 + t111) (-m(4) - m(5) + t102) * (-g(3) * t50 + t105) -m(5) * (g(1) * t69 + g(2) * t74 + g(3) * t68) - m(6) * (g(1) * (t24 + t58) + g(2) * (t21 + t59) + g(3) * (t54 - t80)) + (-g(1) * (t24 + t60) - g(2) * (t21 + t61) - g(3) * (t57 - t80) + t111) * m(7), t102 * (g(2) * t51 - t98) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-rSges(7,1) * t72 + rSges(7,2) * t73) + g(3) * (-t44 * rSges(7,1) - t48 * rSges(7,2)) * t13)];
taug  = t7(:);
