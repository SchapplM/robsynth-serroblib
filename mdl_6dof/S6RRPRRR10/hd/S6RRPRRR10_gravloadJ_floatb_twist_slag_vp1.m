% Calculate Gravitation load on the joints for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:21
% EndTime: 2019-03-09 14:18:25
% DurationCPUTime: 1.30s
% Computational Cost: add. (835->179), mult. (1342->252), div. (0->0), fcn. (1583->14), ass. (0->75)
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t59 = cos(qJ(2));
t80 = cos(pkin(6));
t96 = cos(qJ(1));
t72 = t80 * t96;
t29 = t56 * t72 + t57 * t59;
t49 = pkin(12) + qJ(4);
t44 = sin(t49);
t45 = cos(t49);
t52 = sin(pkin(6));
t78 = t52 * t96;
t13 = t29 * t45 - t44 * t78;
t28 = t56 * t57 - t59 * t72;
t50 = qJ(5) + qJ(6);
t46 = sin(t50);
t47 = cos(t50);
t117 = t13 * t46 - t28 * t47;
t116 = -t13 * t47 - t28 * t46;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t115 = t13 * t55 - t28 * t58;
t91 = t28 * t55;
t114 = -t13 * t58 - t91;
t86 = t52 * t59;
t113 = g(3) * t86;
t112 = rSges(6,1) * t55 + rSges(6,2) * t58;
t101 = g(2) * t28;
t75 = t57 * t80;
t30 = t96 * t56 + t59 * t75;
t111 = g(1) * t30 + t101;
t31 = -t56 * t75 + t96 * t59;
t87 = t52 * t57;
t17 = t31 * t45 + t44 * t87;
t88 = t52 * t56;
t23 = t80 * t44 + t45 * t88;
t110 = g(1) * t17 + g(2) * t13 + g(3) * t23;
t16 = t31 * t44 - t45 * t87;
t22 = -t44 * t88 + t80 * t45;
t76 = -t29 * t44 - t45 * t78;
t109 = -g(1) * t16 + g(2) * t76 + g(3) * t22;
t43 = pkin(5) * t58 + pkin(4);
t66 = rSges(7,1) * t47 - rSges(7,2) * t46 + t43;
t82 = pkin(11) + pkin(10) + rSges(7,3);
t108 = t82 * t44 + t66 * t45;
t68 = rSges(6,1) * t58 - rSges(6,2) * t55 + pkin(4);
t97 = pkin(10) + rSges(6,3);
t107 = t97 * t44 + t68 * t45;
t100 = g(2) * t29;
t53 = cos(pkin(12));
t42 = pkin(3) * t53 + pkin(2);
t99 = t42 * t113;
t98 = g(3) * t52;
t89 = t30 * t55;
t54 = -pkin(9) - qJ(3);
t85 = -t28 * t42 - t29 * t54;
t84 = -t30 * t42 - t31 * t54;
t83 = t96 * pkin(1) + pkin(8) * t87;
t81 = qJ(3) + rSges(4,3);
t51 = sin(pkin(12));
t79 = t51 * t87;
t77 = -t57 * pkin(1) + pkin(8) * t78;
t74 = t51 * t78;
t73 = pkin(3) * t79 - t30 * t54 + t31 * t42 + t83;
t71 = rSges(5,1) * t45 - rSges(5,2) * t44;
t7 = -t17 * t55 + t30 * t58;
t69 = rSges(4,1) * t53 - rSges(4,2) * t51 + pkin(2);
t67 = -t23 * t55 - t58 * t86;
t65 = pkin(3) * t74 + t28 * t54 - t29 * t42 + t77;
t64 = rSges(7,1) * t46 + rSges(7,2) * t47 + pkin(5) * t55;
t5 = -t17 * t46 + t30 * t47;
t6 = t17 * t47 + t30 * t46;
t61 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (-rSges(7,1) * t117 + t116 * rSges(7,2)) + g(3) * ((-t23 * t46 - t47 * t86) * rSges(7,1) + (-t23 * t47 + t46 * t86) * rSges(7,2)));
t8 = t17 * t58 + t89;
t1 = [-m(2) * (g(1) * (-t57 * rSges(2,1) - t96 * rSges(2,2)) + g(2) * (t96 * rSges(2,1) - t57 * rSges(2,2))) - m(3) * (g(1) * (-t29 * rSges(3,1) + t28 * rSges(3,2) + rSges(3,3) * t78 + t77) + g(2) * (rSges(3,1) * t31 - rSges(3,2) * t30 + rSges(3,3) * t87 + t83)) - m(4) * (g(1) * (-t29 * pkin(2) + (-t29 * t53 + t74) * rSges(4,1) + (t29 * t51 + t53 * t78) * rSges(4,2) - t81 * t28 + t77) + g(2) * (t31 * pkin(2) + (t31 * t53 + t79) * rSges(4,1) + (-t31 * t51 + t53 * t87) * rSges(4,2) + t81 * t30 + t83)) - m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t76 - rSges(5,3) * t28 + t65) + g(2) * (rSges(5,1) * t17 - rSges(5,2) * t16 + rSges(5,3) * t30 + t73)) - m(6) * (g(1) * (rSges(6,1) * t114 + rSges(6,2) * t115 - t13 * pkin(4) + t97 * t76 + t65) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 + pkin(4) * t17 + t97 * t16 + t73)) - m(7) * (g(1) * (t116 * rSges(7,1) + rSges(7,2) * t117 - pkin(5) * t91 - t13 * t43 + t82 * t76 + t65) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + pkin(5) * t89 + t82 * t16 + t17 * t43 + t73)) -m(3) * (g(1) * (-rSges(3,1) * t30 - rSges(3,2) * t31) + g(2) * (-rSges(3,1) * t28 - rSges(3,2) * t29) + (rSges(3,1) * t59 - rSges(3,2) * t56) * t98) - m(4) * (g(1) * (-t69 * t30 + t81 * t31) + t81 * t100 - t69 * t101 + (t81 * t56 + t69 * t59) * t98) - m(5) * (g(1) * (rSges(5,3) * t31 - t71 * t30 + t84) + g(2) * (rSges(5,3) * t29 - t71 * t28 + t85) + t99 + (t71 * t59 + (rSges(5,3) - t54) * t56) * t98) - m(6) * (g(1) * (t112 * t31 + t84) + g(2) * (t112 * t29 + t85) + t99 + ((-t54 + t112) * t56 + t107 * t59) * t98 - t111 * t107) - m(7) * (g(2) * t85 + t99 + t64 * t100 + ((-t54 + t64) * t56 + t108 * t59) * t98 - t111 * t108 + (t64 * t31 + t84) * g(1)) (-m(4) - m(5) - m(6) - m(7)) * (t111 - t113) -m(5) * (g(1) * (-rSges(5,1) * t16 - rSges(5,2) * t17) + g(2) * (rSges(5,1) * t76 - rSges(5,2) * t13) + g(3) * (rSges(5,1) * t22 - rSges(5,2) * t23)) - m(6) * (t109 * t68 + t110 * t97) - m(7) * (t109 * t66 + t110 * t82) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t115 + rSges(6,2) * t114) + g(3) * (t67 * rSges(6,1) + (-t23 * t58 + t55 * t86) * rSges(6,2))) - t61 - m(7) * (g(1) * t7 - g(2) * t115 + g(3) * t67) * pkin(5), -t61];
taug  = t1(:);
