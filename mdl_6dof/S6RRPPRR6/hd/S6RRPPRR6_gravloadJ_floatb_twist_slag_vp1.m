% Calculate Gravitation load on the joints for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:06
% EndTime: 2019-03-09 09:13:08
% DurationCPUTime: 0.91s
% Computational Cost: add. (426->148), mult. (648->197), div. (0->0), fcn. (671->10), ass. (0->70)
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t70 = pkin(10) + qJ(5);
t64 = sin(t70);
t65 = cos(t70);
t12 = t42 * t64 + t45 * t65;
t13 = t42 * t65 - t45 * t64;
t43 = sin(qJ(1));
t3 = t13 * t43;
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t46 = cos(qJ(1));
t51 = t46 * t64;
t52 = t46 * t65;
t6 = -t42 * t52 + t45 * t51;
t98 = (g(1) * t6 - g(2) * t3 + g(3) * t12) * (t44 * rSges(7,1) - t41 * rSges(7,2) + pkin(5));
t84 = rSges(7,3) + pkin(9);
t38 = sin(pkin(10));
t81 = t38 * t42;
t25 = pkin(4) * t81;
t39 = cos(pkin(10));
t30 = pkin(4) * t39 + pkin(3);
t97 = t45 * t30 + t25;
t33 = t42 * qJ(3);
t74 = t45 * pkin(2) + t33;
t96 = g(1) * t46 + g(2) * t43;
t95 = -rSges(5,3) - qJ(4);
t94 = t96 * t42;
t91 = g(1) * t43;
t87 = g(3) * t13;
t86 = t45 * pkin(3);
t83 = -pkin(2) - t30;
t82 = rSges(4,1) * t45;
t80 = t38 * t45;
t79 = t42 * t46;
t78 = t45 * t46;
t72 = qJ(3) * t45;
t26 = t43 * t72;
t69 = pkin(4) * t80;
t77 = t43 * t69 + t26;
t28 = t46 * t72;
t76 = t46 * t69 + t28;
t36 = t46 * pkin(7);
t40 = -pkin(8) - qJ(4);
t75 = t46 * t40 + t36;
t73 = t46 * pkin(1) + t43 * pkin(7);
t71 = -m(5) - m(6) - m(7);
t68 = t42 * (-pkin(2) - pkin(3));
t67 = t74 + t97;
t66 = pkin(2) * t78 + t46 * t33 + t73;
t4 = t12 * t43;
t63 = rSges(6,1) * t3 - rSges(6,2) * t4;
t5 = -t42 * t51 - t45 * t52;
t62 = -t6 * rSges(6,1) + t5 * rSges(6,2);
t61 = -t4 * t44 - t41 * t46;
t60 = t4 * t41 - t44 * t46;
t59 = rSges(3,1) * t45 - rSges(3,2) * t42;
t57 = -rSges(6,1) * t12 - rSges(6,2) * t13;
t56 = -t39 * t42 + t80;
t55 = t39 * t45 + t81;
t49 = -pkin(1) - t74;
t48 = t46 * t25 + t30 * t78 + t43 * t40 + t66;
t47 = t83 * t94;
t11 = t55 * t46;
t10 = t56 * t46;
t9 = t55 * t43;
t8 = t56 * t43;
t2 = -t41 * t43 - t44 * t5;
t1 = t41 * t5 - t43 * t44;
t7 = [-m(2) * (g(1) * (-t43 * rSges(2,1) - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - t43 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t46 + t36) + g(2) * (rSges(3,1) * t78 - rSges(3,2) * t79 + t73) + (g(1) * (-pkin(1) - t59) + g(2) * rSges(3,3)) * t43) - m(4) * (g(1) * (rSges(4,2) * t46 + t36) + g(2) * (rSges(4,1) * t78 + rSges(4,3) * t79 + t66) + (g(1) * (-rSges(4,3) * t42 + t49 - t82) + g(2) * rSges(4,2)) * t43) - m(5) * (g(1) * (-t9 * rSges(5,1) + t8 * rSges(5,2) + t95 * t46 + t36) + g(2) * (t11 * rSges(5,1) - t10 * rSges(5,2) + pkin(3) * t78 + t66) + (g(1) * (t49 - t86) + g(2) * t95) * t43) - m(6) * (g(1) * (-t4 * rSges(6,1) - t3 * rSges(6,2) - rSges(6,3) * t46 + t75) + g(2) * (-rSges(6,1) * t5 - rSges(6,2) * t6 + t48) + (g(1) * (t49 - t97) - g(2) * rSges(6,3)) * t43) - m(7) * (g(1) * (t61 * rSges(7,1) + t60 * rSges(7,2) - t4 * pkin(5) + t84 * t3 + t75) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 - pkin(5) * t5 + t84 * t6 + t48) + (-pkin(1) + t83 * t45 + (-pkin(4) * t38 - qJ(3)) * t42) * t91) -m(3) * (g(3) * t59 + t96 * (-rSges(3,1) * t42 - rSges(3,2) * t45)) - m(4) * (g(1) * (rSges(4,3) * t78 + t28) + g(2) * (rSges(4,3) * t43 * t45 + t26) + g(3) * (t74 + t82) + (g(3) * rSges(4,3) + t96 * (-rSges(4,1) - pkin(2))) * t42) - m(5) * (g(1) * (t10 * rSges(5,1) + t11 * rSges(5,2) + t46 * t68 + t28) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t9 + t43 * t68 + t26) + g(3) * (t55 * rSges(5,1) - t56 * rSges(5,2) + t74 + t86)) - m(6) * (g(1) * (-t62 + t76) + g(2) * (-t63 + t77) + g(3) * (-t57 + t67) + t47) + (-g(1) * (t84 * t5 + t76) - g(2) * (-t84 * t4 + t77) - g(3) * (-t84 * t13 + t67) - t47 - t98) * m(7) (-m(4) + t71) * (-g(3) * t45 + t94) t71 * (g(2) * t46 - t91) -m(6) * (g(1) * t62 + g(2) * t63 + g(3) * t57) - m(7) * ((-g(1) * t5 + g(2) * t4 + t87) * t84 - t98) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t60 * rSges(7,1) + t61 * rSges(7,2)) + (-t41 * rSges(7,1) - t44 * rSges(7,2)) * t87)];
taug  = t7(:);
