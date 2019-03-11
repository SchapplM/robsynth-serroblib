% Calculate Gravitation load on the joints for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:51
% EndTime: 2019-03-09 15:31:55
% DurationCPUTime: 1.30s
% Computational Cost: add. (552->177), mult. (774->236), div. (0->0), fcn. (805->10), ass. (0->73)
t37 = qJ(3) + pkin(10);
t32 = sin(t37);
t33 = cos(t37);
t46 = cos(qJ(1));
t42 = sin(qJ(1));
t45 = cos(qJ(2));
t82 = t42 * t45;
t11 = t32 * t82 + t33 * t46;
t80 = t46 * t32;
t12 = t33 * t82 - t80;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t102 = t11 * t43 - t12 * t39;
t58 = t11 * t39 + t12 * t43;
t103 = t102 * rSges(7,1) - t58 * rSges(7,2);
t40 = sin(qJ(3));
t44 = cos(qJ(3));
t81 = t45 * t46;
t18 = -t40 * t81 + t42 * t44;
t31 = pkin(3) * t44 + pkin(2);
t24 = t45 * t31;
t71 = -pkin(1) - t24;
t93 = g(2) * t42;
t100 = g(1) * t46 + t93;
t41 = sin(qJ(2));
t89 = rSges(4,3) + pkin(8);
t99 = t45 * pkin(2) + t41 * t89;
t98 = -m(6) - m(7);
t97 = -pkin(4) - pkin(5);
t96 = pkin(3) * t40;
t95 = g(1) * t42;
t92 = g(3) * t41;
t90 = -rSges(6,1) - pkin(4);
t88 = rSges(7,3) + pkin(9);
t87 = rSges(3,2) * t41;
t85 = t40 * t46;
t84 = t42 * t40;
t79 = t46 * t44;
t38 = -qJ(4) - pkin(8);
t78 = rSges(6,2) - t38;
t77 = rSges(5,3) - t38;
t76 = t46 * pkin(1) + t42 * pkin(7);
t75 = rSges(6,3) + qJ(5);
t74 = -t38 - t88;
t35 = t46 * pkin(7);
t72 = t42 * t41 * t38 + pkin(3) * t85 + t35;
t69 = t78 * t46;
t68 = t77 * t46;
t67 = pkin(3) * t84 + t31 * t81 + t76;
t66 = g(3) * (t24 + (t33 * pkin(4) + t32 * qJ(5)) * t45);
t65 = t74 * t46;
t14 = t42 * t32 + t33 * t81;
t64 = t14 * pkin(4) + t67;
t13 = -t42 * t33 + t45 * t80;
t2 = t13 * t43 - t14 * t39;
t3 = t13 * t39 + t14 * t43;
t63 = rSges(7,1) * t2 - rSges(7,2) * t3;
t56 = t32 * t39 + t33 * t43;
t57 = t32 * t43 - t33 * t39;
t62 = (rSges(7,1) * t57 - rSges(7,2) * t56) * t41;
t61 = rSges(3,1) * t45 - t87;
t59 = rSges(5,1) * t33 - rSges(5,2) * t32;
t55 = t18 * pkin(3);
t54 = rSges(4,1) * t44 - rSges(4,2) * t40 + pkin(2);
t16 = t40 * t82 + t79;
t52 = -t13 * pkin(4) + t55;
t51 = -t12 * pkin(4) - t11 * qJ(5) + t72;
t50 = t16 * pkin(3);
t49 = -t11 * pkin(4) - t50;
t21 = t41 * t33 * qJ(5);
t19 = t45 * t79 + t84;
t17 = -t44 * t82 + t85;
t1 = [-m(2) * (g(1) * (-t42 * rSges(2,1) - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - t42 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t46 + t35) + g(2) * (rSges(3,1) * t81 - t46 * t87 + t76) + (g(1) * (-pkin(1) - t61) + g(2) * rSges(3,3)) * t42) - m(4) * (g(1) * (rSges(4,1) * t17 + rSges(4,2) * t16 + t35) + (-pkin(1) - t99) * t95 + (t19 * rSges(4,1) + t18 * rSges(4,2) + t99 * t46 + t76) * g(2)) - m(5) * (g(1) * (-rSges(5,1) * t12 + rSges(5,2) * t11 + t72) + g(2) * (t14 * rSges(5,1) - t13 * rSges(5,2) + t41 * t68 + t67) + (-rSges(5,3) * t41 + t71) * t95) - m(6) * (g(1) * (-rSges(6,1) * t12 - rSges(6,3) * t11 + t51) + g(2) * (t14 * rSges(6,1) + t13 * t75 + t41 * t69 + t64) + (-rSges(6,2) * t41 + t71) * t95) - m(7) * (g(1) * (-rSges(7,1) * t58 - rSges(7,2) * t102 - t12 * pkin(5) + t71 * t42 + t51) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t14 * pkin(5) + t13 * qJ(5) + t64) + (g(2) * t65 + t88 * t95) * t41) -m(3) * (g(3) * t61 + t100 * (-rSges(3,1) * t41 - rSges(3,2) * t45)) - m(4) * ((g(3) * t54 + t100 * t89) * t45 + (g(3) * t89 - t100 * t54) * t41) - m(5) * (g(3) * t24 + (g(1) * t68 + g(3) * t59 + t77 * t93) * t45 + (g(3) * t77 + t100 * (-t31 - t59)) * t41) - m(6) * (t66 + (g(3) * (rSges(6,1) * t33 + rSges(6,3) * t32) + g(1) * t69 + t78 * t93) * t45 + (g(3) * t78 + t100 * (-t32 * t75 + t33 * t90 - t31)) * t41) - m(7) * (t66 + (g(3) * (rSges(7,1) * t56 + rSges(7,2) * t57 + t33 * pkin(5)) + g(1) * t65 + t74 * t93) * t45 + (g(3) * t74 + t100 * (-t31 + (-t39 * rSges(7,1) - t43 * rSges(7,2) - qJ(5)) * t32 + (-t43 * rSges(7,1) + t39 * rSges(7,2) + t97) * t33)) * t41) -m(4) * (g(1) * (rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t16 + rSges(4,2) * t17) + (-rSges(4,1) * t40 - rSges(4,2) * t44) * t92) - m(5) * (g(1) * (-t13 * rSges(5,1) - t14 * rSges(5,2) + t55) + g(2) * (-t11 * rSges(5,1) - t12 * rSges(5,2) - t50) + (-rSges(5,1) * t32 - rSges(5,2) * t33 - t96) * t92) - m(6) * (g(1) * (-t13 * rSges(6,1) + t14 * t75 + t52) + g(2) * (-t11 * rSges(6,1) + t12 * t75 + t49) + g(3) * t21 + (rSges(6,3) * t33 + t32 * t90 - t96) * t92) - m(7) * (g(1) * (-t13 * pkin(5) + t14 * qJ(5) + t52 - t63) + g(2) * (-t11 * pkin(5) + t12 * qJ(5) - t103 + t49) + g(3) * (t21 - t62) + (t32 * t97 - t96) * t92) (-m(5) + t98) * (-g(3) * t45 + t100 * t41) t98 * (g(1) * t13 + g(2) * t11 + t32 * t92) -m(7) * (g(1) * t63 + g(2) * t103 + g(3) * t62)];
taug  = t1(:);
