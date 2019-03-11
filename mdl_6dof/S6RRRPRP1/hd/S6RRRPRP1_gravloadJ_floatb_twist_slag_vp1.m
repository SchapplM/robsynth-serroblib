% Calculate Gravitation load on the joints for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:44
% EndTime: 2019-03-09 16:30:46
% DurationCPUTime: 0.77s
% Computational Cost: add. (533->143), mult. (483->184), div. (0->0), fcn. (433->10), ass. (0->70)
t98 = rSges(6,3) + pkin(9);
t97 = rSges(7,1) + pkin(5);
t40 = qJ(2) + qJ(3);
t35 = pkin(10) + t40;
t30 = sin(t35);
t31 = cos(t35);
t45 = cos(qJ(5));
t33 = pkin(5) * t45 + pkin(4);
t96 = rSges(7,3) * t30 + t31 * t33;
t60 = rSges(5,1) * t31 - rSges(5,2) * t30;
t36 = sin(t40);
t37 = cos(t40);
t65 = rSges(4,1) * t37 - rSges(4,2) * t36;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t95 = g(1) * t47 + g(2) * t44;
t56 = t31 * pkin(4) + t30 * t98;
t94 = t95 * t30;
t48 = -pkin(8) - pkin(7);
t43 = sin(qJ(2));
t93 = pkin(2) * t43;
t92 = pkin(3) * t36;
t42 = sin(qJ(5));
t91 = pkin(5) * t42;
t88 = rSges(3,3) + pkin(7);
t79 = t42 * t44;
t71 = t30 * t79;
t82 = t31 * t44;
t87 = rSges(7,2) * t71 + rSges(7,3) * t82;
t46 = cos(qJ(2));
t38 = t46 * pkin(2);
t34 = t38 + pkin(1);
t41 = -qJ(6) - pkin(9);
t84 = t30 * t41;
t83 = t31 * t42;
t81 = t31 * t45;
t80 = t31 * t47;
t78 = t42 * t47;
t77 = t44 * t45;
t76 = t45 * t47;
t75 = rSges(4,3) - t48;
t39 = -qJ(4) + t48;
t74 = rSges(5,3) - t39;
t70 = t30 * t78;
t73 = rSges(7,2) * t70 + rSges(7,3) * t80;
t72 = rSges(6,2) * t71 + t82 * t98;
t69 = rSges(6,2) * t70 + t80 * t98;
t67 = -t39 + t91;
t66 = -rSges(7,1) * t45 - t33;
t32 = pkin(3) * t37;
t64 = t32 + t60;
t63 = rSges(3,1) * t46 - rSges(3,2) * t43;
t61 = -rSges(4,1) * t36 - rSges(4,2) * t37;
t59 = -rSges(5,1) * t30 - rSges(5,2) * t31;
t58 = pkin(1) + t63;
t3 = -t31 * t78 + t77;
t1 = t31 * t79 + t76;
t57 = t34 + t65;
t55 = rSges(7,1) * t81 - rSges(7,2) * t83 + t32 + t96;
t52 = -t84 + t96;
t51 = rSges(6,1) * t81 - rSges(6,2) * t83 + t32 + t56;
t49 = (-rSges(6,1) * t45 - pkin(4)) * t94;
t14 = -t92 - t93;
t13 = t32 + t34;
t7 = t47 * t14;
t6 = t44 * t14;
t5 = t47 * t13;
t4 = t31 * t76 + t79;
t2 = -t31 * t77 + t78;
t8 = [-m(2) * (g(1) * (-rSges(2,1) * t44 - rSges(2,2) * t47) + g(2) * (rSges(2,1) * t47 - rSges(2,2) * t44)) - m(3) * ((g(1) * t88 + g(2) * t58) * t47 + (-g(1) * t58 + g(2) * t88) * t44) - m(4) * ((g(1) * t75 + g(2) * t57) * t47 + (-g(1) * t57 + g(2) * t75) * t44) - m(5) * (g(2) * t5 + (g(1) * t74 + g(2) * t60) * t47 + (g(1) * (-t13 - t60) + g(2) * t74) * t44) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t39 + g(2) * t56) * t47 + (g(1) * (-t13 - t56) - g(2) * t39) * t44) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t5) + (g(1) * t67 + g(2) * t52) * t47 + (g(1) * (-t13 - t52) + g(2) * t67) * t44) -m(3) * (g(3) * t63 + t95 * (-rSges(3,1) * t43 - rSges(3,2) * t46)) - m(4) * (g(3) * (t38 + t65) + t95 * (t61 - t93)) - m(5) * (g(1) * (t47 * t59 + t7) + g(2) * (t44 * t59 + t6) + g(3) * (t38 + t64)) - m(6) * (g(1) * (t7 + t69) + g(2) * (t6 + t72) + g(3) * (t38 + t51) + t49) - m(7) * (g(1) * (-t41 * t80 + t7 + t73) + g(2) * (-t41 * t82 + t6 + t87) + g(3) * (t38 + t55) + (-g(3) * t41 + t66 * t95) * t30) -m(4) * (g(3) * t65 + t61 * t95) - m(5) * (g(3) * t64 + t95 * (t59 - t92)) - m(6) * (g(1) * (-t47 * t92 + t69) + g(2) * (-t44 * t92 + t72) + g(3) * t51 + t49) - m(7) * (g(1) * t73 + g(2) * t87 + g(3) * (t55 - t84) + t95 * (t30 * t66 - t31 * t41 - t92)) (-m(5) - m(6) - m(7)) * (g(1) * t44 - g(2) * t47) -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2)) - m(7) * (g(1) * (-rSges(7,2) * t4 + t3 * t97) + g(2) * (rSges(7,2) * t2 - t1 * t97)) + (-m(6) * (-rSges(6,1) * t42 - rSges(6,2) * t45) - m(7) * (-rSges(7,1) * t42 - rSges(7,2) * t45 - t91)) * g(3) * t30, -m(7) * (-g(3) * t31 + t94)];
taug  = t8(:);
