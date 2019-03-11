% Calculate Gravitation load on the joints for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:44
% EndTime: 2019-03-10 03:46:47
% DurationCPUTime: 0.88s
% Computational Cost: add. (731->182), mult. (698->244), div. (0->0), fcn. (676->12), ass. (0->77)
t61 = -pkin(9) - pkin(8);
t79 = rSges(5,3) - t61;
t53 = -pkin(10) + t61;
t78 = rSges(6,3) - t53;
t77 = rSges(7,3) + pkin(11) - t53;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t95 = g(1) * t60 + g(2) * t57;
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t84 = rSges(4,3) + pkin(8);
t94 = t59 * pkin(2) + t84 * t56;
t54 = qJ(3) + qJ(4);
t47 = qJ(5) + t54;
t44 = qJ(6) + t47;
t38 = sin(t44);
t39 = cos(t44);
t81 = t57 * t59;
t5 = t38 * t81 + t39 * t60;
t6 = t38 * t60 - t39 * t81;
t93 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t80 = t59 * t60;
t7 = -t38 * t80 + t39 * t57;
t8 = t38 * t57 + t39 * t80;
t92 = t7 * rSges(7,1) - t8 * rSges(7,2);
t55 = sin(qJ(3));
t91 = pkin(3) * t55;
t45 = sin(t54);
t90 = pkin(4) * t45;
t41 = sin(t47);
t89 = pkin(5) * t41;
t86 = g(3) * t56;
t42 = cos(t47);
t13 = t41 * t81 + t42 * t60;
t14 = t41 * t60 - t42 * t81;
t83 = -t13 * rSges(6,1) + t14 * rSges(6,2);
t82 = rSges(3,2) * t56;
t15 = -t41 * t80 + t42 * t57;
t16 = t41 * t57 + t42 * t80;
t76 = t15 * rSges(6,1) - t16 * rSges(6,2);
t46 = cos(t54);
t40 = pkin(4) * t46;
t58 = cos(qJ(3));
t50 = t58 * pkin(3);
t36 = t40 + t50;
t75 = t60 * pkin(1) + t57 * pkin(7);
t37 = pkin(5) * t42;
t27 = t37 + t36;
t32 = -t89 - t90;
t74 = rSges(3,1) * t59 - t82;
t72 = -rSges(5,1) * t45 - rSges(5,2) * t46;
t71 = -rSges(6,1) * t41 - rSges(6,2) * t42;
t70 = -rSges(7,1) * t38 - rSges(7,2) * t39;
t69 = rSges(4,1) * t58 - rSges(4,2) * t55 + pkin(2);
t24 = -t45 * t80 + t46 * t57;
t22 = t45 * t81 + t46 * t60;
t30 = -t55 * t80 + t57 * t58;
t28 = t55 * t81 + t58 * t60;
t43 = t50 + pkin(2);
t68 = rSges(5,1) * t46 - rSges(5,2) * t45 + t43;
t34 = pkin(2) + t36;
t67 = rSges(6,1) * t42 - rSges(6,2) * t41 + t34;
t21 = pkin(2) + t27;
t66 = rSges(7,1) * t39 - rSges(7,2) * t38 + t21;
t65 = t59 * t43 + t79 * t56;
t64 = t34 * t59 + t78 * t56;
t63 = t21 * t59 + t77 * t56;
t23 = t45 * t60 - t46 * t81;
t25 = t45 * t57 + t46 * t80;
t62 = g(1) * (t24 * rSges(5,1) - t25 * rSges(5,2)) + g(2) * (-t22 * rSges(5,1) + t23 * rSges(5,2));
t51 = t60 * pkin(7);
t35 = t90 + t91;
t33 = t37 + t40;
t31 = t55 * t57 + t58 * t80;
t29 = t55 * t60 - t58 * t81;
t26 = -t32 + t91;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t57 - rSges(2,2) * t60) + g(2) * (rSges(2,1) * t60 - rSges(2,2) * t57)) - m(3) * (g(1) * (rSges(3,3) * t60 + t51) + g(2) * (rSges(3,1) * t80 - t60 * t82 + t75) + (g(1) * (-pkin(1) - t74) + g(2) * rSges(3,3)) * t57) - m(4) * ((rSges(4,1) * t31 + rSges(4,2) * t30 + t94 * t60 + t75) * g(2) + (rSges(4,1) * t29 + rSges(4,2) * t28 + t51 + (-pkin(1) - t94) * t57) * g(1)) - m(5) * (g(1) * (t23 * rSges(5,1) + t22 * rSges(5,2) + t51) + g(2) * (t25 * rSges(5,1) + t24 * rSges(5,2) + t75) + (g(1) * t91 + g(2) * t65) * t60 + (g(1) * (-pkin(1) - t65) + g(2) * t91) * t57) - m(6) * (g(1) * (rSges(6,1) * t14 + rSges(6,2) * t13 + t51) + g(2) * (rSges(6,1) * t16 + rSges(6,2) * t15 + t75) + (g(1) * t35 + g(2) * t64) * t60 + (g(1) * (-pkin(1) - t64) + g(2) * t35) * t57) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t51) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t75) + (g(1) * t26 + g(2) * t63) * t60 + (g(1) * (-pkin(1) - t63) + g(2) * t26) * t57) -m(3) * (g(3) * t74 + t95 * (-rSges(3,1) * t56 - rSges(3,2) * t59)) - m(4) * ((g(3) * t69 + t95 * t84) * t59 + (g(3) * t84 - t95 * t69) * t56) - m(5) * ((g(3) * t68 + t95 * t79) * t59 + (g(3) * t79 - t95 * t68) * t56) - m(6) * ((g(3) * t67 + t95 * t78) * t59 + (g(3) * t78 - t95 * t67) * t56) - m(7) * ((g(3) * t66 + t95 * t77) * t59 + (g(3) * t77 - t95 * t66) * t56) -m(4) * (g(1) * (rSges(4,1) * t30 - rSges(4,2) * t31) + g(2) * (-rSges(4,1) * t28 + rSges(4,2) * t29) + (-rSges(4,1) * t55 - rSges(4,2) * t58) * t86) - m(5) * ((g(1) * t30 - g(2) * t28) * pkin(3) + (t72 - t91) * t86 + t62) - m(6) * (g(1) * (-t35 * t80 + t36 * t57 + t76) + g(2) * (-t35 * t81 - t36 * t60 + t83) + (-t35 + t71) * t86) - m(7) * (g(1) * (-t26 * t80 + t27 * t57 + t92) + g(2) * (-t26 * t81 - t27 * t60 + t93) + (-t26 + t70) * t86) -m(5) * t62 - m(6) * (g(1) * (pkin(4) * t24 + t76) + g(2) * (-pkin(4) * t22 + t83)) - m(7) * (g(1) * (t32 * t80 + t33 * t57 + t92) + g(2) * (t32 * t81 - t33 * t60 + t93)) + (-m(5) * t72 - m(6) * (t71 - t90) - m(7) * (t32 + t70)) * t86, -m(6) * (g(1) * t76 + g(2) * t83) - m(7) * (g(1) * (pkin(5) * t15 + t92) + g(2) * (-pkin(5) * t13 + t93)) + (-m(6) * t71 - m(7) * (t70 - t89)) * t86, -m(7) * (g(1) * t92 + g(2) * t93 + t70 * t86)];
taug  = t1(:);
