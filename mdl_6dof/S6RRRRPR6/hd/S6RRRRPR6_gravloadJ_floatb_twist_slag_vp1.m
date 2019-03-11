% Calculate Gravitation load on the joints for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:42
% EndTime: 2019-03-09 22:18:44
% DurationCPUTime: 0.83s
% Computational Cost: add. (652->176), mult. (650->236), div. (0->0), fcn. (626->12), ass. (0->76)
t62 = -pkin(9) - pkin(8);
t80 = rSges(5,3) - t62;
t54 = -qJ(5) + t62;
t79 = rSges(6,3) - t54;
t78 = rSges(7,3) + pkin(10) - t54;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t95 = g(1) * t61 + g(2) * t58;
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t85 = rSges(4,3) + pkin(8);
t94 = t60 * pkin(2) + t57 * t85;
t55 = qJ(3) + qJ(4);
t46 = pkin(11) + t55;
t44 = qJ(6) + t46;
t39 = sin(t44);
t40 = cos(t44);
t82 = t58 * t60;
t5 = t39 * t82 + t40 * t61;
t6 = t39 * t61 - t40 * t82;
t93 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t81 = t60 * t61;
t7 = -t39 * t81 + t40 * t58;
t8 = t39 * t58 + t40 * t81;
t92 = t7 * rSges(7,1) - t8 * rSges(7,2);
t56 = sin(qJ(3));
t91 = pkin(3) * t56;
t47 = sin(t55);
t90 = pkin(4) * t47;
t87 = g(3) * t57;
t41 = sin(t46);
t42 = cos(t46);
t13 = t41 * t82 + t42 * t61;
t14 = t41 * t61 - t42 * t82;
t84 = -t13 * rSges(6,1) + t14 * rSges(6,2);
t83 = rSges(3,2) * t57;
t15 = -t41 * t81 + t42 * t58;
t16 = t41 * t58 + t42 * t81;
t77 = t15 * rSges(6,1) - t16 * rSges(6,2);
t48 = cos(t55);
t43 = pkin(4) * t48;
t59 = cos(qJ(3));
t51 = t59 * pkin(3);
t37 = t43 + t51;
t76 = t61 * pkin(1) + t58 * pkin(7);
t38 = pkin(5) * t42;
t27 = t38 + t37;
t33 = -pkin(5) * t41 - t90;
t75 = rSges(3,1) * t60 - t83;
t73 = -rSges(5,1) * t47 - rSges(5,2) * t48;
t72 = -rSges(6,1) * t41 - rSges(6,2) * t42;
t71 = -rSges(7,1) * t39 - rSges(7,2) * t40;
t70 = rSges(4,1) * t59 - rSges(4,2) * t56 + pkin(2);
t24 = -t47 * t81 + t48 * t58;
t22 = t47 * t82 + t48 * t61;
t31 = -t56 * t81 + t58 * t59;
t29 = t56 * t82 + t59 * t61;
t45 = t51 + pkin(2);
t69 = rSges(5,1) * t48 - rSges(5,2) * t47 + t45;
t35 = pkin(2) + t37;
t68 = rSges(6,1) * t42 - rSges(6,2) * t41 + t35;
t21 = pkin(2) + t27;
t67 = rSges(7,1) * t40 - rSges(7,2) * t39 + t21;
t66 = t60 * t45 + t80 * t57;
t65 = t35 * t60 + t79 * t57;
t64 = t21 * t60 + t78 * t57;
t23 = t47 * t61 - t48 * t82;
t25 = t47 * t58 + t48 * t81;
t63 = g(1) * (t24 * rSges(5,1) - t25 * rSges(5,2)) + g(2) * (-t22 * rSges(5,1) + t23 * rSges(5,2));
t52 = t61 * pkin(7);
t36 = t90 + t91;
t34 = t38 + t43;
t32 = t56 * t58 + t59 * t81;
t30 = t56 * t61 - t59 * t82;
t26 = -t33 + t91;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t58 - rSges(2,2) * t61) + g(2) * (rSges(2,1) * t61 - rSges(2,2) * t58)) - m(3) * (g(1) * (rSges(3,3) * t61 + t52) + g(2) * (rSges(3,1) * t81 - t61 * t83 + t76) + (g(1) * (-pkin(1) - t75) + g(2) * rSges(3,3)) * t58) - m(4) * ((rSges(4,1) * t32 + rSges(4,2) * t31 + t94 * t61 + t76) * g(2) + (rSges(4,1) * t30 + rSges(4,2) * t29 + t52 + (-pkin(1) - t94) * t58) * g(1)) - m(5) * (g(1) * (t23 * rSges(5,1) + t22 * rSges(5,2) + t52) + g(2) * (t25 * rSges(5,1) + t24 * rSges(5,2) + t76) + (g(1) * t91 + g(2) * t66) * t61 + (g(1) * (-pkin(1) - t66) + g(2) * t91) * t58) - m(6) * (g(1) * (rSges(6,1) * t14 + rSges(6,2) * t13 + t52) + g(2) * (rSges(6,1) * t16 + rSges(6,2) * t15 + t76) + (g(1) * t36 + g(2) * t65) * t61 + (g(1) * (-pkin(1) - t65) + g(2) * t36) * t58) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t52) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t76) + (g(1) * t26 + g(2) * t64) * t61 + (g(1) * (-pkin(1) - t64) + g(2) * t26) * t58) -m(3) * (g(3) * t75 + t95 * (-rSges(3,1) * t57 - rSges(3,2) * t60)) - m(4) * ((g(3) * t70 + t95 * t85) * t60 + (g(3) * t85 - t95 * t70) * t57) - m(5) * ((g(3) * t69 + t95 * t80) * t60 + (g(3) * t80 - t95 * t69) * t57) - m(6) * ((g(3) * t68 + t95 * t79) * t60 + (g(3) * t79 - t95 * t68) * t57) - m(7) * ((g(3) * t67 + t95 * t78) * t60 + (g(3) * t78 - t95 * t67) * t57) -m(4) * (g(1) * (rSges(4,1) * t31 - rSges(4,2) * t32) + g(2) * (-rSges(4,1) * t29 + rSges(4,2) * t30) + (-rSges(4,1) * t56 - rSges(4,2) * t59) * t87) - m(5) * ((g(1) * t31 - g(2) * t29) * pkin(3) + (t73 - t91) * t87 + t63) - m(6) * (g(1) * (-t36 * t81 + t37 * t58 + t77) + g(2) * (-t36 * t82 - t37 * t61 + t84) + (-t36 + t72) * t87) - m(7) * (g(1) * (-t26 * t81 + t27 * t58 + t92) + g(2) * (-t26 * t82 - t27 * t61 + t93) + (-t26 + t71) * t87) -m(5) * t63 - m(6) * (g(1) * (pkin(4) * t24 + t77) + g(2) * (-pkin(4) * t22 + t84)) - m(7) * (g(1) * (t33 * t81 + t34 * t58 + t92) + g(2) * (t33 * t82 - t34 * t61 + t93)) + (-m(5) * t73 - m(6) * (t72 - t90) - m(7) * (t33 + t71)) * t87 (-m(6) - m(7)) * (-g(3) * t60 + t95 * t57) -m(7) * (g(1) * t92 + g(2) * t93 + t71 * t87)];
taug  = t1(:);
