% Calculate Gravitation load on the joints for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:24
% EndTime: 2019-03-08 18:45:25
% DurationCPUTime: 0.76s
% Computational Cost: add. (758->122), mult. (1977->197), div. (0->0), fcn. (2518->16), ass. (0->73)
t65 = sin(pkin(12));
t66 = sin(pkin(11));
t51 = t66 * t65;
t69 = cos(pkin(12));
t70 = cos(pkin(11));
t58 = t70 * t69;
t72 = cos(pkin(6));
t44 = -t58 * t72 + t51;
t67 = sin(pkin(7));
t68 = sin(pkin(6));
t55 = t68 * t67;
t71 = cos(pkin(7));
t99 = t44 * t71 + t70 * t55;
t52 = t66 * t69;
t56 = t70 * t65;
t45 = t52 * t72 + t56;
t54 = t68 * t66;
t98 = t45 * t71 - t67 * t54;
t97 = t69 * t71 * t68 + t72 * t67;
t37 = sin(qJ(3));
t53 = t68 * t65;
t82 = cos(qJ(3));
t18 = t97 * t37 + t82 * t53;
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t43 = -t55 * t69 + t71 * t72;
t12 = t18 * t36 - t38 * t43;
t57 = t70 * t68;
t39 = t44 * t67 - t57 * t71;
t23 = t56 * t72 + t52;
t9 = t23 * t82 - t99 * t37;
t2 = t36 * t9 - t38 * t39;
t24 = -t51 * t72 + t58;
t11 = t24 * t82 - t98 * t37;
t40 = t45 * t67 + t54 * t71;
t4 = t11 * t36 - t38 * t40;
t96 = g(1) * t4 + g(2) * t2 + g(3) * t12;
t13 = t18 * t38 + t36 * t43;
t3 = t36 * t39 + t9 * t38;
t5 = t11 * t38 + t36 * t40;
t95 = g(1) * t5 + g(2) * t3 + g(3) * t13;
t10 = t24 * t37 + t98 * t82;
t17 = t37 * t53 - t97 * t82;
t8 = t23 * t37 + t99 * t82;
t94 = t36 * (-g(1) * t10 - g(2) * t8 - g(3) * t17);
t90 = -m(6) - m(7);
t85 = t38 * pkin(4);
t33 = sin(pkin(13));
t84 = t9 * t33;
t83 = rSges(5,3) + pkin(9);
t81 = t11 * t33;
t80 = t18 * t33;
t32 = pkin(13) + qJ(6);
t30 = sin(t32);
t79 = t30 * t38;
t31 = cos(t32);
t78 = t31 * t38;
t77 = t33 * t38;
t34 = cos(pkin(13));
t76 = t34 * t38;
t29 = pkin(5) * t34 + pkin(4);
t75 = t38 * t29;
t74 = rSges(7,3) + pkin(10) + qJ(5);
t73 = rSges(6,3) + qJ(5);
t6 = t8 * pkin(3);
t64 = t9 * pkin(9) - t6;
t7 = t10 * pkin(3);
t63 = t11 * pkin(9) - t7;
t16 = t17 * pkin(3);
t62 = t18 * pkin(9) - t16;
t61 = -m(3) - m(4) - m(5) + t90;
t60 = -rSges(5,1) * t38 + rSges(5,2) * t36;
t1 = [(-m(2) + t61) * g(3), t61 * (g(1) * t54 - g(2) * t57 + g(3) * t72) -m(4) * (g(1) * (-rSges(4,1) * t10 - rSges(4,2) * t11) + g(2) * (-rSges(4,1) * t8 - rSges(4,2) * t9) + g(3) * (-rSges(4,1) * t17 - rSges(4,2) * t18)) - m(5) * (g(1) * (t10 * t60 + t11 * t83 - t7) + g(2) * (t60 * t8 + t83 * t9 - t6) + g(3) * (t17 * t60 + t18 * t83 - t16)) + (-g(1) * (-t10 * t75 + pkin(5) * t81 + (-t10 * t78 + t11 * t30) * rSges(7,1) + (t10 * t79 + t11 * t31) * rSges(7,2) + t63) - g(2) * (-t8 * t75 + pkin(5) * t84 + (t9 * t30 - t78 * t8) * rSges(7,1) + (t9 * t31 + t79 * t8) * rSges(7,2) + t64) - g(3) * (-t17 * t75 + pkin(5) * t80 + (-t17 * t78 + t18 * t30) * rSges(7,1) + (t17 * t79 + t18 * t31) * rSges(7,2) + t62) - t74 * t94) * m(7) + (-g(1) * (-t10 * t85 + (-t10 * t76 + t81) * rSges(6,1) + (t10 * t77 + t11 * t34) * rSges(6,2) + t63) - g(2) * (-t8 * t85 + (-t76 * t8 + t84) * rSges(6,1) + (t9 * t34 + t77 * t8) * rSges(6,2) + t64) - g(3) * (-t17 * t85 + (-t17 * t76 + t80) * rSges(6,1) + (t17 * t77 + t18 * t34) * rSges(6,2) + t62) - t73 * t94) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(2) * (-rSges(5,1) * t2 - rSges(5,2) * t3) + g(3) * (-rSges(5,1) * t12 - rSges(5,2) * t13)) - m(6) * (t95 * t73 + t96 * (-rSges(6,1) * t34 + rSges(6,2) * t33 - pkin(4))) - m(7) * (t95 * t74 + t96 * (-rSges(7,1) * t31 + rSges(7,2) * t30 - t29)) t90 * t96, -m(7) * (g(1) * ((t10 * t31 - t30 * t5) * rSges(7,1) + (-t10 * t30 - t31 * t5) * rSges(7,2)) + g(2) * ((-t3 * t30 + t31 * t8) * rSges(7,1) + (-t3 * t31 - t30 * t8) * rSges(7,2)) + g(3) * ((-t13 * t30 + t17 * t31) * rSges(7,1) + (-t13 * t31 - t17 * t30) * rSges(7,2)))];
taug  = t1(:);
