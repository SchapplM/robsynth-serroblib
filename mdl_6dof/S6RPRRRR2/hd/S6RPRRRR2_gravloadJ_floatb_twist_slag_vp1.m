% Calculate Gravitation load on the joints for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:01
% EndTime: 2019-03-09 06:57:03
% DurationCPUTime: 0.71s
% Computational Cost: add. (531->131), mult. (432->178), div. (0->0), fcn. (389->12), ass. (0->63)
t96 = rSges(6,3) + pkin(9);
t44 = qJ(3) + qJ(4);
t37 = sin(t44);
t39 = cos(t44);
t95 = t39 * rSges(5,1) - t37 * rSges(5,2);
t42 = qJ(1) + pkin(11);
t34 = sin(t42);
t35 = cos(t42);
t94 = g(1) * t35 + g(2) * t34;
t48 = cos(qJ(5));
t32 = t48 * pkin(5) + pkin(4);
t51 = -pkin(10) - pkin(9);
t59 = t39 * t32 + (rSges(7,3) - t51) * t37;
t62 = t39 * pkin(4) + t96 * t37;
t43 = qJ(5) + qJ(6);
t38 = cos(t43);
t36 = sin(t43);
t82 = t36 * t39;
t5 = t34 * t82 + t35 * t38;
t79 = t38 * t39;
t6 = -t34 * t79 + t35 * t36;
t93 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = t34 * t38 - t35 * t82;
t8 = t34 * t36 + t35 * t79;
t92 = t7 * rSges(7,1) - t8 * rSges(7,2);
t46 = sin(qJ(3));
t91 = pkin(3) * t46;
t45 = sin(qJ(5));
t90 = pkin(5) * t45;
t87 = g(3) * t37;
t47 = sin(qJ(1));
t86 = t47 * pkin(1);
t85 = rSges(4,3) + pkin(7);
t84 = t34 * t39;
t83 = t35 * t39;
t78 = t39 * t45;
t77 = t39 * t48;
t52 = -pkin(8) - pkin(7);
t76 = rSges(5,3) - t52;
t49 = cos(qJ(3));
t40 = t49 * pkin(3);
t33 = t40 + pkin(2);
t50 = cos(qJ(1));
t41 = t50 * pkin(1);
t75 = t35 * t33 + t41;
t74 = g(1) * t86;
t73 = rSges(6,2) * t37 * t45;
t72 = rSges(7,2) * t36 * t37;
t71 = -rSges(6,1) * t48 - pkin(4);
t70 = -t52 + t90;
t69 = -rSges(7,1) * t38 - t32;
t67 = t49 * rSges(4,1) - t46 * rSges(4,2);
t64 = -rSges(7,1) * t36 - rSges(7,2) * t38;
t63 = pkin(2) + t67;
t11 = t34 * t48 - t35 * t78;
t9 = t34 * t78 + t35 * t48;
t61 = rSges(6,1) * t77 - rSges(6,2) * t78 + t62;
t57 = g(1) * (rSges(7,3) * t83 + t35 * t72) + g(2) * (rSges(7,3) * t84 + t34 * t72);
t56 = rSges(7,1) * t79 - rSges(7,2) * t82 + t59;
t55 = g(1) * (t35 * t73 + t96 * t83) + g(2) * (t34 * t73 + t96 * t84);
t12 = t34 * t45 + t35 * t77;
t10 = -t34 * t77 + t35 * t45;
t1 = [-m(2) * (g(1) * (-t47 * rSges(2,1) - t50 * rSges(2,2)) + g(2) * (t50 * rSges(2,1) - t47 * rSges(2,2))) - m(3) * (g(1) * (-t34 * rSges(3,1) - t35 * rSges(3,2) - t86) + g(2) * (t35 * rSges(3,1) - t34 * rSges(3,2) + t41)) - m(4) * (-t74 + g(2) * t41 + (g(1) * t85 + g(2) * t63) * t35 + (-g(1) * t63 + g(2) * t85) * t34) - m(5) * (-t74 + g(2) * t75 + (g(1) * t76 + g(2) * t95) * t35 + (g(1) * (-t33 - t95) + g(2) * t76) * t34) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) - t86) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t75) + (-g(1) * t52 + g(2) * t62) * t35 + (g(1) * (-t33 - t62) - g(2) * t52) * t34) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) - t86) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t75) + (g(1) * t70 + g(2) * t59) * t35 + (g(1) * (-t33 - t59) + g(2) * t70) * t34) (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(4) * (g(3) * t67 + t94 * (-rSges(4,1) * t46 - rSges(4,2) * t49)) - m(5) * (g(3) * (t40 + t95) + t94 * (-rSges(5,1) * t37 - rSges(5,2) * t39 - t91)) - m(6) * (g(3) * (t40 + t61) + t55 + t94 * (t71 * t37 - t91)) - m(7) * (g(3) * (t40 + t56) + t57 + t94 * (t69 * t37 - t39 * t51 - t91)) -m(5) * g(3) * t95 - m(6) * (g(3) * t61 + t55) - m(7) * (g(3) * t56 + t57) + t94 * ((m(5) * rSges(5,2) + m(7) * t51) * t39 + (m(5) * rSges(5,1) - m(6) * t71 - m(7) * t69) * t37) -m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2))) - m(7) * (g(1) * (t11 * pkin(5) + t92) + g(2) * (-t9 * pkin(5) + t93)) + (-m(6) * (-rSges(6,1) * t45 - rSges(6,2) * t48) - m(7) * (t64 - t90)) * t87, -m(7) * (g(1) * t92 + g(2) * t93 + t64 * t87)];
taug  = t1(:);
