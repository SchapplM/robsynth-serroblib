% Calculate Gravitation load on the joints for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:16
% EndTime: 2019-03-08 19:03:18
% DurationCPUTime: 0.86s
% Computational Cost: add. (898->136), mult. (2336->223), div. (0->0), fcn. (2987->16), ass. (0->76)
t71 = sin(pkin(13));
t72 = sin(pkin(12));
t54 = t72 * t71;
t75 = cos(pkin(13));
t76 = cos(pkin(12));
t61 = t76 * t75;
t78 = cos(pkin(6));
t46 = -t78 * t61 + t54;
t73 = sin(pkin(7));
t74 = sin(pkin(6));
t58 = t74 * t73;
t77 = cos(pkin(7));
t104 = t46 * t77 + t76 * t58;
t55 = t72 * t75;
t59 = t76 * t71;
t47 = t78 * t55 + t59;
t57 = t74 * t72;
t103 = t47 * t77 - t73 * t57;
t102 = t75 * t77 * t74 + t78 * t73;
t40 = sin(qJ(3));
t56 = t74 * t71;
t88 = cos(qJ(3));
t21 = t102 * t40 + t88 * t56;
t27 = -t75 * t58 + t78 * t77;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t17 = -t21 * t39 + t27 * t42;
t28 = t78 * t59 + t55;
t14 = -t104 * t40 + t28 * t88;
t60 = t76 * t74;
t22 = t46 * t73 - t77 * t60;
t7 = -t14 * t39 + t22 * t42;
t29 = -t78 * t54 + t61;
t16 = -t103 * t40 + t29 * t88;
t23 = t47 * t73 + t77 * t57;
t9 = -t16 * t39 + t23 * t42;
t101 = g(1) * t9 + g(2) * t7 + g(3) * t17;
t10 = t16 * t42 + t23 * t39;
t18 = t21 * t42 + t27 * t39;
t8 = t14 * t42 + t22 * t39;
t100 = g(1) * t10 + g(2) * t8 + g(3) * t18;
t13 = t104 * t88 + t28 * t40;
t15 = t103 * t88 + t29 * t40;
t20 = -t102 * t88 + t40 * t56;
t99 = (-g(1) * t15 - g(2) * t13 - g(3) * t20) * t39;
t91 = t42 * pkin(4);
t90 = rSges(5,3) + pkin(9);
t89 = rSges(6,3) + pkin(10);
t38 = sin(qJ(5));
t87 = t14 * t38;
t86 = t16 * t38;
t85 = t21 * t38;
t37 = qJ(5) + qJ(6);
t35 = sin(t37);
t84 = t35 * t42;
t36 = cos(t37);
t83 = t36 * t42;
t82 = t38 * t42;
t41 = cos(qJ(5));
t81 = t41 * t42;
t34 = pkin(5) * t41 + pkin(4);
t80 = t42 * t34;
t79 = rSges(7,3) + pkin(11) + pkin(10);
t11 = t13 * pkin(3);
t70 = t14 * pkin(9) - t11;
t12 = t15 * pkin(3);
t69 = t16 * pkin(9) - t12;
t19 = t20 * pkin(3);
t68 = t21 * pkin(9) - t19;
t67 = -m(3) - m(4) - m(5) - m(6) - m(7);
t66 = t13 * t41 - t38 * t8;
t65 = -rSges(5,1) * t42 + rSges(5,2) * t39;
t64 = -t10 * t38 + t15 * t41;
t63 = -t18 * t38 + t20 * t41;
t48 = m(7) * (g(1) * ((-t10 * t35 + t15 * t36) * rSges(7,1) + (-t10 * t36 - t15 * t35) * rSges(7,2)) + g(2) * ((t13 * t36 - t35 * t8) * rSges(7,1) + (-t13 * t35 - t36 * t8) * rSges(7,2)) + g(3) * ((-t18 * t35 + t20 * t36) * rSges(7,1) + (-t18 * t36 - t20 * t35) * rSges(7,2)));
t1 = [(-m(2) + t67) * g(3), t67 * (g(1) * t57 - g(2) * t60 + g(3) * t78) -m(4) * (g(1) * (-rSges(4,1) * t15 - rSges(4,2) * t16) + g(2) * (-rSges(4,1) * t13 - rSges(4,2) * t14) + g(3) * (-rSges(4,1) * t20 - rSges(4,2) * t21)) - m(5) * (g(1) * (t65 * t15 + t90 * t16 - t12) + g(2) * (t65 * t13 + t90 * t14 - t11) + g(3) * (t65 * t20 + t90 * t21 - t19)) + (-g(1) * (-t15 * t80 + pkin(5) * t86 + (-t15 * t83 + t16 * t35) * rSges(7,1) + (t15 * t84 + t16 * t36) * rSges(7,2) + t69) - g(2) * (-t13 * t80 + pkin(5) * t87 + (-t13 * t83 + t14 * t35) * rSges(7,1) + (t13 * t84 + t14 * t36) * rSges(7,2) + t70) - g(3) * (-t20 * t80 + pkin(5) * t85 + (-t20 * t83 + t21 * t35) * rSges(7,1) + (t20 * t84 + t21 * t36) * rSges(7,2) + t68) - t79 * t99) * m(7) + (-g(1) * (-t15 * t91 + (-t15 * t81 + t86) * rSges(6,1) + (t15 * t82 + t16 * t41) * rSges(6,2) + t69) - g(2) * (-t13 * t91 + (-t13 * t81 + t87) * rSges(6,1) + (t13 * t82 + t14 * t41) * rSges(6,2) + t70) - g(3) * (-t20 * t91 + (-t20 * t81 + t85) * rSges(6,1) + (t20 * t82 + t21 * t41) * rSges(6,2) + t68) - t89 * t99) * m(6), -m(5) * (g(1) * (rSges(5,1) * t9 - rSges(5,2) * t10) + g(2) * (rSges(5,1) * t7 - rSges(5,2) * t8) + g(3) * (rSges(5,1) * t17 - rSges(5,2) * t18)) - m(6) * (t100 * t89 + t101 * (rSges(6,1) * t41 - rSges(6,2) * t38 + pkin(4))) - m(7) * (t100 * t79 + t101 * (rSges(7,1) * t36 - rSges(7,2) * t35 + t34)) -m(6) * (g(1) * (t64 * rSges(6,1) + (-t10 * t41 - t15 * t38) * rSges(6,2)) + g(2) * (t66 * rSges(6,1) + (-t13 * t38 - t41 * t8) * rSges(6,2)) + g(3) * (t63 * rSges(6,1) + (-t18 * t41 - t20 * t38) * rSges(6,2))) - t48 - m(7) * (g(1) * t64 + g(2) * t66 + g(3) * t63) * pkin(5), -t48];
taug  = t1(:);
