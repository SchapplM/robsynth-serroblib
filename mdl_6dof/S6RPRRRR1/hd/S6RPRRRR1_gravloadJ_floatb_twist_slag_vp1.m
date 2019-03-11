% Calculate Gravitation load on the joints for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:10
% EndTime: 2019-03-09 06:54:12
% DurationCPUTime: 0.63s
% Computational Cost: add. (520->113), mult. (369->138), div. (0->0), fcn. (314->12), ass. (0->66)
t92 = rSges(7,3) + pkin(10);
t35 = qJ(3) + qJ(4);
t30 = qJ(5) + t35;
t23 = sin(t30);
t24 = cos(t30);
t36 = sin(qJ(6));
t78 = rSges(7,2) * t36;
t91 = t23 * t78 + t24 * t92;
t89 = t24 * rSges(6,1) - t23 * rSges(6,2);
t28 = sin(t35);
t29 = cos(t35);
t62 = t29 * rSges(5,1) - t28 * rSges(5,2);
t33 = qJ(1) + pkin(11);
t26 = sin(t33);
t27 = cos(t33);
t88 = g(1) * t27 + g(2) * t26;
t49 = t24 * pkin(5) + t92 * t23;
t39 = cos(qJ(6));
t79 = rSges(7,1) * t39;
t87 = (-pkin(5) - t79) * t23;
t42 = -pkin(8) - pkin(7);
t86 = pkin(4) * t28;
t37 = sin(qJ(3));
t83 = t37 * pkin(3);
t38 = sin(qJ(1));
t82 = t38 * pkin(1);
t81 = rSges(4,3) + pkin(7);
t40 = cos(qJ(3));
t31 = t40 * pkin(3);
t25 = t31 + pkin(2);
t22 = pkin(4) * t29;
t14 = t22 + t25;
t41 = cos(qJ(1));
t32 = t41 * pkin(1);
t80 = t27 * t14 + t32;
t74 = t26 * t36;
t73 = t26 * t39;
t72 = t27 * t36;
t71 = t27 * t39;
t69 = rSges(5,3) - t42;
t34 = -pkin(9) + t42;
t68 = rSges(6,3) - t34;
t67 = g(1) * t82;
t66 = t91 * t26;
t65 = t91 * t27;
t60 = t22 + t89;
t59 = t40 * rSges(4,1) - t37 * rSges(4,2);
t57 = -rSges(5,1) * t28 - rSges(5,2) * t29;
t55 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t54 = g(2) * t32 - t67;
t53 = pkin(2) + t59;
t52 = t25 + t62;
t51 = t55 * t26;
t50 = t55 * t27;
t48 = t49 + (-t78 + t79) * t24;
t46 = t22 + t48;
t45 = g(1) * t65 + g(2) * t66;
t44 = t88 * t87;
t15 = -t83 - t86;
t7 = t27 * t15;
t6 = t26 * t15;
t4 = t24 * t71 + t74;
t3 = -t24 * t72 + t73;
t2 = -t24 * t73 + t72;
t1 = t24 * t74 + t71;
t5 = [-m(2) * (g(1) * (-t38 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t38 * rSges(2,2))) - m(3) * (g(1) * (-t26 * rSges(3,1) - t27 * rSges(3,2) - t82) + g(2) * (t27 * rSges(3,1) - t26 * rSges(3,2) + t32)) - m(4) * ((g(1) * t81 + g(2) * t53) * t27 + (-g(1) * t53 + g(2) * t81) * t26 + t54) - m(5) * ((g(1) * t69 + g(2) * t52) * t27 + (-g(1) * t52 + g(2) * t69) * t26 + t54) - m(6) * (-t67 + g(2) * t80 + (g(1) * t68 + g(2) * t89) * t27 + (g(1) * (-t14 - t89) + g(2) * t68) * t26) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t82) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t80) + (-g(1) * t34 + g(2) * t49) * t27 + (g(1) * (-t14 - t49) - g(2) * t34) * t26) (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(4) * (g(3) * t59 + t88 * (-rSges(4,1) * t37 - rSges(4,2) * t40)) - m(5) * (g(3) * (t31 + t62) + t88 * (t57 - t83)) - m(6) * (g(1) * (t7 + t50) + g(2) * (t6 + t51) + g(3) * (t31 + t60)) - m(7) * (g(1) * (t7 + t65) + g(2) * (t6 + t66) + g(3) * (t31 + t46) + t44) -m(7) * t45 + (-m(5) * t62 - m(6) * t60 - m(7) * t46) * g(3) + t88 * (-m(5) * t57 - m(6) * (t55 - t86) - m(7) * (-t86 + t87)) -m(6) * (g(1) * t50 + g(2) * t51 + g(3) * t89) - m(7) * (g(3) * t48 + t44 + t45) -m(7) * (g(1) * (t3 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (-t1 * rSges(7,1) + t2 * rSges(7,2)) + g(3) * (-rSges(7,1) * t36 - rSges(7,2) * t39) * t23)];
taug  = t5(:);
