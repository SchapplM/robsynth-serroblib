% Calculate Gravitation load on the joints for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:44
% EndTime: 2019-03-08 20:04:45
% DurationCPUTime: 0.70s
% Computational Cost: add. (519->123), mult. (877->191), div. (0->0), fcn. (1020->12), ass. (0->67)
t45 = sin(qJ(2));
t47 = cos(qJ(2));
t59 = sin(pkin(10));
t61 = cos(pkin(6));
t52 = t61 * t59;
t60 = cos(pkin(10));
t26 = t45 * t60 + t47 * t52;
t53 = t61 * t60;
t24 = t45 * t59 - t47 * t53;
t77 = g(2) * t24;
t88 = g(1) * t26 + t77;
t40 = sin(pkin(6));
t68 = t40 * t47;
t87 = g(3) * t68 - t88;
t46 = cos(qJ(5));
t35 = pkin(5) * t46 + pkin(4);
t38 = pkin(11) + qJ(4);
t36 = sin(t38);
t37 = cos(t38);
t65 = rSges(7,3) + qJ(6) + pkin(9);
t86 = t35 * t37 + t36 * t65;
t25 = t45 * t53 + t47 * t59;
t56 = t40 * t60;
t10 = t25 * t36 + t37 * t56;
t27 = -t45 * t52 + t47 * t60;
t55 = t40 * t59;
t12 = t27 * t36 - t37 * t55;
t69 = t40 * t45;
t20 = t36 * t69 - t37 * t61;
t85 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t11 = t25 * t37 - t36 * t56;
t13 = t27 * t37 + t36 * t55;
t21 = t36 * t61 + t37 * t69;
t84 = g(1) * t13 + g(2) * t11 + g(3) * t21;
t83 = pkin(4) * t37;
t76 = g(3) * t40;
t75 = rSges(6,3) + pkin(9);
t44 = sin(qJ(5));
t74 = t25 * t44;
t73 = t27 * t44;
t71 = t37 * t44;
t70 = t37 * t46;
t67 = t44 * t47;
t66 = t46 * t47;
t41 = cos(pkin(11));
t34 = pkin(3) * t41 + pkin(2);
t43 = -pkin(8) - qJ(3);
t64 = -t24 * t34 - t25 * t43;
t63 = -t26 * t34 - t27 * t43;
t62 = rSges(4,3) + qJ(3);
t57 = -m(4) - m(5) - m(6) - m(7);
t54 = rSges(5,1) * t37 - rSges(5,2) * t36;
t1 = -t11 * t44 + t24 * t46;
t3 = -t13 * t44 + t26 * t46;
t51 = rSges(4,1) * t41 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t14 = -t21 * t44 - t40 * t66;
t28 = t34 * t68;
t17 = (t37 * t66 + t44 * t45) * t40;
t16 = (-t37 * t67 + t45 * t46) * t40;
t15 = -t21 * t46 + t40 * t67;
t8 = -t26 * t70 + t73;
t7 = t26 * t71 + t27 * t46;
t6 = -t24 * t70 + t74;
t5 = t24 * t71 + t25 * t46;
t4 = -t13 * t46 - t26 * t44;
t2 = -t11 * t46 - t24 * t44;
t9 = [(-m(2) - m(3) + t57) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t26 - rSges(3,2) * t27) + g(2) * (-rSges(3,1) * t24 - rSges(3,2) * t25) + (rSges(3,1) * t47 - rSges(3,2) * t45) * t76) - m(4) * (g(1) * (-t26 * t51 + t27 * t62) + g(2) * t62 * t25 - t51 * t77 + (t45 * t62 + t47 * t51) * t76) - m(5) * (g(1) * (rSges(5,3) * t27 - t26 * t54 + t63) + g(2) * (rSges(5,3) * t25 - t24 * t54 + t64) + g(3) * t28 + (t54 * t47 + (rSges(5,3) - t43) * t45) * t76) - m(6) * (g(1) * (rSges(6,1) * t8 + rSges(6,2) * t7 - t26 * t83 + t63) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 - t24 * t83 + t64) + g(3) * (t17 * rSges(6,1) + t16 * rSges(6,2) - t43 * t69 + t68 * t83 + t28) + t87 * t36 * t75) - m(7) * (g(1) * (rSges(7,1) * t8 + rSges(7,2) * t7 + pkin(5) * t73 + t63) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t74 + t64) + g(3) * (t17 * rSges(7,1) + t16 * rSges(7,2) + t28) + ((pkin(5) * t44 - t43) * t45 + t86 * t47) * t76 - t88 * t86) -t57 * t87, -m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (-rSges(5,1) * t10 - rSges(5,2) * t11) + g(3) * (-rSges(5,1) * t20 - rSges(5,2) * t21)) - m(6) * (t84 * t75 + t85 * (-rSges(6,1) * t46 + rSges(6,2) * t44 - pkin(4))) - m(7) * (t84 * t65 + t85 * (-rSges(7,1) * t46 + rSges(7,2) * t44 - t35)) -m(6) * (g(1) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (rSges(6,1) * t14 + rSges(6,2) * t15)) + (-g(1) * (rSges(7,1) * t3 + rSges(7,2) * t4) - g(2) * (rSges(7,1) * t1 + rSges(7,2) * t2) - g(3) * (t14 * rSges(7,1) + t15 * rSges(7,2)) - (g(1) * t3 + g(2) * t1 + g(3) * t14) * pkin(5)) * m(7), -m(7) * t85];
taug  = t9(:);
