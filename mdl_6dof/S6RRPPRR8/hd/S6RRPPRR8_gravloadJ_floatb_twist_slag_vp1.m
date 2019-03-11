% Calculate Gravitation load on the joints for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:53
% EndTime: 2019-03-09 09:22:55
% DurationCPUTime: 0.99s
% Computational Cost: add. (399->172), mult. (807->229), div. (0->0), fcn. (864->10), ass. (0->66)
t39 = sin(qJ(2));
t84 = g(3) * t39;
t43 = cos(qJ(1));
t40 = sin(qJ(1));
t85 = g(2) * t40;
t90 = g(1) * t43 + t85;
t36 = sin(pkin(10));
t37 = cos(pkin(10));
t42 = cos(qJ(2));
t77 = t40 * t42;
t15 = t36 * t77 + t37 * t43;
t75 = t43 * t36;
t16 = t37 * t77 - t75;
t35 = qJ(5) + qJ(6);
t28 = sin(t35);
t29 = cos(t35);
t56 = -t15 * t28 - t16 * t29;
t57 = t15 * t29 - t16 * t28;
t89 = t57 * rSges(7,1) + t56 * rSges(7,2);
t17 = -t40 * t37 + t42 * t75;
t76 = t42 * t43;
t18 = t36 * t40 + t37 * t76;
t5 = t17 * t29 - t18 * t28;
t6 = t17 * t28 + t18 * t29;
t88 = t5 * rSges(7,1) - t6 * rSges(7,2);
t87 = g(1) * t40;
t32 = t42 * pkin(2);
t83 = -pkin(8) - rSges(6,3);
t38 = sin(qJ(5));
t82 = t15 * t38;
t81 = t36 * t38;
t80 = t36 * t42;
t79 = t37 * t42;
t78 = t39 * t43;
t30 = t39 * qJ(3);
t74 = t30 + t32;
t73 = t43 * pkin(1) + t40 * pkin(7);
t72 = -pkin(9) - pkin(8) - rSges(7,3);
t71 = qJ(3) * t42;
t70 = rSges(5,3) + qJ(4);
t69 = -m(5) - m(6) - m(7);
t68 = -pkin(1) - t32;
t52 = t28 * t36 + t29 * t37;
t53 = -t28 * t37 + t29 * t36;
t67 = (rSges(7,1) * t53 - rSges(7,2) * t52) * t84;
t66 = t83 * t43;
t65 = t72 * t43;
t64 = t38 * pkin(5) + qJ(4);
t63 = pkin(3) * t79 + qJ(4) * t80 + t74;
t62 = pkin(2) * t76 + t43 * t30 + t73;
t33 = t43 * pkin(7);
t61 = -t16 * pkin(3) - t15 * qJ(4) + t33;
t60 = t18 * pkin(3) + t62;
t59 = rSges(3,1) * t42 - rSges(3,2) * t39;
t41 = cos(qJ(5));
t55 = t15 * t41 - t16 * t38;
t54 = -t16 * t41 - t82;
t8 = t17 * t41 - t18 * t38;
t51 = t36 * t41 - t37 * t38;
t50 = t37 * t41 + t81;
t20 = t40 * t71;
t23 = t43 * t71;
t47 = g(1) * t23 + g(2) * t20 + g(3) * t63;
t27 = pkin(5) * t41 + pkin(4);
t9 = t17 * t38 + t18 * t41;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t40 - rSges(2,2) * t43) + g(2) * (rSges(2,1) * t43 - rSges(2,2) * t40)) - m(3) * (g(1) * (rSges(3,3) * t43 + t33) + g(2) * (rSges(3,1) * t76 - rSges(3,2) * t78 + t73) + (g(1) * (-pkin(1) - t59) + g(2) * rSges(3,3)) * t40) - m(4) * (g(1) * (-rSges(4,1) * t16 + rSges(4,2) * t15 + t33) + g(2) * (rSges(4,1) * t18 - rSges(4,2) * t17 + rSges(4,3) * t78 + t62) + ((-rSges(4,3) - qJ(3)) * t39 + t68) * t87) - m(5) * (g(1) * (-rSges(5,1) * t16 - rSges(5,3) * t15 + t61) + g(2) * (rSges(5,1) * t18 + rSges(5,2) * t78 + t70 * t17 + t60) + ((-rSges(5,2) - qJ(3)) * t39 + t68) * t87) - m(6) * (g(1) * (t54 * rSges(6,1) - t55 * rSges(6,2) - t16 * pkin(4) + t61) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + pkin(4) * t18 + t17 * qJ(4) + t39 * t66 + t60) + ((-qJ(3) - t83) * t39 + t68) * t87) - m(7) * (g(1) * (t56 * rSges(7,1) - t57 * rSges(7,2) - pkin(5) * t82 - t16 * t27 + t61) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t64 * t17 + t18 * t27 + t39 * t65 + t60) + ((-qJ(3) - t72) * t39 + t68) * t87) -m(3) * (g(3) * t59 + t90 * (-rSges(3,1) * t39 - rSges(3,2) * t42)) - m(4) * (g(1) * (rSges(4,3) * t76 + t23) + g(2) * (rSges(4,3) * t77 + t20) + g(3) * (rSges(4,1) * t79 - rSges(4,2) * t80 + t74) + (g(3) * rSges(4,3) + t90 * (-rSges(4,1) * t37 + rSges(4,2) * t36 - pkin(2))) * t39) - m(5) * (g(1) * (rSges(5,2) * t76 + t23) + g(2) * (rSges(5,2) * t77 + t20) + g(3) * (rSges(5,1) * t79 + rSges(5,3) * t80 + t63) + (g(3) * rSges(5,2) + t90 * (-pkin(2) + (-rSges(5,1) - pkin(3)) * t37 - t70 * t36)) * t39) - m(6) * ((g(3) * (t50 * rSges(6,1) + t51 * rSges(6,2) + t37 * pkin(4)) + g(1) * t66 + t83 * t85) * t42 + (g(3) * t83 + t90 * (-pkin(2) + (-t38 * rSges(6,1) - t41 * rSges(6,2) - qJ(4)) * t36 + (-t41 * rSges(6,1) + t38 * rSges(6,2) - pkin(3) - pkin(4)) * t37)) * t39 + t47) - m(7) * ((g(3) * (t52 * rSges(7,1) + t53 * rSges(7,2) + pkin(5) * t81 + t37 * t27) + g(1) * t65 + t72 * t85) * t42 + (g(3) * t72 + t90 * (-pkin(2) + (-t29 * rSges(7,1) + t28 * rSges(7,2) - pkin(3) - t27) * t37 + (-t28 * rSges(7,1) - t29 * rSges(7,2) - t64) * t36)) * t39 + t47) (-m(4) + t69) * (-g(3) * t42 + t90 * t39) t69 * (g(1) * t17 + g(2) * t15 + t36 * t84) -m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (t55 * rSges(6,1) + t54 * rSges(6,2))) - m(7) * (g(1) * (t8 * pkin(5) + t88) + g(2) * (t55 * pkin(5) + t89) + t67) + (-m(6) * (t51 * rSges(6,1) - t50 * rSges(6,2)) - m(7) * t51 * pkin(5)) * t84, -m(7) * (g(1) * t88 + g(2) * t89 + t67)];
taug  = t1(:);
