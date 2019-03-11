% Calculate Gravitation load on the joints for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:13
% EndTime: 2019-03-08 20:45:14
% DurationCPUTime: 0.68s
% Computational Cost: add. (427->113), mult. (929->171), div. (0->0), fcn. (1088->12), ass. (0->56)
t84 = -rSges(5,3) - pkin(8);
t31 = sin(pkin(11));
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t62 = cos(pkin(11));
t63 = cos(pkin(6));
t53 = t63 * t62;
t16 = t31 * t38 + t35 * t53;
t60 = t31 * t63;
t18 = -t35 * t60 + t38 * t62;
t83 = g(1) * t18 + g(2) * t16;
t15 = t31 * t35 - t38 * t53;
t17 = t35 * t62 + t38 * t60;
t82 = g(1) * t17 + g(2) * t15;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t32 = sin(pkin(6));
t59 = t32 * t62;
t11 = -t15 * t34 + t37 * t59;
t67 = t32 * t38;
t20 = -t34 * t67 + t37 * t63;
t69 = t31 * t32;
t9 = t17 * t34 + t37 * t69;
t81 = g(1) * t9 - g(2) * t11 + g(3) * t20;
t10 = t15 * t37 + t34 * t59;
t19 = -t34 * t63 - t37 * t67;
t8 = t17 * t37 - t34 * t69;
t80 = g(1) * t8 + g(2) * t10 + g(3) * t19;
t71 = g(3) * t32;
t70 = rSges(6,3) + pkin(9);
t68 = t32 * t35;
t66 = rSges(7,3) + pkin(10) + pkin(9);
t65 = pkin(2) * t67 + qJ(3) * t68;
t64 = rSges(4,3) + qJ(3);
t61 = -m(4) - m(5) - m(6) - m(7);
t58 = g(3) * (pkin(8) * t67 + t65);
t33 = sin(qJ(5));
t36 = cos(qJ(5));
t57 = t18 * t36 - t33 * t9;
t56 = rSges(5,1) * t34 + rSges(5,2) * t37;
t55 = t33 * rSges(6,1) + t36 * rSges(6,2);
t54 = t11 * t33 + t16 * t36;
t52 = rSges(6,1) * t36 - rSges(6,2) * t33 + pkin(4);
t50 = -t20 * t33 + t36 * t68;
t30 = qJ(5) + qJ(6);
t28 = sin(t30);
t29 = cos(t30);
t49 = rSges(7,1) * t29 - rSges(7,2) * t28 + pkin(5) * t36 + pkin(4);
t47 = t28 * rSges(7,1) + t29 * rSges(7,2) + t33 * pkin(5);
t13 = t15 * pkin(2);
t14 = t17 * pkin(2);
t45 = -g(1) * t14 - g(2) * t13 + t58;
t44 = t34 * t52 - t37 * t70;
t43 = m(7) * (g(1) * ((t18 * t29 - t28 * t9) * rSges(7,1) + (-t18 * t28 - t29 * t9) * rSges(7,2)) + g(2) * ((t11 * t28 + t16 * t29) * rSges(7,1) + (t11 * t29 - t16 * t28) * rSges(7,2)) + g(3) * ((-t20 * t28 + t29 * t68) * rSges(7,1) + (-t20 * t29 - t28 * t68) * rSges(7,2)));
t42 = t34 * t49 - t37 * t66;
t1 = [(-m(2) - m(3) + t61) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t17 - rSges(3,2) * t18) + g(2) * (-rSges(3,1) * t15 - rSges(3,2) * t16) + (rSges(3,1) * t38 - rSges(3,2) * t35) * t71) - m(4) * (g(1) * (rSges(4,2) * t17 + t18 * t64 - t14) + g(2) * (rSges(4,2) * t15 + t16 * t64 - t13) + g(3) * ((-rSges(4,2) * t38 + rSges(4,3) * t35) * t32 + t65)) - m(5) * (g(1) * (t84 * t17 - t14) + g(2) * (t84 * t15 - t13) + t58 + (rSges(5,3) * t38 + t35 * t56) * t71 + t83 * (qJ(3) + t56)) - m(6) * ((t35 * t44 + t38 * t55) * t71 + t45 + t82 * (-pkin(8) - t55) + t83 * (qJ(3) + t44)) - m(7) * ((t35 * t42 + t38 * t47) * t71 + t45 + t82 * (-pkin(8) - t47) + t83 * (qJ(3) + t42)) t61 * (-g(3) * t67 + t82) -m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (rSges(5,1) * t10 + rSges(5,2) * t11) + g(3) * (rSges(5,1) * t19 - rSges(5,2) * t20)) - m(6) * (t80 * t52 + t81 * t70) - m(7) * (t80 * t49 + t81 * t66) -m(6) * (g(1) * (t57 * rSges(6,1) + (-t18 * t33 - t36 * t9) * rSges(6,2)) + g(2) * (t54 * rSges(6,1) + (t11 * t36 - t16 * t33) * rSges(6,2)) + g(3) * (t50 * rSges(6,1) + (-t20 * t36 - t33 * t68) * rSges(6,2))) - t43 - m(7) * (g(1) * t57 + g(2) * t54 + g(3) * t50) * pkin(5), -t43];
taug  = t1(:);
