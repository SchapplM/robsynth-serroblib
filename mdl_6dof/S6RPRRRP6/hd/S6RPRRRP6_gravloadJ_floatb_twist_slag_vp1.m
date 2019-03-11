% Calculate Gravitation load on the joints for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:15
% EndTime: 2019-03-09 06:14:17
% DurationCPUTime: 0.65s
% Computational Cost: add. (471->134), mult. (480->183), div. (0->0), fcn. (452->10), ass. (0->63)
t67 = rSges(5,3) + pkin(8);
t40 = -pkin(9) - pkin(8);
t56 = rSges(6,3) - t40;
t55 = rSges(7,3) + qJ(6) - t40;
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t77 = g(1) * t39 + g(2) * t37;
t31 = pkin(10) + qJ(3);
t26 = cos(t31);
t32 = qJ(4) + qJ(5);
t28 = cos(t32);
t63 = t28 * t37;
t27 = sin(t32);
t64 = t27 * t39;
t10 = -t26 * t63 + t64;
t62 = t28 * t39;
t65 = t27 * t37;
t9 = t26 * t65 + t62;
t76 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t75 = -t9 * rSges(7,1) + t10 * rSges(7,2);
t11 = -t26 * t64 + t63;
t12 = t26 * t62 + t65;
t74 = t11 * rSges(7,1) - t12 * rSges(7,2);
t73 = t11 * rSges(6,1) - t12 * rSges(6,2);
t36 = sin(qJ(4));
t72 = pkin(4) * t36;
t71 = pkin(5) * t27;
t25 = sin(t31);
t68 = g(3) * t25;
t18 = t71 + t72;
t66 = t18 * t26;
t61 = t36 * t37;
t60 = t36 * t39;
t38 = cos(qJ(4));
t59 = t37 * t38;
t58 = t38 * t39;
t35 = -pkin(7) - qJ(2);
t57 = rSges(4,3) - t35;
t54 = t18 - t35;
t29 = t38 * pkin(4);
t19 = pkin(5) * t28 + t29;
t53 = rSges(3,3) + qJ(2);
t52 = -t35 + t72;
t51 = rSges(4,1) * t26 - rSges(4,2) * t25;
t49 = -rSges(6,1) * t27 - rSges(6,2) * t28;
t48 = -rSges(7,1) * t27 - rSges(7,2) * t28;
t34 = cos(pkin(10));
t47 = rSges(3,1) * t34 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t46 = rSges(5,1) * t38 - rSges(5,2) * t36 + pkin(3);
t15 = -t26 * t60 + t59;
t13 = t26 * t61 + t58;
t24 = t29 + pkin(3);
t45 = rSges(6,1) * t28 - rSges(6,2) * t27 + t24;
t17 = pkin(3) + t19;
t44 = rSges(7,1) * t28 - rSges(7,2) * t27 + t17;
t43 = pkin(3) * t26 + t67 * t25;
t42 = t26 * t24 + t56 * t25;
t41 = t17 * t26 + t55 * t25;
t22 = pkin(2) * t34 + pkin(1);
t20 = t39 * t22;
t16 = t26 * t58 + t61;
t14 = -t26 * t59 + t60;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t37 - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 - rSges(2,2) * t37)) - m(3) * ((g(1) * t53 + g(2) * t47) * t39 + (-g(1) * t47 + g(2) * t53) * t37) - m(4) * (g(2) * t20 + (g(1) * t57 + g(2) * t51) * t39 + (g(1) * (-t22 - t51) + g(2) * t57) * t37) - m(5) * (g(1) * (rSges(5,1) * t14 + rSges(5,2) * t13) + g(2) * (rSges(5,1) * t16 + rSges(5,2) * t15 + t20) + (-g(1) * t35 + g(2) * t43) * t39 + (g(1) * (-t22 - t43) - g(2) * t35) * t37) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2)) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t20) + (g(1) * t52 + g(2) * t42) * t39 + (g(1) * (-t22 - t42) + g(2) * t52) * t37) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9) + g(2) * (rSges(7,1) * t12 + rSges(7,2) * t11 + t20) + (g(1) * t54 + g(2) * t41) * t39 + (g(1) * (-t22 - t41) + g(2) * t54) * t37) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t37 - g(2) * t39) -m(4) * (g(3) * t51 + t77 * (-rSges(4,1) * t25 - rSges(4,2) * t26)) - m(5) * ((g(3) * t46 + t67 * t77) * t26 + (g(3) * t67 - t77 * t46) * t25) - m(6) * ((g(3) * t45 + t56 * t77) * t26 + (g(3) * t56 - t45 * t77) * t25) - m(7) * ((g(3) * t44 + t55 * t77) * t26 + (g(3) * t55 - t44 * t77) * t25) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (-rSges(5,1) * t13 + rSges(5,2) * t14)) - m(6) * (g(1) * (t15 * pkin(4) + t73) + g(2) * (-t13 * pkin(4) + t76)) - m(7) * (g(1) * (t19 * t37 - t39 * t66 + t74) + g(2) * (-t19 * t39 - t37 * t66 + t75)) + (-m(5) * (-rSges(5,1) * t36 - rSges(5,2) * t38) - m(6) * (t49 - t72) - m(7) * (-t18 + t48)) * t68, -m(6) * (g(1) * t73 + g(2) * t76) - m(7) * (g(1) * (t11 * pkin(5) + t74) + g(2) * (-t9 * pkin(5) + t75)) + (-m(6) * t49 - m(7) * (t48 - t71)) * t68, -m(7) * (-g(3) * t26 + t25 * t77)];
taug  = t1(:);
