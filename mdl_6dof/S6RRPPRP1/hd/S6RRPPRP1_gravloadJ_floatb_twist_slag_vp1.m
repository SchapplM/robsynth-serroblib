% Calculate Gravitation load on the joints for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:31
% EndTime: 2019-03-09 08:25:33
% DurationCPUTime: 0.79s
% Computational Cost: add. (450->118), mult. (466->153), div. (0->0), fcn. (442->10), ass. (0->55)
t23 = qJ(2) + pkin(9);
t18 = sin(t23);
t20 = cos(t23);
t24 = sin(pkin(10));
t25 = cos(pkin(10));
t39 = rSges(5,1) * t25 - rSges(5,2) * t24 + pkin(3);
t51 = rSges(5,3) + qJ(4);
t71 = t51 * t18 + t39 * t20;
t62 = rSges(7,1) + pkin(5);
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t45 = g(1) * t31 + g(2) * t29;
t50 = rSges(7,3) + qJ(6);
t22 = pkin(10) + qJ(5);
t17 = sin(t22);
t19 = cos(t22);
t70 = t17 * t50 + t19 * t62;
t28 = sin(qJ(2));
t69 = pkin(2) * t28;
t68 = pkin(4) * t24;
t30 = cos(qJ(2));
t21 = t30 * pkin(2);
t16 = t21 + pkin(1);
t10 = t31 * t16;
t66 = g(2) * t10;
t26 = -qJ(3) - pkin(7);
t65 = g(2) * t26;
t63 = g(3) * t18;
t61 = rSges(3,3) + pkin(7);
t15 = pkin(4) * t25 + pkin(3);
t7 = t20 * t15;
t60 = t21 + t7;
t59 = t20 * t31;
t58 = t29 * t17;
t57 = t29 * t19;
t56 = t31 * t17;
t55 = t31 * t18;
t27 = -pkin(8) - qJ(4);
t54 = rSges(7,2) - t27;
t53 = rSges(4,3) - t26;
t52 = rSges(6,3) - t27;
t49 = -m(5) - m(6) - m(7);
t48 = -t16 - t7;
t46 = t29 * t18 * t27 + (-t26 + t68) * t31;
t44 = rSges(3,1) * t30 - rSges(3,2) * t28;
t42 = rSges(4,1) * t20 - rSges(4,2) * t18;
t41 = rSges(6,1) * t19 - rSges(6,2) * t17;
t40 = pkin(1) + t44;
t38 = t24 * rSges(5,1) + t25 * rSges(5,2) - t26;
t37 = t15 * t59 - t27 * t55 + t29 * t68 + t10;
t5 = t19 * t59 + t58;
t4 = t20 * t56 - t57;
t3 = t20 * t57 - t56;
t2 = t19 * t31 + t20 * t58;
t1 = [-m(2) * (g(1) * (-t29 * rSges(2,1) - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - t29 * rSges(2,2))) - m(3) * ((g(1) * t61 + g(2) * t40) * t31 + (-g(1) * t40 + g(2) * t61) * t29) - m(4) * (t66 + (g(1) * t53 + g(2) * t42) * t31 + (g(1) * (-t16 - t42) + g(2) * t53) * t29) - m(5) * (t66 + (g(1) * t38 + t71 * g(2)) * t31 + (g(2) * t38 + (-t16 - t71) * g(1)) * t29) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t46) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + rSges(6,3) * t55 + t37) + (g(1) * (-t18 * rSges(6,3) + t48) - t65) * t29) - m(7) * (g(1) * (-t50 * t2 - t62 * t3 + t46) + g(2) * (rSges(7,2) * t55 + t50 * t4 + t62 * t5 + t37) + (g(1) * (-rSges(7,2) * t18 + t48) - t65) * t29) -m(3) * (g(3) * t44 + t45 * (-rSges(3,1) * t28 - rSges(3,2) * t30)) - m(4) * (g(3) * (t21 + t42) + t45 * (-rSges(4,1) * t18 - rSges(4,2) * t20 - t69)) - m(5) * (g(3) * (t21 + t71) + t45 * (-t18 * t39 + t20 * t51 - t69)) - m(6) * (g(3) * (t18 * t52 + t20 * t41 + t60) + t45 * (-t69 + t52 * t20 + (-t15 - t41) * t18)) - m(7) * (g(3) * t60 - t45 * t69 + (g(3) * t70 + t45 * t54) * t20 + (g(3) * t54 + t45 * (-t15 - t70)) * t18) (-m(4) + t49) * (g(1) * t29 - g(2) * t31) t49 * (-g(3) * t20 + t18 * t45) -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3)) - m(7) * (g(1) * (-t4 * t62 + t5 * t50) + g(2) * (-t2 * t62 + t3 * t50)) + (-m(6) * (-rSges(6,1) * t17 - rSges(6,2) * t19) - m(7) * (-t62 * t17 + t50 * t19)) * t63, -m(7) * (g(1) * t4 + g(2) * t2 + t17 * t63)];
taug  = t1(:);
