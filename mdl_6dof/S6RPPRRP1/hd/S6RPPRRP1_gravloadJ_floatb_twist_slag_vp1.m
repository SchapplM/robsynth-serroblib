% Calculate Gravitation load on the joints for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:37
% EndTime: 2019-03-09 01:57:37
% DurationCPUTime: 0.42s
% Computational Cost: add. (364->94), mult. (310->126), div. (0->0), fcn. (278->10), ass. (0->42)
t51 = rSges(7,1) + pkin(5);
t50 = rSges(6,3) + pkin(8);
t49 = rSges(7,3) + qJ(6) + pkin(8);
t15 = qJ(1) + pkin(9);
t10 = sin(t15);
t12 = cos(t15);
t48 = g(1) * t12 + g(2) * t10;
t20 = sin(qJ(5));
t22 = cos(qJ(5));
t8 = t22 * pkin(5) + pkin(4);
t47 = m(6) * (rSges(6,1) * t22 - rSges(6,2) * t20 + pkin(4)) + m(7) * (rSges(7,1) * t22 - rSges(7,2) * t20 + t8) + m(5) * rSges(5,1);
t45 = pkin(5) * t20;
t21 = sin(qJ(1));
t42 = t21 * pkin(1);
t23 = cos(qJ(1));
t13 = t23 * pkin(1);
t17 = cos(pkin(10));
t7 = t17 * pkin(3) + pkin(2);
t41 = t12 * t7 + t13;
t40 = t10 * t20;
t39 = t10 * t22;
t38 = t12 * t20;
t37 = t12 * t22;
t19 = -pkin(7) - qJ(3);
t36 = rSges(5,3) - t19;
t35 = rSges(4,3) + qJ(3);
t34 = g(1) * t42;
t33 = -m(4) - m(5) - m(6) - m(7);
t32 = -t19 + t45;
t14 = pkin(10) + qJ(4);
t11 = cos(t14);
t9 = sin(t14);
t31 = t11 * rSges(5,1) - t9 * rSges(5,2);
t30 = rSges(4,1) * t17 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t3 = -t11 * t38 + t39;
t1 = t11 * t40 + t37;
t27 = t11 * pkin(4) + t50 * t9;
t26 = t11 * t8 + t49 * t9;
t25 = m(5) * rSges(5,2) - m(6) * t50 - m(7) * t49;
t4 = t11 * t37 + t40;
t2 = -t11 * t39 + t38;
t5 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - t23 * rSges(2,2)) + g(2) * (t23 * rSges(2,1) - t21 * rSges(2,2))) - m(3) * (g(1) * (-t10 * rSges(3,1) - t12 * rSges(3,2) - t42) + g(2) * (t12 * rSges(3,1) - t10 * rSges(3,2) + t13)) - m(4) * (-t34 + g(2) * t13 + (g(1) * t35 + g(2) * t30) * t12 + (-g(1) * t30 + g(2) * t35) * t10) - m(5) * (-t34 + g(2) * t41 + (g(1) * t36 + g(2) * t31) * t12 + (g(1) * (-t31 - t7) + g(2) * t36) * t10) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t42) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t41) + (-g(1) * t19 + g(2) * t27) * t12 + (g(1) * (-t27 - t7) - g(2) * t19) * t10) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t42) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t41) + (g(1) * t32 + g(2) * t26) * t12 + (g(1) * (-t26 - t7) + g(2) * t32) * t10) (-m(3) + t33) * g(3), t33 * (g(1) * t10 - g(2) * t12) (-t47 * t11 + t25 * t9) * g(3) + t48 * (t25 * t11 + t47 * t9) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2))) - m(7) * (g(1) * (-t4 * rSges(7,2) + t51 * t3) + g(2) * (t2 * rSges(7,2) - t51 * t1)) + (-m(6) * (-rSges(6,1) * t20 - rSges(6,2) * t22) - m(7) * (-rSges(7,1) * t20 - rSges(7,2) * t22 - t45)) * g(3) * t9, -m(7) * (-g(3) * t11 + t48 * t9)];
taug  = t5(:);
