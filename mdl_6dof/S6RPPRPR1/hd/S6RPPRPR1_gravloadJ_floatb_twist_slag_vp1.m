% Calculate Gravitation load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:08
% EndTime: 2019-03-09 01:39:09
% DurationCPUTime: 0.43s
% Computational Cost: add. (360->85), mult. (280->113), div. (0->0), fcn. (247->12), ass. (0->45)
t18 = pkin(10) + qJ(4);
t11 = sin(t18);
t14 = cos(t18);
t20 = sin(pkin(11));
t22 = cos(pkin(11));
t34 = rSges(6,1) * t22 - rSges(6,2) * t20 + pkin(4);
t40 = rSges(6,3) + qJ(5);
t56 = t11 * t40 + t34 * t14;
t53 = rSges(7,3) + pkin(8) + qJ(5);
t19 = qJ(1) + pkin(9);
t12 = sin(t19);
t15 = cos(t19);
t52 = g(1) * t15 + g(2) * t12;
t17 = pkin(11) + qJ(6);
t10 = sin(t17);
t13 = cos(t17);
t8 = t22 * pkin(5) + pkin(4);
t51 = m(6) * t34 + m(7) * (rSges(7,1) * t13 - rSges(7,2) * t10 + t8) + m(5) * rSges(5,1);
t50 = -m(6) - m(7);
t26 = sin(qJ(1));
t46 = t26 * pkin(1);
t27 = cos(qJ(1));
t16 = t27 * pkin(1);
t23 = cos(pkin(10));
t9 = t23 * pkin(3) + pkin(2);
t45 = t15 * t9 + t16;
t44 = t12 * t14;
t43 = t15 * t14;
t25 = -pkin(7) - qJ(3);
t42 = rSges(5,3) - t25;
t41 = rSges(4,3) + qJ(3);
t39 = g(1) * t46;
t38 = -m(4) - m(5) + t50;
t37 = pkin(5) * t20 - t25;
t36 = t14 * rSges(5,1) - t11 * rSges(5,2);
t35 = rSges(4,1) * t23 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t32 = t20 * rSges(6,1) + t22 * rSges(6,2) - t25;
t31 = g(2) * t45 - t39;
t30 = t53 * t11 + t14 * t8;
t29 = m(5) * rSges(5,2) - m(6) * t40 - m(7) * t53;
t5 = t12 * t10 + t13 * t43;
t4 = -t10 * t43 + t12 * t13;
t3 = t15 * t10 - t13 * t44;
t2 = t10 * t44 + t15 * t13;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (-t12 * rSges(3,1) - t15 * rSges(3,2) - t46) + g(2) * (t15 * rSges(3,1) - t12 * rSges(3,2) + t16)) - m(4) * (-t39 + g(2) * t16 + (g(1) * t41 + g(2) * t35) * t15 + (-g(1) * t35 + g(2) * t41) * t12) - m(5) * ((g(1) * t42 + g(2) * t36) * t15 + (g(1) * (-t36 - t9) + g(2) * t42) * t12 + t31) - m(6) * ((g(1) * t32 + t56 * g(2)) * t15 + (g(2) * t32 + (-t9 - t56) * g(1)) * t12 + t31) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2) - t46) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t45) + (g(1) * t37 + g(2) * t30) * t15 + (g(1) * (-t30 - t9) + g(2) * t37) * t12) (-m(3) + t38) * g(3), t38 * (g(1) * t12 - g(2) * t15) (t29 * t11 - t51 * t14) * g(3) + t52 * (t51 * t11 + t29 * t14) t50 * (-g(3) * t14 + t52 * t11) -m(7) * (g(1) * (t4 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (-t2 * rSges(7,1) + t3 * rSges(7,2)) + g(3) * (-rSges(7,1) * t10 - rSges(7,2) * t13) * t11)];
taug  = t1(:);
