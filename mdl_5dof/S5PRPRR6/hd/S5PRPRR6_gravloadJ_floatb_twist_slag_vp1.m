% Calculate Gravitation load on the joints for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:28
% EndTime: 2019-12-05 15:56:30
% DurationCPUTime: 0.45s
% Computational Cost: add. (290->82), mult. (509->129), div. (0->0), fcn. (582->12), ass. (0->46)
t26 = sin(pkin(5));
t32 = cos(qJ(2));
t46 = t26 * t32;
t59 = g(3) * t46;
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t58 = t29 * rSges(6,1) + t31 * rSges(6,2);
t30 = sin(qJ(2));
t25 = sin(pkin(9));
t42 = cos(pkin(5));
t39 = t25 * t42;
t41 = cos(pkin(9));
t14 = t41 * t30 + t32 * t39;
t36 = t42 * t41;
t12 = t25 * t30 - t32 * t36;
t54 = g(2) * t12;
t57 = g(1) * t14 + t54;
t23 = pkin(10) + qJ(4);
t21 = sin(t23);
t22 = cos(t23);
t34 = rSges(6,1) * t31 - rSges(6,2) * t29 + pkin(4);
t51 = rSges(6,3) + pkin(8);
t56 = t51 * t21 + t34 * t22;
t27 = cos(pkin(10));
t20 = t27 * pkin(3) + pkin(2);
t53 = t20 * t59;
t52 = g(3) * t26;
t13 = t25 * t32 + t30 * t36;
t28 = -pkin(7) - qJ(3);
t50 = -t12 * t20 - t13 * t28;
t15 = -t30 * t39 + t41 * t32;
t49 = -t14 * t20 - t15 * t28;
t48 = t25 * t26;
t47 = t26 * t30;
t43 = rSges(4,3) + qJ(3);
t40 = -m(4) - m(5) - m(6);
t38 = t26 * t41;
t37 = rSges(5,1) * t22 - rSges(5,2) * t21;
t35 = rSges(4,1) * t27 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t9 = t42 * t21 + t22 * t47;
t8 = -t21 * t47 + t42 * t22;
t5 = t15 * t22 + t21 * t48;
t4 = -t15 * t21 + t22 * t48;
t3 = t13 * t22 - t21 * t38;
t2 = -t13 * t21 - t22 * t38;
t1 = [(-m(2) - m(3) + t40) * g(3), -m(3) * (g(1) * (-t14 * rSges(3,1) - t15 * rSges(3,2)) + g(2) * (-t12 * rSges(3,1) - t13 * rSges(3,2)) + (rSges(3,1) * t32 - rSges(3,2) * t30) * t52) - m(4) * (g(1) * (-t35 * t14 + t43 * t15) + g(2) * t43 * t13 - t35 * t54 + (t43 * t30 + t35 * t32) * t52) - m(5) * (g(1) * (t15 * rSges(5,3) - t37 * t14 + t49) + g(2) * (t13 * rSges(5,3) - t37 * t12 + t50) + t53 + (t37 * t32 + (rSges(5,3) - t28) * t30) * t52) - m(6) * (g(1) * (t58 * t15 + t49) + g(2) * (t58 * t13 + t50) + t53 + ((-t28 + t58) * t30 + t56 * t32) * t52 - t57 * t56), t40 * (t57 - t59), -m(5) * (g(1) * (t4 * rSges(5,1) - t5 * rSges(5,2)) + g(2) * (t2 * rSges(5,1) - t3 * rSges(5,2)) + g(3) * (t8 * rSges(5,1) - t9 * rSges(5,2))) - m(6) * (g(3) * (t34 * t8 + t51 * t9) + (t34 * t2 + t51 * t3) * g(2) + (t34 * t4 + t51 * t5) * g(1)), -m(6) * (g(1) * ((t14 * t31 - t5 * t29) * rSges(6,1) + (-t14 * t29 - t5 * t31) * rSges(6,2)) + g(2) * ((t12 * t31 - t3 * t29) * rSges(6,1) + (-t12 * t29 - t3 * t31) * rSges(6,2)) + g(3) * ((-t9 * t29 - t31 * t46) * rSges(6,1) + (t29 * t46 - t9 * t31) * rSges(6,2)))];
taug = t1(:);
