% Calculate Gravitation load on the joints for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:27
% EndTime: 2019-12-05 16:02:28
% DurationCPUTime: 0.42s
% Computational Cost: add. (208->81), mult. (507->128), div. (0->0), fcn. (582->10), ass. (0->42)
t56 = -rSges(5,3) - pkin(7);
t21 = sin(pkin(9));
t23 = cos(pkin(9));
t29 = cos(qJ(2));
t26 = sin(qJ(2));
t41 = cos(pkin(5));
t39 = t26 * t41;
t10 = t21 * t29 + t23 * t39;
t12 = -t21 * t39 + t23 * t29;
t55 = g(1) * t12 + g(2) * t10;
t38 = t29 * t41;
t11 = t21 * t38 + t23 * t26;
t9 = t21 * t26 - t23 * t38;
t54 = g(1) * t11 + g(2) * t9;
t22 = sin(pkin(5));
t49 = g(3) * t22;
t48 = rSges(6,3) + pkin(8);
t25 = sin(qJ(4));
t47 = t22 * t25;
t46 = t22 * t26;
t28 = cos(qJ(4));
t45 = t22 * t28;
t44 = t22 * t29;
t43 = pkin(2) * t44 + qJ(3) * t46;
t42 = rSges(4,3) + qJ(3);
t40 = -m(4) - m(5) - m(6);
t37 = g(3) * (pkin(7) * t44 + t43);
t36 = rSges(5,1) * t25 + rSges(5,2) * t28;
t24 = sin(qJ(5));
t27 = cos(qJ(5));
t35 = t24 * rSges(6,1) + t27 * rSges(6,2);
t34 = rSges(6,1) * t27 - rSges(6,2) * t24 + pkin(4);
t31 = t34 * t25 - t48 * t28;
t14 = -t25 * t44 + t41 * t28;
t13 = -t41 * t25 - t28 * t44;
t8 = t11 * pkin(2);
t7 = t9 * pkin(2);
t5 = t23 * t45 - t9 * t25;
t4 = t23 * t47 + t9 * t28;
t3 = t11 * t25 + t21 * t45;
t2 = t11 * t28 - t21 * t47;
t1 = [(-m(2) - m(3) + t40) * g(3), -m(3) * (g(1) * (-t11 * rSges(3,1) - t12 * rSges(3,2)) + g(2) * (-t9 * rSges(3,1) - t10 * rSges(3,2)) + (rSges(3,1) * t29 - rSges(3,2) * t26) * t49) - m(4) * (g(1) * (t11 * rSges(4,2) + t42 * t12 - t8) + g(2) * (t9 * rSges(4,2) + t42 * t10 - t7) + g(3) * ((-rSges(4,2) * t29 + rSges(4,3) * t26) * t22 + t43)) - m(5) * (g(1) * (t56 * t11 - t8) + g(2) * (t56 * t9 - t7) + t37 + (rSges(5,3) * t29 + t36 * t26) * t49 + t55 * (qJ(3) + t36)) - m(6) * (-g(1) * t8 - g(2) * t7 + t37 + (t31 * t26 + t35 * t29) * t49 + t54 * (-pkin(7) - t35) + t55 * (qJ(3) + t31)), t40 * (-g(3) * t44 + t54), -m(5) * (g(1) * (t2 * rSges(5,1) - t3 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t5 * rSges(5,2)) + g(3) * (t13 * rSges(5,1) - t14 * rSges(5,2))) - m(6) * (g(2) * (t34 * t4 - t48 * t5) + (t34 * t13 + t48 * t14) * g(3) + (t34 * t2 + t48 * t3) * g(1)), -m(6) * (g(1) * ((t12 * t27 - t3 * t24) * rSges(6,1) + (-t12 * t24 - t3 * t27) * rSges(6,2)) + g(2) * ((t10 * t27 + t5 * t24) * rSges(6,1) + (-t10 * t24 + t5 * t27) * rSges(6,2)) + g(3) * ((-t14 * t24 + t27 * t46) * rSges(6,1) + (-t14 * t27 - t24 * t46) * rSges(6,2)))];
taug = t1(:);
