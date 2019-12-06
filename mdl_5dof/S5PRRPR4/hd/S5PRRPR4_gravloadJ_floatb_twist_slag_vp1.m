% Calculate Gravitation load on the joints for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:29
% EndTime: 2019-12-05 16:21:31
% DurationCPUTime: 0.39s
% Computational Cost: add. (207->77), mult. (270->114), div. (0->0), fcn. (247->10), ass. (0->36)
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t50 = g(1) * t20 + g(2) * t19;
t49 = -m(5) - m(6);
t18 = qJ(3) + pkin(9);
t15 = qJ(5) + t18;
t10 = sin(t15);
t11 = cos(t15);
t25 = cos(qJ(2));
t41 = t19 * t25;
t48 = (-t10 * t41 - t20 * t11) * rSges(6,1) + (t20 * t10 - t11 * t41) * rSges(6,2);
t40 = t20 * t25;
t47 = (-t10 * t40 + t19 * t11) * rSges(6,1) + (-t19 * t10 - t11 * t40) * rSges(6,2);
t23 = sin(qJ(2));
t44 = g(3) * t23;
t22 = sin(qJ(3));
t43 = t22 * pkin(3);
t42 = rSges(4,3) + pkin(6);
t14 = cos(t18);
t24 = cos(qJ(3));
t16 = t24 * pkin(3);
t8 = pkin(4) * t14 + t16;
t39 = t22 * t25;
t38 = t24 * t25;
t21 = -qJ(4) - pkin(6);
t37 = rSges(5,3) - t21;
t36 = rSges(6,3) + pkin(7) - t21;
t34 = -rSges(6,1) * t10 - rSges(6,2) * t11;
t33 = t24 * rSges(4,1) - t22 * rSges(4,2) + pkin(2);
t32 = rSges(6,1) * t11 - rSges(6,2) * t10 + pkin(2) + t8;
t31 = t19 * t24 - t20 * t39;
t30 = -t19 * t39 - t20 * t24;
t13 = sin(t18);
t29 = t14 * rSges(5,1) - t13 * rSges(5,2) + pkin(2) + t16;
t7 = -pkin(4) * t13 - t43;
t1 = [(-m(2) - m(3) - m(4) + t49) * g(3), -m(3) * (g(3) * (t25 * rSges(3,1) - t23 * rSges(3,2)) + t50 * (-rSges(3,1) * t23 - rSges(3,2) * t25)) - m(4) * (g(3) * (t42 * t23 + t33 * t25) + t50 * (-t33 * t23 + t42 * t25)) - m(5) * (g(3) * (t37 * t23 + t29 * t25) + t50 * (-t29 * t23 + t37 * t25)) - m(6) * (g(3) * (t36 * t23 + t32 * t25) + t50 * (-t32 * t23 + t36 * t25)), -m(4) * (g(1) * (t31 * rSges(4,1) + (-t19 * t22 - t20 * t38) * rSges(4,2)) + g(2) * (t30 * rSges(4,1) + (-t19 * t38 + t20 * t22) * rSges(4,2))) - m(5) * (g(1) * ((-t13 * t40 + t19 * t14) * rSges(5,1) + (-t19 * t13 - t14 * t40) * rSges(5,2) + t31 * pkin(3)) + g(2) * ((-t13 * t41 - t20 * t14) * rSges(5,1) + (t20 * t13 - t14 * t41) * rSges(5,2) + t30 * pkin(3))) - m(6) * (g(1) * (t19 * t8 + t7 * t40 + t47) + g(2) * (-t20 * t8 + t7 * t41 + t48)) + (-m(4) * (-rSges(4,1) * t22 - rSges(4,2) * t24) - m(5) * (-rSges(5,1) * t13 - rSges(5,2) * t14 - t43) - m(6) * (t34 + t7)) * t44, t49 * (-g(3) * t25 + t50 * t23), -m(6) * (g(1) * t47 + g(2) * t48 + t34 * t44)];
taug = t1(:);
