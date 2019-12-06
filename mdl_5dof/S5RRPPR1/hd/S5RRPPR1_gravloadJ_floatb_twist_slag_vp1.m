% Calculate Gravitation load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (243->59), mult. (150->69), div. (0->0), fcn. (114->10), ass. (0->33)
t20 = pkin(9) + qJ(5);
t15 = sin(t20);
t16 = cos(t20);
t51 = rSges(6,1) * t16 - rSges(6,2) * t15;
t50 = rSges(6,3) + pkin(7) + qJ(4);
t23 = cos(pkin(9));
t49 = -pkin(4) * t23 - pkin(3) - t51;
t48 = rSges(5,3) + qJ(4);
t47 = -rSges(5,1) * t23 - pkin(3);
t46 = -m(5) - m(6);
t25 = sin(qJ(1));
t45 = pkin(1) * t25;
t26 = cos(qJ(1));
t44 = pkin(1) * t26;
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t43 = pkin(2) * t18;
t19 = cos(t21);
t42 = pkin(2) * t19;
t39 = rSges(5,2) * sin(pkin(9));
t37 = -t19 * rSges(3,1) + t18 * rSges(3,2);
t17 = pkin(8) + t21;
t12 = sin(t17);
t13 = cos(t17);
t35 = t12 * t39 + t48 * t13 - t43;
t34 = -rSges(3,1) * t18 - rSges(3,2) * t19;
t32 = -t13 * rSges(4,1) + t12 * rSges(4,2) - t42;
t31 = -rSges(4,1) * t12 - rSges(4,2) * t13 - t43;
t30 = -t42 + (t39 + t47) * t13;
t29 = (-g(2) * t48 + g(3) * t47) * t12;
t28 = t49 * t12 + t50 * t13 - t43;
t27 = -t50 * t12 + t49 * t13 - t42;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t26 + t25 * rSges(2,2)) + g(3) * (-t25 * rSges(2,1) - rSges(2,2) * t26)) - m(3) * (g(2) * (t37 - t44) + g(3) * (t34 - t45)) - m(4) * (g(2) * (t32 - t44) + g(3) * (t31 - t45)) - m(5) * (g(2) * (t30 - t44) + g(3) * (t35 - t45) + t29) - m(6) * (g(2) * (t27 - t44) + g(3) * (t28 - t45)), -m(3) * (g(2) * t37 + g(3) * t34) - m(4) * (g(2) * t32 + g(3) * t31) - m(5) * (g(2) * t30 + g(3) * t35 + t29) - m(6) * (g(2) * t27 + g(3) * t28), (-m(4) + t46) * g(1), t46 * (g(2) * t13 + g(3) * t12), -m(6) * (g(1) * t51 + (g(2) * t12 - g(3) * t13) * (rSges(6,1) * t15 + rSges(6,2) * t16))];
taug = t1(:);
