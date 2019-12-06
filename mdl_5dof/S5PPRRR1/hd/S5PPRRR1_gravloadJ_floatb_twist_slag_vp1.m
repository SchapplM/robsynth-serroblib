% Calculate Gravitation load on the joints for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:38
% DurationCPUTime: 0.23s
% Computational Cost: add. (181->39), mult. (163->58), div. (0->0), fcn. (134->8), ass. (0->27)
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t44 = g(1) * t21 + g(2) * t20;
t47 = rSges(6,3) + pkin(7);
t23 = cos(qJ(5));
t46 = -t23 * rSges(6,1) - pkin(4);
t19 = pkin(9) + qJ(3);
t18 = qJ(4) + t19;
t14 = sin(t18);
t43 = t46 * t14;
t16 = sin(t19);
t42 = pkin(3) * t16;
t22 = sin(qJ(5));
t39 = rSges(6,2) * t22;
t36 = t20 * t22;
t35 = t20 * t23;
t34 = t21 * t22;
t33 = t21 * t23;
t30 = -m(3) - m(4) - m(5) - m(6);
t15 = cos(t18);
t28 = t15 * rSges(5,1) - rSges(5,2) * t14;
t27 = -rSges(5,1) * t14 - rSges(5,2) * t15;
t26 = t47 * t14 + (-t39 - t46) * t15;
t25 = t44 * (t14 * t39 + t15 * t47);
t17 = cos(t19);
t13 = pkin(3) * t17;
t1 = [(-m(2) + t30) * g(3), t30 * (g(1) * t20 - g(2) * t21), -m(6) * t25 + (-m(4) * (rSges(4,1) * t17 - rSges(4,2) * t16) - m(5) * (t13 + t28) - m(6) * (t13 + t26)) * g(3) + t44 * (-m(4) * (-rSges(4,1) * t16 - rSges(4,2) * t17) - m(5) * (t27 - t42) - m(6) * (-t42 + t43)), -m(5) * (g(3) * t28 + t44 * t27) - m(6) * (g(3) * t26 + t44 * t43 + t25), -m(6) * (g(1) * ((-t15 * t34 + t35) * rSges(6,1) + (-t15 * t33 - t36) * rSges(6,2)) + g(2) * ((-t15 * t36 - t33) * rSges(6,1) + (-t15 * t35 + t34) * rSges(6,2)) + g(3) * (-rSges(6,1) * t22 - rSges(6,2) * t23) * t14)];
taug = t1(:);
