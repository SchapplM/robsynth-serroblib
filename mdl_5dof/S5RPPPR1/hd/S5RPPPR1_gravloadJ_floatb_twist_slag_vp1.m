% Calculate Gravitation load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:47
% DurationCPUTime: 0.29s
% Computational Cost: add. (178->65), mult. (164->88), div. (0->0), fcn. (148->10), ass. (0->30)
t15 = sin(pkin(9));
t16 = sin(pkin(8));
t17 = cos(pkin(9));
t18 = cos(pkin(8));
t36 = -pkin(2) + (-rSges(5,3) - qJ(4)) * t16 + (-rSges(5,1) * t17 + rSges(5,2) * t15 - pkin(3)) * t18;
t35 = -m(5) - m(6);
t34 = pkin(4) * t15;
t20 = sin(qJ(1));
t33 = t20 * pkin(1);
t21 = cos(qJ(1));
t32 = t21 * pkin(1);
t14 = qJ(1) + pkin(7);
t10 = sin(t14);
t31 = t10 * t18;
t12 = cos(t14);
t30 = t12 * t18;
t29 = -m(4) + t35;
t28 = t12 * qJ(3) - t33;
t27 = t15 * rSges(5,1) + t17 * rSges(5,2);
t26 = -rSges(4,1) * t18 + rSges(4,2) * t16 - pkin(2);
t23 = -t18 * (pkin(4) * t17 + pkin(3)) - pkin(2) + (-rSges(6,3) - pkin(6) - qJ(4)) * t16;
t22 = -g(2) * t32 + g(3) * t28;
t13 = pkin(9) + qJ(5);
t11 = cos(t13);
t9 = sin(t13);
t4 = -t10 * t9 - t11 * t30;
t3 = -t10 * t11 + t9 * t30;
t2 = t11 * t31 - t12 * t9;
t1 = t11 * t12 + t9 * t31;
t5 = [-m(2) * (g(2) * (-rSges(2,1) * t21 + t20 * rSges(2,2)) + g(3) * (-t20 * rSges(2,1) - rSges(2,2) * t21)) - m(3) * (g(2) * (-t12 * rSges(3,1) + t10 * rSges(3,2) - t32) + g(3) * (-rSges(3,1) * t10 - rSges(3,2) * t12 - t33)) - m(4) * ((g(3) * rSges(4,3) + g(2) * t26) * t12 + (g(2) * (-rSges(4,3) - qJ(3)) + g(3) * t26) * t10 + t22) - m(5) * ((t36 * g(2) + g(3) * t27) * t12 + (g(2) * (-qJ(3) - t27) + t36 * g(3)) * t10 + t22) - m(6) * (g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) - t32) + g(3) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t28) + (g(2) * t23 + g(3) * t34) * t12 + (g(2) * (-qJ(3) - t34) + g(3) * t23) * t10), (-m(3) + t29) * g(1), t29 * (g(2) * t12 + g(3) * t10), t35 * (-g(1) * t18 + (-g(2) * t10 + g(3) * t12) * t16), -m(6) * (g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t3 + rSges(6,2) * t4) + g(1) * (-rSges(6,1) * t9 - rSges(6,2) * t11) * t16)];
taug = t5(:);
