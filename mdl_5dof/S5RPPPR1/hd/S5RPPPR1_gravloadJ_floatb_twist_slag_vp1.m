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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:05
% EndTime: 2020-01-03 11:20:07
% DurationCPUTime: 0.28s
% Computational Cost: add. (178->66), mult. (164->92), div. (0->0), fcn. (148->10), ass. (0->30)
t19 = sin(pkin(9));
t20 = sin(pkin(8));
t21 = cos(pkin(9));
t22 = cos(pkin(8));
t38 = (rSges(5,3) + qJ(4)) * t20 + (rSges(5,1) * t21 - rSges(5,2) * t19 + pkin(3)) * t22;
t37 = -m(5) - m(6);
t36 = pkin(4) * t19;
t18 = qJ(1) + pkin(7);
t12 = sin(t18);
t24 = sin(qJ(1));
t15 = t24 * pkin(1);
t35 = t12 * pkin(2) + t15;
t34 = rSges(4,2) * t20;
t33 = t12 * t22;
t14 = cos(t18);
t32 = t14 * t22;
t31 = -m(4) + t37;
t25 = cos(qJ(1));
t16 = t25 * pkin(1);
t30 = t14 * pkin(2) + t12 * qJ(3) + t16;
t28 = t19 * rSges(5,1) + t21 * rSges(5,2);
t26 = (pkin(4) * t21 + pkin(3)) * t22 + (rSges(6,3) + pkin(6) + qJ(4)) * t20;
t17 = pkin(9) + qJ(5);
t13 = cos(t17);
t11 = sin(t17);
t4 = t11 * t12 + t13 * t32;
t3 = t11 * t32 - t12 * t13;
t2 = -t11 * t14 + t13 * t33;
t1 = -t11 * t33 - t13 * t14;
t5 = [-m(2) * (g(2) * (rSges(2,1) * t25 - t24 * rSges(2,2)) + g(3) * (t24 * rSges(2,1) + rSges(2,2) * t25)) - m(3) * (g(2) * (rSges(3,1) * t14 - rSges(3,2) * t12 + t16) + g(3) * (rSges(3,1) * t12 + rSges(3,2) * t14 + t15)) - m(4) * (g(2) * (rSges(4,3) * t12 + t30) + g(3) * (rSges(4,1) * t33 - t12 * t34 + t35) + (g(2) * (rSges(4,1) * t22 - t34) + g(3) * (-rSges(4,3) - qJ(3))) * t14) - m(5) * (g(2) * t30 + g(3) * t35 + (g(2) * t28 + t38 * g(3)) * t12 + (g(3) * (-qJ(3) - t28) + t38 * g(2)) * t14) - m(6) * (g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t30) + g(3) * (rSges(6,1) * t2 + rSges(6,2) * t1 + t35) + (g(2) * t36 + g(3) * t26) * t12 + (g(2) * t26 + g(3) * (-qJ(3) - t36)) * t14), (-m(3) + t31) * g(1), t31 * (-g(2) * t14 - g(3) * t12), t37 * (-g(1) * t22 + (g(2) * t12 - g(3) * t14) * t20), -m(6) * (g(2) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(1) * (-rSges(6,1) * t11 - rSges(6,2) * t13) * t20)];
taug = t5(:);
