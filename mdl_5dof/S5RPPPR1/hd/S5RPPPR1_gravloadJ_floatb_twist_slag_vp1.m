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
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
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
% StartTime: 2022-01-20 09:12:13
% EndTime: 2022-01-20 09:12:14
% DurationCPUTime: 0.34s
% Computational Cost: add. (178->66), mult. (164->92), div. (0->0), fcn. (148->10), ass. (0->30)
t18 = sin(pkin(9));
t19 = sin(pkin(8));
t20 = cos(pkin(9));
t21 = cos(pkin(8));
t39 = (rSges(5,3) + qJ(4)) * t19 + (rSges(5,1) * t20 - rSges(5,2) * t18 + pkin(3)) * t21;
t36 = -m(5) - m(6);
t35 = pkin(4) * t18;
t23 = sin(qJ(1));
t34 = t23 * pkin(1);
t33 = rSges(4,2) * t19;
t17 = qJ(1) + pkin(7);
t12 = sin(t17);
t32 = t12 * t21;
t14 = cos(t17);
t31 = t14 * t21;
t30 = -m(4) + t36;
t24 = cos(qJ(1));
t15 = t24 * pkin(1);
t29 = t14 * pkin(2) + t12 * qJ(3) + t15;
t28 = t14 * qJ(3) - t34;
t27 = t18 * rSges(5,1) + t20 * rSges(5,2);
t25 = (pkin(4) * t20 + pkin(3)) * t21 + (rSges(6,3) + pkin(6) + qJ(4)) * t19;
t16 = pkin(9) + qJ(5);
t13 = cos(t16);
t11 = sin(t16);
t4 = t11 * t12 + t13 * t31;
t3 = -t11 * t31 + t12 * t13;
t2 = t11 * t14 - t13 * t32;
t1 = t11 * t32 + t13 * t14;
t5 = [-m(2) * (g(1) * (-t23 * rSges(2,1) - rSges(2,2) * t24) + g(2) * (rSges(2,1) * t24 - t23 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t12 - rSges(3,2) * t14 - t34) + g(2) * (rSges(3,1) * t14 - rSges(3,2) * t12 + t15)) - m(4) * (g(1) * (rSges(4,3) * t14 + t28) + g(2) * (rSges(4,1) * t31 - t14 * t33 + t29) + (g(1) * (-rSges(4,1) * t21 - pkin(2) + t33) + g(2) * rSges(4,3)) * t12) - m(5) * (g(1) * t28 + g(2) * t29 + (g(1) * t27 + t39 * g(2)) * t14 + (g(2) * t27 + (-pkin(2) - t39) * g(1)) * t12) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1 + t28) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t3 + t29) + (g(1) * t35 + g(2) * t25) * t14 + (g(1) * (-pkin(2) - t25) + g(2) * t35) * t12), (-m(3) + t30) * g(3), t30 * (g(1) * t12 - g(2) * t14), t36 * (-g(3) * t21 + (g(1) * t14 + g(2) * t12) * t19), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t11 - rSges(6,2) * t13) * t19)];
taug = t5(:);
