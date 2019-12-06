% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:04
% EndTime: 2019-12-05 17:44:06
% DurationCPUTime: 0.38s
% Computational Cost: add. (217->86), mult. (259->119), div. (0->0), fcn. (249->10), ass. (0->40)
t27 = sin(pkin(9));
t28 = sin(pkin(8));
t29 = cos(pkin(9));
t30 = cos(pkin(8));
t50 = (-rSges(4,1) * t29 + rSges(4,2) * t27 - pkin(2)) * t30 - pkin(1) + (-rSges(4,3) - qJ(3)) * t28;
t26 = pkin(9) + qJ(4);
t22 = qJ(5) + t26;
t17 = sin(t22);
t18 = cos(t22);
t33 = cos(qJ(1));
t32 = sin(qJ(1));
t43 = t32 * t30;
t5 = t17 * t43 + t33 * t18;
t6 = -t33 * t17 + t18 * t43;
t49 = t5 * rSges(6,1) + t6 * rSges(6,2);
t42 = t33 * t30;
t7 = t17 * t42 - t32 * t18;
t8 = -t32 * t17 - t18 * t42;
t48 = -t7 * rSges(6,1) + t8 * rSges(6,2);
t20 = sin(t26);
t47 = pkin(4) * t20;
t46 = g(1) * t28;
t24 = t33 * qJ(2);
t45 = g(3) * t24;
t44 = t27 * pkin(3);
t19 = t29 * pkin(3) + pkin(2);
t31 = -pkin(6) - qJ(3);
t41 = -m(4) - m(5) - m(6);
t40 = t27 * rSges(4,1) + t29 * rSges(4,2);
t39 = -rSges(6,1) * t17 - rSges(6,2) * t18;
t38 = -rSges(3,1) * t30 + rSges(3,2) * t28 - pkin(1);
t21 = cos(t26);
t11 = t20 * t42 - t32 * t21;
t9 = t20 * t43 + t33 * t21;
t35 = -t19 * t30 - pkin(1) + (-rSges(5,3) + t31) * t28;
t34 = -(pkin(4) * t21 + t19) * t30 - pkin(1) + (-rSges(6,3) - pkin(7) + t31) * t28;
t15 = t44 + t47;
t12 = -t32 * t20 - t21 * t42;
t10 = -t33 * t20 + t21 * t43;
t1 = [-m(2) * (g(2) * (-t33 * rSges(2,1) + t32 * rSges(2,2)) + g(3) * (-t32 * rSges(2,1) - t33 * rSges(2,2))) - m(3) * (t45 + (g(3) * rSges(3,3) + g(2) * t38) * t33 + (g(2) * (-rSges(3,3) - qJ(2)) + g(3) * t38) * t32) - m(4) * (t45 + (t50 * g(2) + g(3) * t40) * t33 + (g(2) * (-qJ(2) - t40) + t50 * g(3)) * t32) - m(5) * (g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2)) + g(3) * (-t10 * rSges(5,1) + t9 * rSges(5,2) + t24) + (g(2) * t35 + g(3) * t44) * t33 + (g(2) * (-qJ(2) - t44) + g(3) * t35) * t32) - m(6) * (g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2)) + g(3) * (-t6 * rSges(6,1) + t5 * rSges(6,2) + t24) + (g(2) * t34 + g(3) * t15) * t33 + (g(2) * (-qJ(2) - t15) + g(3) * t34) * t32), (-m(3) + t41) * (g(2) * t33 + g(3) * t32), t41 * (-g(1) * t30 + (-g(2) * t32 + g(3) * t33) * t28), -m(5) * (g(2) * (t9 * rSges(5,1) + t10 * rSges(5,2)) + g(3) * (-t11 * rSges(5,1) + t12 * rSges(5,2))) - m(6) * (g(2) * (pkin(4) * t9 + t49) + g(3) * (-pkin(4) * t11 + t48)) + (-m(5) * (-rSges(5,1) * t20 - rSges(5,2) * t21) - m(6) * (t39 - t47)) * t46, -m(6) * (g(2) * t49 + g(3) * t48 + t39 * t46)];
taug = t1(:);
