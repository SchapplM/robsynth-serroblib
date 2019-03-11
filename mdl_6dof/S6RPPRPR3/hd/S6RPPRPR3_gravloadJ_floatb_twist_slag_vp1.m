% Calculate Gravitation load on the joints for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:43:59
% EndTime: 2019-03-09 01:44:01
% DurationCPUTime: 0.38s
% Computational Cost: add. (276->90), mult. (238->119), div. (0->0), fcn. (201->10), ass. (0->43)
t19 = qJ(4) + pkin(10);
t16 = cos(t19);
t42 = rSges(7,3) + pkin(8);
t55 = t42 * t16;
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t22 = sin(qJ(6));
t25 = cos(qJ(6));
t29 = rSges(7,1) * t25 - rSges(7,2) * t22 + pkin(5);
t54 = -(m(6) * rSges(6,1) + m(7) * t29) * t16 - m(5) * (rSges(5,1) * t26 - rSges(5,2) * t23);
t20 = qJ(1) + pkin(9);
t15 = sin(t20);
t47 = g(1) * t15;
t48 = pkin(4) * t26;
t51 = t48 * t47;
t50 = -m(6) - m(7);
t17 = cos(t20);
t46 = g(2) * t17;
t45 = t23 * pkin(4);
t24 = sin(qJ(1));
t44 = t24 * pkin(1);
t43 = rSges(5,3) + pkin(7);
t41 = t15 * t22;
t40 = t15 * t25;
t39 = t17 * t22;
t38 = t17 * t25;
t27 = cos(qJ(1));
t18 = t27 * pkin(1);
t37 = t17 * pkin(2) + t15 * qJ(3) + t18;
t36 = -m(4) - m(5) + t50;
t35 = t17 * qJ(3) - t44;
t34 = t15 * t45 + t37;
t32 = t23 * rSges(5,1) + t26 * rSges(5,2);
t14 = sin(t19);
t31 = t14 * rSges(6,1) + t16 * rSges(6,2);
t21 = -qJ(5) - pkin(7);
t30 = t15 * t21 + t17 * t45 + t35;
t28 = t14 * pkin(5) - t55;
t4 = t14 * t38 - t41;
t3 = t14 * t39 + t40;
t2 = t14 * t40 + t39;
t1 = -t14 * t41 + t38;
t5 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(1) * (-t15 * rSges(3,1) - t17 * rSges(3,2) - t44) + g(2) * (t17 * rSges(3,1) - t15 * rSges(3,2) + t18)) - m(4) * (g(1) * (t17 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t15 + t35) + g(2) * (-t17 * rSges(4,2) + t15 * rSges(4,3) + t37)) - m(5) * (g(1) * t35 + g(2) * t37 + (g(1) * t32 + g(2) * t43) * t17 + (g(1) * (-pkin(2) - t43) + g(2) * t32) * t15) - m(6) * (g(1) * t30 + g(2) * t34 + (g(1) * t31 + g(2) * (rSges(6,3) - t21)) * t17 + (g(1) * (-rSges(6,3) - pkin(2)) + g(2) * t31) * t15) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) + t30) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t34) + (g(1) * t28 - g(2) * t21) * t17 + (-g(1) * pkin(2) + g(2) * t28) * t15) (-m(3) + t36) * g(3), t36 * (-t46 + t47) m(5) * g(3) * t32 - m(6) * (t51 + g(3) * (-t31 - t45)) - m(7) * (t51 + g(3) * (-t29 * t14 - t45 + t55)) + ((m(6) * rSges(6,2) - m(7) * t42) * t14 + t54) * t47 + (-m(6) * (rSges(6,2) * t14 - t48) - m(7) * (-t42 * t14 - t48) - t54) * t46, t50 * (g(1) * t17 + g(2) * t15) -m(7) * (g(1) * (t1 * rSges(7,1) - rSges(7,2) * t2) + g(2) * (t3 * rSges(7,1) + t4 * rSges(7,2)) + g(3) * (-rSges(7,1) * t22 - rSges(7,2) * t25) * t16)];
taug  = t5(:);
