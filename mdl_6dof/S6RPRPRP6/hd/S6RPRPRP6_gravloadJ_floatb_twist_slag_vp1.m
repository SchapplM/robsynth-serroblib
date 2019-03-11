% Calculate Gravitation load on the joints for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:03
% EndTime: 2019-03-09 03:18:05
% DurationCPUTime: 0.56s
% Computational Cost: add. (328->114), mult. (389->158), div. (0->0), fcn. (354->8), ass. (0->54)
t63 = rSges(7,1) + pkin(5);
t22 = -pkin(7) - qJ(2);
t56 = pkin(4) - t22;
t47 = rSges(7,3) + qJ(6) + pkin(8);
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t60 = g(2) * t24;
t34 = g(1) * t26 + t60;
t25 = cos(qJ(5));
t62 = pkin(5) * t25;
t18 = pkin(9) + qJ(3);
t17 = cos(t18);
t59 = g(3) * t17;
t13 = t17 * pkin(3);
t57 = rSges(6,3) + pkin(8);
t55 = rSges(5,2) * t17;
t54 = t17 * t26;
t23 = sin(qJ(5));
t53 = t23 * t26;
t52 = t24 * t23;
t51 = t24 * t25;
t50 = t25 * t26;
t49 = rSges(5,1) - t22;
t48 = rSges(4,3) - t22;
t16 = sin(t18);
t12 = t16 * qJ(4);
t46 = t12 + t13;
t45 = t62 + t56;
t44 = qJ(4) * t17;
t43 = rSges(3,3) + qJ(2);
t42 = -m(5) - m(6) - m(7);
t41 = -pkin(3) - t57;
t20 = cos(pkin(9));
t14 = pkin(2) * t20 + pkin(1);
t9 = t26 * t14;
t40 = pkin(3) * t54 + t26 * t12 + t9;
t39 = pkin(5) * t16 * t23;
t38 = -pkin(3) - t47;
t37 = -t14 - t12;
t36 = g(1) * t41;
t35 = g(1) * t38;
t33 = rSges(4,1) * t17 - rSges(4,2) * t16;
t31 = rSges(6,1) * t23 + rSges(6,2) * t25;
t30 = rSges(5,3) * t16 - t55;
t29 = rSges(3,1) * t20 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t2 = t16 * t50 - t52;
t4 = t16 * t51 + t53;
t28 = rSges(7,2) * t25 + t63 * t23;
t6 = t24 * t44;
t8 = t26 * t44;
t27 = g(1) * t8 + g(2) * t6 + g(3) * t46;
t5 = -t16 * t52 + t50;
t3 = t16 * t53 + t51;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2))) - m(3) * ((g(1) * t43 + g(2) * t29) * t26 + (-g(1) * t29 + g(2) * t43) * t24) - m(4) * (g(2) * t9 + (g(1) * t48 + g(2) * t33) * t26 + (g(1) * (-t14 - t33) + g(2) * t48) * t24) - m(5) * (g(2) * t40 + (g(1) * t49 + g(2) * t30) * t26 + (g(1) * (-t30 + t37 - t13) + g(2) * t49) * t24) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t40) + (g(2) * t57 * t17 + g(1) * t56) * t26 + (g(1) * t37 + g(2) * t56 + t17 * t36) * t24) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t40) + (g(1) * t45 + g(2) * (t47 * t17 + t39)) * t26 + (g(1) * (t37 - t39) + g(2) * t45 + t17 * t35) * t24) (-m(3) - m(4) + t42) * (g(1) * t24 - g(2) * t26) -m(4) * (g(3) * t33 + t34 * (-rSges(4,1) * t16 - rSges(4,2) * t17)) - m(5) * (g(1) * (rSges(5,3) * t54 + t8) + g(2) * (t24 * t17 * rSges(5,3) + t6) + g(3) * (t46 - t55) + (g(3) * rSges(5,3) + t34 * (rSges(5,2) - pkin(3))) * t16) - m(6) * ((g(3) * t57 + t34 * t31) * t17 + (g(3) * t31 + t26 * t36 + t41 * t60) * t16 + t27) - m(7) * ((g(3) * t47 + t34 * t28) * t17 + (g(3) * t28 + t26 * t35 + t38 * t60) * t16 + t27) t42 * (t34 * t16 - t59) -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5)) - m(7) * (g(1) * (-t3 * rSges(7,2) + t63 * t2) + g(2) * (t5 * rSges(7,2) + t63 * t4)) + (-m(6) * (-rSges(6,1) * t25 + rSges(6,2) * t23) - m(7) * (-rSges(7,1) * t25 + rSges(7,2) * t23 - t62)) * t59, -m(7) * (g(3) * t16 + t34 * t17)];
taug  = t1(:);
