% Calculate Gravitation load on the joints for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:53
% EndTime: 2019-03-09 01:52:55
% DurationCPUTime: 0.55s
% Computational Cost: add. (261->94), mult. (298->127), div. (0->0), fcn. (265->10), ass. (0->44)
t27 = sin(qJ(1));
t28 = cos(qJ(1));
t65 = -g(1) * t27 + g(2) * t28;
t43 = rSges(7,3) + pkin(8) + qJ(5);
t62 = g(1) * t28 + g(2) * t27;
t23 = cos(pkin(10));
t10 = t23 * pkin(5) + pkin(4);
t19 = pkin(10) + qJ(6);
t12 = sin(t19);
t14 = cos(t19);
t21 = sin(pkin(10));
t29 = m(6) * (rSges(6,1) * t23 - rSges(6,2) * t21 + pkin(4)) + m(7) * (rSges(7,1) * t14 - rSges(7,2) * t12 + t10) + m(5) * rSges(5,1);
t40 = rSges(6,3) + qJ(5);
t61 = -m(5) * rSges(5,2) + m(6) * t40 + m(7) * t43;
t60 = -m(6) - m(7);
t22 = sin(pkin(9));
t57 = pkin(3) * t22;
t56 = pkin(5) * t21;
t20 = pkin(9) + qJ(4);
t13 = sin(t20);
t51 = t27 * t13;
t50 = t27 * t14;
t49 = t27 * t21;
t48 = t27 * t23;
t47 = t28 * t13;
t46 = t28 * t14;
t45 = t28 * t21;
t44 = t28 * t23;
t42 = t28 * pkin(1) + t27 * qJ(2);
t41 = rSges(4,3) + qJ(3);
t17 = t28 * qJ(2);
t26 = -pkin(7) - qJ(3);
t39 = t27 * t26 + t28 * t57 + t17;
t38 = t27 * t57 + t42;
t37 = -m(4) - m(5) + t60;
t35 = rSges(4,1) * t22 + rSges(4,2) * cos(pkin(9));
t15 = cos(t20);
t34 = t13 * rSges(5,1) + t15 * rSges(5,2);
t31 = t13 * t10 - t43 * t15;
t5 = -t27 * t12 + t13 * t46;
t4 = t12 * t47 + t50;
t3 = t28 * t12 + t13 * t50;
t2 = -t12 * t51 + t46;
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t27 * rSges(2,2))) - m(3) * (g(1) * (t28 * rSges(3,3) + t17 + (rSges(3,2) - pkin(1)) * t27) + g(2) * (-t28 * rSges(3,2) + t27 * rSges(3,3) + t42)) - m(4) * (g(1) * t17 + g(2) * t42 + (g(1) * t35 + g(2) * t41) * t28 + (g(1) * (-pkin(1) - t41) + g(2) * t35) * t27) - m(5) * (g(1) * t39 + g(2) * t38 + (g(1) * t34 + g(2) * (rSges(5,3) - t26)) * t28 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t34) * t27) - m(6) * (g(1) * (pkin(4) * t47 - t27 * pkin(1) + (t13 * t44 - t49) * rSges(6,1) + (-t13 * t45 - t48) * rSges(6,2) + t39) + g(2) * (-t28 * t26 + pkin(4) * t51 + (t13 * t48 + t45) * rSges(6,1) + (-t13 * t49 + t44) * rSges(6,2) + t38) - t62 * t15 * t40) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2) + t39) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t38) + (g(1) * t31 + g(2) * (-t26 + t56)) * t28 + (g(1) * (-pkin(1) - t56) + g(2) * t31) * t27) -(-m(3) + t37) * t65, t37 * t62 (t13 * t29 - t15 * t61) * g(3) + t65 * (t61 * t13 + t29 * t15) t60 * (g(3) * t13 + t65 * t15) -m(7) * (g(1) * (t2 * rSges(7,1) - t3 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t5 * rSges(7,2)) + g(3) * (-rSges(7,1) * t12 - rSges(7,2) * t14) * t15)];
taug  = t1(:);
