% Calculate Gravitation load on the joints for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:12
% EndTime: 2019-03-09 02:00:13
% DurationCPUTime: 0.61s
% Computational Cost: add. (402->100), mult. (346->137), div. (0->0), fcn. (324->10), ass. (0->46)
t60 = -m(6) - m(7);
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t55 = g(2) * t17;
t59 = rSges(7,1) + pkin(5);
t21 = pkin(10) + qJ(4);
t16 = sin(t21);
t18 = cos(t21);
t58 = t18 * pkin(4) + t16 * pkin(8);
t19 = cos(t22);
t57 = g(1) * t19 + t55;
t41 = rSges(7,3) + qJ(6);
t25 = -pkin(7) - qJ(3);
t54 = g(2) * t25;
t53 = g(3) * t16;
t27 = sin(qJ(1));
t52 = t27 * pkin(1);
t24 = cos(pkin(10));
t15 = t24 * pkin(3) + pkin(2);
t29 = cos(qJ(1));
t20 = t29 * pkin(1);
t50 = t19 * t15 + t20;
t49 = t16 * t19;
t26 = sin(qJ(5));
t48 = t17 * t26;
t28 = cos(qJ(5));
t47 = t17 * t28;
t46 = t18 * t19;
t45 = t19 * t26;
t44 = t19 * t28;
t43 = rSges(5,3) - t25;
t42 = rSges(4,3) + qJ(3);
t40 = g(1) * t52;
t39 = -m(4) - m(5) + t60;
t38 = pkin(4) * t46 + pkin(8) * t49 + t50;
t37 = -t19 * t25 - t52;
t36 = t18 * rSges(5,1) - t16 * rSges(5,2);
t35 = rSges(6,1) * t28 - rSges(6,2) * t26;
t34 = -t15 - t58;
t33 = rSges(4,1) * t24 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t31 = t41 * t26 + t59 * t28;
t4 = t18 * t44 + t48;
t3 = t18 * t45 - t47;
t2 = t18 * t47 - t45;
t1 = t18 * t48 + t44;
t5 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t27 * rSges(2,2))) - m(3) * (g(1) * (-t17 * rSges(3,1) - t19 * rSges(3,2) - t52) + g(2) * (t19 * rSges(3,1) - t17 * rSges(3,2) + t20)) - m(4) * (-t40 + g(2) * t20 + (g(1) * t42 + g(2) * t33) * t19 + (-g(1) * t33 + g(2) * t42) * t17) - m(5) * (-t40 + g(2) * t50 + (g(1) * t43 + g(2) * t36) * t19 + (g(1) * (-t15 - t36) + g(2) * t43) * t17) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) + t37) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t49 + t38) + (g(1) * (-t16 * rSges(6,3) + t34) - t54) * t17) - m(7) * (g(1) * (-t41 * t1 - t59 * t2 + t37) + g(2) * (rSges(7,2) * t49 + t41 * t3 + t59 * t4 + t38) + (g(1) * (-t16 * rSges(7,2) + t34) - t54) * t17) (-m(3) + t39) * g(3), t39 * (g(1) * t17 - g(2) * t19) (-m(5) * (g(3) * rSges(5,1) - t57 * rSges(5,2)) - m(6) * (t57 * rSges(6,3) + g(3) * t35) - m(7) * (t57 * rSges(7,2) + g(3) * t31)) * t18 + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * g(3) + t57 * (m(5) * rSges(5,1) - m(6) * (-pkin(4) - t35) - m(7) * (-pkin(4) - t31))) * t16 + t60 * (g(3) * t58 + (g(1) * t46 + t18 * t55) * pkin(8)) -m(6) * (g(1) * (-t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) - t2 * rSges(6,2))) - m(7) * (g(1) * (-t3 * t59 + t41 * t4) + g(2) * (-t1 * t59 + t41 * t2)) + (-m(6) * (-rSges(6,1) * t26 - rSges(6,2) * t28) - m(7) * (-t59 * t26 + t41 * t28)) * t53, -m(7) * (g(1) * t3 + g(2) * t1 + t26 * t53)];
taug  = t5(:);
