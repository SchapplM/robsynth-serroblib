% Calculate Gravitation load on the joints for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:43
% EndTime: 2019-03-09 02:02:44
% DurationCPUTime: 0.47s
% Computational Cost: add. (307->100), mult. (340->136), div. (0->0), fcn. (318->8), ass. (0->46)
t57 = -m(6) - m(7);
t47 = rSges(7,1) + pkin(5);
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t40 = t24 * rSges(5,2);
t56 = t21 * rSges(5,1) + t40;
t37 = rSges(7,3) + qJ(6);
t55 = -pkin(2) - pkin(7);
t54 = m(5) * rSges(5,1);
t53 = m(5) * rSges(5,2);
t19 = qJ(1) + pkin(9);
t15 = sin(t19);
t52 = g(1) * t15;
t16 = cos(t19);
t51 = g(2) * t16;
t50 = g(2) * t24;
t49 = g(3) * t24;
t22 = sin(qJ(1));
t48 = t22 * pkin(1);
t17 = t24 * pkin(8);
t46 = -rSges(7,2) - pkin(8);
t45 = -rSges(6,3) - pkin(8);
t44 = t15 * t21;
t20 = sin(qJ(5));
t43 = t20 * t21;
t23 = cos(qJ(5));
t41 = t21 * t23;
t39 = t24 * rSges(7,2);
t38 = t24 * rSges(6,3);
t25 = cos(qJ(1));
t18 = t25 * pkin(1);
t36 = t16 * pkin(2) + t15 * qJ(3) + t18;
t35 = g(1) * t55;
t33 = -m(4) - m(5) + t57;
t32 = t16 * qJ(3) - t48;
t31 = t16 * pkin(7) + t36;
t30 = pkin(4) * t44 + t31;
t29 = rSges(6,1) * t23 - rSges(6,2) * t20;
t28 = t32 + (t21 * pkin(4) - t17) * t16;
t27 = t37 * t20 + t47 * t23;
t26 = t54 - m(6) * (-pkin(4) - t29) - m(7) * (-pkin(4) - t27);
t4 = -t15 * t20 + t16 * t41;
t3 = t15 * t23 + t16 * t43;
t2 = t15 * t41 + t16 * t20;
t1 = t15 * t43 - t16 * t23;
t5 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * (g(1) * (-t15 * rSges(3,1) - t16 * rSges(3,2) - t48) + g(2) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t18)) - m(4) * (g(1) * (t16 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t15 + t32) + g(2) * (-t16 * rSges(4,2) + t15 * rSges(4,3) + t36)) - m(5) * (g(1) * (t56 * t16 + t32) + g(2) * (t16 * rSges(5,3) + t31) + (g(1) * (-rSges(5,3) + t55) + g(2) * t56) * t15) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t16 * t38 + t28) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t30) + (t45 * t50 + t35) * t15) - m(7) * (g(1) * (-t16 * t39 + t37 * t3 + t47 * t4 + t28) + g(2) * (t37 * t1 + t47 * t2 + t30) + (t46 * t50 + t35) * t15) (-m(3) + t33) * g(3), t33 * (-t51 + t52) ((-m(6) * rSges(6,3) - m(7) * rSges(7,2) + t53) * t21 + (-m(6) * t29 - m(7) * t27 - t54) * t24) * t52 + ((-m(6) * t45 - m(7) * t46 - t53) * t21 + t26 * t24) * t51 + t57 * g(1) * (t15 * t24 * pkin(4) + pkin(8) * t44) + (m(5) * t40 - m(6) * (t17 + t38) - m(7) * (t17 + t39) + t26 * t21) * g(3), -m(6) * (g(1) * (-t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t4 * rSges(6,2))) - m(7) * (g(1) * (-t47 * t1 + t37 * t2) + g(2) * (t47 * t3 - t37 * t4)) + (-m(6) * (-rSges(6,1) * t20 - rSges(6,2) * t23) - m(7) * (-t47 * t20 + t37 * t23)) * t49, -m(7) * (g(1) * t1 - g(2) * t3 + t20 * t49)];
taug  = t5(:);
