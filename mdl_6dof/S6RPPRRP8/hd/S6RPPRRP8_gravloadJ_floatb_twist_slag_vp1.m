% Calculate Gravitation load on the joints for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:27
% EndTime: 2019-03-09 02:15:28
% DurationCPUTime: 0.56s
% Computational Cost: add. (277->102), mult. (364->139), div. (0->0), fcn. (342->8), ass. (0->47)
t62 = -m(6) - m(7);
t51 = rSges(7,1) + pkin(5);
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t61 = g(1) * t28 + g(2) * t26;
t40 = rSges(7,3) + qJ(6);
t21 = pkin(9) + qJ(4);
t17 = cos(t21);
t60 = t61 * t17;
t59 = m(5) * rSges(5,1);
t58 = m(5) * rSges(5,2);
t22 = sin(pkin(9));
t57 = pkin(3) * t22;
t56 = g(1) * t26;
t53 = g(2) * t28;
t52 = g(3) * t17;
t50 = -rSges(7,2) - pkin(8);
t49 = -rSges(6,3) - pkin(8);
t16 = sin(t21);
t48 = t16 * t26;
t47 = t17 * rSges(5,2);
t25 = sin(qJ(5));
t46 = t26 * t25;
t27 = cos(qJ(5));
t45 = t26 * t27;
t44 = t28 * t25;
t43 = t28 * t27;
t42 = t28 * pkin(1) + t26 * qJ(2);
t41 = rSges(4,3) + qJ(3);
t39 = t26 * t57 + t42;
t19 = t28 * qJ(2);
t24 = -pkin(7) - qJ(3);
t38 = t26 * t24 + t28 * t57 + t19;
t36 = -m(4) - m(5) + t62;
t35 = rSges(4,1) * t22 + rSges(4,2) * cos(pkin(9));
t34 = t16 * rSges(5,1) + t47;
t33 = rSges(6,1) * t27 - rSges(6,2) * t25;
t32 = t28 * t16 * pkin(4) - t26 * pkin(1) + t38;
t31 = pkin(4) * t48 - t28 * t24 + t39;
t30 = t40 * t25 + t51 * t27;
t29 = t59 - m(6) * (-pkin(4) - t33) - m(7) * (-pkin(4) - t30);
t13 = t17 * pkin(8);
t4 = t16 * t43 - t46;
t3 = t16 * t44 + t45;
t2 = t16 * t45 + t44;
t1 = t16 * t46 - t43;
t5 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (t28 * rSges(3,3) + t19 + (rSges(3,2) - pkin(1)) * t26) + g(2) * (-t28 * rSges(3,2) + t26 * rSges(3,3) + t42)) - m(4) * (g(1) * t19 + g(2) * t42 + (g(1) * t35 + g(2) * t41) * t28 + (g(1) * (-pkin(1) - t41) + g(2) * t35) * t26) - m(5) * (g(1) * t38 + g(2) * t39 + (g(1) * t34 + g(2) * (rSges(5,3) - t24)) * t28 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t34) * t26) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t32) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t31) + t49 * t60) - m(7) * (g(1) * (t40 * t3 + t51 * t4 + t32) + g(2) * (t40 * t1 + t51 * t2 + t31) + t50 * t60) (-m(3) + t36) * (-t53 + t56) t36 * t61 ((-m(6) * rSges(6,3) - m(7) * rSges(7,2) + t58) * t16 + (-m(6) * t33 - m(7) * t30 - t59) * t17) * t56 + ((-m(6) * t49 - m(7) * t50 - t58) * t16 + t29 * t17) * t53 + t62 * g(1) * (t26 * t17 * pkin(4) + pkin(8) * t48) + (m(5) * t47 - m(6) * (t17 * rSges(6,3) + t13) - m(7) * (t17 * rSges(7,2) + t13) + t29 * t16) * g(3), -m(6) * (g(1) * (-t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t4 * rSges(6,2))) - m(7) * (g(1) * (-t51 * t1 + t40 * t2) + g(2) * (t51 * t3 - t40 * t4)) + (-m(6) * (-rSges(6,1) * t25 - rSges(6,2) * t27) - m(7) * (-t51 * t25 + t40 * t27)) * t52, -m(7) * (g(1) * t1 - g(2) * t3 + t25 * t52)];
taug  = t5(:);
