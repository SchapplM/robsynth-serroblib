% Calculate Gravitation load on the joints for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:09
% EndTime: 2019-12-31 20:32:11
% DurationCPUTime: 0.51s
% Computational Cost: add. (299->111), mult. (375->162), div. (0->0), fcn. (353->10), ass. (0->48)
t30 = -pkin(7) - qJ(3);
t47 = rSges(5,3) - t30;
t46 = rSges(6,3) + pkin(8) - t30;
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t53 = g(2) * t32;
t58 = g(1) * t34 + t53;
t27 = pkin(9) + qJ(4);
t21 = qJ(5) + t27;
t16 = sin(t21);
t17 = cos(t21);
t33 = cos(qJ(2));
t49 = t32 * t33;
t5 = t16 * t49 + t34 * t17;
t6 = t34 * t16 - t17 * t49;
t57 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t48 = t34 * t33;
t7 = -t16 * t48 + t32 * t17;
t8 = t32 * t16 + t17 * t48;
t56 = t7 * rSges(6,1) - t8 * rSges(6,2);
t19 = sin(t27);
t55 = pkin(4) * t19;
t31 = sin(qJ(2));
t52 = g(3) * t31;
t28 = sin(pkin(9));
t51 = t28 * pkin(3);
t29 = cos(pkin(9));
t18 = t29 * pkin(3) + pkin(2);
t50 = t31 * rSges(3,2);
t45 = t34 * pkin(1) + t32 * pkin(6);
t44 = rSges(4,3) + qJ(3);
t43 = t44 * t34;
t42 = t33 * rSges(3,1) - t50;
t40 = -rSges(6,1) * t16 - rSges(6,2) * t17;
t39 = rSges(4,1) * t29 - rSges(4,2) * t28 + pkin(2);
t20 = cos(t27);
t11 = -t19 * t48 + t32 * t20;
t9 = t19 * t49 + t34 * t20;
t38 = rSges(5,1) * t20 - rSges(5,2) * t19 + t18;
t14 = pkin(4) * t20 + t18;
t37 = rSges(6,1) * t17 - rSges(6,2) * t16 + t14;
t36 = t33 * t18 + t47 * t31;
t35 = t33 * t14 + t46 * t31;
t24 = t34 * pkin(6);
t15 = t51 + t55;
t12 = t32 * t19 + t20 * t48;
t10 = t34 * t19 - t20 * t49;
t1 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) - t32 * rSges(2,2))) - m(3) * (g(1) * (t34 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t48 - t34 * t50 + t45) + (g(1) * (-pkin(1) - t42) + g(2) * rSges(3,3)) * t32) - m(4) * (g(1) * (-pkin(2) * t49 - t32 * pkin(1) + t24 + (t34 * t28 - t29 * t49) * rSges(4,1) + (t28 * t49 + t34 * t29) * rSges(4,2)) + g(2) * (pkin(2) * t48 + (t32 * t28 + t29 * t48) * rSges(4,1) + (-t28 * t48 + t32 * t29) * rSges(4,2) + t45) + (-g(1) * t44 * t32 + g(2) * t43) * t31) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t24) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t45) + (g(1) * t51 + g(2) * t36) * t34 + (g(1) * (-pkin(1) - t36) + g(2) * t51) * t32) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t24) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t45) + (g(1) * t15 + g(2) * t35) * t34 + (g(1) * (-pkin(1) - t35) + g(2) * t15) * t32), -m(3) * (g(3) * t42 + t58 * (-rSges(3,1) * t31 - rSges(3,2) * t33)) - m(4) * ((g(1) * t43 + g(3) * t39 + t44 * t53) * t33 + (g(3) * t44 - t58 * t39) * t31) - m(5) * ((g(3) * t38 + t58 * t47) * t33 + (g(3) * t47 - t58 * t38) * t31) - m(6) * ((g(3) * t37 + t58 * t46) * t33 + (g(3) * t46 - t58 * t37) * t31), (-m(4) - m(5) - m(6)) * (-g(3) * t33 + t58 * t31), -m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2))) - m(6) * (g(1) * (t11 * pkin(4) + t56) + g(2) * (-t9 * pkin(4) + t57)) + (-m(5) * (-rSges(5,1) * t19 - rSges(5,2) * t20) - m(6) * (t40 - t55)) * t52, -m(6) * (g(1) * t56 + g(2) * t57 + t40 * t52)];
taug = t1(:);
