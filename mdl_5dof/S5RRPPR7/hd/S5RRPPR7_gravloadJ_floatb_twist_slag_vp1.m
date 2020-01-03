% Calculate Gravitation load on the joints for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:34:58
% DurationCPUTime: 0.40s
% Computational Cost: add. (213->88), mult. (265->119), div. (0->0), fcn. (233->8), ass. (0->48)
t18 = qJ(2) + pkin(8);
t15 = sin(t18);
t12 = t15 * qJ(4);
t16 = cos(t18);
t13 = t16 * pkin(3);
t55 = t12 + t13;
t25 = cos(qJ(1));
t22 = sin(qJ(1));
t50 = g(2) * t22;
t54 = g(1) * t25 + t50;
t53 = -m(5) - m(6);
t21 = sin(qJ(2));
t52 = pkin(2) * t21;
t49 = g(3) * t16;
t48 = rSges(3,3) + pkin(6);
t47 = rSges(6,3) + pkin(7);
t19 = -qJ(3) - pkin(6);
t46 = pkin(4) - t19;
t20 = sin(qJ(5));
t45 = t22 * t20;
t23 = cos(qJ(5));
t44 = t22 * t23;
t43 = t25 * t20;
t42 = t25 * t23;
t41 = rSges(5,1) - t19;
t40 = rSges(4,3) - t19;
t39 = qJ(4) * t16;
t38 = -pkin(3) - t47;
t24 = cos(qJ(2));
t17 = t24 * pkin(2);
t14 = t17 + pkin(1);
t11 = t25 * t14;
t37 = t55 * t25 + t11;
t36 = t17 + t55;
t35 = -t14 - t12;
t34 = g(1) * t38;
t33 = t24 * rSges(3,1) - t21 * rSges(3,2);
t31 = t16 * rSges(4,1) - t15 * rSges(4,2);
t30 = rSges(6,1) * t20 + rSges(6,2) * t23;
t29 = t16 * rSges(5,2) - t15 * rSges(5,3);
t28 = pkin(1) + t33;
t8 = t25 * t39;
t6 = t22 * t39;
t5 = -t15 * t45 + t42;
t4 = t15 * t44 + t43;
t3 = t15 * t43 + t44;
t2 = t15 * t42 - t45;
t1 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * ((g(1) * t48 + g(2) * t28) * t25 + (-g(1) * t28 + g(2) * t48) * t22) - m(4) * (g(2) * t11 + (g(1) * t40 + g(2) * t31) * t25 + (g(1) * (-t14 - t31) + g(2) * t40) * t22) - m(5) * (g(2) * t37 + (g(1) * t41 - g(2) * t29) * t25 + (g(1) * (t29 + t35 - t13) + g(2) * t41) * t22) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t37) + (g(2) * t47 * t16 + g(1) * t46) * t25 + (g(1) * t35 + g(2) * t46 + t16 * t34) * t22), -m(3) * (g(3) * t33 + t54 * (-rSges(3,1) * t21 - rSges(3,2) * t24)) - m(4) * (g(3) * (t17 + t31) + t54 * (-rSges(4,1) * t15 - rSges(4,2) * t16 - t52)) - m(5) * (g(1) * t8 + g(2) * t6 + g(3) * (-t29 + t36) + t54 * (rSges(5,3) * t16 - t52 + (rSges(5,2) - pkin(3)) * t15)) - m(6) * (g(1) * (-t25 * t52 + t8) + g(2) * (-t22 * t52 + t6) + g(3) * t36 + (g(3) * t47 + t54 * t30) * t16 + (g(3) * t30 + t25 * t34 + t38 * t50) * t15), (-m(4) + t53) * (g(1) * t22 - g(2) * t25), t53 * (t54 * t15 - t49), -m(6) * (g(1) * (t2 * rSges(6,1) - t3 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t5 * rSges(6,2)) + (-rSges(6,1) * t23 + rSges(6,2) * t20) * t49)];
taug = t1(:);
