% Calculate Gravitation load on the joints for
% S5RPPPR2
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:58
% EndTime: 2019-12-05 17:31:00
% DurationCPUTime: 0.46s
% Computational Cost: add. (167->79), mult. (365->117), div. (0->0), fcn. (411->10), ass. (0->42)
t21 = sin(pkin(8));
t25 = cos(pkin(7));
t29 = cos(qJ(1));
t37 = t29 * t25;
t24 = cos(pkin(8));
t27 = sin(qJ(1));
t39 = t27 * t24;
t15 = t21 * t37 - t39;
t26 = sin(qJ(5));
t28 = cos(qJ(5));
t40 = t27 * t21;
t16 = t24 * t37 + t40;
t20 = sin(pkin(9));
t23 = cos(pkin(9));
t22 = sin(pkin(7));
t38 = t29 * t22;
t6 = t16 * t23 + t20 * t38;
t50 = -t15 * t28 + t6 * t26;
t49 = -t15 * t26 - t6 * t28;
t48 = rSges(6,3) + pkin(6);
t47 = -m(5) - m(6);
t46 = g(2) * t29;
t45 = g(2) * qJ(2);
t42 = t21 * t22;
t41 = t22 * t27;
t36 = -m(4) + t47;
t35 = -pkin(2) * t25 - pkin(1);
t34 = -t16 * pkin(3) - t15 * qJ(4);
t13 = t29 * t24 + t25 * t40;
t14 = t29 * t21 - t25 * t39;
t19 = t29 * qJ(2);
t33 = t14 * pkin(3) - t13 * qJ(4) + t19;
t32 = -rSges(3,1) * t25 + rSges(3,2) * t22 - pkin(1);
t31 = -qJ(3) * t22 + t35;
t30 = (g(3) * t31 - t45) * t27 + t31 * t46;
t12 = t22 * t24 * t23 - t25 * t20;
t7 = -t16 * t20 + t23 * t38;
t5 = t14 * t23 - t20 * t41;
t4 = t14 * t20 + t23 * t41;
t2 = -t13 * t26 + t5 * t28;
t1 = -t13 * t28 - t5 * t26;
t3 = [-m(2) * (g(2) * (-t29 * rSges(2,1) + t27 * rSges(2,2)) + g(3) * (-t27 * rSges(2,1) - t29 * rSges(2,2))) - m(3) * (g(3) * t19 + (g(3) * rSges(3,3) + g(2) * t32) * t29 + (g(2) * (-rSges(3,3) - qJ(2)) + g(3) * t32) * t27) - m(4) * (g(2) * (-t16 * rSges(4,1) + t15 * rSges(4,2)) + g(3) * (t14 * rSges(4,1) + t13 * rSges(4,2) + t19) + ((-rSges(4,3) - qJ(3)) * t22 + t35) * t46 + (-t45 + g(3) * (-rSges(4,3) * t22 + t31)) * t27) - m(5) * (g(2) * (-rSges(5,1) * t6 - t7 * rSges(5,2) - t15 * rSges(5,3) + t34) + g(3) * (t5 * rSges(5,1) - t4 * rSges(5,2) - t13 * rSges(5,3) + t33) + t30) - m(6) * (g(2) * (t49 * rSges(6,1) + t50 * rSges(6,2) - t6 * pkin(4) + t48 * t7 + t34) + g(3) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t5 * pkin(4) + t4 * t48 + t33) + t30), (-m(3) + t36) * (g(3) * t27 + t46), t36 * (-g(1) * t25 + (-g(2) * t27 + g(3) * t29) * t22), t47 * (g(1) * t42 - g(2) * t13 + g(3) * t15), -m(6) * (g(1) * ((-t12 * t26 + t28 * t42) * rSges(6,1) + (-t12 * t28 - t26 * t42) * rSges(6,2)) + g(2) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (-t50 * rSges(6,1) + t49 * rSges(6,2)))];
taug = t3(:);
