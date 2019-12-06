% Calculate Gravitation load on the joints for
% S5RRPRR4
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:49
% EndTime: 2019-12-05 18:31:50
% DurationCPUTime: 0.27s
% Computational Cost: add. (272->70), mult. (177->85), div. (0->0), fcn. (135->10), ass. (0->44)
t57 = rSges(5,3) + pkin(7);
t56 = rSges(6,3) + pkin(8) + pkin(7);
t27 = cos(qJ(4));
t48 = rSges(5,1) * t27;
t55 = -pkin(3) - t48;
t24 = qJ(1) + qJ(2);
t17 = pkin(9) + t24;
t14 = sin(t17);
t23 = qJ(4) + qJ(5);
t18 = sin(t23);
t44 = t14 * t18;
t20 = cos(t23);
t45 = rSges(6,2) * t20;
t54 = rSges(6,1) * t44 + t14 * t45;
t28 = cos(qJ(1));
t53 = pkin(1) * t28;
t19 = sin(t24);
t52 = pkin(2) * t19;
t21 = cos(t24);
t51 = pkin(2) * t21;
t15 = cos(t17);
t50 = g(3) * t15;
t26 = sin(qJ(1));
t49 = t26 * pkin(1);
t25 = sin(qJ(4));
t47 = rSges(5,2) * t25;
t46 = rSges(6,2) * t18;
t43 = t14 * t25;
t13 = t20 * rSges(6,1);
t42 = -rSges(3,1) * t21 + t19 * rSges(3,2);
t22 = t27 * pkin(4);
t41 = -t22 - pkin(3) - t13;
t40 = t13 - t46;
t39 = -rSges(3,1) * t19 - rSges(3,2) * t21;
t38 = rSges(5,1) * t25 + rSges(5,2) * t27;
t37 = -rSges(6,1) * t18 - t45;
t36 = rSges(5,2) * t43 + t57 * t15 - t52;
t35 = -rSges(4,1) * t15 + t14 * rSges(4,2) - t51;
t34 = -rSges(4,1) * t14 - rSges(4,2) * t15 - t52;
t33 = -t51 + (t47 + t55) * t15;
t32 = (-g(2) * t57 + g(3) * t55) * t14;
t31 = -t51 + (t41 + t46) * t15 - t56 * t14;
t30 = rSges(6,2) * t44 + t41 * t14 + t56 * t15 - t52;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t28 + rSges(2,2) * t26) + g(3) * (-rSges(2,1) * t26 - rSges(2,2) * t28)) - m(3) * (g(2) * (t42 - t53) + g(3) * (t39 - t49)) - m(4) * (g(2) * (t35 - t53) + g(3) * (t34 - t49)) - m(5) * (g(2) * (t33 - t53) + g(3) * (t36 - t49) + t32) - m(6) * (g(2) * (t31 - t53) + g(3) * (t30 - t49)), -m(3) * (g(2) * t42 + g(3) * t39) - m(4) * (g(2) * t35 + g(3) * t34) - m(5) * (g(2) * t33 + g(3) * t36 + t32) - m(6) * (g(2) * t31 + g(3) * t30), (-m(4) - m(5) - m(6)) * g(1), -m(5) * (g(1) * (-t47 + t48) + g(2) * t38 * t14) - m(6) * (g(1) * (t22 + t40) + g(2) * (pkin(4) * t43 + t54)) + (m(5) * t38 - m(6) * (-pkin(4) * t25 + t37)) * t50, -m(6) * (g(1) * t40 + g(2) * t54 + t37 * t50)];
taug = t1(:);
