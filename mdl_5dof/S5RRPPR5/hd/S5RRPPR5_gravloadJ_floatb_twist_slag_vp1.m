% Calculate Gravitation load on the joints for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:28:55
% DurationCPUTime: 0.40s
% Computational Cost: add. (229->85), mult. (278->112), div. (0->0), fcn. (254->8), ass. (0->46)
t21 = qJ(2) + pkin(8);
t18 = sin(t21);
t15 = t18 * qJ(4);
t19 = cos(t21);
t58 = -t19 * pkin(3) - t15;
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t57 = g(1) * t28 + g(2) * t25;
t56 = t57 * t18;
t55 = -m(5) - m(6);
t24 = sin(qJ(2));
t53 = pkin(2) * t24;
t50 = t19 * pkin(4);
t49 = rSges(3,3) + pkin(6);
t26 = cos(qJ(5));
t48 = t18 * t26;
t47 = t19 * t28;
t22 = -qJ(3) - pkin(6);
t46 = rSges(5,2) - t22;
t45 = rSges(4,3) - t22;
t44 = qJ(4) * t19;
t43 = -rSges(6,3) - pkin(7) - t22;
t27 = cos(qJ(2));
t20 = t27 * pkin(2);
t17 = t20 + pkin(1);
t14 = t28 * t17;
t42 = pkin(3) * t47 + t28 * t15 + t14;
t41 = t20 - t58;
t23 = sin(qJ(5));
t7 = -t19 * t23 + t48;
t2 = t7 * t25;
t33 = t18 * t23 + t19 * t26;
t3 = t33 * t25;
t40 = t2 * rSges(6,1) - t3 * rSges(6,2);
t4 = t23 * t47 - t28 * t48;
t5 = t33 * t28;
t39 = -t4 * rSges(6,1) - t5 * rSges(6,2);
t38 = -rSges(6,1) * t33 - t7 * rSges(6,2);
t37 = t27 * rSges(3,1) - t24 * rSges(3,2);
t35 = t19 * rSges(4,1) - t18 * rSges(4,2);
t34 = t19 * rSges(5,1) + t18 * rSges(5,3);
t32 = pkin(1) + t37;
t31 = -t17 + t58;
t11 = t28 * t44;
t9 = t25 * t44;
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t25 * rSges(2,2))) - m(3) * ((g(1) * t49 + g(2) * t32) * t28 + (-g(1) * t32 + g(2) * t49) * t25) - m(4) * (g(2) * t14 + (g(1) * t45 + g(2) * t35) * t28 + (g(1) * (-t17 - t35) + g(2) * t45) * t25) - m(5) * (g(2) * t42 + (g(1) * t46 + g(2) * t34) * t28 + (g(1) * (t31 - t34) + g(2) * t46) * t25) - m(6) * (g(1) * (-t3 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t42) + (g(1) * t43 + g(2) * t50) * t28 + (g(1) * (t31 - t50) + g(2) * t43) * t25), -m(3) * (g(3) * t37 + t57 * (-rSges(3,1) * t24 - rSges(3,2) * t27)) - m(4) * (g(3) * (t20 + t35) + t57 * (-rSges(4,1) * t18 - rSges(4,2) * t19 - t53)) - m(5) * (g(1) * t11 + g(2) * t9 + g(3) * (t34 + t41) + t57 * (rSges(5,3) * t19 - t53 + (-rSges(5,1) - pkin(3)) * t18)) - m(6) * (g(1) * (-t28 * t53 + t11 - t39) + g(2) * (-t25 * t53 - t40 + t9) + g(3) * (-t38 + t41 + t50) + (-pkin(3) - pkin(4)) * t56), (-m(4) + t55) * (g(1) * t25 - g(2) * t28), t55 * (-g(3) * t19 + t56), -m(6) * (g(1) * t39 + g(2) * t40 + g(3) * t38)];
taug = t1(:);
