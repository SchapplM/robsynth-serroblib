% Calculate Gravitation load on the joints for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:43
% EndTime: 2019-12-31 19:09:45
% DurationCPUTime: 0.38s
% Computational Cost: add. (274->89), mult. (306->125), div. (0->0), fcn. (283->10), ass. (0->48)
t58 = rSges(5,3) + pkin(7);
t57 = rSges(6,3) + pkin(8) + pkin(7);
t28 = cos(qJ(4));
t16 = t28 * pkin(4) + pkin(3);
t22 = qJ(4) + qJ(5);
t19 = sin(t22);
t20 = cos(t22);
t26 = sin(qJ(4));
t56 = m(5) * (rSges(5,1) * t28 - rSges(5,2) * t26 + pkin(3)) + m(6) * (rSges(6,1) * t20 - rSges(6,2) * t19 + t16) + m(4) * rSges(4,1);
t21 = pkin(9) + qJ(3);
t18 = cos(t21);
t29 = cos(qJ(1));
t45 = t29 * t20;
t27 = sin(qJ(1));
t50 = t27 * t19;
t5 = t18 * t50 + t45;
t46 = t29 * t19;
t49 = t27 * t20;
t6 = -t18 * t49 + t46;
t55 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = -t18 * t46 + t49;
t8 = t18 * t45 + t50;
t54 = t7 * rSges(6,1) - t8 * rSges(6,2);
t52 = pkin(4) * t26;
t17 = sin(t21);
t51 = g(3) * t17;
t48 = t27 * t26;
t47 = t27 * t28;
t44 = t29 * t26;
t43 = t29 * t28;
t25 = -pkin(6) - qJ(2);
t42 = rSges(4,3) - t25;
t41 = rSges(3,3) + qJ(2);
t40 = -t25 + t52;
t39 = t18 * rSges(4,1) - t17 * rSges(4,2);
t38 = -rSges(6,1) * t19 - rSges(6,2) * t20;
t24 = cos(pkin(9));
t37 = rSges(3,1) * t24 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t11 = -t18 * t44 + t47;
t9 = t18 * t48 + t43;
t34 = t18 * pkin(3) + t58 * t17;
t33 = t18 * t16 + t57 * t17;
t32 = m(4) * rSges(4,2) - m(5) * t58 - m(6) * t57;
t15 = t24 * pkin(2) + pkin(1);
t13 = t29 * t15;
t12 = t18 * t43 + t48;
t10 = -t18 * t47 + t44;
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t27 * rSges(2,2))) - m(3) * ((g(1) * t41 + g(2) * t37) * t29 + (-g(1) * t37 + g(2) * t41) * t27) - m(4) * (g(2) * t13 + (g(1) * t42 + g(2) * t39) * t29 + (g(1) * (-t15 - t39) + g(2) * t42) * t27) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2)) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t13) + (-g(1) * t25 + g(2) * t34) * t29 + (g(1) * (-t15 - t34) - g(2) * t25) * t27) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2)) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t13) + (g(1) * t40 + g(2) * t33) * t29 + (g(1) * (-t15 - t33) + g(2) * t40) * t27), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t27 - g(2) * t29), (t32 * t17 - t18 * t56) * g(3) + (g(1) * t29 + g(2) * t27) * (t17 * t56 + t32 * t18), -m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2))) - m(6) * (g(1) * (t11 * pkin(4) + t54) + g(2) * (-t9 * pkin(4) + t55)) + (-m(5) * (-rSges(5,1) * t26 - rSges(5,2) * t28) - m(6) * (t38 - t52)) * t51, -m(6) * (g(1) * t54 + g(2) * t55 + t38 * t51)];
taug = t1(:);
