% Calculate Gravitation load on the joints for
% S5RPRRR2
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:16
% EndTime: 2019-12-05 18:11:18
% DurationCPUTime: 0.28s
% Computational Cost: add. (254->62), mult. (201->76), div. (0->0), fcn. (158->10), ass. (0->35)
t21 = pkin(9) + qJ(3);
t17 = qJ(4) + t21;
t12 = cos(t17);
t14 = qJ(5) + t17;
t8 = sin(t14);
t9 = cos(t14);
t38 = t9 * rSges(6,1) - t8 * rSges(6,2);
t36 = pkin(4) * t12 + t38;
t11 = sin(t17);
t37 = t12 * rSges(5,1) - t11 * rSges(5,2);
t35 = -rSges(6,1) * t8 - rSges(6,2) * t9;
t50 = -pkin(4) * t11 + t35;
t25 = sin(qJ(1));
t26 = cos(qJ(1));
t49 = g(1) * t26 + g(2) * t25;
t15 = sin(t21);
t48 = pkin(3) * t15;
t23 = cos(pkin(9));
t13 = t23 * pkin(2) + pkin(1);
t24 = -pkin(6) - qJ(2);
t42 = rSges(4,3) - t24;
t20 = -pkin(7) + t24;
t41 = rSges(5,3) - t20;
t40 = rSges(6,3) + pkin(8) - t20;
t39 = rSges(3,3) + qJ(2);
t16 = cos(t21);
t10 = pkin(3) * t16;
t3 = t10 + t13;
t34 = t16 * rSges(4,1) - t15 * rSges(4,2);
t33 = -rSges(5,1) * t11 - rSges(5,2) * t12;
t32 = t3 + t36;
t31 = rSges(3,1) * t23 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t30 = t3 + t37;
t29 = t13 + t34;
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - t26 * rSges(2,2)) + g(2) * (t26 * rSges(2,1) - t25 * rSges(2,2))) - m(3) * ((g(1) * t39 + g(2) * t31) * t26 + (-g(1) * t31 + g(2) * t39) * t25) - m(4) * ((g(1) * t42 + g(2) * t29) * t26 + (-g(1) * t29 + g(2) * t42) * t25) - m(5) * ((g(1) * t41 + g(2) * t30) * t26 + (-g(1) * t30 + g(2) * t41) * t25) - m(6) * ((g(1) * t40 + g(2) * t32) * t26 + (-g(1) * t32 + g(2) * t40) * t25), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t25 - g(2) * t26), (-m(4) * t34 - m(5) * (t10 + t37) - m(6) * (t10 + t36)) * g(3) + t49 * (-m(4) * (-rSges(4,1) * t15 - rSges(4,2) * t16) - m(5) * (t33 - t48) - m(6) * (-t48 + t50)), (-m(5) * t37 - m(6) * t36) * g(3) + t49 * (-m(5) * t33 - m(6) * t50), -m(6) * (g(3) * t38 + t49 * t35)];
taug = t1(:);
