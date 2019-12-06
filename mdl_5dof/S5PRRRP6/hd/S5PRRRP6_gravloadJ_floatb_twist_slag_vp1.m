% Calculate Gravitation load on the joints for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:47
% EndTime: 2019-12-05 16:50:49
% DurationCPUTime: 0.43s
% Computational Cost: add. (238->76), mult. (341->110), div. (0->0), fcn. (329->8), ass. (0->43)
t57 = rSges(6,1) + pkin(4);
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t29 = cos(qJ(3));
t27 = sin(qJ(3));
t30 = cos(qJ(2));
t50 = t27 * t30;
t67 = t25 * t29 - t26 * t50;
t66 = rSges(6,3) + qJ(5);
t65 = g(1) * t26 + g(2) * t25;
t24 = qJ(3) + qJ(4);
t22 = sin(t24);
t23 = cos(t24);
t64 = t66 * t22 + t57 * t23;
t52 = t25 * t30;
t11 = t22 * t52 + t26 * t23;
t12 = -t26 * t22 + t23 * t52;
t63 = -t11 * rSges(5,1) - t12 * rSges(5,2);
t51 = t26 * t30;
t13 = t22 * t51 - t25 * t23;
t14 = t25 * t22 + t23 * t51;
t62 = -t13 * rSges(5,1) - t14 * rSges(5,2);
t61 = pkin(3) * t27;
t28 = sin(qJ(2));
t58 = g(3) * t28;
t56 = rSges(4,3) + pkin(6);
t55 = t22 * t28;
t49 = t29 * t30;
t31 = -pkin(7) - pkin(6);
t48 = rSges(6,2) - t31;
t47 = rSges(5,3) - t31;
t46 = t66 * t23 * t28;
t42 = rSges(5,1) * t23 - rSges(5,2) * t22;
t41 = -rSges(5,1) * t22 - rSges(5,2) * t23;
t40 = -t57 * t11 + t66 * t12;
t39 = t67 * pkin(3);
t38 = -t57 * t13 + t66 * t14;
t37 = t29 * rSges(4,1) - t27 * rSges(4,2) + pkin(2);
t36 = -t25 * t50 - t26 * t29;
t35 = t36 * pkin(3);
t21 = t29 * pkin(3) + pkin(2);
t16 = t30 * t21;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(3) * (t30 * rSges(3,1) - t28 * rSges(3,2)) + t65 * (-rSges(3,1) * t28 - rSges(3,2) * t30)) - m(4) * (g(3) * (t56 * t28 + t37 * t30) + t65 * (-t37 * t28 + t56 * t30)) - m(5) * (g(3) * (t47 * t28 + t42 * t30 + t16) + t65 * (t47 * t30 + (-t21 - t42) * t28)) - m(6) * (g(3) * t16 + (g(3) * t64 + t65 * t48) * t30 + (g(3) * t48 + t65 * (-t21 - t64)) * t28), -m(4) * (g(1) * (t67 * rSges(4,1) + (-t25 * t27 - t26 * t49) * rSges(4,2)) + g(2) * (t36 * rSges(4,1) + (-t25 * t49 + t26 * t27) * rSges(4,2))) - m(5) * (g(1) * (t39 + t62) + g(2) * (t35 + t63)) - m(6) * (g(1) * (t38 + t39) + g(2) * (t35 + t40) + g(3) * t46) + (-m(4) * (-rSges(4,1) * t27 - rSges(4,2) * t29) - m(5) * (t41 - t61) - m(6) * (-t57 * t22 - t61)) * t58, -m(5) * (g(1) * t62 + g(2) * t63 + t41 * t58) - m(6) * (g(1) * t38 + g(2) * t40 + g(3) * (-t57 * t55 + t46)), -m(6) * (g(1) * t13 + g(2) * t11 + g(3) * t55)];
taug = t1(:);
