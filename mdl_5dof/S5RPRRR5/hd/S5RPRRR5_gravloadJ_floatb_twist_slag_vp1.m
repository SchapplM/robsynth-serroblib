% Calculate Gravitation load on the joints for
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:06
% EndTime: 2019-12-05 18:16:07
% DurationCPUTime: 0.23s
% Computational Cost: add. (253->65), mult. (165->82), div. (0->0), fcn. (125->10), ass. (0->42)
t54 = rSges(5,3) + pkin(7);
t53 = rSges(6,3) + pkin(8) + pkin(7);
t26 = cos(qJ(4));
t47 = rSges(5,1) * t26;
t52 = -pkin(3) - t47;
t22 = qJ(1) + pkin(9);
t18 = qJ(3) + t22;
t13 = sin(t18);
t23 = qJ(4) + qJ(5);
t19 = sin(t23);
t43 = t13 * t19;
t20 = cos(t23);
t44 = rSges(6,2) * t20;
t51 = rSges(6,1) * t43 + t13 * t44;
t27 = cos(qJ(1));
t50 = pkin(1) * t27;
t14 = cos(t18);
t49 = g(3) * t14;
t25 = sin(qJ(1));
t48 = t25 * pkin(1);
t24 = sin(qJ(4));
t46 = rSges(5,2) * t24;
t45 = rSges(6,2) * t19;
t42 = t13 * t24;
t12 = t20 * rSges(6,1);
t41 = rSges(5,2) * t42 + t54 * t14;
t40 = -rSges(4,1) * t14 + t13 * rSges(4,2);
t21 = t26 * pkin(4);
t39 = -t21 - pkin(3) - t12;
t38 = t12 - t45;
t16 = sin(t22);
t37 = -pkin(2) * t16 - t48;
t17 = cos(t22);
t36 = -pkin(2) * t17 - t50;
t35 = -rSges(4,1) * t13 - rSges(4,2) * t14;
t34 = rSges(5,1) * t24 + rSges(5,2) * t26;
t33 = -rSges(6,1) * t19 - t44;
t32 = (t46 + t52) * t14;
t31 = (t39 + t45) * t14 - t53 * t13;
t30 = rSges(6,2) * t43 + t39 * t13 + t53 * t14;
t29 = (-g(2) * t54 + g(3) * t52) * t13;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t27 + rSges(2,2) * t25) + g(3) * (-rSges(2,1) * t25 - rSges(2,2) * t27)) - m(3) * (g(2) * (-rSges(3,1) * t17 + rSges(3,2) * t16 - t50) + g(3) * (-rSges(3,1) * t16 - rSges(3,2) * t17 - t48)) - m(4) * (g(2) * (t36 + t40) + g(3) * (t35 + t37)) - m(5) * (g(2) * (t32 + t36) + g(3) * (t37 + t41) + t29) - m(6) * (g(2) * (t31 + t36) + g(3) * (t30 + t37)), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t40 + g(3) * t35) - m(5) * (g(2) * t32 + g(3) * t41 + t29) - m(6) * (g(2) * t31 + g(3) * t30), -m(5) * (g(1) * (-t46 + t47) + g(2) * t34 * t13) - m(6) * (g(1) * (t21 + t38) + g(2) * (pkin(4) * t42 + t51)) + (m(5) * t34 - m(6) * (-pkin(4) * t24 + t33)) * t49, -m(6) * (g(1) * t38 + g(2) * t51 + t33 * t49)];
taug = t1(:);
