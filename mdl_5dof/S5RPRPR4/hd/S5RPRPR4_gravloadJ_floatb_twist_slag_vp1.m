% Calculate Gravitation load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:15
% EndTime: 2020-01-03 11:38:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (200->62), mult. (158->76), div. (0->0), fcn. (121->10), ass. (0->37)
t21 = qJ(1) + pkin(8);
t14 = cos(t21);
t52 = g(3) * t14;
t20 = qJ(3) + pkin(9);
t13 = cos(t20);
t25 = cos(qJ(3));
t17 = t25 * pkin(3);
t15 = qJ(5) + t20;
t8 = sin(t15);
t9 = cos(t15);
t36 = t9 * rSges(6,1) - rSges(6,2) * t8;
t51 = pkin(4) * t13 + t17 + t36;
t50 = rSges(6,1) * t8 + rSges(6,2) * t9;
t11 = sin(t20);
t49 = rSges(5,1) * t13 - rSges(5,2) * t11 + t17;
t23 = sin(qJ(3));
t43 = pkin(3) * t23;
t48 = m(4) * (rSges(4,1) * t23 + rSges(4,2) * t25) + m(5) * (rSges(5,1) * t11 + rSges(5,2) * t13 + t43);
t47 = -m(5) - m(6);
t12 = sin(t21);
t42 = g(2) * t12;
t41 = rSges(4,3) + pkin(6);
t22 = -qJ(4) - pkin(6);
t39 = rSges(5,3) - t22;
t38 = rSges(6,3) + pkin(7) - t22;
t37 = t50 * t52;
t24 = sin(qJ(1));
t16 = t24 * pkin(1);
t26 = cos(qJ(1));
t18 = t26 * pkin(1);
t34 = g(2) * t18 + g(3) * t16;
t33 = rSges(4,1) * t25 - rSges(4,2) * t23;
t30 = pkin(2) + t51;
t29 = pkin(2) + t33;
t28 = pkin(2) + t49;
t5 = -pkin(4) * t11 - t43;
t1 = [-m(2) * (g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2)) + g(3) * (t24 * rSges(2,1) + rSges(2,2) * t26)) - m(3) * (g(2) * (rSges(3,1) * t14 - rSges(3,2) * t12 + t18) + g(3) * (rSges(3,1) * t12 + rSges(3,2) * t14 + t16)) - m(4) * ((g(2) * t29 - g(3) * t41) * t14 + (g(2) * t41 + g(3) * t29) * t12 + t34) - m(5) * ((g(2) * t28 - g(3) * t39) * t14 + (g(2) * t39 + g(3) * t28) * t12 + t34) - m(6) * ((g(2) * t30 - g(3) * t38) * t14 + (g(2) * t38 + g(3) * t30) * t12 + t34), (-m(3) - m(4) + t47) * g(1), -m(6) * t37 + (-m(4) * t33 - m(5) * t49 - m(6) * t51) * g(1) + (m(6) * t5 - t48) * t52 + (-m(6) * (-t50 + t5) + t48) * t42, t47 * (-g(2) * t14 - g(3) * t12), -m(6) * (g(1) * t36 - t42 * t50 + t37)];
taug = t1(:);
