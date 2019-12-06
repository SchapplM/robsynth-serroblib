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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:53:09
% EndTime: 2019-12-05 17:53:10
% DurationCPUTime: 0.30s
% Computational Cost: add. (200->62), mult. (158->77), div. (0->0), fcn. (121->10), ass. (0->35)
t19 = qJ(1) + pkin(8);
t12 = sin(t19);
t50 = g(2) * t12;
t18 = qJ(3) + pkin(9);
t13 = cos(t18);
t23 = cos(qJ(3));
t16 = t23 * pkin(3);
t15 = qJ(5) + t18;
t8 = sin(t15);
t9 = cos(t15);
t34 = t9 * rSges(6,1) - t8 * rSges(6,2);
t49 = pkin(4) * t13 + t16 + t34;
t48 = rSges(6,1) * t8 + rSges(6,2) * t9;
t11 = sin(t18);
t47 = t13 * rSges(5,1) - t11 * rSges(5,2) + t16;
t21 = sin(qJ(3));
t41 = t21 * pkin(3);
t46 = m(4) * (rSges(4,1) * t21 + rSges(4,2) * t23) + m(5) * (rSges(5,1) * t11 + rSges(5,2) * t13 + t41);
t45 = -m(5) - m(6);
t14 = cos(t19);
t42 = g(3) * t14;
t39 = rSges(4,3) + pkin(6);
t20 = -qJ(4) - pkin(6);
t37 = rSges(5,3) - t20;
t36 = rSges(6,3) + pkin(7) - t20;
t35 = t48 * t50;
t32 = t23 * rSges(4,1) - t21 * rSges(4,2);
t29 = -pkin(2) - t49;
t28 = -pkin(2) - t32;
t27 = -pkin(2) - t47;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t26 = (-g(2) * t24 - g(3) * t22) * pkin(1);
t5 = -pkin(4) * t11 - t41;
t1 = [-m(2) * (g(2) * (-t24 * rSges(2,1) + t22 * rSges(2,2)) + g(3) * (-t22 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(2) * (-t14 * rSges(3,1) + t12 * rSges(3,2) - t24 * pkin(1)) + g(3) * (-t12 * rSges(3,1) - t14 * rSges(3,2) - t22 * pkin(1))) - m(4) * (t26 + (g(2) * t28 + g(3) * t39) * t14 + (-g(2) * t39 + g(3) * t28) * t12) - m(5) * (t26 + (g(2) * t27 + g(3) * t37) * t14 + (-g(2) * t37 + g(3) * t27) * t12) - m(6) * (t26 + (g(2) * t29 + g(3) * t36) * t14 + (-g(2) * t36 + g(3) * t29) * t12), (-m(3) - m(4) + t45) * g(1), -m(6) * t35 + (-m(4) * t32 - m(5) * t47 - m(6) * t49) * g(1) + (m(6) * t5 - t46) * t50 + (-m(6) * (-t48 + t5) + t46) * t42, t45 * (g(2) * t14 + g(3) * t12), -m(6) * (g(1) * t34 - t42 * t48 + t35)];
taug = t1(:);
