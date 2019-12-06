% Calculate Gravitation load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:48
% EndTime: 2019-12-05 18:37:49
% DurationCPUTime: 0.37s
% Computational Cost: add. (289->71), mult. (242->85), div. (0->0), fcn. (193->10), ass. (0->42)
t24 = qJ(2) + qJ(3);
t18 = pkin(9) + t24;
t14 = cos(t18);
t16 = qJ(5) + t18;
t10 = sin(t16);
t11 = cos(t16);
t46 = t11 * rSges(6,1) - t10 * rSges(6,2);
t62 = pkin(4) * t14 + t46;
t13 = sin(t18);
t39 = -rSges(6,1) * t10 - rSges(6,2) * t11;
t19 = sin(t24);
t57 = pkin(3) * t19;
t61 = -pkin(4) * t13 + t39 - t57;
t60 = t14 * rSges(5,1) - t13 * rSges(5,2);
t20 = cos(t24);
t45 = t20 * rSges(4,1) - t19 * rSges(4,2);
t59 = -rSges(5,1) * t13 - rSges(5,2) * t14 - t57;
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t58 = g(1) * t28 + g(2) * t26;
t29 = -pkin(7) - pkin(6);
t25 = sin(qJ(2));
t54 = t25 * pkin(2);
t53 = rSges(3,3) + pkin(6);
t27 = cos(qJ(2));
t22 = t27 * pkin(2);
t17 = t22 + pkin(1);
t49 = rSges(4,3) - t29;
t23 = -qJ(4) + t29;
t48 = rSges(5,3) - t23;
t47 = rSges(6,3) + pkin(8) - t23;
t15 = pkin(3) * t20;
t4 = t15 + t17;
t44 = t15 + t60;
t43 = t27 * rSges(3,1) - t25 * rSges(3,2);
t41 = -rSges(4,1) * t19 - rSges(4,2) * t20;
t38 = t15 + t62;
t37 = pkin(1) + t43;
t36 = t4 + t60;
t34 = t4 + t62;
t32 = t17 + t45;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * ((g(1) * t53 + g(2) * t37) * t28 + (-g(1) * t37 + g(2) * t53) * t26) - m(4) * ((g(1) * t49 + g(2) * t32) * t28 + (-g(1) * t32 + g(2) * t49) * t26) - m(5) * ((g(1) * t48 + g(2) * t36) * t28 + (-g(1) * t36 + g(2) * t48) * t26) - m(6) * ((g(1) * t47 + g(2) * t34) * t28 + (-g(1) * t34 + g(2) * t47) * t26), -m(3) * (g(3) * t43 + t58 * (-rSges(3,1) * t25 - rSges(3,2) * t27)) - m(4) * (g(3) * (t22 + t45) + t58 * (t41 - t54)) - m(5) * (g(3) * (t22 + t44) + t58 * (-t54 + t59)) - m(6) * (g(3) * (t22 + t38) + t58 * (-t54 + t61)), (-m(4) * t45 - m(5) * t44 - m(6) * t38) * g(3) + t58 * (-m(4) * t41 - m(5) * t59 - m(6) * t61), (-m(5) - m(6)) * (g(1) * t26 - g(2) * t28), -m(6) * (g(3) * t46 + t58 * t39)];
taug = t1(:);
