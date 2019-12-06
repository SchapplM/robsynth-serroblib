% Calculate Gravitation load on the joints for
% S5RPRPR3
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:15
% DurationCPUTime: 0.26s
% Computational Cost: add. (262->65), mult. (180->84), div. (0->0), fcn. (158->10), ass. (0->38)
t20 = qJ(1) + pkin(8);
t19 = qJ(3) + t20;
t15 = sin(t19);
t16 = cos(t19);
t25 = cos(qJ(5));
t22 = cos(pkin(9));
t23 = sin(qJ(5));
t37 = t22 * t23;
t5 = t15 * t37 + t16 * t25;
t36 = t22 * t25;
t6 = t15 * t36 - t16 * t23;
t48 = -t6 * rSges(6,1) + t5 * rSges(6,2);
t47 = -rSges(5,1) * t22 - pkin(3);
t21 = sin(pkin(9));
t38 = rSges(5,2) * t21;
t46 = t16 * rSges(5,3) + t15 * t38;
t45 = -pkin(4) * t22 - pkin(3) + (-rSges(6,3) - pkin(7)) * t21;
t44 = -m(5) - m(6);
t7 = -t15 * t25 + t16 * t37;
t8 = -t15 * t23 - t16 * t36;
t43 = t8 * rSges(6,1) + t7 * rSges(6,2);
t24 = sin(qJ(1));
t42 = pkin(1) * t24;
t26 = cos(qJ(1));
t41 = pkin(1) * t26;
t40 = g(2) * t16;
t34 = -t16 * rSges(4,1) + t15 * rSges(4,2);
t17 = sin(t20);
t33 = -pkin(2) * t17 - t42;
t18 = cos(t20);
t32 = -pkin(2) * t18 - t41;
t31 = -rSges(4,1) * t15 - rSges(4,2) * t16;
t12 = t16 * qJ(4);
t30 = t12 + t33;
t29 = (t38 + t47) * t16;
t28 = (g(2) * (-rSges(5,3) - qJ(4)) + g(3) * t47) * t15;
t27 = (-g(2) * qJ(4) + g(3) * t45) * t15 + t45 * t40;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t26 + t24 * rSges(2,2)) + g(3) * (-t24 * rSges(2,1) - rSges(2,2) * t26)) - m(3) * (g(2) * (-t18 * rSges(3,1) + t17 * rSges(3,2) - t41) + g(3) * (-rSges(3,1) * t17 - rSges(3,2) * t18 - t42)) - m(4) * (g(2) * (t32 + t34) + g(3) * (t31 + t33)) - m(5) * (g(2) * (t29 + t32) + g(3) * (t30 + t46) + t28) - m(6) * (g(2) * (t32 + t43) + g(3) * (t30 + t48) + t27), (-m(3) - m(4) + t44) * g(1), -m(4) * (g(2) * t34 + g(3) * t31) - m(5) * (g(2) * t29 + g(3) * (t12 + t46) + t28) - m(6) * (g(2) * t43 + g(3) * (t12 + t48) + t27), t44 * (g(3) * t15 + t40), -m(6) * (g(2) * (rSges(6,1) * t5 + rSges(6,2) * t6) + g(3) * (-rSges(6,1) * t7 + rSges(6,2) * t8) + g(1) * (-rSges(6,1) * t23 - rSges(6,2) * t25) * t21)];
taug = t1(:);
