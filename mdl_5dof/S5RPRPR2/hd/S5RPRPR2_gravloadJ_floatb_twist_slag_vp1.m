% Calculate Gravitation load on the joints for
% S5RPRPR2
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:20
% EndTime: 2019-12-05 17:49:21
% DurationCPUTime: 0.22s
% Computational Cost: add. (224->54), mult. (138->66), div. (0->0), fcn. (104->10), ass. (0->31)
t19 = pkin(9) + qJ(5);
t14 = sin(t19);
t16 = cos(t19);
t48 = rSges(6,1) * t16 - rSges(6,2) * t14;
t47 = rSges(6,3) + pkin(7) + qJ(4);
t22 = cos(pkin(9));
t46 = -pkin(4) * t22 - pkin(3) - t48;
t45 = rSges(5,3) + qJ(4);
t44 = -rSges(5,1) * t22 - pkin(3);
t43 = -m(5) - m(6);
t24 = sin(qJ(1));
t42 = pkin(1) * t24;
t25 = cos(qJ(1));
t41 = pkin(1) * t25;
t38 = rSges(5,2) * sin(pkin(9));
t20 = qJ(1) + pkin(8);
t18 = qJ(3) + t20;
t11 = sin(t18);
t12 = cos(t18);
t36 = t11 * t38 + t45 * t12;
t35 = -t12 * rSges(4,1) + t11 * rSges(4,2);
t15 = sin(t20);
t33 = -pkin(2) * t15 - t42;
t17 = cos(t20);
t32 = -pkin(2) * t17 - t41;
t31 = -rSges(4,1) * t11 - rSges(4,2) * t12;
t29 = (t38 + t44) * t12;
t28 = -t47 * t11 + t46 * t12;
t27 = t46 * t11 + t47 * t12;
t26 = (-g(2) * t45 + g(3) * t44) * t11;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t25 + t24 * rSges(2,2)) + g(3) * (-t24 * rSges(2,1) - rSges(2,2) * t25)) - m(3) * (g(2) * (-t17 * rSges(3,1) + t15 * rSges(3,2) - t41) + g(3) * (-rSges(3,1) * t15 - rSges(3,2) * t17 - t42)) - m(4) * (g(2) * (t32 + t35) + g(3) * (t31 + t33)) - m(5) * (g(2) * (t29 + t32) + g(3) * (t33 + t36) + t26) - m(6) * (g(2) * (t28 + t32) + g(3) * (t27 + t33)), (-m(3) - m(4) + t43) * g(1), -m(4) * (g(2) * t35 + g(3) * t31) - m(5) * (g(2) * t29 + g(3) * t36 + t26) - m(6) * (g(2) * t28 + g(3) * t27), t43 * (g(2) * t12 + g(3) * t11), -m(6) * (g(1) * t48 + (g(2) * t11 - g(3) * t12) * (rSges(6,1) * t14 + rSges(6,2) * t16))];
taug = t1(:);
