% Calculate Gravitation load on the joints for
% S5RRRPR2
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:44
% DurationCPUTime: 0.25s
% Computational Cost: add. (308->64), mult. (162->74), div. (0->0), fcn. (120->10), ass. (0->39)
t46 = rSges(6,3) + pkin(8);
t19 = cos(qJ(5));
t38 = rSges(6,1) * t19;
t45 = -pkin(4) - t38;
t18 = sin(qJ(1));
t44 = pkin(1) * t18;
t20 = cos(qJ(1));
t43 = pkin(1) * t20;
t16 = qJ(1) + qJ(2);
t13 = sin(t16);
t42 = pkin(2) * t13;
t14 = cos(t16);
t41 = pkin(2) * t14;
t15 = qJ(3) + t16;
t11 = sin(t15);
t40 = pkin(3) * t11;
t12 = cos(t15);
t39 = pkin(3) * t12;
t17 = sin(qJ(5));
t37 = rSges(6,2) * t17;
t36 = -t14 * rSges(3,1) + t13 * rSges(3,2);
t35 = -t12 * rSges(4,1) + t11 * rSges(4,2);
t10 = pkin(9) + t15;
t7 = sin(t10);
t8 = cos(t10);
t34 = t7 * t37 + t46 * t8 - t40;
t33 = -rSges(3,1) * t13 - rSges(3,2) * t14;
t32 = -rSges(4,1) * t11 - rSges(4,2) * t12;
t30 = -t8 * rSges(5,1) + t7 * rSges(5,2) - t39;
t29 = t35 - t41;
t28 = -rSges(5,1) * t7 - rSges(5,2) * t8 - t40;
t27 = t34 - t42;
t26 = t32 - t42;
t25 = t30 - t41;
t24 = -t39 + (t37 + t45) * t8;
t23 = t28 - t42;
t22 = (-g(2) * t46 + g(3) * t45) * t7;
t21 = t24 - t41;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t20 + t18 * rSges(2,2)) + g(3) * (-t18 * rSges(2,1) - rSges(2,2) * t20)) - m(3) * (g(2) * (t36 - t43) + g(3) * (t33 - t44)) - m(4) * (g(2) * (t29 - t43) + g(3) * (t26 - t44)) - m(5) * (g(2) * (t25 - t43) + g(3) * (t23 - t44)) - m(6) * (g(2) * (t21 - t43) + g(3) * (t27 - t44) + t22), -m(3) * (g(2) * t36 + g(3) * t33) - m(4) * (g(2) * t29 + g(3) * t26) - m(5) * (g(2) * t25 + g(3) * t23) - m(6) * (g(2) * t21 + g(3) * t27 + t22), -m(4) * (g(2) * t35 + g(3) * t32) - m(5) * (g(2) * t30 + g(3) * t28) - m(6) * (g(2) * t24 + g(3) * t34 + t22), (-m(5) - m(6)) * g(1), -m(6) * (g(1) * (-t37 + t38) + (g(2) * t7 - g(3) * t8) * (rSges(6,1) * t17 + rSges(6,2) * t19))];
taug = t1(:);
