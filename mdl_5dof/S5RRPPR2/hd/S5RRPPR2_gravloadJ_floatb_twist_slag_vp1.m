% Calculate Gravitation load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:03
% EndTime: 2019-12-05 18:20:04
% DurationCPUTime: 0.25s
% Computational Cost: add. (281->69), mult. (192->87), div. (0->0), fcn. (168->10), ass. (0->40)
t21 = qJ(1) + qJ(2);
t18 = pkin(8) + t21;
t16 = sin(t18);
t17 = cos(t18);
t26 = cos(qJ(5));
t23 = cos(pkin(9));
t24 = sin(qJ(5));
t40 = t23 * t24;
t5 = t16 * t40 + t17 * t26;
t39 = t23 * t26;
t6 = t16 * t39 - t17 * t24;
t52 = -t6 * rSges(6,1) + t5 * rSges(6,2);
t51 = -rSges(5,1) * t23 - pkin(3);
t22 = sin(pkin(9));
t41 = rSges(5,2) * t22;
t50 = t17 * rSges(5,3) + t16 * t41;
t49 = -pkin(4) * t23 - pkin(3) + (-rSges(6,3) - pkin(7)) * t22;
t48 = -m(5) - m(6);
t25 = sin(qJ(1));
t47 = pkin(1) * t25;
t27 = cos(qJ(1));
t46 = pkin(1) * t27;
t19 = sin(t21);
t45 = pkin(2) * t19;
t20 = cos(t21);
t44 = pkin(2) * t20;
t43 = g(2) * t17;
t37 = t17 * qJ(4) - t45;
t36 = -t20 * rSges(3,1) + t19 * rSges(3,2);
t7 = -t16 * t26 + t17 * t40;
t8 = -t16 * t24 - t17 * t39;
t35 = t8 * rSges(6,1) + t7 * rSges(6,2) - t44;
t34 = -rSges(3,1) * t19 - rSges(3,2) * t20;
t33 = t37 - t47;
t32 = -t17 * rSges(4,1) + t16 * rSges(4,2) - t44;
t31 = -rSges(4,1) * t16 - rSges(4,2) * t17 - t45;
t30 = -t44 + (t41 + t51) * t17;
t29 = (g(2) * (-rSges(5,3) - qJ(4)) + g(3) * t51) * t16;
t28 = (-g(2) * qJ(4) + g(3) * t49) * t16 + t49 * t43;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t27 + t25 * rSges(2,2)) + g(3) * (-t25 * rSges(2,1) - rSges(2,2) * t27)) - m(3) * (g(2) * (t36 - t46) + g(3) * (t34 - t47)) - m(4) * (g(2) * (t32 - t46) + g(3) * (t31 - t47)) - m(5) * (g(2) * (t30 - t46) + g(3) * (t33 + t50) + t29) - m(6) * (g(2) * (t35 - t46) + g(3) * (t33 + t52) + t28), -m(3) * (g(2) * t36 + g(3) * t34) - m(4) * (g(2) * t32 + g(3) * t31) - m(5) * (g(2) * t30 + g(3) * (t37 + t50) + t29) - m(6) * (g(2) * t35 + g(3) * (t37 + t52) + t28), (-m(4) + t48) * g(1), t48 * (g(3) * t16 + t43), -m(6) * (g(2) * (rSges(6,1) * t5 + rSges(6,2) * t6) + g(3) * (-rSges(6,1) * t7 + rSges(6,2) * t8) + g(1) * (-rSges(6,1) * t24 - rSges(6,2) * t26) * t22)];
taug = t1(:);
