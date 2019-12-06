% Calculate potential energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:33
% EndTime: 2019-12-05 15:12:33
% DurationCPUTime: 0.36s
% Computational Cost: add. (168->81), mult. (130->88), div. (0->0), fcn. (110->10), ass. (0->34)
t40 = rSges(6,3) + pkin(7);
t15 = pkin(9) + qJ(3);
t11 = qJ(4) + t15;
t6 = sin(t11);
t39 = rSges(5,2) * t6;
t19 = cos(pkin(8));
t7 = cos(t11);
t38 = t19 * t7;
t18 = cos(pkin(9));
t8 = t18 * pkin(2) + pkin(1);
t17 = sin(pkin(8));
t21 = sin(qJ(5));
t36 = t17 * t21;
t22 = cos(qJ(5));
t35 = t17 * t22;
t34 = t19 * t21;
t33 = t19 * t22;
t20 = -pkin(5) - qJ(2);
t32 = rSges(4,3) - t20;
t31 = rSges(3,3) + qJ(2);
t10 = cos(t15);
t3 = pkin(3) * t10 + t8;
t30 = t19 * t3 + r_base(1);
t29 = qJ(1) + r_base(3);
t14 = -pkin(6) + t20;
t28 = t19 * t14 + t17 * t3 + r_base(2);
t16 = sin(pkin(9));
t27 = t16 * pkin(2) + t29;
t9 = sin(t15);
t26 = pkin(3) * t9 + t27;
t25 = rSges(4,1) * t10 - rSges(4,2) * t9 + t8;
t24 = rSges(3,1) * t18 - rSges(3,2) * t16 + pkin(1);
t23 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t29)) - m(3) * (g(3) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t29) + (g(1) * t24 - g(2) * t31) * t19 + (g(1) * t31 + g(2) * t24) * t17 + t23) - m(4) * (g(3) * (t9 * rSges(4,1) + t10 * rSges(4,2) + t27) + (g(1) * t25 - g(2) * t32) * t19 + (g(1) * t32 + g(2) * t25) * t17 + t23) - m(5) * (g(1) * (rSges(5,1) * t38 - t19 * t39 + t30) + g(2) * (-t19 * rSges(5,3) + t28) + g(3) * (t6 * rSges(5,1) + t7 * rSges(5,2) + t26) + (g(1) * (rSges(5,3) - t14) + g(2) * (rSges(5,1) * t7 - t39)) * t17) - m(6) * (g(1) * (pkin(4) * t38 - t17 * t14 + (t33 * t7 + t36) * rSges(6,1) + (-t34 * t7 + t35) * rSges(6,2) + t30) + g(2) * (t17 * t7 * pkin(4) + (t35 * t7 - t34) * rSges(6,1) + (-t36 * t7 - t33) * rSges(6,2) + t28) + g(3) * (-t40 * t7 + t26) + (g(3) * (rSges(6,1) * t22 - rSges(6,2) * t21 + pkin(4)) + (g(1) * t19 + g(2) * t17) * t40) * t6);
U = t1;
