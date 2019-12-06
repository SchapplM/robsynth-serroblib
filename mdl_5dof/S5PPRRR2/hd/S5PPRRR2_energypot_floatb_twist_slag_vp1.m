% Calculate potential energy for
% S5PPRRR2
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:06
% EndTime: 2019-12-05 15:14:06
% DurationCPUTime: 0.46s
% Computational Cost: add. (169->89), mult. (154->101), div. (0->0), fcn. (138->10), ass. (0->37)
t45 = rSges(5,3) + pkin(6);
t44 = rSges(6,3) + pkin(7) + pkin(6);
t14 = sin(pkin(8));
t16 = cos(pkin(8));
t43 = g(1) * t16 + g(2) * t14;
t11 = pkin(9) + qJ(3);
t6 = sin(t11);
t42 = rSges(4,2) * t6;
t7 = cos(t11);
t39 = t14 * t7;
t12 = qJ(4) + qJ(5);
t8 = sin(t12);
t38 = t14 * t8;
t9 = cos(t12);
t37 = t14 * t9;
t36 = t16 * t7;
t35 = t16 * t8;
t34 = t16 * t9;
t18 = sin(qJ(4));
t32 = t14 * t18;
t19 = cos(qJ(4));
t31 = t14 * t19;
t30 = t16 * t18;
t29 = t16 * t19;
t27 = rSges(3,3) + qJ(2);
t15 = cos(pkin(9));
t4 = t15 * pkin(2) + pkin(1);
t26 = t16 * t4 + r_base(1);
t25 = qJ(1) + r_base(3);
t17 = -pkin(5) - qJ(2);
t24 = t14 * t4 + t16 * t17 + r_base(2);
t13 = sin(pkin(9));
t23 = t13 * pkin(2) + t25;
t22 = -t14 * t17 + t26;
t21 = rSges(3,1) * t15 - rSges(3,2) * t13 + pkin(1);
t5 = t19 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t16 * rSges(2,1) - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t13 * rSges(3,1) + t15 * rSges(3,2) + t25) + (g(1) * t21 - g(2) * t27) * t16 + (g(1) * t27 + g(2) * t21) * t14) - m(4) * (g(1) * (rSges(4,1) * t36 - t16 * t42 + t26) + g(2) * (-t16 * rSges(4,3) + t24) + g(3) * (t6 * rSges(4,1) + t7 * rSges(4,2) + t23) + (g(1) * (rSges(4,3) - t17) + g(2) * (rSges(4,1) * t7 - t42)) * t14) - m(5) * (g(1) * (pkin(3) * t36 + (t29 * t7 + t32) * rSges(5,1) + (-t30 * t7 + t31) * rSges(5,2) + t22) + g(2) * (pkin(3) * t39 + (t31 * t7 - t30) * rSges(5,1) + (-t32 * t7 - t29) * rSges(5,2) + t24) + g(3) * (-t45 * t7 + t23) + (g(3) * (rSges(5,1) * t19 - rSges(5,2) * t18 + pkin(3)) + t43 * t45) * t6) - m(6) * (g(1) * (t5 * t36 + pkin(4) * t32 + (t34 * t7 + t38) * rSges(6,1) + (-t35 * t7 + t37) * rSges(6,2) + t22) + g(2) * (t5 * t39 - pkin(4) * t30 + (t37 * t7 - t35) * rSges(6,1) + (-t38 * t7 - t34) * rSges(6,2) + t24) + g(3) * (-t44 * t7 + t23) + (g(3) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t5) + t43 * t44) * t6);
U = t1;
