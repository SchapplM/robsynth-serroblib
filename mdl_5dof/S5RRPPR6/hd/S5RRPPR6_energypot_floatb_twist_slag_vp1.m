% Calculate potential energy for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:30
% EndTime: 2019-12-31 19:31:31
% DurationCPUTime: 0.42s
% Computational Cost: add. (169->89), mult. (154->101), div. (0->0), fcn. (138->10), ass. (0->33)
t41 = rSges(6,3) + pkin(7) + qJ(4);
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t40 = g(1) * t20 + g(2) * t18;
t39 = rSges(5,3) + qJ(4);
t12 = qJ(2) + pkin(8);
t7 = sin(t12);
t38 = rSges(4,2) * t7;
t9 = cos(t12);
t35 = t18 * t9;
t34 = t20 * t9;
t33 = rSges(3,3) + pkin(6);
t13 = sin(pkin(9));
t32 = t18 * t13;
t14 = cos(pkin(9));
t31 = t18 * t14;
t30 = t20 * t13;
t29 = t20 * t14;
t26 = pkin(5) + r_base(3);
t19 = cos(qJ(2));
t5 = t19 * pkin(2) + pkin(1);
t25 = t20 * t5 + r_base(1);
t15 = -qJ(3) - pkin(6);
t24 = t20 * t15 + t18 * t5 + r_base(2);
t17 = sin(qJ(2));
t23 = t17 * pkin(2) + t26;
t22 = -t18 * t15 + t25;
t21 = rSges(3,1) * t19 - rSges(3,2) * t17 + pkin(1);
t11 = pkin(9) + qJ(5);
t8 = cos(t11);
t6 = sin(t11);
t3 = t14 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t17 * rSges(3,1) + t19 * rSges(3,2) + t26) + (g(1) * t21 - g(2) * t33) * t20 + (g(1) * t33 + g(2) * t21) * t18) - m(4) * (g(1) * (rSges(4,1) * t34 - t20 * t38 + t25) + g(2) * (-t20 * rSges(4,3) + t24) + g(3) * (t7 * rSges(4,1) + t9 * rSges(4,2) + t23) + (g(1) * (rSges(4,3) - t15) + g(2) * (rSges(4,1) * t9 - t38)) * t18) - m(5) * (g(1) * (pkin(3) * t34 + (t29 * t9 + t32) * rSges(5,1) + (-t30 * t9 + t31) * rSges(5,2) + t22) + g(2) * (pkin(3) * t35 + (t31 * t9 - t30) * rSges(5,1) + (-t32 * t9 - t29) * rSges(5,2) + t24) + g(3) * (-t39 * t9 + t23) + (g(3) * (rSges(5,1) * t14 - rSges(5,2) * t13 + pkin(3)) + t40 * t39) * t7) - m(6) * (g(1) * (t3 * t34 + pkin(4) * t32 + (t18 * t6 + t34 * t8) * rSges(6,1) + (t18 * t8 - t6 * t34) * rSges(6,2) + t22) + g(2) * (t3 * t35 - pkin(4) * t30 + (-t20 * t6 + t35 * t8) * rSges(6,1) + (-t20 * t8 - t35 * t6) * rSges(6,2) + t24) + g(3) * (-t41 * t9 + t23) + (g(3) * (rSges(6,1) * t8 - rSges(6,2) * t6 + t3) + t40 * t41) * t7);
U = t1;
