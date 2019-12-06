% Calculate potential energy for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:41
% EndTime: 2019-12-05 16:50:41
% DurationCPUTime: 0.44s
% Computational Cost: add. (167->89), mult. (195->103), div. (0->0), fcn. (187->8), ass. (0->34)
t47 = rSges(6,1) + pkin(4);
t46 = rSges(4,3) + pkin(6);
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t45 = g(1) * t19 + g(2) * t18;
t44 = rSges(6,3) + qJ(5);
t21 = sin(qJ(2));
t40 = rSges(3,2) * t21;
t20 = sin(qJ(3));
t39 = t18 * t20;
t23 = cos(qJ(2));
t38 = t18 * t23;
t37 = t19 * t20;
t36 = t19 * t23;
t35 = t20 * t23;
t22 = cos(qJ(3));
t34 = t22 * t23;
t31 = t18 * pkin(1) + r_base(2);
t30 = qJ(1) + r_base(3);
t29 = t19 * pkin(1) + t18 * pkin(5) + r_base(1);
t10 = t22 * pkin(3) + pkin(2);
t24 = -pkin(7) - pkin(6);
t28 = t21 * t10 + t23 * t24 + t30;
t27 = -t19 * pkin(5) + t31;
t26 = pkin(3) * t39 + t10 * t36 + t29;
t25 = -pkin(3) * t37 + t10 * t38 + t27;
t17 = qJ(3) + qJ(4);
t13 = cos(t17);
t12 = sin(t17);
t4 = t18 * t12 + t13 * t36;
t3 = t12 * t36 - t18 * t13;
t2 = -t19 * t12 + t13 * t38;
t1 = t12 * t38 + t19 * t13;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t18 * rSges(3,3) + t29) + g(2) * (rSges(3,1) * t38 - t18 * t40 + t31) + g(3) * (t21 * rSges(3,1) + t23 * rSges(3,2) + t30) + (g(1) * (rSges(3,1) * t23 - t40) + g(2) * (-rSges(3,3) - pkin(5))) * t19) - m(4) * (g(1) * (pkin(2) * t36 + (t19 * t34 + t39) * rSges(4,1) + (t18 * t22 - t19 * t35) * rSges(4,2) + t29) + g(2) * (pkin(2) * t38 + (t18 * t34 - t37) * rSges(4,1) + (-t18 * t35 - t19 * t22) * rSges(4,2) + t27) + g(3) * (-t46 * t23 + t30) + (g(3) * (rSges(4,1) * t22 - rSges(4,2) * t20 + pkin(2)) + t45 * t46) * t21) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t26) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t25) + g(3) * (-t23 * rSges(5,3) + t28) + (g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12) + t45 * (rSges(5,3) - t24)) * t21) - m(6) * (g(1) * (t44 * t3 + t47 * t4 + t26) + g(2) * (t44 * t1 + t47 * t2 + t25) + g(3) * (-t23 * rSges(6,2) + t28) + (g(3) * (t44 * t12 + t47 * t13) + t45 * (rSges(6,2) - t24)) * t21);
U = t5;
