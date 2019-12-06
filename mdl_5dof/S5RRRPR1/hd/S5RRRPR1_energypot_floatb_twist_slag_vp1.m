% Calculate potential energy for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:43
% EndTime: 2019-12-05 18:37:44
% DurationCPUTime: 0.22s
% Computational Cost: add. (156->70), mult. (110->68), div. (0->0), fcn. (86->10), ass. (0->30)
t23 = -pkin(7) - pkin(6);
t35 = rSges(3,3) + pkin(6);
t21 = cos(qJ(2));
t10 = t21 * pkin(2) + pkin(1);
t34 = rSges(4,3) - t23;
t17 = -qJ(4) + t23;
t33 = rSges(5,3) - t17;
t32 = rSges(6,3) + pkin(8) - t17;
t18 = qJ(2) + qJ(3);
t31 = pkin(5) + r_base(3);
t13 = cos(t18);
t2 = pkin(3) * t13 + t10;
t19 = sin(qJ(2));
t30 = t19 * pkin(2) + t31;
t11 = pkin(9) + t18;
t12 = sin(t18);
t29 = pkin(3) * t12 + t30;
t5 = sin(t11);
t6 = cos(t11);
t28 = rSges(5,1) * t6 - rSges(5,2) * t5 + t2;
t9 = qJ(5) + t11;
t3 = sin(t9);
t4 = cos(t9);
t27 = rSges(6,1) * t4 - rSges(6,2) * t3 + pkin(4) * t6 + t2;
t26 = rSges(3,1) * t21 - rSges(3,2) * t19 + pkin(1);
t25 = rSges(4,1) * t13 - rSges(4,2) * t12 + t10;
t24 = g(1) * r_base(1) + g(2) * r_base(2);
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t22 - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + rSges(2,2) * t22 + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(3) * (t19 * rSges(3,1) + rSges(3,2) * t21 + t31) + (g(1) * t26 - g(2) * t35) * t22 + (g(1) * t35 + g(2) * t26) * t20 + t24) - m(4) * (g(3) * (rSges(4,1) * t12 + rSges(4,2) * t13 + t30) + (g(1) * t25 - g(2) * t34) * t22 + (g(1) * t34 + g(2) * t25) * t20 + t24) - m(5) * (g(3) * (rSges(5,1) * t5 + rSges(5,2) * t6 + t29) + (g(1) * t28 - g(2) * t33) * t22 + (g(1) * t33 + g(2) * t28) * t20 + t24) - m(6) * (g(3) * (rSges(6,1) * t3 + rSges(6,2) * t4 + pkin(4) * t5 + t29) + (g(1) * t27 - g(2) * t32) * t22 + (g(1) * t32 + g(2) * t27) * t20 + t24);
U = t1;
