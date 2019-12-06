% Calculate potential energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:01
% EndTime: 2019-12-05 18:55:02
% DurationCPUTime: 0.28s
% Computational Cost: add. (144->82), mult. (144->99), div. (0->0), fcn. (128->10), ass. (0->32)
t12 = qJ(2) + qJ(3);
t7 = sin(t12);
t36 = pkin(5) * t7;
t17 = cos(qJ(2));
t35 = pkin(1) * t17;
t15 = sin(qJ(1));
t9 = cos(t12);
t34 = t15 * t9;
t18 = cos(qJ(1));
t33 = t18 * t9;
t13 = sin(qJ(4));
t32 = t15 * t13;
t16 = cos(qJ(4));
t31 = t15 * t16;
t30 = t18 * t13;
t29 = t18 * t16;
t28 = pkin(4) + r_base(3);
t27 = t15 * t35 + r_base(2);
t26 = t18 * t35 + r_base(1);
t25 = t15 * t36 + t27;
t24 = t18 * t36 + t26;
t14 = sin(qJ(2));
t23 = t14 * pkin(1) + t28;
t22 = rSges(4,1) * t9 - rSges(4,2) * t7;
t21 = g(1) * t18 + g(2) * t15;
t20 = rSges(3,1) * t17 - rSges(3,2) * t14;
t19 = -t9 * pkin(5) + t23;
t11 = qJ(4) + qJ(5);
t8 = cos(t11);
t6 = sin(t11);
t5 = t16 * pkin(3) + pkin(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t18 * rSges(2,1) - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + t18 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(1) * (t15 * rSges(3,3) + t20 * t18 + r_base(1)) + g(2) * (-t18 * rSges(3,3) + t20 * t15 + r_base(2)) + g(3) * (t14 * rSges(3,1) + t17 * rSges(3,2) + t28)) - m(4) * (g(1) * (t15 * rSges(4,3) + t22 * t18 + t26) + g(2) * (-t18 * rSges(4,3) + t22 * t15 + t27) + g(3) * (t7 * rSges(4,1) + t9 * rSges(4,2) + t23)) - m(5) * (g(1) * (pkin(2) * t33 + (t9 * t29 + t32) * rSges(5,1) + (-t9 * t30 + t31) * rSges(5,2) + t24) + g(2) * (pkin(2) * t34 + (t9 * t31 - t30) * rSges(5,1) + (-t9 * t32 - t29) * rSges(5,2) + t25) + g(3) * (-t9 * rSges(5,3) + t19) + (g(3) * (rSges(5,1) * t16 - rSges(5,2) * t13 + pkin(2)) + t21 * rSges(5,3)) * t7) - m(6) * (g(1) * (t5 * t33 + pkin(3) * t32 + (t15 * t6 + t8 * t33) * rSges(6,1) + (t15 * t8 - t6 * t33) * rSges(6,2) + t24) + g(2) * (t5 * t34 - pkin(3) * t30 + (-t18 * t6 + t8 * t34) * rSges(6,1) + (-t18 * t8 - t6 * t34) * rSges(6,2) + t25) + g(3) * (-t9 * rSges(6,3) + t19) + (g(3) * (rSges(6,1) * t8 - rSges(6,2) * t6 + t5) + t21 * rSges(6,3)) * t7);
U = t1;
