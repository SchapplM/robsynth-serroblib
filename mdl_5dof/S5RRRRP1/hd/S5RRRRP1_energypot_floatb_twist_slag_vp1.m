% Calculate potential energy for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:40
% EndTime: 2019-12-05 18:44:41
% DurationCPUTime: 0.27s
% Computational Cost: add. (150->68), mult. (110->66), div. (0->0), fcn. (86->8), ass. (0->28)
t33 = rSges(6,1) + pkin(4);
t20 = -pkin(7) - pkin(6);
t32 = rSges(3,3) + pkin(6);
t18 = cos(qJ(2));
t7 = t18 * pkin(2) + pkin(1);
t31 = rSges(4,3) - t20;
t14 = -pkin(8) + t20;
t30 = rSges(5,3) - t14;
t29 = rSges(6,3) + qJ(5) - t14;
t15 = qJ(2) + qJ(3);
t28 = pkin(5) + r_base(3);
t9 = cos(t15);
t2 = pkin(3) * t9 + t7;
t16 = sin(qJ(2));
t27 = t16 * pkin(2) + t28;
t8 = sin(t15);
t26 = pkin(3) * t8 + t27;
t25 = rSges(4,1) * t9 - rSges(4,2) * t8 + t7;
t11 = qJ(4) + t15;
t5 = sin(t11);
t6 = cos(t11);
t24 = rSges(5,1) * t6 - rSges(5,2) * t5 + t2;
t23 = -rSges(6,2) * t5 + t33 * t6 + t2;
t22 = rSges(3,1) * t18 - rSges(3,2) * t16 + pkin(1);
t21 = g(1) * r_base(1) + g(2) * r_base(2);
t19 = cos(qJ(1));
t17 = sin(qJ(1));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t19 - rSges(2,2) * t17 + r_base(1)) + g(2) * (rSges(2,1) * t17 + rSges(2,2) * t19 + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(3) * (rSges(3,1) * t16 + rSges(3,2) * t18 + t28) + (g(1) * t22 - g(2) * t32) * t19 + (g(1) * t32 + g(2) * t22) * t17 + t21) - m(4) * (g(3) * (rSges(4,1) * t8 + rSges(4,2) * t9 + t27) + (g(1) * t25 - g(2) * t31) * t19 + (g(1) * t31 + g(2) * t25) * t17 + t21) - m(5) * (g(3) * (rSges(5,1) * t5 + rSges(5,2) * t6 + t26) + (g(1) * t24 - g(2) * t30) * t19 + (g(1) * t30 + g(2) * t24) * t17 + t21) - m(6) * (g(3) * (rSges(6,2) * t6 + t33 * t5 + t26) + (g(1) * t23 - g(2) * t29) * t19 + (g(1) * t29 + g(2) * t23) * t17 + t21);
U = t1;
