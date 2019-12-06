% Calculate potential energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:49:41
% EndTime: 2019-12-05 18:49:42
% DurationCPUTime: 0.31s
% Computational Cost: add. (150->73), mult. (122->83), div. (0->0), fcn. (102->10), ass. (0->30)
t33 = rSges(6,3) + pkin(6);
t11 = qJ(2) + qJ(3);
t9 = qJ(4) + t11;
t5 = cos(t9);
t32 = t5 * pkin(4);
t16 = cos(qJ(2));
t6 = t16 * pkin(2) + pkin(1);
t12 = sin(qJ(5));
t14 = sin(qJ(1));
t30 = t14 * t12;
t15 = cos(qJ(5));
t29 = t14 * t15;
t17 = cos(qJ(1));
t28 = t17 * t12;
t27 = t17 * t15;
t26 = pkin(5) + r_base(3);
t8 = cos(t11);
t3 = pkin(3) * t8 + t6;
t25 = t14 * t3 + r_base(2);
t24 = t17 * t3 + r_base(1);
t4 = sin(t9);
t23 = rSges(5,1) * t5 - rSges(5,2) * t4;
t13 = sin(qJ(2));
t22 = -t13 * pkin(2) + t26;
t7 = sin(t11);
t21 = rSges(4,1) * t8 - rSges(4,2) * t7 + t6;
t20 = rSges(3,1) * t16 - rSges(3,2) * t13 + pkin(1);
t19 = g(1) * r_base(1) + g(2) * r_base(2);
t18 = -pkin(3) * t7 + t22;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t17 * rSges(2,1) - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + t17 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(3) * (-t13 * rSges(3,1) - t16 * rSges(3,2) + t26) + (g(2) * rSges(3,3) + g(1) * t20) * t17 + (-g(1) * rSges(3,3) + g(2) * t20) * t14 + t19) - m(4) * (g(3) * (-t7 * rSges(4,1) - t8 * rSges(4,2) + t22) + (g(2) * rSges(4,3) + g(1) * t21) * t17 + (-g(1) * rSges(4,3) + g(2) * t21) * t14 + t19) - m(5) * (g(1) * (-t14 * rSges(5,3) + t17 * t23 + t24) + g(2) * (t17 * rSges(5,3) + t14 * t23 + t25) + g(3) * (-t4 * rSges(5,1) - t5 * rSges(5,2) + t18)) - m(6) * (g(1) * (t17 * t32 + (t27 * t5 - t30) * rSges(6,1) + (-t28 * t5 - t29) * rSges(6,2) + t24) + g(2) * (t14 * t32 + (t29 * t5 + t28) * rSges(6,1) + (-t30 * t5 + t27) * rSges(6,2) + t25) + g(3) * (t33 * t5 + t18) + (g(3) * (-rSges(6,1) * t15 + rSges(6,2) * t12 - pkin(4)) + (g(1) * t17 + g(2) * t14) * t33) * t4);
U = t1;
