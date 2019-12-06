% Calculate potential energy for
% S5RRPPR2
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:00
% EndTime: 2019-12-05 18:20:00
% DurationCPUTime: 0.32s
% Computational Cost: add. (170->73), mult. (105->76), div. (0->0), fcn. (85->10), ass. (0->27)
t34 = rSges(6,3) + pkin(7);
t12 = sin(pkin(9));
t13 = cos(pkin(9));
t33 = rSges(5,1) * t13 - rSges(5,2) * t12;
t32 = t13 * pkin(4);
t14 = sin(qJ(5));
t28 = t13 * t14;
t16 = cos(qJ(5));
t27 = t13 * t16;
t11 = qJ(1) + qJ(2);
t26 = pkin(5) + r_base(1);
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + r_base(3);
t24 = pkin(6) + t26;
t9 = cos(t11);
t23 = pkin(2) * t9 + t25;
t15 = sin(qJ(1));
t22 = -t15 * pkin(1) + r_base(2);
t21 = qJ(3) + t24;
t7 = pkin(8) + t11;
t4 = sin(t7);
t5 = cos(t7);
t20 = t5 * pkin(3) + t4 * qJ(4) + t23;
t8 = sin(t11);
t19 = -pkin(2) * t8 + t22;
t18 = t5 * qJ(4) + t19;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t26) + g(2) * (-t15 * rSges(2,1) - rSges(2,2) * t17 + r_base(2)) + g(3) * (rSges(2,1) * t17 - t15 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t24) + g(2) * (-rSges(3,1) * t8 - rSges(3,2) * t9 + t22) + g(3) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t25)) - m(4) * (g(1) * (rSges(4,3) + t21) + g(2) * (-rSges(4,1) * t4 - rSges(4,2) * t5 + t19) + g(3) * (rSges(4,1) * t5 - rSges(4,2) * t4 + t23)) - m(5) * (g(1) * (rSges(5,1) * t12 + rSges(5,2) * t13 + t21) + g(2) * (rSges(5,3) * t5 + t18) + g(3) * (t33 * t5 + t20) + (g(2) * (-pkin(3) - t33) + g(3) * rSges(5,3)) * t4) - m(6) * (g(1) * (-t34 * t13 + t21) + g(2) * (t18 + (t14 * rSges(6,1) + t16 * rSges(6,2)) * t5 + (-t27 * rSges(6,1) + t28 * rSges(6,2) - pkin(3) - t32) * t4) + g(3) * (t5 * t32 + (t14 * t4 + t5 * t27) * rSges(6,1) + (t16 * t4 - t5 * t28) * rSges(6,2) + t20) + (g(1) * (rSges(6,1) * t16 - rSges(6,2) * t14 + pkin(4)) + (-g(2) * t4 + g(3) * t5) * t34) * t12);
U = t1;
