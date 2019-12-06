% Calculate potential energy for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:37
% EndTime: 2019-12-05 17:07:37
% DurationCPUTime: 0.21s
% Computational Cost: add. (154->64), mult. (85->60), div. (0->0), fcn. (61->10), ass. (0->26)
t31 = rSges(5,3) + pkin(7);
t30 = rSges(6,3) + pkin(8) + pkin(7);
t13 = pkin(9) + qJ(2);
t15 = sin(pkin(9));
t29 = t15 * pkin(1) + r_base(2);
t16 = cos(pkin(9));
t28 = t16 * pkin(1) + r_base(1);
t27 = qJ(1) + r_base(3);
t6 = sin(t13);
t26 = pkin(2) * t6 + t29;
t7 = cos(t13);
t25 = pkin(2) * t7 + t28;
t24 = pkin(5) + t27;
t23 = pkin(6) + t24;
t14 = qJ(4) + qJ(5);
t10 = cos(t14);
t18 = cos(qJ(4));
t9 = sin(t14);
t22 = rSges(6,1) * t10 - rSges(6,2) * t9 + pkin(4) * t18 + pkin(3);
t17 = sin(qJ(4));
t21 = rSges(5,1) * t18 - rSges(5,2) * t17 + pkin(3);
t20 = g(1) * t25 + g(2) * t26;
t8 = qJ(3) + t13;
t4 = cos(t8);
t3 = sin(t8);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t16 - rSges(2,2) * t15 + r_base(1)) + g(2) * (rSges(2,1) * t15 + rSges(2,2) * t16 + r_base(2)) + g(3) * (rSges(2,3) + t27)) - m(3) * (g(1) * (rSges(3,1) * t7 - rSges(3,2) * t6 + t28) + g(2) * (rSges(3,1) * t6 + rSges(3,2) * t7 + t29) + g(3) * (rSges(3,3) + t24)) - m(4) * (g(1) * (rSges(4,1) * t4 - rSges(4,2) * t3 + t25) + g(2) * (rSges(4,1) * t3 + rSges(4,2) * t4 + t26) + g(3) * (rSges(4,3) + t23)) - m(5) * (g(3) * (rSges(5,1) * t17 + rSges(5,2) * t18 + t23) + (g(1) * t21 - g(2) * t31) * t4 + (g(1) * t31 + g(2) * t21) * t3 + t20) - m(6) * (g(3) * (rSges(6,1) * t9 + rSges(6,2) * t10 + t17 * pkin(4) + t23) + (g(1) * t22 - g(2) * t30) * t4 + (g(1) * t30 + g(2) * t22) * t3 + t20);
U = t1;
