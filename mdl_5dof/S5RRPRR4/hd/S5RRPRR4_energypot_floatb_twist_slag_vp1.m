% Calculate potential energy for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:46
% EndTime: 2020-01-03 12:01:46
% DurationCPUTime: 0.21s
% Computational Cost: add. (154->64), mult. (85->60), div. (0->0), fcn. (61->10), ass. (0->26)
t29 = -rSges(5,3) - pkin(7);
t28 = -rSges(6,3) - pkin(8) - pkin(7);
t12 = qJ(1) + qJ(2);
t27 = pkin(5) + r_base(1);
t14 = sin(qJ(1));
t26 = t14 * pkin(1) + r_base(2);
t25 = pkin(6) + t27;
t7 = sin(t12);
t24 = pkin(2) * t7 + t26;
t16 = cos(qJ(1));
t23 = -t16 * pkin(1) + r_base(3);
t22 = qJ(3) + t25;
t15 = cos(qJ(4));
t11 = qJ(4) + qJ(5);
t6 = sin(t11);
t8 = cos(t11);
t21 = rSges(6,1) * t8 - rSges(6,2) * t6 + pkin(4) * t15 + pkin(3);
t13 = sin(qJ(4));
t20 = rSges(5,1) * t15 - rSges(5,2) * t13 + pkin(3);
t9 = cos(t12);
t19 = -pkin(2) * t9 + t23;
t18 = g(2) * t24 + g(3) * t19;
t5 = pkin(9) + t12;
t2 = cos(t5);
t1 = sin(t5);
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t27) + g(2) * (rSges(2,1) * t14 + rSges(2,2) * t16 + r_base(2)) + g(3) * (-rSges(2,1) * t16 + rSges(2,2) * t14 + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t25) + g(2) * (rSges(3,1) * t7 + rSges(3,2) * t9 + t26) + g(3) * (-t9 * rSges(3,1) + rSges(3,2) * t7 + t23)) - m(4) * (g(1) * (rSges(4,3) + t22) + g(2) * (rSges(4,1) * t1 + rSges(4,2) * t2 + t24) + g(3) * (-rSges(4,1) * t2 + rSges(4,2) * t1 + t19)) - m(5) * (g(1) * (rSges(5,1) * t13 + rSges(5,2) * t15 + t22) + (g(2) * t29 - g(3) * t20) * t2 + (g(2) * t20 + g(3) * t29) * t1 + t18) - m(6) * (g(1) * (rSges(6,1) * t6 + rSges(6,2) * t8 + pkin(4) * t13 + t22) + (g(2) * t28 - g(3) * t21) * t2 + (g(2) * t21 + g(3) * t28) * t1 + t18);
U = t3;
