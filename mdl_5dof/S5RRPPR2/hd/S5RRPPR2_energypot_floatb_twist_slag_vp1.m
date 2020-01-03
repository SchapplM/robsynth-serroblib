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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:57:17
% EndTime: 2020-01-03 11:57:18
% DurationCPUTime: 0.33s
% Computational Cost: add. (170->74), mult. (105->76), div. (0->0), fcn. (85->10), ass. (0->27)
t30 = rSges(6,3) + pkin(7);
t11 = cos(pkin(9));
t29 = t11 * pkin(4);
t12 = sin(qJ(5));
t27 = t11 * t12;
t14 = cos(qJ(5));
t26 = t11 * t14;
t25 = -rSges(5,3) - qJ(4);
t9 = qJ(1) + qJ(2);
t24 = pkin(5) + r_base(1);
t13 = sin(qJ(1));
t23 = t13 * pkin(1) + r_base(2);
t22 = pkin(6) + t24;
t6 = sin(t9);
t21 = pkin(2) * t6 + t23;
t5 = pkin(8) + t9;
t2 = sin(t5);
t20 = t2 * pkin(3) + t21;
t15 = cos(qJ(1));
t19 = -t15 * pkin(1) + r_base(3);
t18 = qJ(3) + t22;
t10 = sin(pkin(9));
t17 = rSges(5,1) * t11 - rSges(5,2) * t10;
t7 = cos(t9);
t16 = -pkin(2) * t7 + t19;
t3 = cos(t5);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t24) + g(2) * (t13 * rSges(2,1) + rSges(2,2) * t15 + r_base(2)) + g(3) * (-rSges(2,1) * t15 + t13 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t22) + g(2) * (rSges(3,1) * t6 + rSges(3,2) * t7 + t23) + g(3) * (-t7 * rSges(3,1) + t6 * rSges(3,2) + t19)) - m(4) * (g(1) * (rSges(4,3) + t18) + g(2) * (rSges(4,1) * t2 + rSges(4,2) * t3 + t21) + g(3) * (-t3 * rSges(4,1) + t2 * rSges(4,2) + t16)) - m(5) * (g(1) * (rSges(5,1) * t10 + rSges(5,2) * t11 + t18) + g(2) * t20 + g(3) * t16 + (g(2) * t17 + g(3) * t25) * t2 + (g(2) * t25 + g(3) * (-pkin(3) - t17)) * t3) - m(6) * (g(1) * (-t30 * t11 + t18) + g(2) * (t2 * t29 - t3 * qJ(4) + (-t12 * t3 + t2 * t26) * rSges(6,1) + (-t14 * t3 - t2 * t27) * rSges(6,2) + t20) + g(3) * (t16 + (-t12 * rSges(6,1) - t14 * rSges(6,2) - qJ(4)) * t2 + (-t26 * rSges(6,1) + t27 * rSges(6,2) - pkin(3) - t29) * t3) + (g(1) * (rSges(6,1) * t14 - rSges(6,2) * t12 + pkin(4)) + (g(2) * t2 - g(3) * t3) * t30) * t10);
U = t1;
