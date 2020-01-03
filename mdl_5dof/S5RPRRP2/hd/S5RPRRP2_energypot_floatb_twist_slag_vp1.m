% Calculate potential energy for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:02
% EndTime: 2020-01-03 11:45:02
% DurationCPUTime: 0.27s
% Computational Cost: add. (148->62), mult. (85->58), div. (0->0), fcn. (61->8), ass. (0->24)
t27 = rSges(6,1) + pkin(4);
t26 = -rSges(5,3) - pkin(7);
t25 = -rSges(6,3) - qJ(5) - pkin(7);
t24 = pkin(5) + r_base(1);
t9 = qJ(1) + pkin(8);
t12 = sin(qJ(1));
t23 = t12 * pkin(1) + r_base(2);
t5 = sin(t9);
t22 = pkin(2) * t5 + t23;
t21 = qJ(2) + t24;
t14 = cos(qJ(1));
t20 = -t14 * pkin(1) + r_base(3);
t19 = pkin(6) + t21;
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t18 = rSges(5,1) * t13 - rSges(5,2) * t11 + pkin(3);
t17 = -rSges(6,2) * t11 + t27 * t13 + pkin(3);
t6 = cos(t9);
t16 = -pkin(2) * t6 + t20;
t15 = g(2) * t22 + g(3) * t16;
t7 = qJ(3) + t9;
t3 = cos(t7);
t2 = sin(t7);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t24) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t14 + r_base(2)) + g(3) * (-rSges(2,1) * t14 + t12 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t21) + g(2) * (rSges(3,1) * t5 + rSges(3,2) * t6 + t23) + g(3) * (-t6 * rSges(3,1) + t5 * rSges(3,2) + t20)) - m(4) * (g(1) * (rSges(4,3) + t19) + g(2) * (rSges(4,1) * t2 + rSges(4,2) * t3 + t22) + g(3) * (-t3 * rSges(4,1) + t2 * rSges(4,2) + t16)) - m(5) * (g(1) * (rSges(5,1) * t11 + rSges(5,2) * t13 + t19) + (g(2) * t26 - g(3) * t18) * t3 + (g(2) * t18 + g(3) * t26) * t2 + t15) - m(6) * (g(1) * (rSges(6,2) * t13 + t27 * t11 + t19) + (g(2) * t25 - g(3) * t17) * t3 + (g(2) * t17 + g(3) * t25) * t2 + t15);
U = t1;
