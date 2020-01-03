% Calculate potential energy for
% S5RRRRP2
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:00
% EndTime: 2020-01-03 12:11:01
% DurationCPUTime: 0.27s
% Computational Cost: add. (147->65), mult. (97->62), div. (0->0), fcn. (73->8), ass. (0->26)
t30 = rSges(6,1) + pkin(4);
t17 = -pkin(8) - pkin(7);
t15 = cos(qJ(3));
t2 = t15 * pkin(3) + pkin(2);
t29 = -rSges(4,3) - pkin(7);
t28 = -rSges(5,3) + t17;
t27 = -rSges(6,3) - qJ(5) + t17;
t26 = pkin(5) + r_base(1);
t14 = sin(qJ(1));
t25 = t14 * pkin(1) + r_base(2);
t24 = pkin(6) + t26;
t13 = sin(qJ(3));
t23 = t13 * pkin(3) + t24;
t16 = cos(qJ(1));
t22 = -t16 * pkin(1) + r_base(3);
t11 = qJ(3) + qJ(4);
t3 = sin(t11);
t5 = cos(t11);
t21 = rSges(5,1) * t5 - rSges(5,2) * t3 + t2;
t20 = -rSges(6,2) * t3 + t30 * t5 + t2;
t19 = rSges(4,1) * t15 - rSges(4,2) * t13 + pkin(2);
t18 = g(2) * t25 + g(3) * t22;
t12 = qJ(1) + qJ(2);
t6 = cos(t12);
t4 = sin(t12);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t26) + g(2) * (rSges(2,1) * t14 + rSges(2,2) * t16 + r_base(2)) + g(3) * (-rSges(2,1) * t16 + rSges(2,2) * t14 + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t24) + g(2) * (rSges(3,1) * t4 + rSges(3,2) * t6 + t25) + g(3) * (-rSges(3,1) * t6 + rSges(3,2) * t4 + t22)) - m(4) * (g(1) * (rSges(4,1) * t13 + rSges(4,2) * t15 + t24) + (g(2) * t29 - g(3) * t19) * t6 + (g(2) * t19 + g(3) * t29) * t4 + t18) - m(5) * (g(1) * (rSges(5,1) * t3 + rSges(5,2) * t5 + t23) + (g(2) * t28 - g(3) * t21) * t6 + (g(2) * t21 + g(3) * t28) * t4 + t18) - m(6) * (g(1) * (rSges(6,2) * t5 + t30 * t3 + t23) + (g(2) * t27 - g(3) * t20) * t6 + (g(2) * t20 + g(3) * t27) * t4 + t18);
U = t1;
