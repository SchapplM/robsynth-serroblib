% Calculate potential energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:10
% EndTime: 2019-12-31 16:30:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (94->67), mult. (138->79), div. (0->0), fcn. (130->6), ass. (0->24)
t32 = rSges(5,1) + pkin(3);
t31 = rSges(5,3) + qJ(4);
t14 = sin(pkin(6));
t17 = sin(qJ(2));
t30 = t14 * t17;
t19 = cos(qJ(2));
t29 = t14 * t19;
t15 = cos(pkin(6));
t28 = t15 * t17;
t16 = sin(qJ(3));
t27 = t16 * t19;
t18 = cos(qJ(3));
t26 = t18 * t19;
t25 = t14 * pkin(1) + r_base(2);
t24 = qJ(1) + r_base(3);
t23 = t15 * pkin(1) + t14 * pkin(4) + r_base(1);
t22 = t17 * pkin(2) + t24;
t21 = t15 * t19 * pkin(2) + pkin(5) * t28 + t23;
t20 = pkin(2) * t29 - t15 * pkin(4) + pkin(5) * t30 + t25;
t4 = t14 * t16 + t15 * t26;
t3 = -t14 * t18 + t15 * t27;
t2 = t14 * t26 - t15 * t16;
t1 = t14 * t27 + t15 * t18;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t15 - rSges(2,2) * t14 + r_base(1)) + g(2) * (rSges(2,1) * t14 + rSges(2,2) * t15 + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (t14 * rSges(3,3) + t23) + g(2) * (rSges(3,1) * t29 - rSges(3,2) * t30 + t25) + g(3) * (t17 * rSges(3,1) + rSges(3,2) * t19 + t24) + (g(1) * (rSges(3,1) * t19 - rSges(3,2) * t17) + g(2) * (-rSges(3,3) - pkin(4))) * t15) - m(4) * (g(1) * (rSges(4,1) * t4 - rSges(4,2) * t3 + rSges(4,3) * t28 + t21) + g(2) * (rSges(4,1) * t2 - rSges(4,2) * t1 + rSges(4,3) * t30 + t20) + g(3) * ((-rSges(4,3) - pkin(5)) * t19 + (rSges(4,1) * t18 - rSges(4,2) * t16) * t17 + t22)) - m(5) * (g(1) * (t31 * t3 + t32 * t4 + t21) + g(2) * (t31 * t1 + t32 * t2 + t20) + g(3) * (t22 + (-rSges(5,2) - pkin(5)) * t19) + (g(3) * (t31 * t16 + t32 * t18) + (g(1) * t15 + g(2) * t14) * rSges(5,2)) * t17);
U = t5;
