% Calculate potential energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:54
% EndTime: 2019-12-31 17:54:54
% DurationCPUTime: 0.28s
% Computational Cost: add. (118->69), mult. (110->65), div. (0->0), fcn. (86->6), ass. (0->26)
t30 = rSges(6,1) + pkin(4);
t29 = rSges(6,3) + qJ(5);
t9 = sin(pkin(7));
t28 = pkin(3) * t9;
t11 = -pkin(6) - qJ(3);
t27 = rSges(6,2) - t11;
t26 = rSges(5,3) - t11;
t25 = rSges(4,3) + qJ(3);
t24 = pkin(5) + r_base(3);
t12 = sin(qJ(1));
t23 = t12 * pkin(1) + r_base(2);
t22 = pkin(2) + t24;
t13 = cos(qJ(1));
t21 = t13 * pkin(1) + t12 * qJ(2) + r_base(1);
t20 = -qJ(2) - t28;
t19 = g(2) * t23;
t10 = cos(pkin(7));
t18 = t10 * pkin(3) + t22;
t8 = pkin(7) + qJ(4);
t2 = sin(t8);
t3 = cos(t8);
t17 = rSges(5,1) * t2 + rSges(5,2) * t3;
t16 = rSges(4,1) * t9 + rSges(4,2) * t10;
t15 = t30 * t2 - t29 * t3;
t14 = g(1) * (t12 * t28 + t21) + t19;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t13 - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t13 + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (-rSges(3,2) * t13 + t12 * rSges(3,3) + t21) + g(2) * (-t12 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t13 + t23) + g(3) * (rSges(3,1) + t24)) - m(4) * (g(1) * t21 + t19 + g(3) * (rSges(4,1) * t10 - rSges(4,2) * t9 + t22) + (g(1) * t16 + g(2) * t25) * t12 + (g(1) * t25 + g(2) * (-qJ(2) - t16)) * t13) - m(5) * (g(3) * (rSges(5,1) * t3 - rSges(5,2) * t2 + t18) + (g(1) * t17 + g(2) * t26) * t12 + (g(1) * t26 + g(2) * (-t17 + t20)) * t13 + t14) - m(6) * (g(3) * (t29 * t2 + t30 * t3 + t18) + (g(1) * t15 + g(2) * t27) * t12 + (g(1) * t27 + g(2) * (-t15 + t20)) * t13 + t14);
U = t1;
