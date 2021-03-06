% Calculate potential energy for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:36
% EndTime: 2019-12-31 16:34:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (102->72), mult. (125->86), div. (0->0), fcn. (113->8), ass. (0->25)
t34 = rSges(4,3) + pkin(5);
t33 = rSges(5,3) + pkin(6) + pkin(5);
t11 = sin(qJ(2));
t13 = cos(qJ(2));
t32 = rSges(3,1) * t13 - rSges(3,2) * t11;
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t31 = g(1) * t9 + g(2) * t8;
t10 = sin(qJ(3));
t28 = t10 * t8;
t27 = t10 * t9;
t26 = t8 * t13;
t25 = t9 * t13;
t21 = t10 * t13;
t12 = cos(qJ(3));
t20 = t12 * t13;
t18 = t8 * pkin(1) + r_base(2);
t17 = qJ(1) + r_base(3);
t16 = t9 * pkin(1) + t8 * pkin(4) + r_base(1);
t15 = -t9 * pkin(4) + t18;
t7 = qJ(3) + qJ(4);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(3) * t12 + pkin(2);
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t9 - rSges(2,2) * t8 + r_base(1)) + g(2) * (rSges(2,1) * t8 + rSges(2,2) * t9 + r_base(2)) + g(3) * (rSges(2,3) + t17)) - m(3) * (g(1) * (rSges(3,3) * t8 + t16) + g(2) * (t32 * t8 + t18) + g(3) * (rSges(3,1) * t11 + rSges(3,2) * t13 + t17) + (g(1) * t32 + g(2) * (-rSges(3,3) - pkin(4))) * t9) - m(4) * (g(1) * (pkin(2) * t25 + (t9 * t20 + t28) * rSges(4,1) + (t12 * t8 - t9 * t21) * rSges(4,2) + t16) + g(2) * (pkin(2) * t26 + (t8 * t20 - t27) * rSges(4,1) + (-t12 * t9 - t8 * t21) * rSges(4,2) + t15) + g(3) * (-t34 * t13 + t17) + (g(3) * (rSges(4,1) * t12 - rSges(4,2) * t10 + pkin(2)) + t31 * t34) * t11) - m(5) * (g(1) * (t1 * t25 + pkin(3) * t28 + (t2 * t8 + t3 * t25) * rSges(5,1) + (-t2 * t25 + t3 * t8) * rSges(5,2) + t16) + g(2) * (t1 * t26 - pkin(3) * t27 + (-t2 * t9 + t3 * t26) * rSges(5,1) + (-t2 * t26 - t3 * t9) * rSges(5,2) + t15) + g(3) * (-t33 * t13 + t17) + (g(3) * (rSges(5,1) * t3 - rSges(5,2) * t2 + t1) + t31 * t33) * t11);
U = t4;
