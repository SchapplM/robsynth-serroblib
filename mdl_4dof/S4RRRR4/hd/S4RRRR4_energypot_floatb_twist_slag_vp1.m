% Calculate potential energy for
% S4RRRR4
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:38
% EndTime: 2019-12-31 17:25:38
% DurationCPUTime: 0.29s
% Computational Cost: add. (104->66), mult. (101->74), div. (0->0), fcn. (85->8), ass. (0->25)
t29 = rSges(5,3) + pkin(7);
t8 = qJ(2) + qJ(3);
t5 = sin(t8);
t28 = rSges(4,2) * t5;
t11 = sin(qJ(1));
t9 = sin(qJ(4));
t27 = t11 * t9;
t14 = cos(qJ(1));
t6 = cos(t8);
t26 = t14 * t6;
t25 = t14 * t9;
t24 = rSges(3,3) + pkin(5);
t12 = cos(qJ(4));
t22 = t11 * t12;
t21 = t12 * t14;
t20 = pkin(4) + r_base(3);
t13 = cos(qJ(2));
t3 = pkin(2) * t13 + pkin(1);
t19 = t14 * t3 + r_base(1);
t10 = sin(qJ(2));
t18 = t10 * pkin(2) + t20;
t15 = -pkin(6) - pkin(5);
t17 = t11 * t3 + t14 * t15 + r_base(2);
t16 = rSges(3,1) * t13 - rSges(3,2) * t10 + pkin(1);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t14 - rSges(2,2) * t11 + r_base(1)) + g(2) * (rSges(2,1) * t11 + rSges(2,2) * t14 + r_base(2)) + g(3) * (rSges(2,3) + t20)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t10 + rSges(3,2) * t13 + t20) + (g(1) * t16 - g(2) * t24) * t14 + (g(1) * t24 + g(2) * t16) * t11) - m(4) * (g(1) * (rSges(4,1) * t26 - t14 * t28 + t19) + g(2) * (-rSges(4,3) * t14 + t17) + g(3) * (rSges(4,1) * t5 + rSges(4,2) * t6 + t18) + (g(1) * (rSges(4,3) - t15) + g(2) * (rSges(4,1) * t6 - t28)) * t11) - m(5) * (g(1) * (pkin(3) * t26 - t11 * t15 + (t6 * t21 + t27) * rSges(5,1) + (-t6 * t25 + t22) * rSges(5,2) + t19) + g(2) * (t11 * t6 * pkin(3) + (t6 * t22 - t25) * rSges(5,1) + (-t6 * t27 - t21) * rSges(5,2) + t17) + g(3) * (-t29 * t6 + t18) + (g(3) * (rSges(5,1) * t12 - rSges(5,2) * t9 + pkin(3)) + (g(1) * t14 + g(2) * t11) * t29) * t5);
U = t1;
