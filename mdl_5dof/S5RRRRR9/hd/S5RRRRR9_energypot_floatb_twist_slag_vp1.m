% Calculate potential energy for
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR9_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:38
% EndTime: 2019-12-31 22:27:38
% DurationCPUTime: 0.50s
% Computational Cost: add. (170->98), mult. (180->114), div. (0->0), fcn. (168->10), ass. (0->30)
t39 = rSges(4,3) + pkin(7);
t21 = -pkin(8) - pkin(7);
t38 = rSges(5,3) - t21;
t37 = rSges(6,3) + pkin(9) - t21;
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t36 = g(1) * t20 + g(2) * t17;
t18 = cos(qJ(3));
t5 = t18 * pkin(3) + pkin(2);
t16 = sin(qJ(2));
t32 = rSges(3,2) * t16;
t15 = sin(qJ(3));
t31 = t17 * t15;
t19 = cos(qJ(2));
t30 = t17 * t19;
t29 = t20 * t15;
t28 = t20 * t19;
t14 = qJ(3) + qJ(4);
t25 = pkin(5) + r_base(3);
t24 = t17 * pkin(1) + r_base(2);
t23 = t20 * pkin(1) + t17 * pkin(6) + r_base(1);
t22 = -t20 * pkin(6) + t24;
t8 = qJ(5) + t14;
t7 = cos(t14);
t6 = sin(t14);
t4 = cos(t8);
t3 = sin(t8);
t2 = t15 * pkin(3) + pkin(4) * t6;
t1 = pkin(4) * t7 + t5;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (t17 * rSges(3,3) + t23) + g(2) * (rSges(3,1) * t30 - t17 * t32 + t24) + g(3) * (t16 * rSges(3,1) + t19 * rSges(3,2) + t25) + (g(1) * (rSges(3,1) * t19 - t32) + g(2) * (-rSges(3,3) - pkin(6))) * t20) - m(4) * (g(1) * (pkin(2) * t28 + (t18 * t28 + t31) * rSges(4,1) + (-t15 * t28 + t17 * t18) * rSges(4,2) + t23) + g(2) * (pkin(2) * t30 + (t18 * t30 - t29) * rSges(4,1) + (-t15 * t30 - t20 * t18) * rSges(4,2) + t22) + g(3) * (-t39 * t19 + t25) + (g(3) * (rSges(4,1) * t18 - rSges(4,2) * t15 + pkin(2)) + t36 * t39) * t16) - m(5) * (g(1) * (t5 * t28 + pkin(3) * t31 + (t17 * t6 + t28 * t7) * rSges(5,1) + (t17 * t7 - t28 * t6) * rSges(5,2) + t23) + g(2) * (t5 * t30 - pkin(3) * t29 + (-t20 * t6 + t30 * t7) * rSges(5,1) + (-t20 * t7 - t30 * t6) * rSges(5,2) + t22) + g(3) * (-t38 * t19 + t25) + (g(3) * (rSges(5,1) * t7 - rSges(5,2) * t6 + t5) + t36 * t38) * t16) - m(6) * (g(1) * (t1 * t28 + t17 * t2 + (t17 * t3 + t28 * t4) * rSges(6,1) + (t17 * t4 - t28 * t3) * rSges(6,2) + t23) + g(2) * (t1 * t30 - t20 * t2 + (-t20 * t3 + t30 * t4) * rSges(6,1) + (-t20 * t4 - t3 * t30) * rSges(6,2) + t22) + g(3) * (-t37 * t19 + t25) + (g(3) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t1) + t36 * t37) * t16);
U = t9;
