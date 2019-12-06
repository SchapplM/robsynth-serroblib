% Calculate potential energy for
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:22
% EndTime: 2019-12-05 16:21:23
% DurationCPUTime: 0.52s
% Computational Cost: add. (170->98), mult. (180->116), div. (0->0), fcn. (168->10), ass. (0->32)
t41 = rSges(4,3) + pkin(6);
t17 = -qJ(4) - pkin(6);
t40 = rSges(5,3) - t17;
t39 = rSges(6,3) + pkin(7) - t17;
t15 = sin(pkin(8));
t16 = cos(pkin(8));
t38 = g(1) * t16 + g(2) * t15;
t20 = cos(qJ(3));
t5 = t20 * pkin(3) + pkin(2);
t19 = sin(qJ(2));
t34 = rSges(3,2) * t19;
t18 = sin(qJ(3));
t33 = t15 * t18;
t21 = cos(qJ(2));
t32 = t15 * t21;
t31 = t16 * t18;
t30 = t16 * t21;
t29 = t18 * t21;
t28 = t20 * t21;
t14 = qJ(3) + pkin(9);
t25 = t15 * pkin(1) + r_base(2);
t24 = qJ(1) + r_base(3);
t23 = t16 * pkin(1) + t15 * pkin(5) + r_base(1);
t22 = -t16 * pkin(5) + t25;
t8 = qJ(5) + t14;
t7 = cos(t14);
t6 = sin(t14);
t4 = cos(t8);
t3 = sin(t8);
t2 = t18 * pkin(3) + pkin(4) * t6;
t1 = pkin(4) * t7 + t5;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t16 * rSges(2,1) - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (t15 * rSges(3,3) + t23) + g(2) * (rSges(3,1) * t32 - t15 * t34 + t25) + g(3) * (t19 * rSges(3,1) + t21 * rSges(3,2) + t24) + (g(1) * (rSges(3,1) * t21 - t34) + g(2) * (-rSges(3,3) - pkin(5))) * t16) - m(4) * (g(1) * (pkin(2) * t30 + (t16 * t28 + t33) * rSges(4,1) + (t15 * t20 - t16 * t29) * rSges(4,2) + t23) + g(2) * (pkin(2) * t32 + (t15 * t28 - t31) * rSges(4,1) + (-t15 * t29 - t16 * t20) * rSges(4,2) + t22) + g(3) * (-t21 * t41 + t24) + (g(3) * (rSges(4,1) * t20 - rSges(4,2) * t18 + pkin(2)) + t38 * t41) * t19) - m(5) * (g(1) * (t5 * t30 + pkin(3) * t33 + (t15 * t6 + t30 * t7) * rSges(5,1) + (t15 * t7 - t30 * t6) * rSges(5,2) + t23) + g(2) * (t5 * t32 - pkin(3) * t31 + (-t16 * t6 + t32 * t7) * rSges(5,1) + (-t16 * t7 - t32 * t6) * rSges(5,2) + t22) + g(3) * (-t21 * t40 + t24) + (g(3) * (rSges(5,1) * t7 - rSges(5,2) * t6 + t5) + t38 * t40) * t19) - m(6) * (g(1) * (t1 * t30 + t15 * t2 + (t15 * t3 + t30 * t4) * rSges(6,1) + (t15 * t4 - t3 * t30) * rSges(6,2) + t23) + g(2) * (t1 * t32 - t16 * t2 + (-t16 * t3 + t32 * t4) * rSges(6,1) + (-t16 * t4 - t3 * t32) * rSges(6,2) + t22) + g(3) * (-t21 * t39 + t24) + (g(3) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t1) + t38 * t39) * t19);
U = t9;
