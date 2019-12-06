% Calculate potential energy for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:20
% EndTime: 2019-12-05 15:53:20
% DurationCPUTime: 0.51s
% Computational Cost: add. (170->98), mult. (180->114), div. (0->0), fcn. (168->10), ass. (0->30)
t19 = -pkin(6) - qJ(3);
t39 = rSges(5,3) - t19;
t38 = rSges(6,3) + pkin(7) - t19;
t16 = sin(pkin(8));
t18 = cos(pkin(8));
t37 = g(1) * t18 + g(2) * t16;
t36 = rSges(4,3) + qJ(3);
t17 = cos(pkin(9));
t5 = t17 * pkin(3) + pkin(2);
t20 = sin(qJ(2));
t33 = rSges(3,2) * t20;
t15 = sin(pkin(9));
t32 = t16 * t15;
t21 = cos(qJ(2));
t31 = t16 * t21;
t30 = t18 * t15;
t29 = t18 * t21;
t14 = pkin(9) + qJ(4);
t25 = t16 * pkin(1) + r_base(2);
t24 = qJ(1) + r_base(3);
t23 = t18 * pkin(1) + t16 * pkin(5) + r_base(1);
t22 = -t18 * pkin(5) + t25;
t8 = qJ(5) + t14;
t7 = cos(t14);
t6 = sin(t14);
t4 = cos(t8);
t3 = sin(t8);
t2 = t15 * pkin(3) + pkin(4) * t6;
t1 = pkin(4) * t7 + t5;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t18 * rSges(2,1) - t16 * rSges(2,2) + r_base(1)) + g(2) * (t16 * rSges(2,1) + t18 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (t16 * rSges(3,3) + t23) + g(2) * (rSges(3,1) * t31 - t16 * t33 + t25) + g(3) * (t20 * rSges(3,1) + t21 * rSges(3,2) + t24) + (g(1) * (rSges(3,1) * t21 - t33) + g(2) * (-rSges(3,3) - pkin(5))) * t18) - m(4) * (g(1) * (pkin(2) * t29 + (t17 * t29 + t32) * rSges(4,1) + (-t15 * t29 + t16 * t17) * rSges(4,2) + t23) + g(2) * (pkin(2) * t31 + (t17 * t31 - t30) * rSges(4,1) + (-t15 * t31 - t18 * t17) * rSges(4,2) + t22) + g(3) * (-t36 * t21 + t24) + (g(3) * (rSges(4,1) * t17 - rSges(4,2) * t15 + pkin(2)) + t37 * t36) * t20) - m(5) * (g(1) * (t5 * t29 + pkin(3) * t32 + (t16 * t6 + t29 * t7) * rSges(5,1) + (t16 * t7 - t29 * t6) * rSges(5,2) + t23) + g(2) * (t5 * t31 - pkin(3) * t30 + (-t18 * t6 + t31 * t7) * rSges(5,1) + (-t18 * t7 - t31 * t6) * rSges(5,2) + t22) + g(3) * (-t39 * t21 + t24) + (g(3) * (rSges(5,1) * t7 - rSges(5,2) * t6 + t5) + t37 * t39) * t20) - m(6) * (g(1) * (t1 * t29 + t16 * t2 + (t16 * t3 + t29 * t4) * rSges(6,1) + (t16 * t4 - t29 * t3) * rSges(6,2) + t23) + g(2) * (t1 * t31 - t18 * t2 + (-t18 * t3 + t31 * t4) * rSges(6,1) + (-t18 * t4 - t3 * t31) * rSges(6,2) + t22) + g(3) * (-t38 * t21 + t24) + (g(3) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t1) + t37 * t38) * t20);
U = t9;
