% Calculate potential energy for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:43:57
% EndTime: 2019-12-05 17:43:58
% DurationCPUTime: 0.52s
% Computational Cost: add. (170->98), mult. (180->114), div. (0->0), fcn. (168->10), ass. (0->32)
t19 = -pkin(6) - qJ(3);
t41 = rSges(5,3) - t19;
t40 = rSges(6,3) + pkin(7) - t19;
t20 = sin(qJ(1));
t21 = cos(qJ(1));
t39 = -g(2) * t20 + g(3) * t21;
t38 = rSges(4,3) + qJ(3);
t17 = cos(pkin(9));
t5 = t17 * pkin(3) + pkin(2);
t16 = sin(pkin(8));
t35 = rSges(3,2) * t16;
t18 = cos(pkin(8));
t34 = t18 * t20;
t33 = t18 * t21;
t15 = sin(pkin(9));
t32 = t20 * t15;
t31 = t20 * t17;
t30 = t21 * t15;
t29 = t21 * t17;
t25 = pkin(5) + r_base(1);
t14 = pkin(9) + qJ(4);
t24 = t21 * qJ(2) + r_base(2);
t23 = t21 * pkin(1) + t20 * qJ(2) + r_base(3);
t22 = -t20 * pkin(1) + t24;
t8 = qJ(5) + t14;
t7 = cos(t14);
t6 = sin(t14);
t4 = cos(t8);
t3 = sin(t8);
t2 = t15 * pkin(3) + pkin(4) * t6;
t1 = pkin(4) * t7 + t5;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t25) + g(2) * (-t20 * rSges(2,1) - t21 * rSges(2,2) + r_base(2)) + g(3) * (t21 * rSges(2,1) - t20 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t25) + g(2) * (t21 * rSges(3,3) + t24) + g(3) * (rSges(3,1) * t33 - t21 * t35 + t23) + (g(2) * (-rSges(3,1) * t18 - pkin(1) + t35) + g(3) * rSges(3,3)) * t20) - m(4) * (g(1) * (-t38 * t18 + t25) + g(2) * (-pkin(2) * t34 + (-t18 * t31 + t30) * rSges(4,1) + (t18 * t32 + t29) * rSges(4,2) + t22) + g(3) * (pkin(2) * t33 + (t18 * t29 + t32) * rSges(4,1) + (-t18 * t30 + t31) * rSges(4,2) + t23) + (g(1) * (rSges(4,1) * t17 - rSges(4,2) * t15 + pkin(2)) + t39 * t38) * t16) - m(5) * (g(1) * (-t41 * t18 + t25) + g(2) * (-t5 * t34 + pkin(3) * t30 + (t21 * t6 - t34 * t7) * rSges(5,1) + (t21 * t7 + t34 * t6) * rSges(5,2) + t22) + g(3) * (t5 * t33 + pkin(3) * t32 + (t20 * t6 + t33 * t7) * rSges(5,1) + (t20 * t7 - t33 * t6) * rSges(5,2) + t23) + (g(1) * (rSges(5,1) * t7 - rSges(5,2) * t6 + t5) + t39 * t41) * t16) - m(6) * (g(1) * (-t40 * t18 + t25) + g(2) * (-t1 * t34 + t21 * t2 + (t21 * t3 - t34 * t4) * rSges(6,1) + (t21 * t4 + t3 * t34) * rSges(6,2) + t22) + g(3) * (t1 * t33 + t20 * t2 + (t20 * t3 + t33 * t4) * rSges(6,1) + (t20 * t4 - t3 * t33) * rSges(6,2) + t23) + (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t1) + t39 * t40) * t16);
U = t9;
