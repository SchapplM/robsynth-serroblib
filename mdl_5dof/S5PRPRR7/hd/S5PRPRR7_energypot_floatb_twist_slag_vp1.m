% Calculate potential energy for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:38
% EndTime: 2019-12-05 15:59:39
% DurationCPUTime: 0.48s
% Computational Cost: add. (134->92), mult. (172->109), div. (0->0), fcn. (156->8), ass. (0->29)
t41 = rSges(5,3) + pkin(6);
t40 = rSges(6,3) + pkin(7) + pkin(6);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t39 = g(1) * t14 + g(2) * t13;
t16 = sin(qJ(2));
t35 = t13 * t16;
t18 = cos(qJ(2));
t34 = t13 * t18;
t33 = t14 * t16;
t15 = sin(qJ(4));
t32 = t15 * t16;
t17 = cos(qJ(4));
t31 = t16 * t17;
t29 = qJ(3) * t16;
t28 = t13 * pkin(1) + r_base(2);
t27 = t13 * t32;
t26 = t14 * t32;
t25 = qJ(1) + r_base(3);
t24 = t14 * pkin(1) + t13 * pkin(5) + r_base(1);
t23 = t16 * pkin(2) + t25;
t22 = pkin(2) * t34 + t13 * t29 + t28;
t21 = t24 + (pkin(2) * t18 + t29) * t14;
t20 = -t14 * pkin(5) + t22;
t12 = qJ(4) + qJ(5);
t7 = cos(t12);
t6 = sin(t12);
t5 = t17 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t14 * rSges(2,1) - t13 * rSges(2,2) + r_base(1)) + g(2) * (t13 * rSges(2,1) + t14 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (t13 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t34 - rSges(3,2) * t35 + t28) + g(3) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t25) + (g(1) * (rSges(3,1) * t18 - rSges(3,2) * t16) + g(2) * (-rSges(3,3) - pkin(5))) * t14) - m(4) * (g(1) * (t13 * rSges(4,1) + t21) + g(2) * (-rSges(4,2) * t34 + rSges(4,3) * t35 + t22) + g(3) * (-t16 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t18 + t23) + (g(1) * (-rSges(4,2) * t18 + rSges(4,3) * t16) + g(2) * (-rSges(4,1) - pkin(5))) * t14) - m(5) * (g(1) * (t13 * pkin(3) + (t13 * t17 + t26) * rSges(5,1) + (-t13 * t15 + t14 * t31) * rSges(5,2) + t21) + g(2) * (-t14 * pkin(3) + (-t14 * t17 + t27) * rSges(5,1) + (t13 * t31 + t14 * t15) * rSges(5,2) + t20) + g(3) * (t41 * t16 + t23) + (g(3) * (-rSges(5,1) * t15 - rSges(5,2) * t17 - qJ(3)) + t39 * t41) * t18) - m(6) * (g(1) * (t13 * t5 + pkin(4) * t26 + (t13 * t7 + t6 * t33) * rSges(6,1) + (-t13 * t6 + t7 * t33) * rSges(6,2) + t21) + g(2) * (-t14 * t5 + pkin(4) * t27 + (-t14 * t7 + t6 * t35) * rSges(6,1) + (t14 * t6 + t7 * t35) * rSges(6,2) + t20) + g(3) * (t40 * t16 + t23) + (g(3) * (-rSges(6,1) * t6 - rSges(6,2) * t7 - pkin(4) * t15 - qJ(3)) + t39 * t40) * t18);
U = t1;
