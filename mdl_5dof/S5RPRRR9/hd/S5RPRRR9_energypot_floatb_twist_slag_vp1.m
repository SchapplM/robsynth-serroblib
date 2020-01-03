% Calculate potential energy for
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR9_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:14
% EndTime: 2019-12-31 19:07:15
% DurationCPUTime: 0.36s
% Computational Cost: add. (168->81), mult. (130->88), div. (0->0), fcn. (110->10), ass. (0->34)
t40 = rSges(6,3) + pkin(8);
t15 = pkin(9) + qJ(3);
t11 = qJ(4) + t15;
t6 = sin(t11);
t39 = rSges(5,2) * t6;
t22 = cos(qJ(1));
t7 = cos(t11);
t38 = t22 * t7;
t17 = cos(pkin(9));
t8 = t17 * pkin(2) + pkin(1);
t19 = sin(qJ(5));
t20 = sin(qJ(1));
t36 = t20 * t19;
t21 = cos(qJ(5));
t35 = t20 * t21;
t34 = t22 * t19;
t33 = t22 * t21;
t18 = -pkin(6) - qJ(2);
t32 = rSges(4,3) - t18;
t31 = rSges(3,3) + qJ(2);
t30 = pkin(5) + r_base(3);
t10 = cos(t15);
t3 = pkin(3) * t10 + t8;
t29 = t22 * t3 + r_base(1);
t14 = -pkin(7) + t18;
t28 = t22 * t14 + t20 * t3 + r_base(2);
t16 = sin(pkin(9));
t27 = t16 * pkin(2) + t30;
t9 = sin(t15);
t26 = pkin(3) * t9 + t27;
t25 = rSges(4,1) * t10 - rSges(4,2) * t9 + t8;
t24 = rSges(3,1) * t17 - rSges(3,2) * t16 + pkin(1);
t23 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t22 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t22 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(3) * (t16 * rSges(3,1) + rSges(3,2) * t17 + t30) + (g(1) * t24 - g(2) * t31) * t22 + (g(1) * t31 + g(2) * t24) * t20 + t23) - m(4) * (g(3) * (rSges(4,1) * t9 + rSges(4,2) * t10 + t27) + (g(1) * t25 - g(2) * t32) * t22 + (g(1) * t32 + g(2) * t25) * t20 + t23) - m(5) * (g(1) * (rSges(5,1) * t38 - t22 * t39 + t29) + g(2) * (-t22 * rSges(5,3) + t28) + g(3) * (rSges(5,1) * t6 + rSges(5,2) * t7 + t26) + (g(1) * (rSges(5,3) - t14) + g(2) * (rSges(5,1) * t7 - t39)) * t20) - m(6) * (g(1) * (pkin(4) * t38 - t20 * t14 + (t7 * t33 + t36) * rSges(6,1) + (-t7 * t34 + t35) * rSges(6,2) + t29) + g(2) * (t20 * t7 * pkin(4) + (t7 * t35 - t34) * rSges(6,1) + (-t7 * t36 - t33) * rSges(6,2) + t28) + g(3) * (-t40 * t7 + t26) + (g(3) * (rSges(6,1) * t21 - rSges(6,2) * t19 + pkin(4)) + (g(1) * t22 + g(2) * t20) * t40) * t6);
U = t1;
