% Calculate potential energy for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:00
% EndTime: 2020-01-03 11:20:01
% DurationCPUTime: 0.40s
% Computational Cost: add. (175->80), mult. (141->85), div. (0->0), fcn. (125->10), ass. (0->30)
t13 = cos(pkin(8));
t12 = cos(pkin(9));
t8 = pkin(9) + qJ(5);
t3 = sin(t8);
t5 = cos(t8);
t19 = rSges(6,1) * t5 - rSges(6,2) * t3 + t12 * pkin(4) + pkin(3);
t37 = t13 * t19;
t36 = rSges(6,3) + pkin(6) + qJ(4);
t9 = qJ(1) + pkin(7);
t4 = sin(t9);
t6 = cos(t9);
t35 = g(2) * t4 - g(3) * t6;
t34 = rSges(5,3) + qJ(4);
t31 = t13 * pkin(3);
t10 = sin(pkin(9));
t30 = t10 * t13;
t29 = t12 * t13;
t27 = -rSges(4,3) - qJ(3);
t25 = pkin(5) + r_base(1);
t15 = sin(qJ(1));
t24 = t15 * pkin(1) + r_base(2);
t23 = t4 * pkin(2) + t24;
t22 = qJ(2) + t25;
t16 = cos(qJ(1));
t21 = -t16 * pkin(1) + r_base(3);
t11 = sin(pkin(8));
t20 = rSges(4,1) * t13 - rSges(4,2) * t11;
t18 = -t3 * rSges(6,1) - t5 * rSges(6,2) - t10 * pkin(4) - qJ(3);
t17 = g(2) * t23 + g(3) * t21;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t25) + g(2) * (t15 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (-t16 * rSges(2,1) + t15 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t22) + g(2) * (t4 * rSges(3,1) + t6 * rSges(3,2) + t24) + g(3) * (-t6 * rSges(3,1) + t4 * rSges(3,2) + t21)) - m(4) * (g(1) * (t11 * rSges(4,1) + t13 * rSges(4,2) + t22) + (g(2) * t20 + g(3) * t27) * t4 + (g(2) * t27 + g(3) * (-pkin(2) - t20)) * t6 + t17) - m(5) * (g(1) * (-t34 * t13 + t22) + g(2) * (t4 * t31 - t6 * qJ(3) + (-t6 * t10 + t29 * t4) * rSges(5,1) + (-t6 * t12 - t30 * t4) * rSges(5,2) + t23) + g(3) * (t21 + (-t10 * rSges(5,1) - t12 * rSges(5,2) - qJ(3)) * t4 + (-t29 * rSges(5,1) + t30 * rSges(5,2) - pkin(2) - t31) * t6) + (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t10 + pkin(3)) + t35 * t34) * t11) - m(6) * (g(1) * (-t36 * t13 + t22) + (g(2) * t37 + g(3) * t18) * t4 + (g(2) * t18 + (-pkin(2) - t37) * g(3)) * t6 + (g(1) * t19 + t35 * t36) * t11 + t17);
U = t1;
