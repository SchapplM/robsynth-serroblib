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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:28:41
% EndTime: 2019-12-05 17:28:42
% DurationCPUTime: 0.42s
% Computational Cost: add. (175->84), mult. (141->96), div. (0->0), fcn. (125->10), ass. (0->31)
t41 = rSges(6,3) + pkin(6) + qJ(4);
t13 = sin(pkin(8));
t15 = cos(pkin(8));
t40 = rSges(4,1) * t15 - rSges(4,2) * t13;
t11 = qJ(1) + pkin(7);
t6 = sin(t11);
t8 = cos(t11);
t39 = -g(2) * t6 + g(3) * t8;
t38 = rSges(5,3) + qJ(4);
t35 = t15 * t6;
t34 = t15 * t8;
t12 = sin(pkin(9));
t33 = t6 * t12;
t32 = t8 * t12;
t29 = t12 * t15;
t14 = cos(pkin(9));
t28 = t14 * t15;
t25 = pkin(5) + r_base(1);
t18 = cos(qJ(1));
t24 = t18 * pkin(1) + r_base(3);
t23 = qJ(2) + t25;
t22 = t8 * pkin(2) + t6 * qJ(3) + t24;
t17 = sin(qJ(1));
t21 = -t17 * pkin(1) + r_base(2);
t20 = t8 * qJ(3) + t21;
t19 = -t6 * pkin(2) + t20;
t10 = pkin(9) + qJ(5);
t7 = cos(t10);
t5 = sin(t10);
t4 = t14 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t25) + g(2) * (-t17 * rSges(2,1) - t18 * rSges(2,2) + r_base(2)) + g(3) * (t18 * rSges(2,1) - t17 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t23) + g(2) * (-t6 * rSges(3,1) - t8 * rSges(3,2) + t21) + g(3) * (t8 * rSges(3,1) - t6 * rSges(3,2) + t24)) - m(4) * (g(1) * (t13 * rSges(4,1) + t15 * rSges(4,2) + t23) + g(2) * (t8 * rSges(4,3) + t20) + g(3) * (t40 * t8 + t22) + (g(2) * (-pkin(2) - t40) + g(3) * rSges(4,3)) * t6) - m(5) * (g(1) * (-t38 * t15 + t23) + g(2) * (-pkin(3) * t35 + (-t28 * t6 + t32) * rSges(5,1) + (t8 * t14 + t29 * t6) * rSges(5,2) + t19) + g(3) * (pkin(3) * t34 + (t28 * t8 + t33) * rSges(5,1) + (t6 * t14 - t29 * t8) * rSges(5,2) + t22) + (g(1) * (rSges(5,1) * t14 - rSges(5,2) * t12 + pkin(3)) + t39 * t38) * t13) - m(6) * (g(1) * (-t41 * t15 + t23) + g(2) * (-t4 * t35 + pkin(4) * t32 + (-t35 * t7 + t8 * t5) * rSges(6,1) + (t35 * t5 + t8 * t7) * rSges(6,2) + t19) + g(3) * (t4 * t34 + pkin(4) * t33 + (t34 * t7 + t6 * t5) * rSges(6,1) + (-t34 * t5 + t6 * t7) * rSges(6,2) + t22) + (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t5 + t4) + t39 * t41) * t13);
U = t1;
