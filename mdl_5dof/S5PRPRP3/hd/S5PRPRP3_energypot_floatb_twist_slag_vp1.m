% Calculate potential energy for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:34
% EndTime: 2019-12-05 15:32:34
% DurationCPUTime: 0.37s
% Computational Cost: add. (159->84), mult. (154->93), div. (0->0), fcn. (138->8), ass. (0->34)
t42 = rSges(5,3) + pkin(6);
t41 = rSges(6,3) + qJ(5) + pkin(6);
t14 = sin(pkin(7));
t15 = cos(pkin(7));
t40 = g(1) * t15 + g(2) * t14;
t37 = rSges(3,3) + pkin(5);
t13 = qJ(2) + pkin(8);
t10 = sin(t13);
t35 = rSges(4,2) * t10;
t11 = cos(t13);
t34 = t11 * t14;
t33 = t11 * t15;
t18 = sin(qJ(4));
t32 = t14 * t18;
t20 = cos(qJ(4));
t31 = t14 * t20;
t30 = t15 * t18;
t29 = t15 * t20;
t21 = cos(qJ(2));
t9 = t21 * pkin(2) + pkin(1);
t27 = t15 * t9 + r_base(1);
t26 = qJ(1) + r_base(3);
t17 = -qJ(3) - pkin(5);
t25 = t14 * t9 + t15 * t17 + r_base(2);
t19 = sin(qJ(2));
t24 = t19 * pkin(2) + t26;
t23 = -t14 * t17 + t27;
t22 = rSges(3,1) * t21 - rSges(3,2) * t19 + pkin(1);
t8 = t20 * pkin(4) + pkin(3);
t4 = t11 * t29 + t32;
t3 = -t11 * t30 + t31;
t2 = t11 * t31 - t30;
t1 = -t11 * t32 - t29;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t15 * rSges(2,1) - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + t15 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t19 * rSges(3,1) + t21 * rSges(3,2) + t26) + (g(1) * t22 - g(2) * t37) * t15 + (g(1) * t37 + g(2) * t22) * t14) - m(4) * (g(1) * (rSges(4,1) * t33 - t15 * t35 + t27) + g(2) * (-t15 * rSges(4,3) + t25) + g(3) * (t10 * rSges(4,1) + t11 * rSges(4,2) + t24) + (g(1) * (rSges(4,3) - t17) + g(2) * (rSges(4,1) * t11 - t35)) * t14) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + pkin(3) * t33 + t23) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + pkin(3) * t34 + t25) + g(3) * (-t42 * t11 + t24) + (g(3) * (rSges(5,1) * t20 - rSges(5,2) * t18 + pkin(3)) + t40 * t42) * t10) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t32 + t8 * t33 + t23) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - pkin(4) * t30 + t8 * t34 + t25) + g(3) * (-t41 * t11 + t24) + (g(3) * (rSges(6,1) * t20 - rSges(6,2) * t18 + t8) + t40 * t41) * t10);
U = t5;
