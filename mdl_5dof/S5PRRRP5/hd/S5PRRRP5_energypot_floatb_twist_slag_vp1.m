% Calculate potential energy for
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:39
% EndTime: 2019-12-05 16:47:40
% DurationCPUTime: 0.45s
% Computational Cost: add. (160->93), mult. (180->108), div. (0->0), fcn. (168->8), ass. (0->33)
t42 = rSges(4,3) + pkin(6);
t22 = -pkin(7) - pkin(6);
t41 = rSges(5,3) - t22;
t40 = rSges(6,3) + qJ(5) - t22;
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t39 = g(1) * t17 + g(2) * t16;
t20 = cos(qJ(3));
t7 = t20 * pkin(3) + pkin(2);
t19 = sin(qJ(2));
t35 = rSges(3,2) * t19;
t18 = sin(qJ(3));
t34 = t16 * t18;
t21 = cos(qJ(2));
t33 = t16 * t21;
t32 = t17 * t18;
t31 = t17 * t21;
t30 = t18 * t21;
t29 = t20 * t21;
t26 = t16 * pkin(1) + r_base(2);
t25 = qJ(1) + r_base(3);
t24 = t17 * pkin(1) + t16 * pkin(5) + r_base(1);
t23 = -t17 * pkin(5) + t26;
t15 = qJ(3) + qJ(4);
t9 = cos(t15);
t8 = sin(t15);
t6 = pkin(3) * t18 + pkin(4) * t8;
t5 = pkin(4) * t9 + t7;
t4 = t16 * t8 + t9 * t31;
t3 = t16 * t9 - t8 * t31;
t2 = -t17 * t8 + t33 * t9;
t1 = -t17 * t9 - t33 * t8;
t10 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t17 - rSges(2,2) * t16 + r_base(1)) + g(2) * (rSges(2,1) * t16 + rSges(2,2) * t17 + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (t16 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t33 - t16 * t35 + t26) + g(3) * (rSges(3,1) * t19 + rSges(3,2) * t21 + t25) + (g(1) * (rSges(3,1) * t21 - t35) + g(2) * (-rSges(3,3) - pkin(5))) * t17) - m(4) * (g(1) * (pkin(2) * t31 + (t17 * t29 + t34) * rSges(4,1) + (t16 * t20 - t17 * t30) * rSges(4,2) + t24) + g(2) * (pkin(2) * t33 + (t16 * t29 - t32) * rSges(4,1) + (-t16 * t30 - t17 * t20) * rSges(4,2) + t23) + g(3) * (-t42 * t21 + t25) + (g(3) * (rSges(4,1) * t20 - rSges(4,2) * t18 + pkin(2)) + t39 * t42) * t19) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + pkin(3) * t34 + t7 * t31 + t24) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) - pkin(3) * t32 + t33 * t7 + t23) + g(3) * (-t41 * t21 + t25) + (g(3) * (rSges(5,1) * t9 - rSges(5,2) * t8 + t7) + t39 * t41) * t19) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t16 * t6 + t5 * t31 + t24) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t17 * t6 + t33 * t5 + t23) + g(3) * (-t40 * t21 + t25) + (g(3) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t5) + t39 * t40) * t19);
U = t10;
