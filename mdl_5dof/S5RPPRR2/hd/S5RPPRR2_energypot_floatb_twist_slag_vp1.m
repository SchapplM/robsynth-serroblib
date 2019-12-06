% Calculate potential energy for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:30
% EndTime: 2019-12-05 17:39:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (120->69), mult. (103->65), div. (0->0), fcn. (79->8), ass. (0->25)
t13 = sin(pkin(8));
t30 = pkin(3) * t13;
t15 = -pkin(6) - qJ(3);
t29 = rSges(5,3) - t15;
t28 = rSges(6,3) + pkin(7) - t15;
t27 = rSges(4,3) + qJ(3);
t26 = pkin(5) + r_base(3);
t16 = sin(qJ(1));
t25 = t16 * pkin(1) + r_base(2);
t12 = pkin(8) + qJ(4);
t24 = pkin(2) + t26;
t17 = cos(qJ(1));
t23 = t17 * pkin(1) + t16 * qJ(2) + r_base(1);
t14 = cos(pkin(8));
t22 = t14 * pkin(3) + t24;
t21 = rSges(4,1) * t13 + rSges(4,2) * t14;
t6 = qJ(5) + t12;
t2 = sin(t6);
t3 = cos(t6);
t4 = sin(t12);
t20 = rSges(6,1) * t2 + rSges(6,2) * t3 + pkin(4) * t4 + t30;
t5 = cos(t12);
t19 = rSges(5,1) * t4 + rSges(5,2) * t5 + t30;
t18 = g(1) * t23 + g(2) * t25;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t17 - t16 * rSges(2,2) + r_base(1)) + g(2) * (t16 * rSges(2,1) + rSges(2,2) * t17 + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (-rSges(3,2) * t17 + t16 * rSges(3,3) + t23) + g(2) * (-t16 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t17 + t25) + g(3) * (rSges(3,1) + t26)) - m(4) * (g(3) * (rSges(4,1) * t14 - rSges(4,2) * t13 + t24) + (g(1) * t21 + g(2) * t27) * t16 + (g(1) * t27 + g(2) * (-qJ(2) - t21)) * t17 + t18) - m(5) * (g(3) * (rSges(5,1) * t5 - rSges(5,2) * t4 + t22) + (g(1) * t19 + g(2) * t29) * t16 + (g(1) * t29 + g(2) * (-qJ(2) - t19)) * t17 + t18) - m(6) * (g(3) * (rSges(6,1) * t3 - rSges(6,2) * t2 + pkin(4) * t5 + t22) + (g(1) * t20 + g(2) * t28) * t16 + (g(1) * t28 + g(2) * (-qJ(2) - t20)) * t17 + t18);
U = t1;
