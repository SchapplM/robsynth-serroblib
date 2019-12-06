% Calculate potential energy for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:33
% EndTime: 2019-12-05 16:41:34
% DurationCPUTime: 0.24s
% Computational Cost: add. (157->65), mult. (92->62), div. (0->0), fcn. (68->8), ass. (0->24)
t32 = rSges(6,1) + pkin(4);
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t31 = rSges(5,1) * t17 - rSges(5,2) * t16;
t30 = rSges(6,3) + qJ(5);
t13 = pkin(8) + qJ(2);
t14 = sin(pkin(8));
t27 = t14 * pkin(1) + r_base(2);
t15 = cos(pkin(8));
t26 = t15 * pkin(1) + r_base(1);
t25 = qJ(1) + r_base(3);
t8 = sin(t13);
t24 = pkin(2) * t8 + t27;
t9 = cos(t13);
t23 = pkin(2) * t9 + t26;
t22 = pkin(5) + t25;
t10 = qJ(3) + t13;
t6 = sin(t10);
t21 = t6 * pkin(3) + t24;
t20 = pkin(6) + t22;
t7 = cos(t10);
t19 = t7 * pkin(3) + t6 * pkin(7) + t23;
t18 = t30 * t16 + t32 * t17;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t15 - rSges(2,2) * t14 + r_base(1)) + g(2) * (rSges(2,1) * t14 + rSges(2,2) * t15 + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t26) + g(2) * (rSges(3,1) * t8 + rSges(3,2) * t9 + t27) + g(3) * (rSges(3,3) + t22)) - m(4) * (g(1) * (rSges(4,1) * t7 - rSges(4,2) * t6 + t23) + g(2) * (rSges(4,1) * t6 + rSges(4,2) * t7 + t24) + g(3) * (rSges(4,3) + t20)) - m(5) * (g(1) * (t6 * rSges(5,3) + t19) + g(2) * (t31 * t6 + t21) + g(3) * (t16 * rSges(5,1) + rSges(5,2) * t17 + t20) + (g(1) * t31 + g(2) * (-rSges(5,3) - pkin(7))) * t7) - m(6) * (g(1) * t19 + g(2) * t21 + g(3) * (t32 * t16 - t30 * t17 + t20) + (g(1) * rSges(6,2) + g(2) * t18) * t6 + (g(1) * t18 + g(2) * (-rSges(6,2) - pkin(7))) * t7);
U = t1;
