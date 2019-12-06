% Calculate potential energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:40
% EndTime: 2019-12-05 17:04:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (104->56), mult. (60->50), div. (0->0), fcn. (36->8), ass. (0->21)
t24 = rSges(6,3) + pkin(6);
t10 = qJ(2) + qJ(3);
t23 = pkin(1) + r_base(1);
t12 = sin(qJ(2));
t22 = t12 * pkin(2) + r_base(2);
t21 = qJ(1) + r_base(3);
t14 = cos(qJ(2));
t20 = t14 * pkin(2) + t23;
t5 = sin(t10);
t19 = pkin(3) * t5 + t22;
t18 = pkin(4) + t21;
t6 = cos(t10);
t17 = pkin(3) * t6 + t20;
t16 = pkin(5) + t18;
t11 = sin(qJ(5));
t13 = cos(qJ(5));
t15 = rSges(6,1) * t13 - rSges(6,2) * t11;
t7 = qJ(4) + t10;
t4 = cos(t7);
t3 = sin(t7);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (r_base(1) + rSges(2,1)) + g(2) * (r_base(2) + rSges(2,2)) + g(3) * (rSges(2,3) + t21)) - m(3) * (g(1) * (rSges(3,1) * t14 - t12 * rSges(3,2) + t23) + g(2) * (t12 * rSges(3,1) + rSges(3,2) * t14 + r_base(2)) + g(3) * (rSges(3,3) + t21)) - m(4) * (g(1) * (rSges(4,1) * t6 - rSges(4,2) * t5 + t20) + g(2) * (rSges(4,1) * t5 + rSges(4,2) * t6 + t22) + g(3) * (rSges(4,3) + t18)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + t17) + g(2) * (rSges(5,1) * t3 + rSges(5,2) * t4 + t19) + g(3) * (rSges(5,3) + t16)) - m(6) * (g(1) * t17 + g(2) * t19 + g(3) * (rSges(6,1) * t11 + t13 * rSges(6,2) + t16) + (g(1) * t15 - g(2) * t24) * t4 + (g(1) * t24 + g(2) * t15) * t3);
U = t1;
