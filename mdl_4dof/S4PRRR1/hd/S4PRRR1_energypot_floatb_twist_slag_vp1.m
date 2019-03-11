% Calculate potential energy for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:03
% EndTime: 2019-03-08 18:25:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (91->50), mult. (48->42), div. (0->0), fcn. (28->8), ass. (0->19)
t13 = pkin(7) + qJ(2);
t14 = sin(pkin(7));
t22 = t14 * pkin(1) + r_base(2);
t15 = cos(pkin(7));
t21 = t15 * pkin(1) + r_base(1);
t20 = qJ(1) + r_base(3);
t8 = sin(t13);
t19 = pkin(2) * t8 + t22;
t9 = cos(t13);
t18 = pkin(2) * t9 + t21;
t17 = pkin(4) + t20;
t10 = qJ(3) + t13;
t16 = pkin(5) + t17;
t7 = qJ(4) + t10;
t6 = cos(t10);
t5 = sin(t10);
t2 = cos(t7);
t1 = sin(t7);
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t15 - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + rSges(2,2) * t15 + r_base(2)) + g(3) * (rSges(2,3) + t20)) - m(3) * (g(1) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t21) + g(2) * (rSges(3,1) * t8 + rSges(3,2) * t9 + t22) + g(3) * (rSges(3,3) + t17)) - m(4) * (g(1) * (rSges(4,1) * t6 - rSges(4,2) * t5 + t18) + g(2) * (rSges(4,1) * t5 + rSges(4,2) * t6 + t19) + g(3) * (rSges(4,3) + t16)) - m(5) * (g(1) * (rSges(5,1) * t2 - rSges(5,2) * t1 + pkin(3) * t6 + t18) + g(2) * (rSges(5,1) * t1 + rSges(5,2) * t2 + pkin(3) * t5 + t19) + g(3) * (pkin(6) + rSges(5,3) + t16));
U  = t3;
