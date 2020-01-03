% Calculate potential energy for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:50
% EndTime: 2019-12-31 17:12:51
% DurationCPUTime: 0.19s
% Computational Cost: add. (93->50), mult. (69->48), div. (0->0), fcn. (49->6), ass. (0->18)
t21 = rSges(5,1) + pkin(3);
t20 = rSges(4,3) + pkin(6);
t19 = rSges(5,3) + qJ(4) + pkin(6);
t18 = pkin(4) + r_base(3);
t9 = sin(qJ(1));
t17 = t9 * pkin(1) + r_base(2);
t11 = cos(qJ(1));
t16 = t11 * pkin(1) + r_base(1);
t15 = pkin(5) + t18;
t10 = cos(qJ(3));
t8 = sin(qJ(3));
t14 = rSges(4,1) * t10 - rSges(4,2) * t8 + pkin(2);
t13 = -rSges(5,2) * t8 + t21 * t10 + pkin(2);
t12 = g(1) * t16 + g(2) * t17;
t6 = qJ(1) + qJ(2);
t3 = cos(t6);
t2 = sin(t6);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t11 - t9 * rSges(2,2) + r_base(1)) + g(2) * (t9 * rSges(2,1) + rSges(2,2) * t11 + r_base(2)) + g(3) * (rSges(2,3) + t18)) - m(3) * (g(1) * (rSges(3,1) * t3 - rSges(3,2) * t2 + t16) + g(2) * (rSges(3,1) * t2 + rSges(3,2) * t3 + t17) + g(3) * (rSges(3,3) + t15)) - m(4) * (g(3) * (rSges(4,1) * t8 + rSges(4,2) * t10 + t15) + (g(1) * t14 - g(2) * t20) * t3 + (g(1) * t20 + g(2) * t14) * t2 + t12) - m(5) * (g(3) * (rSges(5,2) * t10 + t21 * t8 + t15) + (g(1) * t13 - g(2) * t19) * t3 + (g(1) * t19 + g(2) * t13) * t2 + t12);
U = t1;
