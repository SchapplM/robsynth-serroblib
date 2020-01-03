% Calculate potential energy for
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:41
% EndTime: 2019-12-31 16:26:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (93->50), mult. (69->48), div. (0->0), fcn. (49->6), ass. (0->18)
t21 = rSges(5,1) + pkin(3);
t20 = rSges(4,3) + pkin(5);
t19 = rSges(5,3) + qJ(4) + pkin(5);
t7 = sin(pkin(6));
t18 = t7 * pkin(1) + r_base(2);
t8 = cos(pkin(6));
t17 = t8 * pkin(1) + r_base(1);
t16 = qJ(1) + r_base(3);
t15 = pkin(4) + t16;
t10 = sin(qJ(3));
t11 = cos(qJ(3));
t14 = rSges(4,1) * t11 - rSges(4,2) * t10 + pkin(2);
t13 = -rSges(5,2) * t10 + t21 * t11 + pkin(2);
t12 = g(1) * t17 + g(2) * t18;
t6 = pkin(6) + qJ(2);
t3 = cos(t6);
t2 = sin(t6);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t8 - rSges(2,2) * t7 + r_base(1)) + g(2) * (rSges(2,1) * t7 + rSges(2,2) * t8 + r_base(2)) + g(3) * (rSges(2,3) + t16)) - m(3) * (g(1) * (rSges(3,1) * t3 - rSges(3,2) * t2 + t17) + g(2) * (rSges(3,1) * t2 + rSges(3,2) * t3 + t18) + g(3) * (rSges(3,3) + t15)) - m(4) * (g(3) * (t10 * rSges(4,1) + rSges(4,2) * t11 + t15) + (g(1) * t14 - g(2) * t20) * t3 + (g(1) * t20 + g(2) * t14) * t2 + t12) - m(5) * (g(3) * (rSges(5,2) * t11 + t21 * t10 + t15) + (g(1) * t13 - g(2) * t19) * t3 + (g(1) * t19 + g(2) * t13) * t2 + t12);
U = t1;
