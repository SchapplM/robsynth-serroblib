% Calculate potential energy for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:07
% EndTime: 2019-12-31 17:15:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (92->53), mult. (81->52), div. (0->0), fcn. (61->6), ass. (0->20)
t23 = rSges(5,1) + pkin(3);
t13 = -pkin(6) - pkin(5);
t11 = cos(qJ(2));
t2 = t11 * pkin(2) + pkin(1);
t22 = rSges(3,3) + pkin(5);
t21 = rSges(5,3) + qJ(4) - t13;
t20 = rSges(4,3) - t13;
t19 = pkin(4) + r_base(3);
t9 = sin(qJ(2));
t18 = t9 * pkin(2) + t19;
t8 = qJ(2) + qJ(3);
t3 = sin(t8);
t4 = cos(t8);
t17 = rSges(4,1) * t4 - rSges(4,2) * t3 + t2;
t16 = -rSges(5,2) * t3 + t23 * t4 + t2;
t15 = rSges(3,1) * t11 - rSges(3,2) * t9 + pkin(1);
t14 = g(1) * r_base(1) + g(2) * r_base(2);
t12 = cos(qJ(1));
t10 = sin(qJ(1));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t12 - rSges(2,2) * t10 + r_base(1)) + g(2) * (rSges(2,1) * t10 + rSges(2,2) * t12 + r_base(2)) + g(3) * (rSges(2,3) + t19)) - m(3) * (g(3) * (rSges(3,1) * t9 + rSges(3,2) * t11 + t19) + (g(1) * t15 - g(2) * t22) * t12 + (g(1) * t22 + g(2) * t15) * t10 + t14) - m(4) * (g(3) * (rSges(4,1) * t3 + rSges(4,2) * t4 + t18) + (g(1) * t17 - g(2) * t20) * t12 + (g(1) * t20 + g(2) * t17) * t10 + t14) - m(5) * (g(3) * (rSges(5,2) * t4 + t23 * t3 + t18) + (g(1) * t16 - g(2) * t21) * t12 + (g(1) * t21 + g(2) * t16) * t10 + t14);
U = t1;
