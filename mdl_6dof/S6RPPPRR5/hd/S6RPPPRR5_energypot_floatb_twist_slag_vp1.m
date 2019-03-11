% Calculate potential energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:09
% EndTime: 2019-03-09 01:37:09
% DurationCPUTime: 0.35s
% Computational Cost: add. (157->87), mult. (197->94), div. (0->0), fcn. (205->8), ass. (0->28)
t37 = rSges(7,3) + pkin(8);
t17 = cos(qJ(1));
t32 = sin(qJ(1));
t27 = t32 * pkin(1) + r_base(2);
t25 = t32 * qJ(3) + t27;
t36 = t25 + (-pkin(3) - qJ(2)) * t17;
t16 = cos(qJ(5));
t35 = t16 * pkin(5);
t34 = -rSges(6,3) - pkin(7);
t13 = sin(qJ(6));
t31 = t13 * t16;
t15 = cos(qJ(6));
t30 = t15 * t16;
t29 = sin(pkin(9));
t28 = pkin(6) + r_base(3);
t26 = pkin(2) + t28;
t24 = t17 * pkin(1) + t32 * qJ(2) + r_base(1);
t23 = t17 * qJ(3) + t24;
t22 = qJ(4) + t26;
t21 = t32 * pkin(3) + t23;
t14 = sin(qJ(5));
t20 = rSges(6,1) * t16 - rSges(6,2) * t14;
t12 = cos(pkin(9));
t4 = t12 * t32 + t17 * t29;
t19 = t4 * pkin(4) + t21;
t3 = t17 * t12 - t29 * t32;
t18 = -t3 * pkin(4) + t36;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t17 * rSges(2,1) - rSges(2,2) * t32 + r_base(1)) + g(2) * (rSges(2,1) * t32 + t17 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(1) * (-t17 * rSges(3,2) + rSges(3,3) * t32 + t24) + g(2) * (-t32 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t17 + t27) + g(3) * (rSges(3,1) + t28)) - m(4) * (g(1) * (rSges(4,1) * t32 + t17 * rSges(4,3) + t23) + g(2) * (t32 * rSges(4,3) + (-rSges(4,1) - qJ(2)) * t17 + t25) + g(3) * (-rSges(4,2) + t26)) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t21) + g(2) * (-t3 * rSges(5,1) + t4 * rSges(5,2) + t36) + g(3) * (rSges(5,3) + t22)) - m(6) * (g(1) * t19 + g(2) * t18 + g(3) * (t14 * rSges(6,1) + t16 * rSges(6,2) + t22) + (g(1) * t20 + g(2) * t34) * t4 + (g(1) * t34 - g(2) * t20) * t3) - m(7) * (g(1) * (t4 * t35 - t3 * pkin(7) + (-t3 * t13 + t30 * t4) * rSges(7,1) + (-t3 * t15 - t31 * t4) * rSges(7,2) + t19) + g(2) * (-t3 * t35 - t4 * pkin(7) + (-t4 * t13 - t3 * t30) * rSges(7,1) + (-t4 * t15 + t3 * t31) * rSges(7,2) + t18) + g(3) * (-t37 * t16 + t22) + (g(3) * (rSges(7,1) * t15 - rSges(7,2) * t13 + pkin(5)) + (g(1) * t4 - g(2) * t3) * t37) * t14);
U  = t1;
