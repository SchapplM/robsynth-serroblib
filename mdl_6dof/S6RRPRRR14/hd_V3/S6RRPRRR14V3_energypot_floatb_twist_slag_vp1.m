% Calculate potential energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14V3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:17
% EndTime: 2019-04-12 15:03:17
% DurationCPUTime: 0.27s
% Computational Cost: add. (124->90), mult. (227->117), div. (0->0), fcn. (240->10), ass. (0->34)
t15 = sin(qJ(4));
t16 = sin(qJ(2));
t35 = t15 * t16;
t17 = sin(qJ(1));
t34 = t16 * t17;
t20 = cos(qJ(4));
t33 = t16 * t20;
t22 = cos(qJ(1));
t32 = t16 * t22;
t21 = cos(qJ(2));
t31 = t17 * t21;
t30 = t22 * t15;
t29 = t22 * t20;
t28 = qJ(3) * t16;
t27 = t17 * t28 + r_base(2);
t26 = t22 * t28 + r_base(1);
t25 = -t21 * qJ(3) + r_base(3);
t24 = rSges(3,1) * t21 - rSges(3,2) * t16;
t23 = rSges(4,1) * t21 + rSges(4,3) * t16;
t19 = cos(qJ(5));
t18 = cos(qJ(6));
t14 = sin(qJ(5));
t13 = sin(qJ(6));
t10 = t17 * t15 + t21 * t29;
t9 = -t17 * t20 + t21 * t30;
t8 = t20 * t31 - t30;
t7 = t15 * t31 + t29;
t6 = -t21 * t14 + t19 * t33;
t5 = t14 * t33 + t21 * t19;
t4 = t10 * t19 + t14 * t32;
t3 = t10 * t14 - t19 * t32;
t2 = t14 * t34 + t8 * t19;
t1 = t8 * t14 - t19 * t34;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t22 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t22 * rSges(2,2) + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (t17 * rSges(3,3) + t24 * t22 + r_base(1)) + g(2) * (-t22 * rSges(3,3) + t24 * t17 + r_base(2)) + g(3) * (t16 * rSges(3,1) + t21 * rSges(3,2) + r_base(3))) - m(4) * (g(1) * (t17 * rSges(4,2) + t23 * t22 + t26) + g(2) * (-t22 * rSges(4,2) + t23 * t17 + t27) + g(3) * (t16 * rSges(4,1) + r_base(3) + (-rSges(4,3) - qJ(3)) * t21)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + rSges(5,3) * t32 + t26) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + rSges(5,3) * t34 + t27) + g(3) * (r_base(3) + (-rSges(5,3) - qJ(3)) * t21 + (rSges(5,1) * t20 - rSges(5,2) * t15) * t16)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t9 * rSges(6,3) + t26) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t7 * rSges(6,3) + t27) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + rSges(6,3) * t35 + t25)) - m(7) * (g(1) * ((t9 * t13 + t4 * t18) * rSges(7,1) + (-t4 * t13 + t9 * t18) * rSges(7,2) + t3 * rSges(7,3) + t26) + g(2) * ((t7 * t13 + t2 * t18) * rSges(7,1) + (-t2 * t13 + t7 * t18) * rSges(7,2) + t1 * rSges(7,3) + t27) + g(3) * ((t13 * t35 + t6 * t18) * rSges(7,1) + (-t6 * t13 + t18 * t35) * rSges(7,2) + t5 * rSges(7,3) + t25));
U  = t11;
