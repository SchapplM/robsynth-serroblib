% Calculate potential energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:28
% EndTime: 2019-03-09 08:33:29
% DurationCPUTime: 0.45s
% Computational Cost: add. (163->106), mult. (226->118), div. (0->0), fcn. (206->6), ass. (0->36)
t52 = rSges(6,3) + pkin(8);
t51 = rSges(7,3) + qJ(6) + pkin(8);
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t50 = g(1) * t24 + g(2) * t21;
t20 = sin(qJ(2));
t46 = t20 * t21;
t45 = t20 * t24;
t19 = sin(qJ(5));
t44 = t21 * t19;
t22 = cos(qJ(5));
t43 = t21 * t22;
t23 = cos(qJ(2));
t42 = t21 * t23;
t41 = t23 * t24;
t40 = t24 * t19;
t39 = t24 * t22;
t37 = qJ(3) * t20;
t36 = pkin(6) + r_base(3);
t35 = t21 * pkin(1) + r_base(2);
t34 = t20 * pkin(2) + t36;
t33 = t24 * pkin(1) + t21 * pkin(7) + r_base(1);
t32 = pkin(2) * t42 + t21 * t37 + t35;
t31 = t20 * pkin(3) + t34;
t30 = rSges(5,1) * t20 - rSges(5,2) * t23;
t29 = pkin(2) * t41 + t24 * t37 + t33;
t28 = pkin(3) * t42 + t24 * qJ(4) + t32;
t27 = pkin(3) * t41 + t29;
t26 = -t24 * pkin(7) + t28;
t25 = -t21 * qJ(4) + t27;
t11 = t22 * pkin(5) + pkin(4);
t4 = t20 * t39 - t44;
t3 = -t20 * t40 - t43;
t2 = t20 * t43 + t40;
t1 = -t20 * t44 + t39;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t21 * rSges(3,3) + t33) + g(2) * (rSges(3,1) * t42 - rSges(3,2) * t46 + t35) + g(3) * (t20 * rSges(3,1) + t23 * rSges(3,2) + t36) + (g(1) * (rSges(3,1) * t23 - rSges(3,2) * t20) + g(2) * (-rSges(3,3) - pkin(7))) * t24) - m(4) * (g(1) * (t21 * rSges(4,2) + t29) + g(2) * (rSges(4,1) * t42 + rSges(4,3) * t46 + t32) + g(3) * (t20 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t23 + t34) + (g(1) * (rSges(4,1) * t23 + rSges(4,3) * t20) + g(2) * (-rSges(4,2) - pkin(7))) * t24) - m(5) * (g(1) * t27 + g(2) * t28 + g(3) * (-t20 * rSges(5,2) + (-rSges(5,1) - qJ(3)) * t23 + t31) + (g(1) * t30 + g(2) * (rSges(5,3) - pkin(7))) * t24 + (g(1) * (-rSges(5,3) - qJ(4)) + g(2) * t30) * t21) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t45 + t25) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t46 + t26) + g(3) * (t52 * t20 + t31) + (g(3) * (-rSges(6,1) * t22 + rSges(6,2) * t19 - pkin(4) - qJ(3)) + t50 * t52) * t23) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) - pkin(5) * t44 + t11 * t45 + t25) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t40 + t11 * t46 + t26) + g(3) * (t51 * t20 + t31) + (g(3) * (-rSges(7,1) * t22 + rSges(7,2) * t19 - qJ(3) - t11) + t50 * t51) * t23);
U  = t5;
