% Calculate potential energy for
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:09
% EndTime: 2019-03-09 17:14:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (215->115), mult. (362->137), div. (0->0), fcn. (388->8), ass. (0->44)
t29 = sin(qJ(1));
t33 = cos(qJ(1));
t59 = g(1) * t33 + g(2) * t29;
t32 = cos(qJ(2));
t56 = g(3) * t32;
t26 = sin(qJ(5));
t27 = sin(qJ(3));
t31 = cos(qJ(3));
t47 = t33 * t31;
t49 = t29 * t32;
t9 = t27 * t49 + t47;
t55 = t9 * t26;
t54 = -rSges(6,3) - pkin(9);
t48 = t33 * t27;
t11 = -t29 * t31 + t32 * t48;
t53 = t11 * t26;
t52 = t26 * t27;
t28 = sin(qJ(2));
t51 = t28 * t29;
t50 = t28 * t33;
t46 = -rSges(7,3) - qJ(6) - pkin(9);
t45 = rSges(5,3) + qJ(4);
t44 = pkin(6) + r_base(3);
t43 = t29 * pkin(1) + r_base(2);
t42 = t28 * pkin(2) + t44;
t41 = t33 * pkin(1) + t29 * pkin(7) + r_base(1);
t40 = t42 + (pkin(3) * t31 + qJ(4) * t27) * t28;
t39 = t33 * t32 * pkin(2) + pkin(8) * t50 + t41;
t12 = t29 * t27 + t32 * t47;
t38 = t12 * pkin(3) + t39;
t37 = pkin(2) * t49 - t33 * pkin(7) + pkin(8) * t51 + t43;
t10 = t31 * t49 - t48;
t36 = t10 * pkin(3) + t37;
t35 = t11 * qJ(4) + t38;
t34 = t9 * qJ(4) + t36;
t30 = cos(qJ(5));
t20 = t30 * pkin(5) + pkin(4);
t6 = (t30 * t31 + t52) * t28;
t5 = (-t26 * t31 + t27 * t30) * t28;
t4 = t12 * t30 + t53;
t3 = t11 * t30 - t12 * t26;
t2 = t10 * t30 + t55;
t1 = -t10 * t26 + t9 * t30;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t33 * rSges(2,1) - t29 * rSges(2,2) + r_base(1)) + g(2) * (t29 * rSges(2,1) + t33 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * (t29 * rSges(3,3) + t41) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t51 + t43) + g(3) * (t28 * rSges(3,1) + t32 * rSges(3,2) + t44) + (g(1) * (rSges(3,1) * t32 - rSges(3,2) * t28) + g(2) * (-rSges(3,3) - pkin(7))) * t33) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + rSges(4,3) * t50 + t39) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + rSges(4,3) * t51 + t37) + g(3) * ((-rSges(4,3) - pkin(8)) * t32 + (rSges(4,1) * t31 - rSges(4,2) * t27) * t28 + t42)) - m(5) * (g(1) * (t12 * rSges(5,1) + rSges(5,2) * t50 + t45 * t11 + t38) + g(2) * (t10 * rSges(5,1) + rSges(5,2) * t51 + t45 * t9 + t36) + g(3) * ((-rSges(5,2) - pkin(8)) * t32 + (rSges(5,1) * t31 + rSges(5,3) * t27) * t28 + t40)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t12 * pkin(4) + t35) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t10 * pkin(4) + t34) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t40) + (-pkin(8) - t54) * t56 + (g(3) * pkin(4) * t31 + t59 * t54) * t28) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t53 + t12 * t20 + t35) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t55 + t10 * t20 + t34) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t40) + (-pkin(8) - t46) * t56 + (g(3) * (pkin(5) * t52 + t20 * t31) + t59 * t46) * t28);
U  = t7;
