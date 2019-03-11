% Calculate potential energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:17
% EndTime: 2019-03-09 02:25:18
% DurationCPUTime: 0.46s
% Computational Cost: add. (204->97), mult. (273->114), div. (0->0), fcn. (305->10), ass. (0->37)
t51 = rSges(6,3) + pkin(8);
t50 = rSges(7,3) + pkin(9) + pkin(8);
t33 = sin(pkin(10));
t34 = cos(pkin(10));
t39 = sin(qJ(1));
t40 = cos(qJ(1));
t3 = -t39 * t33 - t40 * t34;
t4 = t40 * t33 - t39 * t34;
t49 = g(1) * t3 + g(2) * t4;
t17 = sin(qJ(5));
t46 = t3 * t17;
t20 = cos(qJ(4));
t45 = t3 * t20;
t44 = t4 * t17;
t43 = t4 * t20;
t42 = rSges(5,3) + pkin(7);
t16 = qJ(5) + qJ(6);
t10 = cos(t16);
t38 = t10 * t20;
t37 = t17 * t20;
t19 = cos(qJ(5));
t36 = t19 * t20;
t32 = pkin(6) + r_base(3);
t31 = -qJ(3) + t32;
t30 = t40 * pkin(1) + t39 * qJ(2) + r_base(1);
t29 = t40 * pkin(2) + t30;
t18 = sin(qJ(4));
t28 = -rSges(5,1) * t20 + rSges(5,2) * t18;
t27 = -t3 * pkin(3) + t29;
t26 = t39 * pkin(1) - t40 * qJ(2) + r_base(2);
t25 = t4 * pkin(7) + t27;
t24 = t39 * pkin(2) + t26;
t23 = -t4 * pkin(3) + t24;
t22 = -t3 * pkin(7) + t23;
t9 = sin(t16);
t8 = t19 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t39 * rSges(2,2) + r_base(1)) + g(2) * (t39 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t32)) - m(3) * (g(1) * (t40 * rSges(3,1) + t39 * rSges(3,3) + t30) + g(2) * (t39 * rSges(3,1) - t40 * rSges(3,3) + t26) + g(3) * (rSges(3,2) + t32)) - m(4) * (g(1) * (-t3 * rSges(4,1) - t4 * rSges(4,2) + t29) + g(2) * (-t4 * rSges(4,1) + t3 * rSges(4,2) + t24) + g(3) * (-rSges(4,3) + t31)) - m(5) * (g(1) * t27 + g(2) * t23 + g(3) * (-t18 * rSges(5,1) - t20 * rSges(5,2) + t31) + (g(1) * t42 + g(2) * t28) * t4 + (g(1) * t28 - g(2) * t42) * t3) - m(6) * (g(1) * (-pkin(4) * t45 + (-t3 * t36 + t44) * rSges(6,1) + (t4 * t19 + t3 * t37) * rSges(6,2) + t25) + g(2) * (-pkin(4) * t43 + (-t4 * t36 - t46) * rSges(6,1) + (-t3 * t19 + t4 * t37) * rSges(6,2) + t22) + g(3) * (t51 * t20 + t31) + (g(3) * (-rSges(6,1) * t19 + rSges(6,2) * t17 - pkin(4)) - t49 * t51) * t18) - m(7) * (g(1) * (-t8 * t45 + pkin(5) * t44 + (-t3 * t38 + t4 * t9) * rSges(7,1) + (t4 * t10 + t9 * t45) * rSges(7,2) + t25) + g(2) * (-t8 * t43 - pkin(5) * t46 + (-t3 * t9 - t4 * t38) * rSges(7,1) + (-t3 * t10 + t9 * t43) * rSges(7,2) + t22) + g(3) * (t50 * t20 + t31) + (g(3) * (-rSges(7,1) * t10 + rSges(7,2) * t9 - t8) - t49 * t50) * t18);
U  = t1;
