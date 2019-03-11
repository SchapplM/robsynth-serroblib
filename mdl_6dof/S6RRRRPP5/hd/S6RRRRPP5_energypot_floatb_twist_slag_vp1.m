% Calculate potential energy for
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:10
% EndTime: 2019-03-09 21:05:11
% DurationCPUTime: 0.47s
% Computational Cost: add. (245->107), mult. (280->122), div. (0->0), fcn. (278->8), ass. (0->35)
t53 = rSges(7,1) + pkin(5);
t52 = rSges(4,3) + pkin(8);
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t51 = g(1) * t27 + g(2) * t24;
t50 = rSges(7,3) + qJ(6);
t23 = sin(qJ(2));
t46 = rSges(3,2) * t23;
t22 = sin(qJ(3));
t45 = t22 * t24;
t44 = t22 * t27;
t26 = cos(qJ(2));
t43 = t24 * t26;
t42 = t26 * t27;
t39 = pkin(6) + r_base(3);
t38 = t24 * pkin(1) + r_base(2);
t36 = t27 * pkin(1) + t24 * pkin(7) + r_base(1);
t25 = cos(qJ(3));
t14 = pkin(3) * t25 + pkin(2);
t28 = -pkin(9) - pkin(8);
t35 = t23 * t14 + t26 * t28 + t39;
t34 = -t27 * pkin(7) + t38;
t33 = pkin(3) * t45 + t14 * t42 + t36;
t21 = qJ(3) + qJ(4);
t16 = sin(t21);
t17 = cos(t21);
t32 = t35 + (pkin(4) * t17 + qJ(5) * t16) * t23;
t5 = t16 * t42 - t24 * t17;
t6 = t16 * t24 + t17 * t42;
t31 = t6 * pkin(4) + t5 * qJ(5) + t33;
t30 = -pkin(3) * t44 + t14 * t43 + t34;
t3 = t16 * t43 + t17 * t27;
t4 = -t16 * t27 + t17 * t43;
t29 = t4 * pkin(4) + t3 * qJ(5) + t30;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t27 - rSges(2,2) * t24 + r_base(1)) + g(2) * (rSges(2,1) * t24 + rSges(2,2) * t27 + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * (t24 * rSges(3,3) + t36) + g(2) * (rSges(3,1) * t43 - t24 * t46 + t38) + g(3) * (rSges(3,1) * t23 + rSges(3,2) * t26 + t39) + (g(1) * (rSges(3,1) * t26 - t46) + g(2) * (-rSges(3,3) - pkin(7))) * t27) - m(4) * (g(1) * (pkin(2) * t42 + (t25 * t42 + t45) * rSges(4,1) + (-t22 * t42 + t24 * t25) * rSges(4,2) + t36) + g(2) * (pkin(2) * t43 + (t25 * t43 - t44) * rSges(4,1) + (-t22 * t43 - t25 * t27) * rSges(4,2) + t34) + g(3) * (-t52 * t26 + t39) + (g(3) * (rSges(4,1) * t25 - rSges(4,2) * t22 + pkin(2)) + t51 * t52) * t23) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t33) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t30) + g(3) * (-t26 * rSges(5,3) + t35) + (g(3) * (rSges(5,1) * t17 - rSges(5,2) * t16) + t51 * (rSges(5,3) - t28)) * t23) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,3) + t31) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,3) + t29) + g(3) * (-t26 * rSges(6,2) + t32) + (g(3) * (rSges(6,1) * t17 + rSges(6,3) * t16) + t51 * (rSges(6,2) - t28)) * t23) - m(7) * (g(1) * (t5 * rSges(7,2) + t53 * t6 + t31) + g(2) * (t3 * rSges(7,2) + t53 * t4 + t29) + g(3) * (t50 * t26 + t32) + (g(3) * (rSges(7,2) * t16 + t53 * t17) + t51 * (-t28 - t50)) * t23);
U  = t1;
