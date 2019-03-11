% Calculate potential energy for
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:55
% EndTime: 2019-03-09 22:11:56
% DurationCPUTime: 0.57s
% Computational Cost: add. (260->111), mult. (274->129), div. (0->0), fcn. (278->10), ass. (0->39)
t52 = rSges(3,3) + pkin(7);
t51 = -rSges(7,3) - pkin(10);
t22 = qJ(2) + qJ(3);
t20 = cos(t22);
t30 = cos(qJ(1));
t50 = t20 * t30;
t19 = sin(t22);
t26 = sin(qJ(1));
t49 = t26 * t19;
t24 = sin(qJ(4));
t48 = t26 * t24;
t28 = cos(qJ(4));
t47 = t26 * t28;
t46 = t30 * t19;
t45 = t30 * t24;
t44 = t30 * t28;
t43 = rSges(6,3) + qJ(5);
t42 = pkin(6) + r_base(3);
t29 = cos(qJ(2));
t17 = pkin(2) * t29 + pkin(1);
t41 = t17 * t30 + r_base(1);
t25 = sin(qJ(2));
t40 = pkin(2) * t25 + t42;
t31 = -pkin(8) - pkin(7);
t39 = t17 * t26 + t30 * t31 + r_base(2);
t38 = pkin(3) * t19 + t40;
t37 = pkin(3) * t20 * t26 + pkin(9) * t49 + t39;
t36 = rSges(3,1) * t29 - rSges(3,2) * t25 + pkin(1);
t35 = t38 + (pkin(4) * t28 + qJ(5) * t24) * t19;
t4 = t20 * t47 - t45;
t34 = pkin(4) * t4 + t37;
t33 = pkin(3) * t50 + pkin(9) * t46 - t26 * t31 + t41;
t6 = t20 * t44 + t48;
t32 = pkin(4) * t6 + t33;
t27 = cos(qJ(6));
t23 = sin(qJ(6));
t5 = t20 * t45 - t47;
t3 = t20 * t48 + t44;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - rSges(2,2) * t26 + r_base(1)) + g(2) * (rSges(2,1) * t26 + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t42)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t25 + rSges(3,2) * t29 + t42) + (g(1) * t36 - g(2) * t52) * t30 + (g(1) * t52 + g(2) * t36) * t26) - m(4) * (g(1) * (rSges(4,1) * t50 - rSges(4,2) * t46 + t41) + g(2) * (-t30 * rSges(4,3) + t39) + g(3) * (rSges(4,1) * t19 + rSges(4,2) * t20 + t40) + (g(1) * (rSges(4,3) - t31) + g(2) * (rSges(4,1) * t20 - rSges(4,2) * t19)) * t26) - m(5) * (g(1) * (rSges(5,1) * t6 - rSges(5,2) * t5 + rSges(5,3) * t46 + t33) + g(2) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t49 + t37) + g(3) * ((-rSges(5,3) - pkin(9)) * t20 + (rSges(5,1) * t28 - rSges(5,2) * t24) * t19 + t38)) - m(6) * (g(1) * (t6 * rSges(6,1) + rSges(6,2) * t46 + t43 * t5 + t32) + g(2) * (t4 * rSges(6,1) + rSges(6,2) * t49 + t3 * t43 + t34) + g(3) * ((-rSges(6,2) - pkin(9)) * t20 + (rSges(6,1) * t28 + rSges(6,3) * t24) * t19 + t35)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t23 * t5 + t27 * t6) * rSges(7,1) + (-t23 * t6 + t27 * t5) * rSges(7,2) + t32) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t23 * t3 + t27 * t4) * rSges(7,1) + (-t23 * t4 + t27 * t3) * rSges(7,2) + t34) + (g(1) * t30 + g(2) * t26) * t19 * t51 + (t35 + (-pkin(9) - t51) * t20 + (t28 * pkin(5) + (t23 * t24 + t27 * t28) * rSges(7,1) + (-t23 * t28 + t24 * t27) * rSges(7,2)) * t19) * g(3));
U  = t1;
