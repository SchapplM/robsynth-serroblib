% Calculate potential energy for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:20
% EndTime: 2019-03-09 05:04:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (271->107), mult. (260->125), div. (0->0), fcn. (264->10), ass. (0->35)
t49 = -rSges(7,3) - pkin(9);
t22 = qJ(1) + pkin(10);
t17 = sin(t22);
t25 = sin(qJ(3));
t48 = t17 * t25;
t29 = cos(qJ(3));
t47 = t17 * t29;
t18 = cos(t22);
t46 = t18 * t25;
t24 = sin(qJ(4));
t45 = t24 * t29;
t28 = cos(qJ(4));
t44 = t28 * t29;
t43 = rSges(6,3) + qJ(5);
t42 = pkin(6) + r_base(3);
t26 = sin(qJ(1));
t41 = t26 * pkin(1) + r_base(2);
t30 = cos(qJ(1));
t40 = t30 * pkin(1) + r_base(1);
t39 = qJ(2) + t42;
t38 = t17 * pkin(2) + t41;
t37 = t25 * pkin(3) + t39;
t36 = t18 * pkin(2) + t17 * pkin(7) + t40;
t35 = t18 * t29 * pkin(3) + pkin(8) * t46 + t36;
t34 = t37 + (pkin(4) * t28 + qJ(5) * t24) * t25;
t6 = t17 * t24 + t18 * t44;
t33 = t6 * pkin(4) + t35;
t32 = pkin(3) * t47 - t18 * pkin(7) + pkin(8) * t48 + t38;
t4 = t17 * t44 - t18 * t24;
t31 = t4 * pkin(4) + t32;
t27 = cos(qJ(6));
t23 = sin(qJ(6));
t5 = -t17 * t28 + t18 * t45;
t3 = t17 * t45 + t18 * t28;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t30 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t30 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t42)) - m(3) * (g(1) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t40) + g(2) * (t17 * rSges(3,1) + t18 * rSges(3,2) + t41) + g(3) * (rSges(3,3) + t39)) - m(4) * (g(1) * (t17 * rSges(4,3) + t36) + g(2) * (rSges(4,1) * t47 - rSges(4,2) * t48 + t38) + g(3) * (t25 * rSges(4,1) + t29 * rSges(4,2) + t39) + (g(1) * (rSges(4,1) * t29 - rSges(4,2) * t25) + g(2) * (-rSges(4,3) - pkin(7))) * t18) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + rSges(5,3) * t46 + t35) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t48 + t32) + g(3) * ((-rSges(5,3) - pkin(8)) * t29 + (rSges(5,1) * t28 - rSges(5,2) * t24) * t25 + t37)) - m(6) * (g(1) * (t6 * rSges(6,1) + rSges(6,2) * t46 + t43 * t5 + t33) + g(2) * (t4 * rSges(6,1) + rSges(6,2) * t48 + t43 * t3 + t31) + g(3) * ((-rSges(6,2) - pkin(8)) * t29 + (rSges(6,1) * t28 + rSges(6,3) * t24) * t25 + t34)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t5 * t23 + t6 * t27) * rSges(7,1) + (-t6 * t23 + t5 * t27) * rSges(7,2) + t33) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t3 * t23 + t4 * t27) * rSges(7,1) + (-t4 * t23 + t3 * t27) * rSges(7,2) + t31) + (g(1) * t18 + g(2) * t17) * t25 * t49 + (t34 + (-pkin(8) - t49) * t29 + (t28 * pkin(5) + (t23 * t24 + t27 * t28) * rSges(7,1) + (-t23 * t28 + t24 * t27) * rSges(7,2)) * t25) * g(3));
U  = t1;
