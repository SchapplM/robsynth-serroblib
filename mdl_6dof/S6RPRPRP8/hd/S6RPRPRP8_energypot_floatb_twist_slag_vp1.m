% Calculate potential energy for
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:14
% EndTime: 2019-03-09 03:24:15
% DurationCPUTime: 0.43s
% Computational Cost: add. (192->98), mult. (193->105), div. (0->0), fcn. (177->8), ass. (0->41)
t50 = rSges(7,1) + pkin(5);
t49 = rSges(7,3) + qJ(6);
t20 = sin(qJ(3));
t48 = pkin(3) * t20;
t17 = qJ(3) + pkin(9);
t11 = sin(t17);
t47 = pkin(4) * t11;
t18 = -qJ(4) - pkin(7);
t46 = g(1) * t18;
t21 = sin(qJ(1));
t45 = g(1) * t21;
t44 = rSges(4,3) + pkin(7);
t19 = sin(qJ(5));
t43 = t21 * t19;
t22 = cos(qJ(5));
t42 = t21 * t22;
t24 = cos(qJ(1));
t41 = t24 * t19;
t40 = t24 * t22;
t39 = rSges(5,3) - t18;
t38 = pkin(6) + r_base(3);
t37 = t21 * pkin(1) + r_base(2);
t36 = pkin(2) + t38;
t35 = -qJ(2) - t48;
t34 = t24 * pkin(1) + t21 * qJ(2) + r_base(1);
t33 = g(2) * t37;
t23 = cos(qJ(3));
t32 = t23 * pkin(3) + t36;
t31 = t21 * t48 + t34;
t30 = rSges(4,1) * t20 + rSges(4,2) * t23;
t12 = cos(t17);
t29 = rSges(5,1) * t11 + rSges(5,2) * t12;
t28 = t21 * t47 + t31;
t27 = t12 * pkin(4) + t11 * pkin(8) + t32;
t26 = t35 - t47;
t25 = t24 * t12 * pkin(8) - t21 * t18 + t37;
t4 = -t11 * t40 + t43;
t3 = t11 * t41 + t42;
t2 = t11 * t42 + t41;
t1 = t11 * t43 - t40;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (-t24 * rSges(3,2) + t21 * rSges(3,3) + t34) + g(2) * (-t21 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t24 + t37) + g(3) * (rSges(3,1) + t38)) - m(4) * (g(1) * t34 + t33 + g(3) * (t23 * rSges(4,1) - t20 * rSges(4,2) + t36) + (g(1) * t30 + g(2) * t44) * t21 + (g(1) * t44 + g(2) * (-qJ(2) - t30)) * t24) - m(5) * (g(1) * t31 + t33 + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t32) + (g(1) * t29 + g(2) * t39) * t21 + (g(1) * t39 + g(2) * (-t29 + t35)) * t24) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t28) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t25) + g(3) * (t11 * rSges(6,3) + t27) + (g(3) * (rSges(6,1) * t22 - rSges(6,2) * t19) + (-rSges(6,3) - pkin(8)) * t45) * t12 + (-t46 + g(2) * (rSges(6,3) * t12 + t26)) * t24) - m(7) * (g(1) * (t49 * t1 + t50 * t2 + t28) + g(2) * (-t49 * t3 + t50 * t4 + t25) + g(3) * (t11 * rSges(7,2) + t27) + (g(2) * t26 - t46) * t24 + (g(2) * rSges(7,2) * t24 + g(3) * (t49 * t19 + t50 * t22) + (-rSges(7,2) - pkin(8)) * t45) * t12);
U  = t5;
