% Calculate potential energy for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:18
% EndTime: 2019-03-09 15:28:18
% DurationCPUTime: 0.45s
% Computational Cost: add. (210->103), mult. (191->114), div. (0->0), fcn. (167->8), ass. (0->34)
t47 = rSges(7,3) + pkin(9);
t24 = -pkin(8) - pkin(7);
t46 = -qJ(5) - t24;
t45 = rSges(3,3) + pkin(7);
t17 = qJ(2) + qJ(3);
t14 = cos(t17);
t20 = sin(qJ(1));
t43 = t20 * t14;
t18 = sin(qJ(6));
t42 = t20 * t18;
t21 = cos(qJ(6));
t41 = t20 * t21;
t13 = sin(t17);
t23 = cos(qJ(1));
t40 = t23 * t13;
t39 = t23 * t14;
t38 = t23 * t18;
t37 = t23 * t21;
t36 = qJ(4) * t13;
t35 = pkin(6) + r_base(3);
t22 = cos(qJ(2));
t11 = t22 * pkin(2) + pkin(1);
t34 = t23 * t11 + r_base(1);
t19 = sin(qJ(2));
t33 = t19 * pkin(2) + t35;
t32 = t20 * t11 + t23 * t24 + r_base(2);
t31 = pkin(3) * t39 + t23 * t36 + t34;
t30 = t13 * pkin(3) + t33;
t29 = pkin(4) * t39 + t31;
t28 = pkin(3) * t43 + t20 * t36 + t32;
t27 = t13 * pkin(4) + t30;
t26 = rSges(3,1) * t22 - rSges(3,2) * t19 + pkin(1);
t25 = pkin(4) * t43 + t23 * qJ(5) + t28;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t19 * rSges(3,1) + t22 * rSges(3,2) + t35) + (g(1) * t26 - g(2) * t45) * t23 + (g(1) * t45 + g(2) * t26) * t20) - m(4) * (g(1) * (rSges(4,1) * t39 - rSges(4,2) * t40 + t34) + g(2) * (-t23 * rSges(4,3) + t32) + g(3) * (t13 * rSges(4,1) + t14 * rSges(4,2) + t33) + (g(1) * (rSges(4,3) - t24) + g(2) * (rSges(4,1) * t14 - rSges(4,2) * t13)) * t20) - m(5) * (g(1) * (rSges(5,1) * t39 + rSges(5,3) * t40 + t31) + g(2) * (-t23 * rSges(5,2) + t28) + g(3) * (t13 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t14 + t30) + (g(1) * (rSges(5,2) - t24) + g(2) * (rSges(5,1) * t14 + rSges(5,3) * t13)) * t20) - m(6) * (g(1) * (rSges(6,1) * t40 - rSges(6,2) * t39 + t29) + g(2) * (t23 * rSges(6,3) + t25) + g(3) * (-t13 * rSges(6,2) + (-rSges(6,1) - qJ(4)) * t14 + t27) + (g(1) * (-rSges(6,3) + t46) + g(2) * (rSges(6,1) * t13 - rSges(6,2) * t14)) * t20) - m(7) * (g(1) * (pkin(5) * t40 + (t13 * t37 - t42) * rSges(7,1) + (-t13 * t38 - t41) * rSges(7,2) + t29 + t46 * t20) + g(2) * (t20 * t13 * pkin(5) + (t13 * t41 + t38) * rSges(7,1) + (-t13 * t42 + t37) * rSges(7,2) + t25) + g(3) * (t47 * t13 + t27) + (g(3) * (-rSges(7,1) * t21 + rSges(7,2) * t18 - pkin(5) - qJ(4)) + (g(1) * t23 + g(2) * t20) * t47) * t14);
U  = t1;
