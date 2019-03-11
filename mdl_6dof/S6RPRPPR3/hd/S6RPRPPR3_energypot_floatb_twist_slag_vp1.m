% Calculate potential energy for
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:01
% EndTime: 2019-03-09 02:44:01
% DurationCPUTime: 0.41s
% Computational Cost: add. (213->99), mult. (177->111), div. (0->0), fcn. (153->8), ass. (0->30)
t44 = rSges(7,3) + pkin(8);
t17 = qJ(1) + pkin(9);
t11 = sin(t17);
t19 = sin(qJ(3));
t42 = t11 * t19;
t22 = cos(qJ(3));
t41 = t11 * t22;
t12 = cos(t17);
t40 = t12 * t22;
t18 = sin(qJ(6));
t39 = t18 * t19;
t21 = cos(qJ(6));
t38 = t19 * t21;
t37 = qJ(4) * t19;
t36 = pkin(6) + r_base(3);
t20 = sin(qJ(1));
t35 = t20 * pkin(1) + r_base(2);
t23 = cos(qJ(1));
t34 = t23 * pkin(1) + r_base(1);
t33 = t11 * pkin(2) + t35;
t32 = qJ(2) + t36;
t31 = t12 * pkin(2) + t11 * pkin(7) + t34;
t30 = t19 * pkin(3) + t32;
t29 = pkin(3) * t41 + t11 * t37 + t33;
t28 = rSges(6,1) * t19 - rSges(6,2) * t22;
t27 = t19 * pkin(4) + t30;
t26 = pkin(3) * t40 + t12 * t37 + t31;
t25 = pkin(4) * t41 + t12 * qJ(5) + t29;
t24 = pkin(4) * t40 + t26;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t34) + g(2) * (t11 * rSges(3,1) + t12 * rSges(3,2) + t35) + g(3) * (rSges(3,3) + t32)) - m(4) * (g(1) * (t11 * rSges(4,3) + t31) + g(2) * (rSges(4,1) * t41 - rSges(4,2) * t42 + t33) + g(3) * (t19 * rSges(4,1) + t22 * rSges(4,2) + t32) + (g(1) * (rSges(4,1) * t22 - rSges(4,2) * t19) + g(2) * (-rSges(4,3) - pkin(7))) * t12) - m(5) * (g(1) * (t11 * rSges(5,2) + t26) + g(2) * (rSges(5,1) * t41 + rSges(5,3) * t42 + t29) + g(3) * (t19 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t22 + t30) + (g(1) * (rSges(5,1) * t22 + rSges(5,3) * t19) + g(2) * (-rSges(5,2) - pkin(7))) * t12) - m(6) * (g(1) * t24 + g(2) * t25 + g(3) * (-t19 * rSges(6,2) + (-rSges(6,1) - qJ(4)) * t22 + t27) + (g(1) * t28 + g(2) * (rSges(6,3) - pkin(7))) * t12 + (g(1) * (-rSges(6,3) - qJ(5)) + g(2) * t28) * t11) - m(7) * (g(1) * (t12 * t19 * pkin(5) - t11 * qJ(5) + (-t11 * t18 + t12 * t38) * rSges(7,1) + (-t11 * t21 - t12 * t39) * rSges(7,2) + t24) + g(2) * (pkin(5) * t42 - t12 * pkin(7) + (t11 * t38 + t12 * t18) * rSges(7,1) + (-t11 * t39 + t12 * t21) * rSges(7,2) + t25) + g(3) * (t44 * t19 + t27) + (g(3) * (-rSges(7,1) * t21 + rSges(7,2) * t18 - pkin(5) - qJ(4)) + (g(1) * t12 + g(2) * t11) * t44) * t22);
U  = t1;
