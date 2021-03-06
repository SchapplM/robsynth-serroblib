% Calculate potential energy for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:32
% EndTime: 2019-03-09 04:42:32
% DurationCPUTime: 0.45s
% Computational Cost: add. (244->102), mult. (248->112), div. (0->0), fcn. (242->8), ass. (0->39)
t52 = rSges(7,1) + pkin(5);
t51 = rSges(7,2) + qJ(5);
t22 = pkin(9) + qJ(3);
t19 = sin(t22);
t27 = sin(qJ(1));
t50 = t19 * t27;
t29 = cos(qJ(1));
t49 = t19 * t29;
t20 = cos(t22);
t48 = t20 * t29;
t26 = sin(qJ(4));
t47 = t27 * t26;
t28 = cos(qJ(4));
t46 = t27 * t28;
t45 = t29 * t26;
t44 = t29 * t28;
t43 = rSges(3,3) + qJ(2);
t42 = rSges(6,3) + qJ(5);
t41 = -rSges(7,3) - qJ(6);
t40 = pkin(6) + r_base(3);
t24 = cos(pkin(9));
t16 = t24 * pkin(2) + pkin(1);
t39 = t29 * t16 + r_base(1);
t23 = sin(pkin(9));
t38 = t23 * pkin(2) + t40;
t25 = -pkin(7) - qJ(2);
t37 = t27 * t16 + t29 * t25 + r_base(2);
t36 = t19 * pkin(3) + t38;
t35 = t27 * t20 * pkin(3) + pkin(8) * t50 + t37;
t34 = rSges(3,1) * t24 - rSges(3,2) * t23 + pkin(1);
t33 = t36 + (pkin(4) * t28 + qJ(5) * t26) * t19;
t4 = t20 * t46 - t45;
t32 = t4 * pkin(4) + t35;
t31 = pkin(3) * t48 + pkin(8) * t49 - t27 * t25 + t39;
t6 = t20 * t44 + t47;
t30 = t6 * pkin(4) + t31;
t5 = t20 * t45 - t46;
t3 = t20 * t47 + t44;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t23 * rSges(3,1) + t24 * rSges(3,2) + t40) + (g(1) * t34 - g(2) * t43) * t29 + (g(1) * t43 + g(2) * t34) * t27) - m(4) * (g(1) * (rSges(4,1) * t48 - rSges(4,2) * t49 + t39) + g(2) * (-t29 * rSges(4,3) + t37) + g(3) * (t19 * rSges(4,1) + t20 * rSges(4,2) + t38) + (g(1) * (rSges(4,3) - t25) + g(2) * (rSges(4,1) * t20 - rSges(4,2) * t19)) * t27) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + rSges(5,3) * t49 + t31) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t50 + t35) + g(3) * ((-rSges(5,3) - pkin(8)) * t20 + (rSges(5,1) * t28 - rSges(5,2) * t26) * t19 + t36)) - m(6) * (g(1) * (t6 * rSges(6,1) + rSges(6,2) * t49 + t42 * t5 + t30) + g(2) * (t4 * rSges(6,1) + rSges(6,2) * t50 + t3 * t42 + t32) + g(3) * ((-rSges(6,2) - pkin(8)) * t20 + (rSges(6,1) * t28 + rSges(6,3) * t26) * t19 + t33)) - m(7) * (g(1) * (t51 * t5 + t52 * t6 + t30) + g(2) * (t51 * t3 + t52 * t4 + t32) + (g(1) * t29 + g(2) * t27) * t19 * t41 + (t33 + (-pkin(8) - t41) * t20 + (rSges(7,2) * t26 + t52 * t28) * t19) * g(3));
U  = t1;
