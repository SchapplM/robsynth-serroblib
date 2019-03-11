% Calculate potential energy for
% S6RRRPRP4
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:28
% EndTime: 2019-03-09 16:43:29
% DurationCPUTime: 0.42s
% Computational Cost: add. (222->100), mult. (214->111), div. (0->0), fcn. (198->8), ass. (0->39)
t52 = rSges(7,1) + pkin(5);
t51 = rSges(7,3) + qJ(6);
t50 = rSges(3,3) + pkin(7);
t22 = qJ(2) + qJ(3);
t18 = sin(t22);
t28 = cos(qJ(1));
t49 = t18 * t28;
t19 = cos(t22);
t25 = sin(qJ(1));
t48 = t19 * t25;
t47 = t19 * t28;
t23 = sin(qJ(5));
t46 = t25 * t23;
t26 = cos(qJ(5));
t45 = t25 * t26;
t44 = t28 * t23;
t43 = t28 * t26;
t42 = qJ(4) * t18;
t41 = pkin(6) + r_base(3);
t27 = cos(qJ(2));
t16 = t27 * pkin(2) + pkin(1);
t40 = t28 * t16 + r_base(1);
t24 = sin(qJ(2));
t39 = t24 * pkin(2) + t41;
t29 = -pkin(8) - pkin(7);
t38 = t25 * t16 + t28 * t29 + r_base(2);
t37 = pkin(3) * t47 + t28 * t42 + t40;
t36 = t18 * pkin(3) + t39;
t35 = g(1) * t28 + g(2) * t25;
t34 = pkin(3) * t48 + t25 * t42 + t38;
t33 = t18 * pkin(9) + t36;
t32 = rSges(3,1) * t27 - rSges(3,2) * t24 + pkin(1);
t31 = -t28 * pkin(4) + pkin(9) * t48 + t34;
t30 = pkin(9) * t47 + t37 + (pkin(4) - t29) * t25;
t4 = t18 * t46 - t43;
t3 = t18 * t45 + t44;
t2 = t18 * t44 + t45;
t1 = -t18 * t43 + t46;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t24 * rSges(3,1) + t27 * rSges(3,2) + t41) + (g(1) * t32 - g(2) * t50) * t28 + (g(1) * t50 + g(2) * t32) * t25) - m(4) * (g(1) * (rSges(4,1) * t47 - rSges(4,2) * t49 + t40) + g(2) * (-t28 * rSges(4,3) + t38) + g(3) * (t18 * rSges(4,1) + t19 * rSges(4,2) + t39) + (g(1) * (rSges(4,3) - t29) + g(2) * (rSges(4,1) * t19 - rSges(4,2) * t18)) * t25) - m(5) * (g(1) * (-rSges(5,2) * t47 + rSges(5,3) * t49 + t37) + g(2) * (-t28 * rSges(5,1) + t34) + g(3) * (-t18 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t19 + t36) + (g(1) * (rSges(5,1) - t29) + g(2) * (-rSges(5,2) * t19 + rSges(5,3) * t18)) * t25) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t30) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t31) + g(3) * (t18 * rSges(6,3) + t33) + (g(3) * (-rSges(6,1) * t23 - rSges(6,2) * t26 - qJ(4)) + t35 * rSges(6,3)) * t19) - m(7) * (g(1) * (t51 * t1 + t52 * t2 + t30) + g(2) * (-t51 * t3 + t52 * t4 + t31) + g(3) * (t18 * rSges(7,2) + t33) + (g(3) * (-t52 * t23 + t51 * t26 - qJ(4)) + t35 * rSges(7,2)) * t19);
U  = t5;
