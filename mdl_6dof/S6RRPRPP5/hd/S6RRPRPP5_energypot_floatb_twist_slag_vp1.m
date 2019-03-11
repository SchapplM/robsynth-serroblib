% Calculate potential energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:24
% EndTime: 2019-03-09 10:03:24
% DurationCPUTime: 0.44s
% Computational Cost: add. (176->103), mult. (267->115), div. (0->0), fcn. (261->6), ass. (0->34)
t52 = rSges(7,1) + pkin(5);
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t33 = g(1) * t26 + g(2) * t23;
t51 = -rSges(7,3) - qJ(6);
t22 = sin(qJ(2));
t48 = t22 * t23;
t21 = sin(qJ(4));
t47 = t23 * t21;
t24 = cos(qJ(4));
t46 = t23 * t24;
t25 = cos(qJ(2));
t45 = t23 * t25;
t44 = t25 * t26;
t43 = t26 * t21;
t42 = t26 * t24;
t41 = qJ(3) * t22;
t39 = pkin(6) + r_base(3);
t38 = t23 * pkin(1) + r_base(2);
t37 = t22 * pkin(2) + t39;
t36 = t26 * pkin(1) + t23 * pkin(7) + r_base(1);
t35 = t22 * pkin(8) + t37;
t34 = pkin(2) * t45 + t23 * t41 + t38;
t32 = t25 * t24 * qJ(5) + t35;
t31 = pkin(2) * t44 + t26 * t41 + t36;
t30 = t23 * pkin(3) + pkin(8) * t44 + t31;
t29 = pkin(8) * t45 + t34 + (-pkin(3) - pkin(7)) * t26;
t3 = -t22 * t42 + t47;
t4 = t22 * t43 + t46;
t28 = t4 * pkin(4) + t3 * qJ(5) + t30;
t5 = t22 * t46 + t43;
t6 = t22 * t47 - t42;
t27 = t6 * pkin(4) - t5 * qJ(5) + t29;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * (t23 * rSges(3,3) + t36) + g(2) * (rSges(3,1) * t45 - rSges(3,2) * t48 + t38) + g(3) * (t22 * rSges(3,1) + t25 * rSges(3,2) + t39) + (g(1) * (rSges(3,1) * t25 - rSges(3,2) * t22) + g(2) * (-rSges(3,3) - pkin(7))) * t26) - m(4) * (g(1) * (t23 * rSges(4,1) + t31) + g(2) * (-rSges(4,2) * t45 + rSges(4,3) * t48 + t34) + g(3) * (-t22 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t25 + t37) + (g(1) * (-rSges(4,2) * t25 + rSges(4,3) * t22) + g(2) * (-rSges(4,1) - pkin(7))) * t26) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t30) + g(2) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t29) + g(3) * (t22 * rSges(5,3) + t35) + (g(3) * (-rSges(5,1) * t21 - rSges(5,2) * t24 - qJ(3)) + t33 * rSges(5,3)) * t25) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,3) + t28) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,3) + t27) + g(3) * (t22 * rSges(6,2) + t32) + (g(3) * (rSges(6,3) * t24 - qJ(3) + (-rSges(6,1) - pkin(4)) * t21) + t33 * rSges(6,2)) * t25) - m(7) * (g(1) * (t3 * rSges(7,2) + t52 * t4 + t28) + g(2) * (-t5 * rSges(7,2) + t52 * t6 + t27) + t33 * t25 * t51 + (t32 + (rSges(7,2) * t24 - qJ(3) + (-pkin(4) - t52) * t21) * t25 + t51 * t22) * g(3));
U  = t1;
