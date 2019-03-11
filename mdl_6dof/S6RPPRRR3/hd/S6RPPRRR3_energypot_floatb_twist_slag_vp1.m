% Calculate potential energy for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:43
% EndTime: 2019-03-09 02:22:44
% DurationCPUTime: 0.41s
% Computational Cost: add. (214->96), mult. (161->107), div. (0->0), fcn. (141->10), ass. (0->35)
t48 = rSges(6,3) + pkin(8);
t47 = rSges(7,3) + pkin(9) + pkin(8);
t16 = sin(qJ(4));
t19 = cos(qJ(4));
t46 = rSges(5,1) * t16 + rSges(5,2) * t19;
t13 = qJ(1) + pkin(10);
t7 = sin(t13);
t8 = cos(t13);
t45 = -g(1) * t7 + g(2) * t8;
t15 = sin(qJ(5));
t42 = t7 * t15;
t41 = t7 * t16;
t40 = t8 * t15;
t39 = t8 * t16;
t14 = qJ(5) + qJ(6);
t10 = cos(t14);
t35 = t10 * t16;
t34 = t15 * t16;
t18 = cos(qJ(5));
t33 = t16 * t18;
t31 = pkin(6) + r_base(3);
t17 = sin(qJ(1));
t30 = t17 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t29 = t20 * pkin(1) + r_base(1);
t28 = t7 * pkin(2) + t30;
t27 = qJ(2) + t31;
t26 = t8 * pkin(2) + t7 * qJ(3) + t29;
t25 = t7 * pkin(7) + t28;
t24 = pkin(3) + t27;
t23 = t8 * pkin(7) + t26;
t22 = -t8 * qJ(3) + t25;
t9 = sin(t14);
t6 = t18 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(1) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t29) + g(2) * (t7 * rSges(3,1) + t8 * rSges(3,2) + t30) + g(3) * (rSges(3,3) + t27)) - m(4) * (g(1) * (-t8 * rSges(4,2) + t7 * rSges(4,3) + t26) + g(2) * (-t7 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t8 + t28) + g(3) * (rSges(4,1) + t27)) - m(5) * (g(1) * (t46 * t7 + t23) + g(2) * (t7 * rSges(5,3) + t25) + g(3) * (t19 * rSges(5,1) - t16 * rSges(5,2) + t24) + (g(1) * rSges(5,3) + g(2) * (-qJ(3) - t46)) * t8) - m(6) * (g(1) * (pkin(4) * t41 + (t33 * t7 + t40) * rSges(6,1) + (t8 * t18 - t34 * t7) * rSges(6,2) + t23) + g(2) * (-pkin(4) * t39 + (-t33 * t8 + t42) * rSges(6,1) + (t7 * t18 + t34 * t8) * rSges(6,2) + t22) + g(3) * (t48 * t16 + t24) + (g(3) * (rSges(6,1) * t18 - rSges(6,2) * t15 + pkin(4)) + t45 * t48) * t19) - m(7) * (g(1) * (t6 * t41 + pkin(5) * t40 + (t35 * t7 + t8 * t9) * rSges(7,1) + (t8 * t10 - t41 * t9) * rSges(7,2) + t23) + g(2) * (-t6 * t39 + pkin(5) * t42 + (-t35 * t8 + t7 * t9) * rSges(7,1) + (t7 * t10 + t39 * t9) * rSges(7,2) + t22) + g(3) * (t47 * t16 + t24) + (g(3) * (rSges(7,1) * t10 - rSges(7,2) * t9 + t6) + t45 * t47) * t19);
U  = t1;
