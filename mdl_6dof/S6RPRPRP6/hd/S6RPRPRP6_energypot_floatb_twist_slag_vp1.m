% Calculate potential energy for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:17:52
% EndTime: 2019-03-09 03:17:53
% DurationCPUTime: 0.42s
% Computational Cost: add. (217->103), mult. (204->113), div. (0->0), fcn. (184->8), ass. (0->39)
t52 = rSges(6,3) + pkin(8);
t51 = rSges(7,3) + qJ(6) + pkin(8);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t50 = g(1) * t26 + g(2) * t24;
t18 = pkin(9) + qJ(3);
t15 = sin(t18);
t46 = t15 * t26;
t16 = cos(t18);
t45 = t16 * t26;
t23 = sin(qJ(5));
t44 = t24 * t23;
t25 = cos(qJ(5));
t43 = t24 * t25;
t42 = t26 * t23;
t41 = t26 * t25;
t39 = qJ(4) * t15;
t38 = rSges(3,3) + qJ(2);
t37 = pkin(6) + r_base(3);
t20 = cos(pkin(9));
t12 = t20 * pkin(2) + pkin(1);
t36 = t26 * t12 + r_base(1);
t35 = t15 * t44;
t34 = t15 * t42;
t19 = sin(pkin(9));
t33 = t19 * pkin(2) + t37;
t22 = -pkin(7) - qJ(2);
t32 = t24 * t12 + t26 * t22 + r_base(2);
t31 = pkin(3) * t45 + t26 * t39 + t36;
t30 = t15 * pkin(3) + t33;
t29 = t32 + (pkin(3) * t16 + t39) * t24;
t28 = rSges(3,1) * t20 - rSges(3,2) * t19 + pkin(1);
t27 = -t24 * t22 + t31;
t14 = t25 * pkin(5) + pkin(4);
t4 = t35 - t41;
t3 = t15 * t43 + t42;
t2 = t34 + t43;
t1 = t15 * t41 - t44;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t19 * rSges(3,1) + t20 * rSges(3,2) + t37) + (g(1) * t28 - g(2) * t38) * t26 + (g(1) * t38 + g(2) * t28) * t24) - m(4) * (g(1) * (rSges(4,1) * t45 - rSges(4,2) * t46 + t36) + g(2) * (-t26 * rSges(4,3) + t32) + g(3) * (t15 * rSges(4,1) + t16 * rSges(4,2) + t33) + (g(1) * (rSges(4,3) - t22) + g(2) * (rSges(4,1) * t16 - rSges(4,2) * t15)) * t24) - m(5) * (g(1) * (-rSges(5,2) * t45 + rSges(5,3) * t46 + t31) + g(2) * (-t26 * rSges(5,1) + t29) + g(3) * (-t15 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t16 + t30) + (g(1) * (rSges(5,1) - t22) + g(2) * (-rSges(5,2) * t16 + rSges(5,3) * t15)) * t24) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t24 * pkin(4) + t27) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) - t26 * pkin(4) + t29) + g(3) * (t52 * t15 + t30) + (g(3) * (-rSges(6,1) * t23 - rSges(6,2) * t25 - qJ(4)) + t50 * t52) * t16) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t34 + t24 * t14 + t27) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t35 - t26 * t14 + t29) + g(3) * (t51 * t15 + t30) + (g(3) * (-rSges(7,2) * t25 - qJ(4) + (-rSges(7,1) - pkin(5)) * t23) + t50 * t51) * t16);
U  = t5;
