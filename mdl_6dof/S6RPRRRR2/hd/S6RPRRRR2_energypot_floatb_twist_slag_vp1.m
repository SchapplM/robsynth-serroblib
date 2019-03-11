% Calculate potential energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:55
% EndTime: 2019-03-09 06:56:56
% DurationCPUTime: 0.48s
% Computational Cost: add. (255->100), mult. (172->113), div. (0->0), fcn. (152->12), ass. (0->40)
t52 = rSges(6,3) + pkin(9);
t51 = rSges(7,3) + pkin(10) + pkin(9);
t17 = qJ(3) + qJ(4);
t11 = cos(t17);
t9 = sin(t17);
t50 = rSges(5,1) * t11 - rSges(5,2) * t9;
t15 = qJ(1) + pkin(11);
t6 = sin(t15);
t7 = cos(t15);
t49 = g(1) * t7 + g(2) * t6;
t45 = t6 * t11;
t18 = sin(qJ(5));
t44 = t6 * t18;
t43 = t7 * t11;
t42 = t7 * t18;
t41 = rSges(4,3) + pkin(7);
t16 = qJ(5) + qJ(6);
t10 = cos(t16);
t38 = t10 * t11;
t37 = t11 * t18;
t21 = cos(qJ(5));
t36 = t11 * t21;
t34 = pkin(6) + r_base(3);
t20 = sin(qJ(1));
t33 = t20 * pkin(1) + r_base(2);
t23 = cos(qJ(1));
t32 = t23 * pkin(1) + r_base(1);
t22 = cos(qJ(3));
t5 = t22 * pkin(3) + pkin(2);
t31 = t7 * t5 + t32;
t30 = qJ(2) + t34;
t25 = -pkin(8) - pkin(7);
t29 = t7 * t25 + t6 * t5 + t33;
t19 = sin(qJ(3));
t28 = t19 * pkin(3) + t30;
t27 = rSges(4,1) * t22 - rSges(4,2) * t19 + pkin(2);
t26 = -t6 * t25 + t31;
t8 = sin(t16);
t4 = t21 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t7 * rSges(3,1) - t6 * rSges(3,2) + t32) + g(2) * (t6 * rSges(3,1) + t7 * rSges(3,2) + t33) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * t32 + g(2) * t33 + g(3) * (t19 * rSges(4,1) + t22 * rSges(4,2) + t30) + (g(1) * t27 - g(2) * t41) * t7 + (g(1) * t41 + g(2) * t27) * t6) - m(5) * (g(1) * (t50 * t7 + t31) + g(2) * (-t7 * rSges(5,3) + t29) + g(3) * (t9 * rSges(5,1) + t11 * rSges(5,2) + t28) + (g(1) * (rSges(5,3) - t25) + g(2) * t50) * t6) - m(6) * (g(1) * (pkin(4) * t43 + (t7 * t36 + t44) * rSges(6,1) + (t6 * t21 - t7 * t37) * rSges(6,2) + t26) + g(2) * (pkin(4) * t45 + (t6 * t36 - t42) * rSges(6,1) + (-t7 * t21 - t6 * t37) * rSges(6,2) + t29) + g(3) * (-t52 * t11 + t28) + (g(3) * (rSges(6,1) * t21 - rSges(6,2) * t18 + pkin(4)) + t49 * t52) * t9) - m(7) * (g(1) * (t4 * t43 + pkin(5) * t44 + (t7 * t38 + t6 * t8) * rSges(7,1) + (t6 * t10 - t8 * t43) * rSges(7,2) + t26) + g(2) * (t4 * t45 - pkin(5) * t42 + (t6 * t38 - t7 * t8) * rSges(7,1) + (-t7 * t10 - t8 * t45) * rSges(7,2) + t29) + g(3) * (-t51 * t11 + t28) + (g(3) * (rSges(7,1) * t10 - rSges(7,2) * t8 + t4) + t49 * t51) * t9);
U  = t1;
