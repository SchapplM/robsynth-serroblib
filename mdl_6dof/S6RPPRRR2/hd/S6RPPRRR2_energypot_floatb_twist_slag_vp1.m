% Calculate potential energy for
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:14
% EndTime: 2019-03-09 02:20:14
% DurationCPUTime: 0.43s
% Computational Cost: add. (255->100), mult. (172->110), div. (0->0), fcn. (152->12), ass. (0->43)
t55 = rSges(6,3) + pkin(8);
t54 = rSges(7,3) + pkin(9) + pkin(8);
t15 = pkin(11) + qJ(4);
t6 = sin(t15);
t8 = cos(t15);
t53 = rSges(5,1) * t8 - rSges(5,2) * t6;
t16 = qJ(1) + pkin(10);
t7 = sin(t16);
t9 = cos(t16);
t52 = g(1) * t9 + g(2) * t7;
t49 = t7 * t8;
t48 = t9 * t8;
t17 = qJ(5) + qJ(6);
t10 = sin(t17);
t45 = t7 * t10;
t11 = cos(t17);
t44 = t7 * t11;
t21 = sin(qJ(5));
t43 = t7 * t21;
t23 = cos(qJ(5));
t42 = t7 * t23;
t41 = t9 * t10;
t40 = t9 * t11;
t39 = t9 * t21;
t38 = t9 * t23;
t35 = rSges(4,3) + qJ(3);
t34 = pkin(6) + r_base(3);
t22 = sin(qJ(1));
t33 = pkin(1) * t22 + r_base(2);
t24 = cos(qJ(1));
t32 = pkin(1) * t24 + r_base(1);
t19 = cos(pkin(11));
t4 = pkin(3) * t19 + pkin(2);
t31 = t4 * t9 + t32;
t30 = qJ(2) + t34;
t20 = -pkin(7) - qJ(3);
t29 = t20 * t9 + t4 * t7 + t33;
t18 = sin(pkin(11));
t28 = pkin(3) * t18 + t30;
t27 = rSges(4,1) * t19 - rSges(4,2) * t18 + pkin(2);
t26 = -t7 * t20 + t31;
t5 = pkin(5) * t23 + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t24 - rSges(2,2) * t22 + r_base(1)) + g(2) * (rSges(2,1) * t22 + rSges(2,2) * t24 + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (rSges(3,1) * t9 - rSges(3,2) * t7 + t32) + g(2) * (rSges(3,1) * t7 + rSges(3,2) * t9 + t33) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * t32 + g(2) * t33 + g(3) * (rSges(4,1) * t18 + rSges(4,2) * t19 + t30) + (g(1) * t27 - g(2) * t35) * t9 + (g(1) * t35 + g(2) * t27) * t7) - m(5) * (g(1) * (t53 * t9 + t31) + g(2) * (-t9 * rSges(5,3) + t29) + g(3) * (rSges(5,1) * t6 + rSges(5,2) * t8 + t28) + (g(1) * (rSges(5,3) - t20) + g(2) * t53) * t7) - m(6) * (g(1) * (pkin(4) * t48 + (t38 * t8 + t43) * rSges(6,1) + (-t39 * t8 + t42) * rSges(6,2) + t26) + g(2) * (pkin(4) * t49 + (t42 * t8 - t39) * rSges(6,1) + (-t43 * t8 - t38) * rSges(6,2) + t29) + g(3) * (-t55 * t8 + t28) + (g(3) * (rSges(6,1) * t23 - rSges(6,2) * t21 + pkin(4)) + t52 * t55) * t6) - m(7) * (g(1) * (t5 * t48 + pkin(5) * t43 + (t40 * t8 + t45) * rSges(7,1) + (-t41 * t8 + t44) * rSges(7,2) + t26) + g(2) * (t5 * t49 - pkin(5) * t39 + (t44 * t8 - t41) * rSges(7,1) + (-t45 * t8 - t40) * rSges(7,2) + t29) + g(3) * (-t54 * t8 + t28) + (g(3) * (rSges(7,1) * t11 - rSges(7,2) * t10 + t5) + t52 * t54) * t6);
U  = t1;
