% Calculate potential energy for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:38:56
% EndTime: 2019-03-09 01:38:56
% DurationCPUTime: 0.45s
% Computational Cost: add. (255->101), mult. (172->111), div. (0->0), fcn. (152->12), ass. (0->39)
t49 = rSges(7,3) + pkin(8) + qJ(5);
t17 = qJ(1) + pkin(9);
t11 = cos(t17);
t8 = sin(t17);
t48 = g(1) * t11 + g(2) * t8;
t47 = rSges(6,3) + qJ(5);
t16 = pkin(10) + qJ(4);
t7 = sin(t16);
t45 = rSges(5,2) * t7;
t10 = cos(t16);
t43 = t8 * t10;
t18 = sin(pkin(11));
t42 = t8 * t18;
t20 = cos(pkin(11));
t41 = t8 * t20;
t40 = t11 * t10;
t39 = t11 * t18;
t38 = t11 * t20;
t36 = rSges(4,3) + qJ(3);
t34 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t33 = t24 * pkin(1) + r_base(2);
t25 = cos(qJ(1));
t32 = t25 * pkin(1) + r_base(1);
t21 = cos(pkin(10));
t5 = t21 * pkin(3) + pkin(2);
t31 = t11 * t5 + t32;
t30 = qJ(2) + t34;
t23 = -pkin(7) - qJ(3);
t29 = t11 * t23 + t8 * t5 + t33;
t19 = sin(pkin(10));
t28 = t19 * pkin(3) + t30;
t27 = rSges(4,1) * t21 - rSges(4,2) * t19 + pkin(2);
t26 = -t8 * t23 + t31;
t15 = pkin(11) + qJ(6);
t9 = cos(t15);
t6 = sin(t15);
t4 = t20 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t11 * rSges(3,1) - t8 * rSges(3,2) + t32) + g(2) * (t8 * rSges(3,1) + t11 * rSges(3,2) + t33) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * t32 + g(2) * t33 + g(3) * (t19 * rSges(4,1) + t21 * rSges(4,2) + t30) + (g(1) * t36 + g(2) * t27) * t8 + (g(1) * t27 - g(2) * t36) * t11) - m(5) * (g(1) * (rSges(5,1) * t40 - t11 * t45 + t31) + g(2) * (-t11 * rSges(5,3) + t29) + g(3) * (t7 * rSges(5,1) + t10 * rSges(5,2) + t28) + (g(1) * (rSges(5,3) - t23) + g(2) * (rSges(5,1) * t10 - t45)) * t8) - m(6) * (g(1) * (pkin(4) * t40 + (t10 * t38 + t42) * rSges(6,1) + (-t10 * t39 + t41) * rSges(6,2) + t26) + g(2) * (pkin(4) * t43 + (t10 * t41 - t39) * rSges(6,1) + (-t10 * t42 - t38) * rSges(6,2) + t29) + g(3) * (-t47 * t10 + t28) + (g(3) * (rSges(6,1) * t20 - rSges(6,2) * t18 + pkin(4)) + t48 * t47) * t7) - m(7) * (g(1) * (t4 * t40 + pkin(5) * t42 + (t40 * t9 + t8 * t6) * rSges(7,1) + (-t40 * t6 + t8 * t9) * rSges(7,2) + t26) + g(2) * (t4 * t43 - pkin(5) * t39 + (-t11 * t6 + t43 * t9) * rSges(7,1) + (-t11 * t9 - t43 * t6) * rSges(7,2) + t29) + g(3) * (-t49 * t10 + t28) + (g(3) * (rSges(7,1) * t9 - rSges(7,2) * t6 + t4) + t48 * t49) * t7);
U  = t1;
