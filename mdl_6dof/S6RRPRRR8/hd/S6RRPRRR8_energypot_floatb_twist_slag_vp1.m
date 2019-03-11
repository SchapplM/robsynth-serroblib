% Calculate potential energy for
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:12
% EndTime: 2019-03-09 14:02:12
% DurationCPUTime: 0.59s
% Computational Cost: add. (266->123), mult. (240->142), div. (0->0), fcn. (228->12), ass. (0->37)
t26 = -pkin(8) - qJ(3);
t50 = rSges(5,3) - t26;
t22 = -pkin(9) + t26;
t49 = rSges(6,3) - t22;
t48 = rSges(7,3) + pkin(10) - t22;
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t47 = g(1) * t30 + g(2) * t28;
t46 = rSges(4,3) + qJ(3);
t25 = cos(pkin(11));
t11 = t25 * pkin(3) + pkin(2);
t23 = pkin(11) + qJ(4);
t13 = sin(t23);
t24 = sin(pkin(11));
t4 = t24 * pkin(3) + pkin(4) * t13;
t27 = sin(qJ(2));
t43 = rSges(3,2) * t27;
t42 = t24 * t30;
t41 = t28 * t24;
t29 = cos(qJ(2));
t40 = t28 * t29;
t39 = t29 * t30;
t34 = pkin(6) + r_base(3);
t14 = cos(t23);
t3 = pkin(4) * t14 + t11;
t33 = t28 * pkin(1) + r_base(2);
t32 = t30 * pkin(1) + t28 * pkin(7) + r_base(1);
t15 = qJ(5) + t23;
t31 = -t30 * pkin(7) + t33;
t12 = qJ(6) + t15;
t10 = cos(t15);
t9 = sin(t15);
t6 = cos(t12);
t5 = sin(t12);
t2 = pkin(5) * t9 + t4;
t1 = pkin(5) * t10 + t3;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t28 * rSges(3,3) + t32) + g(2) * (rSges(3,1) * t40 - t28 * t43 + t33) + g(3) * (rSges(3,1) * t27 + rSges(3,2) * t29 + t34) + (g(1) * (rSges(3,1) * t29 - t43) + g(2) * (-rSges(3,3) - pkin(7))) * t30) - m(4) * (g(1) * (pkin(2) * t39 + (t25 * t39 + t41) * rSges(4,1) + (-t24 * t39 + t28 * t25) * rSges(4,2) + t32) + g(2) * (pkin(2) * t40 + (t25 * t40 - t42) * rSges(4,1) + (-t24 * t40 - t25 * t30) * rSges(4,2) + t31) + g(3) * (-t46 * t29 + t34) + (g(3) * (rSges(4,1) * t25 - rSges(4,2) * t24 + pkin(2)) + t47 * t46) * t27) - m(5) * (g(1) * (t11 * t39 + pkin(3) * t41 + (t28 * t13 + t14 * t39) * rSges(5,1) + (-t13 * t39 + t28 * t14) * rSges(5,2) + t32) + g(2) * (t11 * t40 - pkin(3) * t42 + (-t13 * t30 + t14 * t40) * rSges(5,1) + (-t13 * t40 - t14 * t30) * rSges(5,2) + t31) + g(3) * (-t50 * t29 + t34) + (g(3) * (rSges(5,1) * t14 - rSges(5,2) * t13 + t11) + t47 * t50) * t27) - m(6) * (g(1) * (t3 * t39 + t28 * t4 + (t10 * t39 + t28 * t9) * rSges(6,1) + (t28 * t10 - t9 * t39) * rSges(6,2) + t32) + g(2) * (t3 * t40 - t30 * t4 + (t10 * t40 - t30 * t9) * rSges(6,1) + (-t10 * t30 - t9 * t40) * rSges(6,2) + t31) + g(3) * (-t49 * t29 + t34) + (g(3) * (rSges(6,1) * t10 - rSges(6,2) * t9 + t3) + t47 * t49) * t27) - m(7) * (g(1) * (t1 * t39 + t28 * t2 + (t28 * t5 + t6 * t39) * rSges(7,1) + (t28 * t6 - t5 * t39) * rSges(7,2) + t32) + g(2) * (t1 * t40 - t30 * t2 + (-t30 * t5 + t6 * t40) * rSges(7,1) + (-t30 * t6 - t5 * t40) * rSges(7,2) + t31) + g(3) * (-t48 * t29 + t34) + (g(3) * (rSges(7,1) * t6 - rSges(7,2) * t5 + t1) + t47 * t48) * t27);
U  = t7;
