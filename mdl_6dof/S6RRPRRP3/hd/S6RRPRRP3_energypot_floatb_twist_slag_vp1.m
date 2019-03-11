% Calculate potential energy for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:13
% EndTime: 2019-03-09 11:47:14
% DurationCPUTime: 0.48s
% Computational Cost: add. (247->109), mult. (212->122), div. (0->0), fcn. (196->10), ass. (0->45)
t55 = rSges(5,3) + pkin(8);
t28 = -pkin(9) - pkin(8);
t54 = rSges(6,3) - t28;
t53 = rSges(7,3) + qJ(6) - t28;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t52 = g(1) * t27 + g(2) * t24;
t49 = rSges(3,3) + pkin(7);
t25 = cos(qJ(4));
t10 = t25 * pkin(4) + pkin(3);
t19 = qJ(2) + pkin(10);
t12 = sin(t19);
t47 = rSges(4,2) * t12;
t13 = cos(t19);
t46 = t24 * t13;
t20 = qJ(4) + qJ(5);
t14 = sin(t20);
t45 = t24 * t14;
t15 = cos(t20);
t44 = t24 * t15;
t22 = sin(qJ(4));
t43 = t24 * t22;
t42 = t24 * t25;
t41 = t27 * t13;
t40 = t27 * t14;
t39 = t27 * t15;
t38 = t27 * t22;
t37 = t27 * t25;
t34 = pkin(6) + r_base(3);
t26 = cos(qJ(2));
t11 = t26 * pkin(2) + pkin(1);
t33 = t27 * t11 + r_base(1);
t21 = -qJ(3) - pkin(7);
t32 = t24 * t11 + t27 * t21 + r_base(2);
t23 = sin(qJ(2));
t31 = t23 * pkin(2) + t34;
t30 = -t24 * t21 + t33;
t29 = rSges(3,1) * t26 - rSges(3,2) * t23 + pkin(1);
t6 = t22 * pkin(4) + pkin(5) * t14;
t5 = pkin(5) * t15 + t10;
t4 = t13 * t39 + t45;
t3 = -t13 * t40 + t44;
t2 = t13 * t44 - t40;
t1 = -t13 * t45 - t39;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t27 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t27 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t23 * rSges(3,1) + t26 * rSges(3,2) + t34) + (g(1) * t29 - g(2) * t49) * t27 + (g(1) * t49 + g(2) * t29) * t24) - m(4) * (g(1) * (rSges(4,1) * t41 - t27 * t47 + t33) + g(2) * (-t27 * rSges(4,3) + t32) + g(3) * (t12 * rSges(4,1) + t13 * rSges(4,2) + t31) + (g(1) * (rSges(4,3) - t21) + g(2) * (rSges(4,1) * t13 - t47)) * t24) - m(5) * (g(1) * (pkin(3) * t41 + (t13 * t37 + t43) * rSges(5,1) + (-t13 * t38 + t42) * rSges(5,2) + t30) + g(2) * (pkin(3) * t46 + (t13 * t42 - t38) * rSges(5,1) + (-t13 * t43 - t37) * rSges(5,2) + t32) + g(3) * (-t55 * t13 + t31) + (g(3) * (rSges(5,1) * t25 - rSges(5,2) * t22 + pkin(3)) + t52 * t55) * t12) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t43 + t10 * t41 + t30) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - pkin(4) * t38 + t10 * t46 + t32) + g(3) * (-t54 * t13 + t31) + (g(3) * (rSges(6,1) * t15 - rSges(6,2) * t14 + t10) + t52 * t54) * t12) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t24 * t6 + t5 * t41 + t30) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t27 * t6 + t5 * t46 + t32) + g(3) * (-t53 * t13 + t31) + (g(3) * (rSges(7,1) * t15 - rSges(7,2) * t14 + t5) + t52 * t53) * t12);
U  = t7;
