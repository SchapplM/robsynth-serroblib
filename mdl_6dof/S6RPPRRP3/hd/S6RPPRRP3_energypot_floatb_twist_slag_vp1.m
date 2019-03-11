% Calculate potential energy for
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:38
% EndTime: 2019-03-09 02:02:38
% DurationCPUTime: 0.39s
% Computational Cost: add. (212->92), mult. (174->100), div. (0->0), fcn. (158->8), ass. (0->35)
t47 = rSges(7,1) + pkin(5);
t46 = rSges(7,3) + qJ(6);
t21 = sin(qJ(4));
t45 = -pkin(4) * t21 - qJ(3);
t19 = qJ(1) + pkin(9);
t13 = sin(t19);
t43 = g(1) * t13;
t14 = cos(t19);
t42 = g(2) * t14;
t24 = cos(qJ(4));
t41 = rSges(5,2) * t24;
t40 = t13 * t21;
t20 = sin(qJ(5));
t39 = t20 * t21;
t23 = cos(qJ(5));
t38 = t21 * t23;
t37 = pkin(6) + r_base(3);
t22 = sin(qJ(1));
t36 = t22 * pkin(1) + r_base(2);
t25 = cos(qJ(1));
t35 = t25 * pkin(1) + r_base(1);
t34 = qJ(2) + t37;
t33 = t13 * pkin(2) + t36;
t32 = pkin(3) + t34;
t31 = t13 * pkin(7) + t33;
t30 = t14 * pkin(2) + t13 * qJ(3) + t35;
t29 = t14 * t24 * pkin(8) + t31;
t28 = t14 * pkin(7) + t30;
t27 = pkin(4) * t40 + t28;
t26 = t24 * pkin(4) + t21 * pkin(8) + t32;
t4 = t13 * t20 - t14 * t38;
t3 = t13 * t23 + t14 * t39;
t2 = t13 * t38 + t14 * t20;
t1 = t13 * t39 - t14 * t23;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t35) + g(2) * (t13 * rSges(3,1) + t14 * rSges(3,2) + t36) + g(3) * (rSges(3,3) + t34)) - m(4) * (g(1) * (-t14 * rSges(4,2) + t13 * rSges(4,3) + t30) + g(2) * (-t13 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t14 + t33) + g(3) * (rSges(4,1) + t34)) - m(5) * (g(1) * (rSges(5,1) * t40 + t13 * t41 + t28) + g(2) * (t13 * rSges(5,3) + t31) + g(3) * (t24 * rSges(5,1) - t21 * rSges(5,2) + t32) + (g(1) * rSges(5,3) + g(2) * (-rSges(5,1) * t21 - qJ(3) - t41)) * t14) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t27) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t29) + g(3) * (t21 * rSges(6,3) + t26) + (g(3) * (rSges(6,1) * t23 - rSges(6,2) * t20) + (-rSges(6,3) - pkin(8)) * t43) * t24 + (rSges(6,3) * t24 + t45) * t42) - m(7) * (g(1) * (t46 * t1 + t47 * t2 + t27) + g(2) * (t45 * t14 - t46 * t3 + t47 * t4 + t29) + g(3) * (t21 * rSges(7,2) + t26) + (rSges(7,2) * t42 + g(3) * (t46 * t20 + t47 * t23) + (-rSges(7,2) - pkin(8)) * t43) * t24);
U  = t5;
