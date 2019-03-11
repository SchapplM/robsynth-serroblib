% Calculate potential energy for
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:17
% EndTime: 2019-03-09 05:34:17
% DurationCPUTime: 0.46s
% Computational Cost: add. (174->106), mult. (262->123), div. (0->0), fcn. (266->8), ass. (0->34)
t45 = -pkin(9) - rSges(7,3);
t22 = sin(qJ(1));
t47 = g(1) * t22;
t26 = cos(qJ(1));
t46 = g(2) * t26;
t25 = cos(qJ(3));
t44 = rSges(4,2) * t25;
t21 = sin(qJ(3));
t43 = t21 * t22;
t20 = sin(qJ(4));
t42 = t22 * t20;
t24 = cos(qJ(4));
t41 = t22 * t24;
t40 = t26 * t20;
t39 = t26 * t24;
t38 = pkin(6) + r_base(3);
t37 = t22 * pkin(1) + r_base(2);
t36 = pkin(2) + t38;
t35 = t26 * pkin(1) + t22 * qJ(2) + r_base(1);
t34 = t22 * pkin(7) + t37;
t33 = t26 * pkin(7) + t35;
t32 = t25 * pkin(3) + t21 * pkin(8) + t36;
t31 = pkin(3) * t43 + t33;
t30 = t32 + (pkin(4) * t24 + qJ(5) * t20) * t25;
t3 = t21 * t42 - t39;
t4 = t21 * t41 + t40;
t29 = t4 * pkin(4) + t3 * qJ(5) + t31;
t28 = t34 + (-t21 * pkin(3) + t25 * pkin(8) - qJ(2)) * t26;
t5 = t21 * t40 + t41;
t6 = -t21 * t39 + t42;
t27 = t6 * pkin(4) - t5 * qJ(5) + t28;
t23 = cos(qJ(6));
t19 = sin(qJ(6));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (-t26 * rSges(3,2) + t22 * rSges(3,3) + t35) + g(2) * (-t22 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t26 + t37) + g(3) * (rSges(3,1) + t38)) - m(4) * (g(1) * (rSges(4,1) * t43 + t22 * t44 + t33) + g(2) * (t22 * rSges(4,3) + t34) + g(3) * (t25 * rSges(4,1) - t21 * rSges(4,2) + t36) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t21 - qJ(2) - t44)) * t26) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t31) + g(2) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t28) + g(3) * (t21 * rSges(5,3) + t32) + (rSges(5,3) * t46 + g(3) * (rSges(5,1) * t24 - rSges(5,2) * t20) + (-rSges(5,3) - pkin(8)) * t47) * t25) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,3) + t29) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,3) + t27) + g(3) * (t21 * rSges(6,2) + t30) + (rSges(6,2) * t46 + g(3) * (rSges(6,1) * t24 + rSges(6,3) * t20) + (-rSges(6,2) - pkin(8)) * t47) * t25) - m(7) * (g(1) * (t4 * pkin(5) + (t3 * t19 + t4 * t23) * rSges(7,1) + (-t4 * t19 + t3 * t23) * rSges(7,2) + t29) + g(2) * (t6 * pkin(5) + (-t5 * t19 + t6 * t23) * rSges(7,1) + (-t6 * t19 - t5 * t23) * rSges(7,2) + t27) + g(3) * (t45 * t21 + t30) + (g(3) * (t24 * pkin(5) + (t19 * t20 + t23 * t24) * rSges(7,1) + (-t19 * t24 + t20 * t23) * rSges(7,2)) + t45 * t46 + (-pkin(8) - t45) * t47) * t25);
U  = t1;
