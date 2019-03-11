% Calculate potential energy for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:04
% EndTime: 2019-03-09 02:05:05
% DurationCPUTime: 0.38s
% Computational Cost: add. (202->90), mult. (298->101), div. (0->0), fcn. (340->8), ass. (0->34)
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t52 = -pkin(4) * t27 - pkin(8) * t25;
t51 = rSges(7,1) + pkin(5);
t50 = rSges(7,3) + qJ(6);
t47 = rSges(5,3) + pkin(7);
t46 = sin(qJ(1));
t24 = sin(qJ(5));
t45 = t24 * t27;
t26 = cos(qJ(5));
t44 = t26 * t27;
t43 = cos(pkin(9));
t42 = sin(pkin(9));
t41 = pkin(6) + r_base(3);
t40 = t46 * pkin(1) + r_base(2);
t39 = -qJ(3) + t41;
t28 = cos(qJ(1));
t38 = t28 * pkin(1) + t46 * qJ(2) + r_base(1);
t37 = t27 * pkin(8) + t39;
t36 = t28 * pkin(2) + t38;
t13 = -t28 * t43 - t42 * t46;
t14 = t28 * t42 - t43 * t46;
t35 = -g(1) * t13 - g(2) * t14;
t34 = -rSges(5,1) * t27 + rSges(5,2) * t25;
t33 = -t13 * pkin(3) + t36;
t32 = t46 * pkin(2) - t28 * qJ(2) + t40;
t31 = -t14 * pkin(3) + t32;
t30 = t14 * pkin(7) + t52 * t13 + t33;
t29 = -t13 * pkin(7) + t52 * t14 + t31;
t4 = -t13 * t44 + t14 * t24;
t3 = -t13 * t45 - t14 * t26;
t2 = -t13 * t24 - t14 * t44;
t1 = t13 * t26 - t14 * t45;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t46 * rSges(2,2) + r_base(1)) + g(2) * (t46 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t28 * rSges(3,1) + t46 * rSges(3,3) + t38) + g(2) * (t46 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t28 + t40) + g(3) * (rSges(3,2) + t41)) - m(4) * (g(1) * (-t13 * rSges(4,1) - t14 * rSges(4,2) + t36) + g(2) * (-t14 * rSges(4,1) + t13 * rSges(4,2) + t32) + g(3) * (-rSges(4,3) + t39)) - m(5) * (g(1) * t33 + g(2) * t31 + g(3) * (-t25 * rSges(5,1) - t27 * rSges(5,2) + t39) + (g(1) * t47 + g(2) * t34) * t14 + (g(1) * t34 - g(2) * t47) * t13) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t30) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t29) + g(3) * (t27 * rSges(6,3) + t37) + (g(3) * (-rSges(6,1) * t26 + rSges(6,2) * t24 - pkin(4)) + t35 * rSges(6,3)) * t25) - m(7) * (g(1) * (t50 * t3 + t51 * t4 + t30) + g(2) * (t50 * t1 + t51 * t2 + t29) + g(3) * (t27 * rSges(7,2) + t37) + (g(3) * (-t50 * t24 - t51 * t26 - pkin(4)) + t35 * rSges(7,2)) * t25);
U  = t5;
