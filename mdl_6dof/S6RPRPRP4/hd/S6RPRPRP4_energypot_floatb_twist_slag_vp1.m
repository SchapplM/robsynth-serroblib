% Calculate potential energy for
% S6RPRPRP4
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:07
% EndTime: 2019-03-09 03:11:07
% DurationCPUTime: 0.36s
% Computational Cost: add. (230->96), mult. (200->107), div. (0->0), fcn. (184->8), ass. (0->35)
t49 = rSges(7,1) + pkin(5);
t48 = rSges(7,3) + qJ(6);
t22 = qJ(1) + pkin(9);
t16 = sin(t22);
t24 = sin(qJ(3));
t47 = t16 * t24;
t27 = cos(qJ(3));
t46 = t16 * t27;
t17 = cos(t22);
t45 = t17 * t27;
t23 = sin(qJ(5));
t44 = t23 * t24;
t26 = cos(qJ(5));
t43 = t24 * t26;
t42 = qJ(4) * t24;
t41 = pkin(6) + r_base(3);
t25 = sin(qJ(1));
t40 = t25 * pkin(1) + r_base(2);
t28 = cos(qJ(1));
t39 = t28 * pkin(1) + r_base(1);
t38 = qJ(2) + t41;
t37 = t16 * pkin(2) + t40;
t36 = t24 * pkin(3) + t38;
t35 = t17 * pkin(2) + t16 * pkin(7) + t39;
t34 = g(1) * t17 + g(2) * t16;
t33 = pkin(3) * t46 + t16 * t42 + t37;
t32 = t24 * pkin(8) + t36;
t31 = pkin(3) * t45 + t17 * t42 + t35;
t30 = t16 * pkin(4) + pkin(8) * t45 + t31;
t29 = t33 + pkin(8) * t46 + (-pkin(4) - pkin(7)) * t17;
t4 = t16 * t44 - t17 * t26;
t3 = t16 * t43 + t17 * t23;
t2 = t16 * t26 + t17 * t44;
t1 = t16 * t23 - t17 * t43;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t17 * rSges(3,1) - t16 * rSges(3,2) + t39) + g(2) * (t16 * rSges(3,1) + t17 * rSges(3,2) + t40) + g(3) * (rSges(3,3) + t38)) - m(4) * (g(1) * (t16 * rSges(4,3) + t35) + g(2) * (rSges(4,1) * t46 - rSges(4,2) * t47 + t37) + g(3) * (t24 * rSges(4,1) + t27 * rSges(4,2) + t38) + (g(1) * (rSges(4,1) * t27 - rSges(4,2) * t24) + g(2) * (-rSges(4,3) - pkin(7))) * t17) - m(5) * (g(1) * (t16 * rSges(5,1) + t31) + g(2) * (-rSges(5,2) * t46 + rSges(5,3) * t47 + t33) + g(3) * (-t24 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t27 + t36) + (g(1) * (-rSges(5,2) * t27 + rSges(5,3) * t24) + g(2) * (-rSges(5,1) - pkin(7))) * t17) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t30) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t29) + g(3) * (t24 * rSges(6,3) + t32) + (g(3) * (-rSges(6,1) * t23 - rSges(6,2) * t26 - qJ(4)) + t34 * rSges(6,3)) * t27) - m(7) * (g(1) * (t48 * t1 + t49 * t2 + t30) + g(2) * (-t48 * t3 + t49 * t4 + t29) + g(3) * (t24 * rSges(7,2) + t32) + (g(3) * (-t49 * t23 + t48 * t26 - qJ(4)) + t34 * rSges(7,2)) * t27);
U  = t5;
