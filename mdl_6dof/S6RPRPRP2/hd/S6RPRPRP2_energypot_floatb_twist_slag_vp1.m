% Calculate potential energy for
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:25
% EndTime: 2019-03-09 03:04:26
% DurationCPUTime: 0.32s
% Computational Cost: add. (259->95), mult. (185->103), div. (0->0), fcn. (169->10), ass. (0->40)
t51 = rSges(7,1) + pkin(5);
t50 = rSges(7,3) + qJ(6);
t49 = rSges(4,3) + pkin(7);
t22 = qJ(3) + pkin(10);
t15 = sin(t22);
t23 = qJ(1) + pkin(9);
t16 = sin(t23);
t48 = t15 * t16;
t18 = cos(t23);
t47 = t15 * t18;
t25 = sin(qJ(5));
t46 = t16 * t25;
t28 = cos(qJ(5));
t45 = t16 * t28;
t17 = cos(t22);
t44 = t17 * t18;
t43 = t18 * t25;
t42 = t18 * t28;
t41 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t40 = t27 * pkin(1) + r_base(2);
t30 = cos(qJ(1));
t39 = t30 * pkin(1) + r_base(1);
t38 = qJ(2) + t41;
t29 = cos(qJ(3));
t14 = t29 * pkin(3) + pkin(2);
t37 = t18 * t14 + t39;
t24 = -qJ(4) - pkin(7);
t36 = t16 * t14 + t18 * t24 + t40;
t26 = sin(qJ(3));
t35 = t26 * pkin(3) + t38;
t34 = t15 * pkin(4) + t35;
t33 = rSges(4,1) * t29 - rSges(4,2) * t26 + pkin(2);
t32 = t16 * t17 * pkin(4) + pkin(8) * t48 + t36;
t31 = pkin(4) * t44 + pkin(8) * t47 - t16 * t24 + t37;
t4 = t17 * t42 + t46;
t3 = t17 * t43 - t45;
t2 = t17 * t45 - t43;
t1 = t17 * t46 + t42;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t30 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t30 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t18 * rSges(3,1) - t16 * rSges(3,2) + t39) + g(2) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t40) + g(3) * (rSges(3,3) + t38)) - m(4) * (g(1) * t39 + g(2) * t40 + g(3) * (t26 * rSges(4,1) + t29 * rSges(4,2) + t38) + (g(1) * t33 - g(2) * t49) * t18 + (g(1) * t49 + g(2) * t33) * t16) - m(5) * (g(1) * (rSges(5,1) * t44 - rSges(5,2) * t47 + t37) + g(2) * (-t18 * rSges(5,3) + t36) + g(3) * (t15 * rSges(5,1) + t17 * rSges(5,2) + t35) + (g(1) * (rSges(5,3) - t24) + g(2) * (rSges(5,1) * t17 - rSges(5,2) * t15)) * t16) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t47 + t31) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + rSges(6,3) * t48 + t32) + g(3) * ((-rSges(6,3) - pkin(8)) * t17 + (rSges(6,1) * t28 - rSges(6,2) * t25) * t15 + t34)) - m(7) * (g(1) * (t50 * t3 + t4 * t51 + t31) + g(2) * (t50 * t1 + t51 * t2 + t32) + g(3) * (t34 + (-rSges(7,2) - pkin(8)) * t17) + (g(3) * (t50 * t25 + t28 * t51) + (g(1) * t18 + g(2) * t16) * rSges(7,2)) * t15);
U  = t5;
