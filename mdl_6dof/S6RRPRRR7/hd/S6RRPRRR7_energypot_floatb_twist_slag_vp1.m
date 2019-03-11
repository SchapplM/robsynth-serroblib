% Calculate potential energy for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:35
% EndTime: 2019-03-09 13:55:35
% DurationCPUTime: 0.54s
% Computational Cost: add. (209->110), mult. (322->123), div. (0->0), fcn. (338->10), ass. (0->37)
t52 = pkin(9) + rSges(6,3);
t46 = pkin(10) + pkin(9) + rSges(7,3);
t26 = sin(qJ(4));
t31 = cos(qJ(2));
t27 = sin(qJ(2));
t30 = cos(qJ(4));
t50 = t27 * t30;
t54 = t31 * t26 - t50;
t53 = g(3) * t54;
t28 = sin(qJ(1));
t51 = t27 * t28;
t49 = t28 * t31;
t32 = cos(qJ(1));
t47 = t31 * t32;
t45 = qJ(3) * t27;
t44 = pkin(6) + r_base(3);
t43 = t28 * pkin(1) + r_base(2);
t42 = t27 * pkin(2) + t44;
t41 = t32 * pkin(1) + t28 * pkin(7) + r_base(1);
t40 = pkin(2) * t49 + t28 * t45 + t43;
t5 = t27 * t26 + t31 * t30;
t39 = pkin(2) * t47 + t32 * t45 + t41;
t24 = qJ(5) + qJ(6);
t16 = sin(t24);
t17 = cos(t24);
t29 = cos(qJ(5));
t38 = t17 * rSges(7,1) - t16 * rSges(7,2) + t29 * pkin(5) + pkin(4);
t37 = pkin(3) * t49 + t32 * pkin(8) + t40;
t36 = pkin(3) * t47 + t39;
t35 = t27 * pkin(3) - t31 * qJ(3) + t42;
t25 = sin(qJ(5));
t34 = t16 * rSges(7,1) + t17 * rSges(7,2) + t25 * pkin(5);
t4 = t5 * t32;
t3 = t26 * t47 - t32 * t50;
t2 = t5 * t28;
t1 = t54 * t28;
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t32 * rSges(2,1) - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + t32 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * (t28 * rSges(3,3) + t41) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t51 + t43) + g(3) * (t27 * rSges(3,1) + t31 * rSges(3,2) + t44) + (g(1) * (rSges(3,1) * t31 - rSges(3,2) * t27) + g(2) * (-rSges(3,3) - pkin(7))) * t32) - m(4) * (g(1) * (t28 * rSges(4,2) + t39) + g(2) * (rSges(4,1) * t49 + rSges(4,3) * t51 + t40) + g(3) * (t27 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t31 + t42) + (g(1) * (rSges(4,1) * t31 + rSges(4,3) * t27) + g(2) * (-rSges(4,2) - pkin(7))) * t32) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + (-rSges(5,3) - pkin(8)) * t28 + t36) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + (rSges(5,3) - pkin(7)) * t32 + t37) + g(3) * (-rSges(5,1) * t54 - t5 * rSges(5,2) + t35)) - m(6) * (g(1) * (t4 * pkin(4) - t28 * pkin(8) + (-t28 * t25 + t4 * t29) * rSges(6,1) + (-t4 * t25 - t28 * t29) * rSges(6,2) + t52 * t3 + t36) + g(2) * (t2 * pkin(4) - t32 * pkin(7) + (t2 * t29 + t32 * t25) * rSges(6,1) + (-t2 * t25 + t32 * t29) * rSges(6,2) + t52 * t1 + t37) + g(3) * (t52 * t5 + t35) - (t29 * rSges(6,1) - t25 * rSges(6,2) + pkin(4)) * t53) - m(7) * (g(1) * (t46 * t3 + t38 * t4 + (-pkin(8) - t34) * t28 + t36) + g(2) * (t38 * t2 + (-pkin(7) + t34) * t32 + t46 * t1 + t37) + g(3) * (t46 * t5 + t35) - t38 * t53);
U  = t6;
