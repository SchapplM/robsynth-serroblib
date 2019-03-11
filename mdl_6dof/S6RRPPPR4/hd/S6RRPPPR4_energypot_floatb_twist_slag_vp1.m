% Calculate potential energy for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:10
% EndTime: 2019-03-09 08:17:10
% DurationCPUTime: 0.50s
% Computational Cost: add. (186->111), mult. (293->130), div. (0->0), fcn. (297->8), ass. (0->33)
t51 = -pkin(8) - rSges(7,3);
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t35 = g(1) * t28 + g(2) * t25;
t24 = sin(qJ(2));
t47 = t24 * t25;
t46 = t24 * t28;
t22 = cos(pkin(9));
t45 = t25 * t22;
t27 = cos(qJ(2));
t44 = t25 * t27;
t43 = qJ(3) * t24;
t42 = qJ(4) * t27;
t41 = pkin(6) + r_base(3);
t40 = t25 * pkin(1) + r_base(2);
t39 = t24 * pkin(2) + t41;
t38 = t28 * pkin(1) + t25 * pkin(7) + r_base(1);
t37 = t24 * qJ(4) + t39;
t36 = pkin(2) * t44 + t25 * t43 + t40;
t34 = t27 * t22 * qJ(5) + t37;
t33 = t38 + (pkin(2) * t27 + t43) * t28;
t32 = t25 * pkin(3) + t28 * t42 + t33;
t31 = t25 * t42 + t36 + (-pkin(3) - pkin(7)) * t28;
t21 = sin(pkin(9));
t3 = t21 * t25 - t22 * t46;
t4 = t21 * t46 + t45;
t30 = t4 * pkin(4) + t3 * qJ(5) + t32;
t5 = t21 * t28 + t24 * t45;
t6 = t21 * t47 - t22 * t28;
t29 = t6 * pkin(4) - t5 * qJ(5) + t31;
t26 = cos(qJ(6));
t23 = sin(qJ(6));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t28 - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + rSges(2,2) * t28 + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t25 * rSges(3,3) + t38) + g(2) * (rSges(3,1) * t44 - rSges(3,2) * t47 + t40) + g(3) * (rSges(3,1) * t24 + rSges(3,2) * t27 + t41) + (g(1) * (rSges(3,1) * t27 - rSges(3,2) * t24) + g(2) * (-rSges(3,3) - pkin(7))) * t28) - m(4) * (g(1) * (t25 * rSges(4,1) + t33) + g(2) * (-rSges(4,2) * t44 + rSges(4,3) * t47 + t36) + g(3) * (-t24 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t27 + t39) + (g(1) * (-rSges(4,2) * t27 + rSges(4,3) * t24) + g(2) * (-rSges(4,1) - pkin(7))) * t28) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t32) + g(2) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t31) + g(3) * (t24 * rSges(5,3) + t37) + (g(3) * (-rSges(5,1) * t21 - rSges(5,2) * t22 - qJ(3)) + t35 * rSges(5,3)) * t27) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,3) + t30) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,3) + t29) + g(3) * (t24 * rSges(6,2) + t34) + (g(3) * (rSges(6,3) * t22 - qJ(3) + (-rSges(6,1) - pkin(4)) * t21) + t35 * rSges(6,2)) * t27) - m(7) * (g(1) * (t4 * pkin(5) + (t23 * t3 + t26 * t4) * rSges(7,1) + (-t23 * t4 + t26 * t3) * rSges(7,2) + t30) + g(2) * (t6 * pkin(5) + (-t23 * t5 + t26 * t6) * rSges(7,1) + (-t23 * t6 - t26 * t5) * rSges(7,2) + t29) + g(3) * (t51 * t24 + t34) + (g(3) * (-qJ(3) + (rSges(7,1) * t23 + rSges(7,2) * t26) * t22 + (-rSges(7,1) * t26 + rSges(7,2) * t23 - pkin(4) - pkin(5)) * t21) + t35 * t51) * t27);
U  = t1;
