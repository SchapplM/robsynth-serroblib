% Calculate potential energy for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:05
% EndTime: 2019-03-09 04:18:06
% DurationCPUTime: 0.45s
% Computational Cost: add. (159->100), mult. (192->111), div. (0->0), fcn. (172->8), ass. (0->32)
t40 = rSges(6,3) + pkin(8);
t16 = sin(qJ(1));
t42 = g(1) * t16;
t19 = cos(qJ(1));
t41 = g(2) * t19;
t15 = sin(qJ(3));
t39 = t15 * t16;
t18 = cos(qJ(3));
t38 = t16 * t18;
t14 = sin(qJ(5));
t37 = t19 * t14;
t17 = cos(qJ(5));
t36 = t19 * t17;
t35 = rSges(7,3) + pkin(9) + pkin(8);
t34 = t18 * qJ(4);
t33 = pkin(6) + r_base(3);
t32 = t16 * pkin(1) + r_base(2);
t31 = pkin(2) + t33;
t30 = t16 * pkin(7) + t32;
t29 = t19 * pkin(1) + t16 * qJ(2) + r_base(1);
t28 = t19 * t34 + t30;
t27 = t19 * pkin(7) + t29;
t26 = t18 * pkin(3) + t15 * qJ(4) + t31;
t25 = -rSges(5,2) * t15 - rSges(5,3) * t18;
t24 = pkin(3) * t39 + t27;
t13 = qJ(5) + qJ(6);
t4 = sin(t13);
t5 = cos(t13);
t23 = t5 * rSges(7,1) - t4 * rSges(7,2) + t17 * pkin(5) + pkin(4);
t22 = rSges(7,1) * t4 + rSges(7,2) * t5 + pkin(5) * t14;
t21 = g(1) * t24 + g(2) * t28;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t16 * rSges(2,2) + r_base(1)) + g(2) * (t16 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (-t19 * rSges(3,2) + t16 * rSges(3,3) + t29) + g(2) * (-t16 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t19 + t32) + g(3) * (rSges(3,1) + t33)) - m(4) * (g(1) * (rSges(4,1) * t39 + rSges(4,2) * t38 + t27) + g(2) * (t16 * rSges(4,3) + t30) + g(3) * (t18 * rSges(4,1) - t15 * rSges(4,2) + t31) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t15 - rSges(4,2) * t18 - qJ(2))) * t19) - m(5) * (g(3) * (-t18 * rSges(5,2) + t15 * rSges(5,3) + t26) + (g(1) * (t25 - t34) + g(2) * rSges(5,1)) * t16 + (g(1) * rSges(5,1) + g(2) * (-t15 * pkin(3) - qJ(2) - t25)) * t19 + t21) - m(6) * (g(1) * (t19 * pkin(4) - t16 * t34 + (-t14 * t38 + t36) * rSges(6,1) + (-t17 * t38 - t37) * rSges(6,2) + t24) + g(2) * (t16 * pkin(4) - t19 * qJ(2) + (t16 * t17 + t18 * t37) * rSges(6,1) + (-t16 * t14 + t18 * t36) * rSges(6,2) + t28) + g(3) * (t40 * t18 + t26) + (g(3) * (rSges(6,1) * t14 + rSges(6,2) * t17) + t40 * t42 + (-pkin(3) - t40) * t41) * t15) - m(7) * (g(3) * t26 + (g(1) * t23 - g(2) * qJ(2)) * t19 + g(2) * t23 * t16 + (g(3) * t35 + t22 * t41 + (-qJ(4) - t22) * t42) * t18 + (g(3) * t22 + t35 * t42 + (-pkin(3) - t35) * t41) * t15 + t21);
U  = t1;
