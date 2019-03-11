% Calculate potential energy for
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:27
% EndTime: 2019-03-09 02:30:27
% DurationCPUTime: 0.41s
% Computational Cost: add. (149->97), mult. (163->105), div. (0->0), fcn. (143->8), ass. (0->31)
t41 = rSges(6,3) + pkin(8);
t40 = rSges(7,3) + pkin(9) + pkin(8);
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t39 = g(1) * t16 + g(2) * t13;
t11 = sin(qJ(5));
t35 = t13 * t11;
t12 = sin(qJ(4));
t34 = t13 * t12;
t14 = cos(qJ(5));
t33 = t13 * t14;
t32 = t16 * t11;
t31 = t16 * t12;
t30 = t16 * t14;
t28 = pkin(6) + r_base(3);
t27 = t13 * pkin(1) + r_base(2);
t26 = pkin(2) + t28;
t25 = t13 * qJ(3) + t27;
t24 = t16 * pkin(1) + t13 * qJ(2) + r_base(1);
t23 = pkin(3) + t26;
t22 = t16 * pkin(7) + t25;
t21 = t16 * qJ(3) + t24;
t15 = cos(qJ(4));
t20 = rSges(5,1) * t12 + rSges(5,2) * t15;
t19 = -t13 * pkin(7) + t21;
t18 = -t16 * qJ(2) + t22;
t10 = qJ(5) + qJ(6);
t3 = cos(t10);
t2 = sin(t10);
t1 = t14 * pkin(5) + pkin(4);
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t16 * rSges(2,1) - t13 * rSges(2,2) + r_base(1)) + g(2) * (t13 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(1) * (-t16 * rSges(3,2) + t13 * rSges(3,3) + t24) + g(2) * (-t13 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t16 + t27) + g(3) * (rSges(3,1) + t28)) - m(4) * (g(1) * (t13 * rSges(4,2) + t16 * rSges(4,3) + t21) + g(2) * (t13 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t16 + t25) + g(3) * (rSges(4,1) + t26)) - m(5) * (g(1) * t21 + g(2) * t22 + g(3) * (t15 * rSges(5,1) - t12 * rSges(5,2) + t23) + (g(1) * t20 + g(2) * (rSges(5,3) - qJ(2))) * t16 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t20) * t13) - m(6) * (g(1) * (pkin(4) * t31 + (t12 * t30 - t35) * rSges(6,1) + (-t11 * t31 - t33) * rSges(6,2) + t19) + g(2) * (pkin(4) * t34 + (t12 * t33 + t32) * rSges(6,1) + (-t11 * t34 + t30) * rSges(6,2) + t18) + g(3) * (t41 * t12 + t23) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t11 + pkin(4)) - t39 * t41) * t15) - m(7) * (g(1) * (t1 * t31 - pkin(5) * t35 + (-t13 * t2 + t3 * t31) * rSges(7,1) + (-t13 * t3 - t2 * t31) * rSges(7,2) + t19) + g(2) * (t1 * t34 + pkin(5) * t32 + (t16 * t2 + t3 * t34) * rSges(7,1) + (t16 * t3 - t2 * t34) * rSges(7,2) + t18) + g(3) * (t40 * t12 + t23) + (g(3) * (rSges(7,1) * t3 - rSges(7,2) * t2 + t1) - t39 * t40) * t15);
U  = t4;
