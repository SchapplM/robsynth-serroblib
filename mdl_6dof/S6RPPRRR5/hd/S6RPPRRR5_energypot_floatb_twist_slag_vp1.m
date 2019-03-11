% Calculate potential energy for
% S6RPPRRR5
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:11
% EndTime: 2019-03-09 02:28:11
% DurationCPUTime: 0.36s
% Computational Cost: add. (151->91), mult. (143->95), div. (0->0), fcn. (119->8), ass. (0->32)
t41 = rSges(7,3) + pkin(9);
t19 = -pkin(8) - pkin(7);
t40 = -t19 - qJ(2);
t12 = qJ(4) + qJ(5);
t5 = cos(t12);
t39 = rSges(6,2) * t5;
t14 = sin(qJ(4));
t38 = pkin(4) * t14;
t15 = sin(qJ(1));
t4 = sin(t12);
t37 = t15 * t4;
t36 = -rSges(5,3) - pkin(7);
t13 = sin(qJ(6));
t34 = t15 * t13;
t16 = cos(qJ(6));
t33 = t15 * t16;
t18 = cos(qJ(1));
t32 = t18 * t13;
t31 = t18 * t16;
t30 = pkin(6) + r_base(3);
t29 = t15 * pkin(1) + r_base(2);
t28 = pkin(2) + t30;
t27 = t15 * qJ(3) + t29;
t26 = t18 * pkin(1) + t15 * qJ(2) + r_base(1);
t25 = pkin(3) + t28;
t24 = t15 * t38 + t27;
t23 = t18 * qJ(3) + t26;
t17 = cos(qJ(4));
t22 = t17 * pkin(4) + t25;
t21 = rSges(5,1) * t14 + rSges(5,2) * t17;
t20 = t15 * t19 + t18 * t38 + t23;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t18 * rSges(2,1) - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + t18 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (-t18 * rSges(3,2) + t15 * rSges(3,3) + t26) + g(2) * (-t15 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t18 + t29) + g(3) * (rSges(3,1) + t30)) - m(4) * (g(1) * (t15 * rSges(4,2) + t18 * rSges(4,3) + t23) + g(2) * (t15 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t18 + t27) + g(3) * (rSges(4,1) + t28)) - m(5) * (g(1) * t23 + g(2) * t27 + g(3) * (t17 * rSges(5,1) - t14 * rSges(5,2) + t25) + (g(1) * t36 + g(2) * t21) * t15 + (g(1) * t21 + g(2) * (-qJ(2) - t36)) * t18) - m(6) * (g(1) * (-t15 * rSges(6,3) + t20) + g(2) * (rSges(6,1) * t37 + t15 * t39 + t24) + g(3) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t22) + (g(1) * (rSges(6,1) * t4 + t39) + g(2) * (rSges(6,3) + t40)) * t18) - m(7) * (g(1) * (t18 * t4 * pkin(5) + (t31 * t4 - t34) * rSges(7,1) + (-t32 * t4 - t33) * rSges(7,2) + t20) + g(2) * (pkin(5) * t37 + (t33 * t4 + t32) * rSges(7,1) + (-t34 * t4 + t31) * rSges(7,2) + t24 + t40 * t18) + g(3) * (t41 * t4 + t22) + (g(3) * (rSges(7,1) * t16 - rSges(7,2) * t13 + pkin(5)) - (g(1) * t18 + g(2) * t15) * t41) * t5);
U  = t1;
