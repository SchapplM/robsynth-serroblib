% Calculate potential energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:15
% EndTime: 2019-12-31 20:03:16
% DurationCPUTime: 0.39s
% Computational Cost: add. (126->86), mult. (179->99), div. (0->0), fcn. (167->6), ass. (0->27)
t36 = -rSges(5,3) - pkin(7);
t17 = sin(qJ(4));
t18 = sin(qJ(2));
t35 = t18 * t17;
t19 = sin(qJ(1));
t34 = t18 * t19;
t21 = cos(qJ(2));
t33 = t19 * t21;
t32 = -rSges(6,3) - qJ(5) - pkin(7);
t31 = qJ(3) * t18;
t30 = pkin(5) + r_base(3);
t29 = t19 * pkin(1) + r_base(2);
t28 = t18 * pkin(2) + t30;
t22 = cos(qJ(1));
t27 = t22 * pkin(1) + t19 * pkin(6) + r_base(1);
t26 = pkin(2) * t33 + t19 * t31 + t29;
t20 = cos(qJ(4));
t6 = -t21 * t17 + t18 * t20;
t25 = t21 * t20 + t35;
t24 = t27 + (pkin(2) * t21 + t31) * t22;
t11 = t20 * pkin(4) + pkin(3);
t23 = pkin(4) * t35 + t11 * t21;
t4 = t25 * t22;
t3 = t6 * t22;
t2 = t25 * t19;
t1 = t6 * t19;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t22 * rSges(2,1) - t19 * rSges(2,2) + r_base(1)) + g(2) * (t19 * rSges(2,1) + t22 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t19 * rSges(3,3) + t27) + g(2) * (rSges(3,1) * t33 - rSges(3,2) * t34 + t29) + g(3) * (t18 * rSges(3,1) + t21 * rSges(3,2) + t30) + (g(1) * (rSges(3,1) * t21 - rSges(3,2) * t18) + g(2) * (-rSges(3,3) - pkin(6))) * t22) - m(4) * (g(1) * (t19 * rSges(4,2) + t24) + g(2) * (rSges(4,1) * t33 + rSges(4,3) * t34 + t26) + g(3) * (t18 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t21 + t28) + (g(1) * (rSges(4,1) * t21 + rSges(4,3) * t18) + g(2) * (-rSges(4,2) - pkin(6))) * t22) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t19 * t36 + t24) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + pkin(3) * t33 + t26) + g(3) * (t6 * rSges(5,1) - rSges(5,2) * t25 + t18 * pkin(3) - t21 * qJ(3) + t28) + (g(1) * pkin(3) * t21 + g(2) * (-pkin(6) - t36)) * t22) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t24) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t26) + g(3) * (t6 * rSges(6,1) - t25 * rSges(6,2) + t18 * t11 + (-pkin(4) * t17 - qJ(3)) * t21 + t28) + (g(1) * t32 + g(2) * t23) * t19 + (g(1) * t23 + g(2) * (-pkin(6) - t32)) * t22);
U = t5;
