% Calculate potential energy for
% S5RRPRP10
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:19
% EndTime: 2019-12-31 20:09:19
% DurationCPUTime: 0.38s
% Computational Cost: add. (124->87), mult. (172->97), div. (0->0), fcn. (156->6), ass. (0->31)
t43 = rSges(5,3) + pkin(7);
t42 = rSges(6,3) + qJ(5) + pkin(7);
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t41 = g(1) * t20 + g(2) * t17;
t16 = sin(qJ(2));
t37 = t16 * t17;
t15 = sin(qJ(4));
t36 = t17 * t15;
t18 = cos(qJ(4));
t35 = t17 * t18;
t19 = cos(qJ(2));
t34 = t17 * t19;
t33 = t20 * t15;
t32 = t20 * t18;
t30 = qJ(3) * t16;
t29 = pkin(5) + r_base(3);
t28 = t17 * pkin(1) + r_base(2);
t27 = t16 * t36;
t26 = t16 * t33;
t25 = t16 * pkin(2) + t29;
t24 = t20 * pkin(1) + t17 * pkin(6) + r_base(1);
t23 = pkin(2) * t34 + t17 * t30 + t28;
t22 = t24 + (pkin(2) * t19 + t30) * t20;
t21 = -t20 * pkin(6) + t23;
t9 = t18 * pkin(4) + pkin(3);
t4 = t27 - t32;
t3 = t16 * t35 + t33;
t2 = t26 + t35;
t1 = t16 * t32 - t36;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t29)) - m(3) * (g(1) * (t17 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t34 - rSges(3,2) * t37 + t28) + g(3) * (t16 * rSges(3,1) + t19 * rSges(3,2) + t29) + (g(1) * (rSges(3,1) * t19 - rSges(3,2) * t16) + g(2) * (-rSges(3,3) - pkin(6))) * t20) - m(4) * (g(1) * (t17 * rSges(4,1) + t22) + g(2) * (-rSges(4,2) * t34 + rSges(4,3) * t37 + t23) + g(3) * (-t16 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t19 + t25) + (g(1) * (-rSges(4,2) * t19 + rSges(4,3) * t16) + g(2) * (-rSges(4,1) - pkin(6))) * t20) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t17 * pkin(3) + t22) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) - t20 * pkin(3) + t21) + g(3) * (t43 * t16 + t25) + (g(3) * (-rSges(5,1) * t15 - rSges(5,2) * t18 - qJ(3)) + t41 * t43) * t19) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t26 + t17 * t9 + t22) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t27 - t20 * t9 + t21) + g(3) * (t42 * t16 + t25) + (g(3) * (-rSges(6,2) * t18 - qJ(3) + (-rSges(6,1) - pkin(4)) * t15) + t41 * t42) * t19);
U = t5;
