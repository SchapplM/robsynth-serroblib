% Calculate potential energy for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:01
% EndTime: 2019-12-31 18:01:01
% DurationCPUTime: 0.20s
% Computational Cost: add. (138->72), mult. (126->71), div. (0->0), fcn. (120->8), ass. (0->27)
t34 = rSges(6,3) + pkin(7);
t15 = sin(pkin(8));
t18 = sin(qJ(1));
t33 = t18 * t15;
t20 = cos(qJ(1));
t32 = t20 * t15;
t31 = pkin(5) + r_base(3);
t30 = pkin(8) + qJ(4);
t29 = t18 * pkin(1) + r_base(2);
t28 = -qJ(3) + t31;
t27 = t20 * pkin(1) + t18 * qJ(2) + r_base(1);
t26 = cos(t30);
t25 = sin(t30);
t24 = -pkin(6) + t28;
t16 = cos(pkin(8));
t11 = t16 * pkin(3) + pkin(2);
t23 = pkin(3) * t33 + t20 * t11 + t27;
t22 = -t20 * qJ(2) + t29;
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t21 = -rSges(6,1) * t19 + rSges(6,2) * t17 - pkin(4);
t5 = t18 * t11;
t4 = t18 * t16 - t32;
t3 = -t20 * t16 - t33;
t2 = -t18 * t26 + t20 * t25;
t1 = -t18 * t25 - t20 * t26;
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(1) * (t20 * rSges(3,1) + t18 * rSges(3,3) + t27) + g(2) * (t18 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t20 + t29) + g(3) * (rSges(3,2) + t31)) - m(4) * (g(1) * (-t3 * rSges(4,1) + t4 * rSges(4,2) + t20 * pkin(2) + t27) + g(2) * (t4 * rSges(4,1) + t3 * rSges(4,2) + t18 * pkin(2) + t22) + g(3) * (-rSges(4,3) + t28)) - m(5) * (g(1) * (-t1 * rSges(5,1) - t2 * rSges(5,2) + t23) + g(2) * (-t2 * rSges(5,1) + t1 * rSges(5,2) + t5 + (-pkin(3) * t15 - qJ(2)) * t20 + t29) + g(3) * (-rSges(5,3) + t24)) - m(6) * (g(1) * t23 + g(2) * (-pkin(3) * t32 + t22 + t5) + g(3) * (-t17 * rSges(6,1) - t19 * rSges(6,2) + t24) + (g(1) * t34 + g(2) * t21) * t2 + (g(1) * t21 - g(2) * t34) * t1);
U = t6;
