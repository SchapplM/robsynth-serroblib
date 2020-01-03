% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:20
% EndTime: 2020-01-03 11:30:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (170->87), mult. (180->97), div. (0->0), fcn. (168->10), ass. (0->35)
t17 = -pkin(6) - qJ(3);
t45 = rSges(5,3) - t17;
t44 = rSges(6,3) + pkin(7) - t17;
t18 = sin(qJ(1));
t19 = cos(qJ(1));
t43 = g(2) * t18 - g(3) * t19;
t42 = rSges(4,3) + qJ(3);
t13 = sin(pkin(9));
t41 = -pkin(3) * t13 - qJ(2);
t40 = g(3) * pkin(1);
t15 = cos(pkin(9));
t5 = t15 * pkin(3) + pkin(2);
t16 = cos(pkin(8));
t38 = g(2) * t16;
t36 = g(3) * t16;
t34 = t16 * t19;
t33 = t18 * t13;
t32 = t18 * t15;
t29 = -rSges(3,3) - qJ(2);
t27 = pkin(5) + r_base(1);
t12 = pkin(9) + qJ(4);
t26 = t18 * pkin(1) + r_base(2);
t14 = sin(pkin(8));
t25 = rSges(3,1) * t16 - rSges(3,2) * t14;
t6 = sin(t12);
t7 = cos(t12);
t24 = rSges(5,1) * t7 - rSges(5,2) * t6 + t5;
t8 = qJ(5) + t12;
t3 = sin(t8);
t4 = cos(t8);
t23 = rSges(6,1) * t4 - rSges(6,2) * t3 + pkin(4) * t7 + t5;
t22 = -rSges(6,1) * t3 - rSges(6,2) * t4 - pkin(4) * t6 + t41;
t21 = g(2) * t26 + g(3) * r_base(3);
t20 = -rSges(5,1) * t6 - rSges(5,2) * t7 + t41;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t27) + g(2) * (t18 * rSges(2,1) + rSges(2,2) * t19 + r_base(2)) + g(3) * (-rSges(2,1) * t19 + t18 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,1) * t14 + rSges(3,2) * t16 + t27) + (g(2) * t25 + g(3) * t29) * t18 + (g(2) * t29 + g(3) * (-pkin(1) - t25)) * t19 + t21) - m(4) * (g(1) * (-t16 * t42 + t27) + g(2) * (t18 * t16 * pkin(2) - t19 * qJ(2) + (-t13 * t19 + t16 * t32) * rSges(4,1) + (-t15 * t19 - t16 * t33) * rSges(4,2) + t26) + g(3) * (-pkin(2) * t34 - t19 * pkin(1) - t18 * qJ(2) + r_base(3) + (-t15 * t34 - t33) * rSges(4,1) + (t13 * t34 - t32) * rSges(4,2)) + (g(1) * (rSges(4,1) * t15 - rSges(4,2) * t13 + pkin(2)) + t43 * t42) * t14) - m(5) * (g(1) * (-t16 * t45 + t27) + (g(3) * t20 + t24 * t38) * t18 + (g(2) * t20 - t24 * t36 - t40) * t19 + (g(1) * t24 + t43 * t45) * t14 + t21) - m(6) * (g(1) * (-t16 * t44 + t27) + (g(3) * t22 + t23 * t38) * t18 + (g(2) * t22 - t23 * t36 - t40) * t19 + (g(1) * t23 + t43 * t44) * t14 + t21);
U = t1;
