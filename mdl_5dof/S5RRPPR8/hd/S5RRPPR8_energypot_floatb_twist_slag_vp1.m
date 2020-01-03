% Calculate potential energy for
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:02
% EndTime: 2019-12-31 19:38:02
% DurationCPUTime: 0.36s
% Computational Cost: add. (138->86), mult. (179->96), div. (0->0), fcn. (167->8), ass. (0->28)
t13 = sin(pkin(8));
t16 = sin(qJ(2));
t36 = t16 * t13;
t17 = sin(qJ(1));
t35 = t16 * t17;
t18 = cos(qJ(2));
t34 = t17 * t18;
t33 = -pkin(7) - qJ(4) - rSges(6,3);
t32 = qJ(3) * t16;
t31 = -qJ(4) - rSges(5,3);
t30 = pkin(5) + r_base(3);
t29 = t17 * pkin(1) + r_base(2);
t28 = t16 * pkin(2) + t30;
t19 = cos(qJ(1));
t27 = t19 * pkin(1) + t17 * pkin(6) + r_base(1);
t26 = pkin(2) * t34 + t17 * t32 + t29;
t25 = t27 + (pkin(2) * t18 + t32) * t19;
t14 = cos(pkin(8));
t24 = -t18 * t13 + t16 * t14;
t23 = t18 * t14 + t36;
t22 = g(1) * t25 + g(2) * t26;
t21 = t23 * rSges(5,1) + t24 * rSges(5,2) + t18 * pkin(3);
t5 = t14 * pkin(4) + pkin(3);
t12 = pkin(8) + qJ(5);
t6 = sin(t12);
t7 = cos(t12);
t20 = t18 * t5 + pkin(4) * t36 + (t16 * t6 + t18 * t7) * rSges(6,1) + (t16 * t7 - t18 * t6) * rSges(6,2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t17 * rSges(3,3) + t27) + g(2) * (rSges(3,1) * t34 - rSges(3,2) * t35 + t29) + g(3) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t30) + (g(1) * (rSges(3,1) * t18 - rSges(3,2) * t16) + g(2) * (-rSges(3,3) - pkin(6))) * t19) - m(4) * (g(1) * (t17 * rSges(4,2) + t25) + g(2) * (rSges(4,1) * t34 + rSges(4,3) * t35 + t26) + g(3) * (t16 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t18 + t28) + (g(1) * (rSges(4,1) * t18 + rSges(4,3) * t16) + g(2) * (-rSges(4,2) - pkin(6))) * t19) - m(5) * (g(3) * (t24 * rSges(5,1) - t23 * rSges(5,2) + t16 * pkin(3) - t18 * qJ(3) + t28) + (g(1) * t31 + g(2) * t21) * t17 + (g(1) * t21 + g(2) * (-pkin(6) - t31)) * t19 + t22) - m(6) * (g(3) * ((t7 * rSges(6,1) - t6 * rSges(6,2) + t5) * t16 + (-t6 * rSges(6,1) - t7 * rSges(6,2) - t13 * pkin(4) - qJ(3)) * t18 + t28) + (g(1) * t33 + g(2) * t20) * t17 + (g(1) * t20 + g(2) * (-pkin(6) - t33)) * t19 + t22);
U = t1;
