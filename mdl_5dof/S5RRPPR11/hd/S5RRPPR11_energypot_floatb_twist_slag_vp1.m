% Calculate potential energy for
% S5RRPPR11
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:19
% EndTime: 2019-12-31 19:46:20
% DurationCPUTime: 0.44s
% Computational Cost: add. (134->92), mult. (172->107), div. (0->0), fcn. (156->8), ass. (0->31)
t43 = rSges(6,3) + pkin(7) + qJ(4);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t42 = g(1) * t19 + g(2) * t17;
t41 = rSges(5,3) + qJ(4);
t16 = sin(qJ(2));
t38 = t16 * t17;
t37 = t16 * t19;
t13 = sin(pkin(8));
t36 = t17 * t13;
t14 = cos(pkin(8));
t35 = t17 * t14;
t18 = cos(qJ(2));
t34 = t17 * t18;
t33 = t19 * t13;
t32 = t19 * t14;
t30 = qJ(3) * t16;
t28 = pkin(5) + r_base(3);
t27 = t17 * pkin(1) + r_base(2);
t26 = t16 * t36;
t25 = t16 * t33;
t24 = t16 * pkin(2) + t28;
t23 = t19 * pkin(1) + t17 * pkin(6) + r_base(1);
t22 = pkin(2) * t34 + t17 * t30 + t27;
t21 = t23 + (pkin(2) * t18 + t30) * t19;
t20 = -t19 * pkin(6) + t22;
t12 = pkin(8) + qJ(5);
t7 = cos(t12);
t6 = sin(t12);
t5 = t14 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(1) * (t17 * rSges(3,3) + t23) + g(2) * (rSges(3,1) * t34 - rSges(3,2) * t38 + t27) + g(3) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t28) + (g(1) * (rSges(3,1) * t18 - rSges(3,2) * t16) + g(2) * (-rSges(3,3) - pkin(6))) * t19) - m(4) * (g(1) * (t17 * rSges(4,1) + t21) + g(2) * (-rSges(4,2) * t34 + rSges(4,3) * t38 + t22) + g(3) * (-t16 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t18 + t24) + (g(1) * (-rSges(4,2) * t18 + rSges(4,3) * t16) + g(2) * (-rSges(4,1) - pkin(6))) * t19) - m(5) * (g(1) * (t17 * pkin(3) + (t25 + t35) * rSges(5,1) + (t16 * t32 - t36) * rSges(5,2) + t21) + g(2) * (-t19 * pkin(3) + (t26 - t32) * rSges(5,1) + (t16 * t35 + t33) * rSges(5,2) + t20) + g(3) * (t41 * t16 + t24) + (g(3) * (-rSges(5,1) * t13 - rSges(5,2) * t14 - qJ(3)) + t42 * t41) * t18) - m(6) * (g(1) * (t17 * t5 + pkin(4) * t25 + (t17 * t7 + t37 * t6) * rSges(6,1) + (-t17 * t6 + t37 * t7) * rSges(6,2) + t21) + g(2) * (-t19 * t5 + pkin(4) * t26 + (-t19 * t7 + t38 * t6) * rSges(6,1) + (t19 * t6 + t38 * t7) * rSges(6,2) + t20) + g(3) * (t43 * t16 + t24) + (g(3) * (-rSges(6,1) * t6 - rSges(6,2) * t7 - pkin(4) * t13 - qJ(3)) + t42 * t43) * t18);
U = t1;
