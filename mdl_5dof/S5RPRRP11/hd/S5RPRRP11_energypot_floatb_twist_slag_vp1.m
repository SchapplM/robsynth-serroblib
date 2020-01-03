% Calculate potential energy for
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:22
% EndTime: 2019-12-31 18:53:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (167->83), mult. (167->93), div. (0->0), fcn. (155->8), ass. (0->34)
t43 = rSges(6,1) + pkin(4);
t42 = rSges(6,3) + qJ(5);
t18 = pkin(8) + qJ(3);
t15 = sin(t18);
t23 = sin(qJ(1));
t41 = t15 * t23;
t25 = cos(qJ(1));
t40 = t15 * t25;
t16 = cos(t18);
t39 = t16 * t25;
t22 = sin(qJ(4));
t38 = t23 * t22;
t24 = cos(qJ(4));
t37 = t23 * t24;
t36 = t25 * t22;
t35 = t25 * t24;
t34 = rSges(3,3) + qJ(2);
t33 = pkin(5) + r_base(3);
t20 = cos(pkin(8));
t12 = pkin(2) * t20 + pkin(1);
t32 = t25 * t12 + r_base(1);
t19 = sin(pkin(8));
t31 = t19 * pkin(2) + t33;
t21 = -pkin(6) - qJ(2);
t30 = t23 * t12 + t25 * t21 + r_base(2);
t29 = t15 * pkin(3) + t31;
t28 = t23 * t16 * pkin(3) + pkin(7) * t41 + t30;
t27 = rSges(3,1) * t20 - rSges(3,2) * t19 + pkin(1);
t26 = pkin(3) * t39 + pkin(7) * t40 - t23 * t21 + t32;
t4 = t16 * t35 + t38;
t3 = t16 * t36 - t37;
t2 = t16 * t37 - t36;
t1 = t16 * t38 + t35;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t25 - rSges(2,2) * t23 + r_base(1)) + g(2) * (rSges(2,1) * t23 + rSges(2,2) * t25 + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t19 + rSges(3,2) * t20 + t33) + (g(1) * t27 - g(2) * t34) * t25 + (g(1) * t34 + g(2) * t27) * t23) - m(4) * (g(1) * (rSges(4,1) * t39 - rSges(4,2) * t40 + t32) + g(2) * (-t25 * rSges(4,3) + t30) + g(3) * (rSges(4,1) * t15 + rSges(4,2) * t16 + t31) + (g(1) * (rSges(4,3) - t21) + g(2) * (rSges(4,1) * t16 - rSges(4,2) * t15)) * t23) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t40 + t26) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t1 + rSges(5,3) * t41 + t28) + g(3) * ((-rSges(5,3) - pkin(7)) * t16 + (rSges(5,1) * t24 - rSges(5,2) * t22) * t15 + t29)) - m(6) * (g(1) * (t3 * t42 + t4 * t43 + t26) + g(2) * (t1 * t42 + t2 * t43 + t28) + g(3) * (t29 + (-rSges(6,2) - pkin(7)) * t16) + (g(3) * (t22 * t42 + t24 * t43) + (g(1) * t25 + g(2) * t23) * rSges(6,2)) * t15);
U = t5;
