% Calculate potential energy for
% S5RRPRP11
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:17
% EndTime: 2019-12-31 20:12:17
% DurationCPUTime: 0.36s
% Computational Cost: add. (126->84), mult. (182->95), div. (0->0), fcn. (170->6), ass. (0->31)
t43 = rSges(6,1) + pkin(4);
t42 = rSges(6,3) + qJ(5);
t19 = sin(qJ(2));
t20 = sin(qJ(1));
t41 = t19 * t20;
t18 = sin(qJ(4));
t40 = t20 * t18;
t21 = cos(qJ(4));
t39 = t20 * t21;
t22 = cos(qJ(2));
t38 = t20 * t22;
t23 = cos(qJ(1));
t37 = t22 * t23;
t36 = t23 * t18;
t35 = t23 * t21;
t34 = qJ(3) * t19;
t33 = pkin(5) + r_base(3);
t32 = t20 * pkin(1) + r_base(2);
t31 = t19 * pkin(2) + t33;
t30 = t23 * pkin(1) + t20 * pkin(6) + r_base(1);
t29 = pkin(2) * t38 + t20 * t34 + t32;
t28 = t19 * pkin(7) + t31;
t27 = g(1) * t23 + g(2) * t20;
t26 = pkin(2) * t37 + t23 * t34 + t30;
t25 = t20 * pkin(3) + pkin(7) * t37 + t26;
t24 = t29 + pkin(7) * t38 + (-pkin(3) - pkin(6)) * t23;
t4 = t19 * t40 - t35;
t3 = t19 * t39 + t36;
t2 = t19 * t36 + t39;
t1 = -t19 * t35 + t40;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t20 * rSges(3,3) + t30) + g(2) * (rSges(3,1) * t38 - rSges(3,2) * t41 + t32) + g(3) * (t19 * rSges(3,1) + t22 * rSges(3,2) + t33) + (g(1) * (rSges(3,1) * t22 - rSges(3,2) * t19) + g(2) * (-rSges(3,3) - pkin(6))) * t23) - m(4) * (g(1) * (t20 * rSges(4,1) + t26) + g(2) * (-rSges(4,2) * t38 + rSges(4,3) * t41 + t29) + g(3) * (-t19 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t22 + t31) + (g(1) * (-rSges(4,2) * t22 + rSges(4,3) * t19) + g(2) * (-rSges(4,1) - pkin(6))) * t23) - m(5) * (g(1) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t25) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t24) + g(3) * (t19 * rSges(5,3) + t28) + (g(3) * (-rSges(5,1) * t18 - rSges(5,2) * t21 - qJ(3)) + t27 * rSges(5,3)) * t22) - m(6) * (g(1) * (t42 * t1 + t2 * t43 + t25) + g(2) * (-t42 * t3 + t4 * t43 + t24) + g(3) * (t19 * rSges(6,2) + t28) + (g(3) * (-t18 * t43 + t42 * t21 - qJ(3)) + t27 * rSges(6,2)) * t22);
U = t5;
