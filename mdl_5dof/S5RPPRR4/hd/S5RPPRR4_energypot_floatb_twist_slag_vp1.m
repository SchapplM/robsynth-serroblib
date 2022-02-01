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
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:16:11
% EndTime: 2022-01-23 09:16:11
% DurationCPUTime: 0.49s
% Computational Cost: add. (170->104), mult. (170->118), div. (0->0), fcn. (158->10), ass. (0->37)
t22 = qJ(3) + pkin(6);
t43 = rSges(6,3) + pkin(7) + t22;
t23 = sin(qJ(1));
t24 = cos(qJ(1));
t25 = g(1) * t24 + g(2) * t23;
t18 = sin(pkin(9));
t40 = t18 * pkin(3);
t20 = cos(pkin(9));
t8 = t20 * pkin(3) + pkin(2);
t19 = sin(pkin(8));
t39 = rSges(3,2) * t19;
t21 = cos(pkin(8));
t38 = t21 * t23;
t37 = t21 * t24;
t17 = pkin(9) + qJ(4);
t10 = cos(t17);
t36 = t23 * t10;
t35 = t23 * t18;
t34 = t23 * t20;
t33 = t24 * t10;
t32 = t24 * t18;
t31 = t24 * t20;
t29 = pkin(5) + r_base(3);
t28 = t23 * qJ(2) + r_base(1);
t27 = t24 * pkin(1) + t28;
t26 = -t24 * qJ(2) + r_base(2);
t14 = t23 * pkin(1);
t11 = qJ(5) + t17;
t9 = sin(t17);
t7 = cos(t11);
t6 = sin(t11);
t5 = qJ(2) + t40;
t4 = pkin(2) * t21 + t19 * qJ(3) + pkin(1);
t3 = pkin(4) * t9 + t40;
t2 = pkin(4) * t10 + t8;
t1 = t22 * t19 + t8 * t21 + pkin(1);
t12 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t29)) - m(3) * (g(1) * (t23 * rSges(3,3) + t27) + g(2) * (rSges(3,1) * t38 - t23 * t39 + t14 + r_base(2)) + g(3) * (t19 * rSges(3,1) + t21 * rSges(3,2) + t29) + (g(1) * (rSges(3,1) * t21 - t39) + g(2) * (-rSges(3,3) - qJ(2))) * t24) - m(4) * (g(1) * (t4 * t24 + (t21 * t31 + t35) * rSges(4,1) + (-t21 * t32 + t34) * rSges(4,2) + t28) + g(2) * (t4 * t23 + (t21 * t34 - t32) * rSges(4,1) + (-t21 * t35 - t31) * rSges(4,2) + t26) + g(3) * (t29 + (-rSges(4,3) - qJ(3)) * t21) + (g(3) * (rSges(4,1) * t20 - rSges(4,2) * t18 + pkin(2)) + t25 * rSges(4,3)) * t19) - m(5) * (g(1) * (t1 * t24 + t5 * t23 + r_base(1) + (t21 * t33 + t23 * t9) * rSges(5,1) + (-t9 * t37 + t36) * rSges(5,2)) + g(2) * (t1 * t23 - t5 * t24 + r_base(2) + (t21 * t36 - t24 * t9) * rSges(5,1) + (-t9 * t38 - t33) * rSges(5,2)) + g(3) * (t29 + (-rSges(5,3) - t22) * t21) + (g(3) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t8) + t25 * rSges(5,3)) * t19) - m(6) * (g(1) * (t2 * t37 + t23 * t3 + (t23 * t6 + t7 * t37) * rSges(6,1) + (t23 * t7 - t6 * t37) * rSges(6,2) + t27) + g(2) * (t2 * t38 - t24 * t3 + t14 + (-t24 * t6 + t7 * t38) * rSges(6,1) + (-t24 * t7 - t6 * t38) * rSges(6,2) + t26) + g(3) * (-t43 * t21 + t29) + (g(3) * (rSges(6,1) * t7 - rSges(6,2) * t6 + t2) + t25 * t43) * t19);
U = t12;
