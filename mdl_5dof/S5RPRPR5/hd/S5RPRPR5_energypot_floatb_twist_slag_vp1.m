% Calculate potential energy for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:48
% EndTime: 2022-01-23 09:24:48
% DurationCPUTime: 0.51s
% Computational Cost: add. (170->97), mult. (173->105), div. (0->0), fcn. (161->10), ass. (0->36)
t45 = rSges(4,3) + pkin(6);
t20 = qJ(4) + pkin(6);
t44 = rSges(6,3) + pkin(7) + t20;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t43 = g(1) * t24 + g(2) * t22;
t23 = cos(qJ(3));
t38 = t23 * pkin(3);
t17 = qJ(3) + pkin(9);
t9 = cos(t17);
t2 = pkin(4) * t9 + pkin(2) + t38;
t10 = qJ(5) + t17;
t5 = sin(t10);
t6 = cos(t10);
t42 = rSges(6,1) * t6 - rSges(6,2) * t5 + t2;
t21 = sin(qJ(3));
t39 = t21 * pkin(3);
t19 = cos(pkin(8));
t37 = pkin(2) * t19 + pkin(1);
t18 = sin(pkin(8));
t36 = rSges(3,2) * t18;
t35 = t19 * t22;
t34 = t19 * t24;
t32 = pkin(5) + r_base(3);
t31 = t22 * qJ(2) + r_base(1);
t30 = t22 * pkin(1) + r_base(2);
t29 = t18 * pkin(2) + t32;
t28 = t24 * pkin(1) + t31;
t27 = rSges(4,1) * t23 - rSges(4,2) * t21;
t26 = t21 * rSges(4,1) + t23 * rSges(4,2);
t25 = t45 * t18 + t27 * t19 + t37;
t8 = sin(t17);
t7 = qJ(2) + t39;
t3 = pkin(4) * t8 + t39;
t1 = t20 * t18 + t19 * t38 + t37;
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t32)) - m(3) * (g(1) * (t22 * rSges(3,3) + t28) + g(2) * (rSges(3,1) * t35 - t22 * t36 + t30) + g(3) * (t18 * rSges(3,1) + t19 * rSges(3,2) + t32) + (g(1) * (rSges(3,1) * t19 - t36) + g(2) * (-rSges(3,3) - qJ(2))) * t24) - m(4) * (g(1) * (t26 * t22 + t25 * t24 + t31) + g(2) * (r_base(2) + (-qJ(2) - t26) * t24 + t25 * t22) + g(3) * (t27 * t18 - t45 * t19 + t29)) - m(5) * (g(1) * (t1 * t24 + t7 * t22 + r_base(1) + (t22 * t8 + t9 * t34) * rSges(5,1) + (t22 * t9 - t8 * t34) * rSges(5,2)) + g(2) * (t1 * t22 - t7 * t24 + r_base(2) + (-t24 * t8 + t9 * t35) * rSges(5,1) + (-t24 * t9 - t8 * t35) * rSges(5,2)) + g(3) * (t29 + (-rSges(5,3) - t20) * t19) + (g(3) * (rSges(5,1) * t9 - rSges(5,2) * t8 + t38) + t43 * rSges(5,3)) * t18) - m(6) * (g(1) * (t2 * t34 + t22 * t3 + (t22 * t5 + t6 * t34) * rSges(6,1) + (t22 * t6 - t5 * t34) * rSges(6,2) + t28) + g(2) * (t30 + t42 * t35 + (-t5 * rSges(6,1) - t6 * rSges(6,2) - qJ(2) - t3) * t24) + g(3) * (-t44 * t19 + t32) + (g(3) * t42 + t43 * t44) * t18);
U = t4;
