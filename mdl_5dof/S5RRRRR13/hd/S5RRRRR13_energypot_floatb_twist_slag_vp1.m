% Calculate potential energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:17
% EndTime: 2024-09-27 17:30:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (213->82), mult. (128->86), div. (0->0), fcn. (106->16), ass. (0->33)
t44 = pkin(9) + pkin(10) + rSges(6,3);
t43 = rSges(5,3) + pkin(9);
t25 = cos(pkin(5));
t26 = sin(qJ(4));
t41 = t25 * t26;
t28 = cos(qJ(4));
t40 = t25 * t28;
t23 = qJ(1) + qJ(2);
t22 = qJ(4) + qJ(5);
t39 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t38 = t27 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t37 = t29 * pkin(1) + r_base(1);
t36 = pkin(7) + t39;
t16 = sin(t23);
t35 = pkin(2) * t16 + t38;
t18 = cos(t23);
t34 = pkin(2) * t18 + t37;
t33 = pkin(8) + t36;
t32 = cos(t22) * rSges(6,1) - sin(t22) * rSges(6,2) + t28 * pkin(4) + pkin(3);
t24 = sin(pkin(5));
t13 = pkin(5) + t22;
t4 = sin(t13) / 0.2e1;
t14 = pkin(5) - t22;
t5 = cos(t14) / 0.2e1;
t6 = sin(t14);
t7 = cos(t13);
t31 = (t4 - t6 / 0.2e1) * rSges(6,1) + (t5 + t7 / 0.2e1) * rSges(6,2) + pkin(4) * t41 - t44 * t24;
t19 = qJ(3) + t23;
t11 = cos(t19);
t10 = sin(t19);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * (t18 * rSges(3,1) - t16 * rSges(3,2) + t37) + g(2) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t38) + g(3) * (rSges(3,3) + t36)) - m(4) * (g(1) * (t11 * rSges(4,1) - t10 * rSges(4,2) + t34) + g(2) * (t10 * rSges(4,1) + t11 * rSges(4,2) + t35) + g(3) * (rSges(4,3) + t33)) - m(5) * (g(1) * (t11 * pkin(3) + (-t10 * t41 + t11 * t28) * rSges(5,1) + (-t10 * t40 - t11 * t26) * rSges(5,2) + t34) + g(2) * (t10 * pkin(3) + (t10 * t28 + t11 * t41) * rSges(5,1) + (-t10 * t26 + t11 * t40) * rSges(5,2) + t35) + g(3) * (t43 * t25 + t33) + (g(3) * (rSges(5,1) * t26 + rSges(5,2) * t28) + (g(1) * t10 - g(2) * t11) * t43) * t24) - m(6) * (g(1) * (-t31 * t10 + t32 * t11 + t34) + g(2) * (t32 * t10 + t31 * t11 + t35) + g(3) * (t24 * t26 * pkin(4) + (t5 - t7 / 0.2e1) * rSges(6,1) + (t4 + t6 / 0.2e1) * rSges(6,2) + t44 * t25 + t33));
U = t1;
