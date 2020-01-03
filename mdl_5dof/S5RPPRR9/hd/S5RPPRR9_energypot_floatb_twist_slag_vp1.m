% Calculate potential energy for
% S5RPPRR9
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:15
% EndTime: 2019-12-31 18:02:15
% DurationCPUTime: 0.33s
% Computational Cost: add. (132->74), mult. (177->85), div. (0->0), fcn. (189->8), ass. (0->25)
t35 = rSges(6,3) + pkin(7);
t16 = cos(qJ(4));
t34 = t16 * pkin(4);
t33 = rSges(5,3) + pkin(6);
t31 = cos(qJ(1));
t30 = sin(qJ(1));
t13 = sin(qJ(5));
t29 = t13 * t16;
t15 = cos(qJ(5));
t28 = t15 * t16;
t27 = cos(pkin(8));
t26 = sin(pkin(8));
t25 = pkin(5) + r_base(3);
t24 = t31 * pkin(1) + t30 * qJ(2) + r_base(1);
t23 = -qJ(3) + t25;
t22 = t31 * pkin(2) + t24;
t14 = sin(qJ(4));
t21 = -rSges(5,1) * t16 + rSges(5,2) * t14;
t3 = -t30 * t26 - t31 * t27;
t20 = -t3 * pkin(3) + t22;
t19 = t30 * pkin(1) - t31 * qJ(2) + r_base(2);
t18 = t30 * pkin(2) + t19;
t4 = t31 * t26 - t30 * t27;
t17 = -t4 * pkin(3) + t18;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t31 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t31 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (t31 * rSges(3,1) + t30 * rSges(3,3) + t24) + g(2) * (t30 * rSges(3,1) - t31 * rSges(3,3) + t19) + g(3) * (rSges(3,2) + t25)) - m(4) * (g(1) * (-t3 * rSges(4,1) - t4 * rSges(4,2) + t22) + g(2) * (-t4 * rSges(4,1) + t3 * rSges(4,2) + t18) + g(3) * (-rSges(4,3) + t23)) - m(5) * (g(1) * t20 + g(2) * t17 + g(3) * (-t14 * rSges(5,1) - t16 * rSges(5,2) + t23) + (g(1) * t33 + g(2) * t21) * t4 + (g(1) * t21 - g(2) * t33) * t3) - m(6) * (g(1) * (-t3 * t34 + t4 * pkin(6) + (t4 * t13 - t3 * t28) * rSges(6,1) + (t4 * t15 + t3 * t29) * rSges(6,2) + t20) + g(2) * (-t4 * t34 - t3 * pkin(6) + (-t3 * t13 - t4 * t28) * rSges(6,1) + (-t3 * t15 + t4 * t29) * rSges(6,2) + t17) + g(3) * (t35 * t16 + t23) + (g(3) * (-rSges(6,1) * t15 + rSges(6,2) * t13 - pkin(4)) - (g(1) * t3 + g(2) * t4) * t35) * t14);
U = t1;
