% Calculate potential energy for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:50
% EndTime: 2019-12-31 18:18:50
% DurationCPUTime: 0.34s
% Computational Cost: add. (167->77), mult. (117->82), div. (0->0), fcn. (97->10), ass. (0->31)
t39 = rSges(6,3) + pkin(7);
t12 = qJ(3) + pkin(9);
t5 = sin(t12);
t7 = cos(t12);
t38 = rSges(5,1) * t7 - rSges(5,2) * t5;
t37 = t7 * pkin(4);
t15 = sin(qJ(5));
t13 = qJ(1) + pkin(8);
t6 = sin(t13);
t34 = t6 * t15;
t18 = cos(qJ(5));
t33 = t6 * t18;
t8 = cos(t13);
t32 = t8 * t15;
t31 = t8 * t18;
t30 = rSges(4,3) + pkin(6);
t28 = pkin(5) + r_base(3);
t17 = sin(qJ(1));
t27 = t17 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t26 = t20 * pkin(1) + r_base(1);
t19 = cos(qJ(3));
t4 = t19 * pkin(3) + pkin(2);
t25 = t8 * t4 + t26;
t24 = qJ(2) + t28;
t14 = -qJ(4) - pkin(6);
t23 = t8 * t14 + t6 * t4 + t27;
t16 = sin(qJ(3));
t22 = t16 * pkin(3) + t24;
t21 = rSges(4,1) * t19 - rSges(4,2) * t16 + pkin(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(1) * (t8 * rSges(3,1) - t6 * rSges(3,2) + t26) + g(2) * (t6 * rSges(3,1) + t8 * rSges(3,2) + t27) + g(3) * (rSges(3,3) + t24)) - m(4) * (g(1) * t26 + g(2) * t27 + g(3) * (t16 * rSges(4,1) + t19 * rSges(4,2) + t24) + (g(1) * t21 - g(2) * t30) * t8 + (g(1) * t30 + g(2) * t21) * t6) - m(5) * (g(1) * (t38 * t8 + t25) + g(2) * (-t8 * rSges(5,3) + t23) + g(3) * (t5 * rSges(5,1) + t7 * rSges(5,2) + t22) + (g(1) * (rSges(5,3) - t14) + g(2) * t38) * t6) - m(6) * (g(1) * (t8 * t37 - t6 * t14 + (t7 * t31 + t34) * rSges(6,1) + (-t7 * t32 + t33) * rSges(6,2) + t25) + g(2) * (t6 * t37 + (t7 * t33 - t32) * rSges(6,1) + (-t7 * t34 - t31) * rSges(6,2) + t23) + g(3) * (-t39 * t7 + t22) + (g(3) * (rSges(6,1) * t18 - rSges(6,2) * t15 + pkin(4)) + (g(1) * t8 + g(2) * t6) * t39) * t5);
U = t1;
