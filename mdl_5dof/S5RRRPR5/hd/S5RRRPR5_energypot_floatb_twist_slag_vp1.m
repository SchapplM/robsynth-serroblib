% Calculate potential energy for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:12:56
% EndTime: 2019-12-31 21:12:57
% DurationCPUTime: 0.35s
% Computational Cost: add. (168->81), mult. (130->88), div. (0->0), fcn. (110->10), ass. (0->34)
t40 = rSges(6,3) + pkin(8);
t22 = -pkin(7) - pkin(6);
t15 = qJ(2) + qJ(3);
t9 = pkin(9) + t15;
t5 = sin(t9);
t39 = rSges(5,2) * t5;
t21 = cos(qJ(1));
t6 = cos(t9);
t38 = t21 * t6;
t37 = rSges(3,3) + pkin(6);
t20 = cos(qJ(2));
t8 = t20 * pkin(2) + pkin(1);
t16 = sin(qJ(5));
t18 = sin(qJ(1));
t35 = t18 * t16;
t19 = cos(qJ(5));
t34 = t18 * t19;
t33 = t21 * t16;
t32 = t21 * t19;
t31 = rSges(4,3) - t22;
t30 = pkin(5) + r_base(3);
t11 = cos(t15);
t3 = pkin(3) * t11 + t8;
t29 = t21 * t3 + r_base(1);
t14 = -qJ(4) + t22;
t28 = t21 * t14 + t18 * t3 + r_base(2);
t17 = sin(qJ(2));
t27 = t17 * pkin(2) + t30;
t10 = sin(t15);
t26 = pkin(3) * t10 + t27;
t25 = rSges(3,1) * t20 - rSges(3,2) * t17 + pkin(1);
t24 = rSges(4,1) * t11 - rSges(4,2) * t10 + t8;
t23 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(3) * (t17 * rSges(3,1) + t20 * rSges(3,2) + t30) + (g(1) * t25 - g(2) * t37) * t21 + (g(1) * t37 + g(2) * t25) * t18 + t23) - m(4) * (g(3) * (rSges(4,1) * t10 + rSges(4,2) * t11 + t27) + (g(1) * t24 - g(2) * t31) * t21 + (g(1) * t31 + g(2) * t24) * t18 + t23) - m(5) * (g(1) * (rSges(5,1) * t38 - t21 * t39 + t29) + g(2) * (-t21 * rSges(5,3) + t28) + g(3) * (t5 * rSges(5,1) + t6 * rSges(5,2) + t26) + (g(1) * (rSges(5,3) - t14) + g(2) * (rSges(5,1) * t6 - t39)) * t18) - m(6) * (g(1) * (pkin(4) * t38 - t18 * t14 + (t6 * t32 + t35) * rSges(6,1) + (-t6 * t33 + t34) * rSges(6,2) + t29) + g(2) * (t18 * t6 * pkin(4) + (t6 * t34 - t33) * rSges(6,1) + (-t6 * t35 - t32) * rSges(6,2) + t28) + g(3) * (-t40 * t6 + t26) + (g(3) * (rSges(6,1) * t19 - rSges(6,2) * t16 + pkin(4)) + (g(1) * t21 + g(2) * t18) * t40) * t5);
U = t1;
