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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:18
% EndTime: 2020-01-03 11:41:19
% DurationCPUTime: 0.52s
% Computational Cost: add. (170->87), mult. (180->95), div. (0->0), fcn. (168->10), ass. (0->37)
t47 = rSges(4,3) + pkin(6);
t15 = -qJ(4) - pkin(6);
t46 = rSges(5,3) - t15;
t45 = rSges(6,3) + pkin(7) - t15;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t44 = g(2) * t17 - g(3) * t19;
t16 = sin(qJ(3));
t43 = -pkin(3) * t16 - qJ(2);
t42 = g(3) * pkin(1);
t14 = cos(pkin(8));
t40 = g(2) * t14;
t38 = g(3) * t14;
t36 = t14 * pkin(2);
t18 = cos(qJ(3));
t5 = t18 * pkin(3) + pkin(2);
t34 = t16 * t19;
t33 = t17 * t16;
t32 = t17 * t18;
t31 = t18 * t19;
t28 = -rSges(3,3) - qJ(2);
t27 = pkin(5) + r_base(1);
t12 = qJ(3) + pkin(9);
t26 = t17 * pkin(1) + r_base(2);
t13 = sin(pkin(8));
t25 = rSges(3,1) * t14 - rSges(3,2) * t13;
t6 = sin(t12);
t7 = cos(t12);
t24 = rSges(5,1) * t7 - rSges(5,2) * t6 + t5;
t8 = qJ(5) + t12;
t3 = sin(t8);
t4 = cos(t8);
t23 = rSges(6,1) * t4 - rSges(6,2) * t3 + pkin(4) * t7 + t5;
t22 = -rSges(6,1) * t3 - rSges(6,2) * t4 - pkin(4) * t6 + t43;
t21 = g(2) * t26 + g(3) * r_base(3);
t20 = -rSges(5,1) * t6 - rSges(5,2) * t7 + t43;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t27) + g(2) * (t17 * rSges(2,1) + rSges(2,2) * t19 + r_base(2)) + g(3) * (-rSges(2,1) * t19 + t17 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,1) * t13 + rSges(3,2) * t14 + t27) + (g(2) * t25 + g(3) * t28) * t17 + (g(2) * t28 + g(3) * (-pkin(1) - t25)) * t19 + t21) - m(4) * (g(1) * (-t47 * t14 + t27) + g(2) * (t17 * t36 - t19 * qJ(2) + (t14 * t32 - t34) * rSges(4,1) + (-t14 * t33 - t31) * rSges(4,2) + t26) + g(3) * (-t17 * qJ(2) + r_base(3) + (-t14 * t31 - t33) * rSges(4,1) + (t14 * t34 - t32) * rSges(4,2) + (-t36 - pkin(1)) * t19) + (g(1) * (rSges(4,1) * t18 - rSges(4,2) * t16 + pkin(2)) + t44 * t47) * t13) - m(5) * (g(1) * (-t46 * t14 + t27) + (g(3) * t20 + t24 * t40) * t17 + (g(2) * t20 - t24 * t38 - t42) * t19 + (g(1) * t24 + t44 * t46) * t13 + t21) - m(6) * (g(1) * (-t45 * t14 + t27) + (g(3) * t22 + t23 * t40) * t17 + (g(2) * t22 - t23 * t38 - t42) * t19 + (g(1) * t23 + t44 * t45) * t13 + t21);
U = t1;
