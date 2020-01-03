% Calculate potential energy for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:48
% EndTime: 2019-12-31 19:28:49
% DurationCPUTime: 0.31s
% Computational Cost: add. (155->80), mult. (144->88), div. (0->0), fcn. (126->8), ass. (0->28)
t36 = rSges(3,3) + pkin(6);
t35 = pkin(7) + rSges(6,3);
t13 = qJ(2) + pkin(8);
t10 = sin(t13);
t20 = cos(qJ(1));
t34 = t10 * t20;
t11 = cos(t13);
t33 = t20 * t11;
t32 = qJ(4) * t10;
t31 = pkin(5) + r_base(3);
t19 = cos(qJ(2));
t9 = pkin(2) * t19 + pkin(1);
t30 = t20 * t9 + r_base(1);
t14 = -qJ(3) - pkin(6);
t17 = sin(qJ(1));
t29 = t20 * t14 + t17 * t9 + r_base(2);
t16 = sin(qJ(2));
t28 = t16 * pkin(2) + t31;
t27 = pkin(3) * t33 + t20 * t32 + t30;
t26 = t10 * pkin(3) + t28;
t25 = t29 + (pkin(3) * t11 + t32) * t17;
t15 = sin(qJ(5));
t18 = cos(qJ(5));
t24 = t10 * t18 - t11 * t15;
t23 = t10 * t15 + t11 * t18;
t22 = rSges(3,1) * t19 - rSges(3,2) * t16 + pkin(1);
t21 = rSges(6,1) * t23 + rSges(6,2) * t24 + t11 * pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t20 - rSges(2,2) * t17 + r_base(1)) + g(2) * (rSges(2,1) * t17 + rSges(2,2) * t20 + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t16 + rSges(3,2) * t19 + t31) + (g(1) * t22 - g(2) * t36) * t20 + (g(1) * t36 + g(2) * t22) * t17) - m(4) * (g(1) * (rSges(4,1) * t33 - rSges(4,2) * t34 + t30) + g(2) * (-t20 * rSges(4,3) + t29) + g(3) * (rSges(4,1) * t10 + rSges(4,2) * t11 + t28) + (g(1) * (rSges(4,3) - t14) + g(2) * (rSges(4,1) * t11 - rSges(4,2) * t10)) * t17) - m(5) * (g(1) * (rSges(5,1) * t33 + rSges(5,3) * t34 + t27) + g(2) * (-t20 * rSges(5,2) + t25) + g(3) * (t10 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t11 + t26) + (g(1) * (rSges(5,2) - t14) + g(2) * (rSges(5,1) * t11 + rSges(5,3) * t10)) * t17) - m(6) * (g(1) * t27 + g(2) * t25 + g(3) * (rSges(6,1) * t24 - rSges(6,2) * t23 + t10 * pkin(4) - t11 * qJ(4) + t26) + (g(1) * t21 + g(2) * t35) * t20 + (g(1) * (-t14 - t35) + g(2) * t21) * t17);
U = t1;
