% Calculate potential energy for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:16
% DurationCPUTime: 0.31s
% Computational Cost: add. (143->77), mult. (128->81), div. (0->0), fcn. (104->6), ass. (0->25)
t33 = rSges(6,3) + qJ(5);
t32 = -rSges(6,1) - pkin(4);
t13 = pkin(7) + qJ(3);
t10 = sin(t13);
t18 = cos(qJ(1));
t31 = t10 * t18;
t11 = cos(t13);
t30 = t11 * t18;
t29 = qJ(4) * t10;
t28 = rSges(3,3) + qJ(2);
t27 = pkin(5) + r_base(3);
t15 = cos(pkin(7));
t8 = pkin(2) * t15 + pkin(1);
t26 = t18 * t8 + r_base(1);
t16 = -pkin(6) - qJ(2);
t17 = sin(qJ(1));
t25 = t18 * t16 + t17 * t8 + r_base(2);
t14 = sin(pkin(7));
t24 = t14 * pkin(2) + t27;
t23 = pkin(3) * t30 + t18 * t29 + t26;
t22 = t10 * pkin(3) + t24;
t21 = t25 + (pkin(3) * t11 + t29) * t17;
t20 = rSges(3,1) * t15 - rSges(3,2) * t14 + pkin(1);
t19 = rSges(6,2) * t10 + t11 * t33;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t18 - rSges(2,2) * t17 + r_base(1)) + g(2) * (rSges(2,1) * t17 + rSges(2,2) * t18 + r_base(2)) + g(3) * (rSges(2,3) + t27)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t14 + rSges(3,2) * t15 + t27) + (g(1) * t20 - g(2) * t28) * t18 + (g(1) * t28 + g(2) * t20) * t17) - m(4) * (g(1) * (rSges(4,1) * t30 - rSges(4,2) * t31 + t26) + g(2) * (-t18 * rSges(4,3) + t25) + g(3) * (rSges(4,1) * t10 + rSges(4,2) * t11 + t24) + (g(1) * (rSges(4,3) - t16) + g(2) * (rSges(4,1) * t11 - rSges(4,2) * t10)) * t17) - m(5) * (g(1) * (-rSges(5,2) * t30 + rSges(5,3) * t31 + t23) + g(2) * (-t18 * rSges(5,1) + t21) + g(3) * (-t10 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t11 + t22) + (g(1) * (rSges(5,1) - t16) + g(2) * (-rSges(5,2) * t11 + rSges(5,3) * t10)) * t17) - m(6) * (g(1) * t23 + g(2) * t21 + g(3) * ((-rSges(6,2) - qJ(4)) * t11 + t33 * t10 + t22) + (g(1) * t19 + g(2) * t32) * t18 + (g(1) * (-t16 - t32) + g(2) * t19) * t17);
U = t1;
