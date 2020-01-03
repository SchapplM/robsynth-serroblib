% Calculate potential energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:04
% EndTime: 2019-12-31 18:16:04
% DurationCPUTime: 0.29s
% Computational Cost: add. (101->73), mult. (117->73), div. (0->0), fcn. (93->4), ass. (0->20)
t28 = rSges(6,1) + pkin(4);
t12 = cos(qJ(3));
t27 = rSges(4,2) * t12;
t26 = rSges(6,2) * t12;
t10 = sin(qJ(3));
t11 = sin(qJ(1));
t25 = t10 * t11;
t24 = qJ(4) * t12;
t23 = -rSges(6,3) - qJ(5);
t22 = pkin(5) + r_base(3);
t21 = t11 * pkin(1) + r_base(2);
t20 = pkin(2) + t22;
t13 = cos(qJ(1));
t19 = t13 * pkin(1) + t11 * qJ(2) + r_base(1);
t18 = t11 * pkin(6) + t21;
t17 = t13 * pkin(6) + t19;
t16 = t12 * pkin(3) + t10 * qJ(4) + t20;
t15 = rSges(5,1) * t10 - rSges(5,3) * t12;
t14 = g(1) * (pkin(3) * t25 + t17) + g(2) * (t13 * t24 + t18);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t13 - t11 * rSges(2,2) + r_base(1)) + g(2) * (t11 * rSges(2,1) + rSges(2,2) * t13 + r_base(2)) + g(3) * (rSges(2,3) + t22)) - m(3) * (g(1) * (-rSges(3,2) * t13 + t11 * rSges(3,3) + t19) + g(2) * (-t11 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t13 + t21) + g(3) * (rSges(3,1) + t22)) - m(4) * (g(1) * (rSges(4,1) * t25 + t11 * t27 + t17) + g(2) * (t11 * rSges(4,3) + t18) + g(3) * (rSges(4,1) * t12 - rSges(4,2) * t10 + t20) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t10 - qJ(2) - t27)) * t13) - m(5) * (g(3) * (rSges(5,1) * t12 + rSges(5,3) * t10 + t16) + (g(1) * (t15 - t24) + g(2) * rSges(5,2)) * t11 + (g(1) * rSges(5,2) + g(2) * (-pkin(3) * t10 - qJ(2) - t15)) * t13 + t14) - m(6) * (g(3) * (rSges(6,2) * t10 + t28 * t12 + t16) + (g(1) * (t28 * t10 - t24 - t26) + g(2) * t23) * t11 + (g(1) * t23 + (-qJ(2) + t26 + (-pkin(3) - t28) * t10) * g(2)) * t13 + t14);
U = t1;
