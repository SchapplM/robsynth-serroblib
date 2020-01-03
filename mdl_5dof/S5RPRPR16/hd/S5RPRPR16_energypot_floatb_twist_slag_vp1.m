% Calculate potential energy for
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR16_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:39
% EndTime: 2019-12-31 18:38:40
% DurationCPUTime: 0.34s
% Computational Cost: add. (106->80), mult. (130->91), div. (0->0), fcn. (110->6), ass. (0->22)
t30 = rSges(6,3) + pkin(7);
t11 = sin(qJ(3));
t12 = sin(qJ(1));
t29 = t12 * t11;
t14 = cos(qJ(3));
t28 = t12 * t14;
t15 = cos(qJ(1));
t27 = t14 * t15;
t26 = t14 * qJ(4);
t25 = pkin(5) + r_base(3);
t24 = t12 * pkin(1) + r_base(2);
t23 = pkin(2) + t25;
t22 = t15 * pkin(1) + t12 * qJ(2) + r_base(1);
t21 = t12 * pkin(6) + t24;
t20 = t15 * t26 + t21;
t19 = t15 * pkin(6) + t22;
t18 = t14 * pkin(3) + t11 * qJ(4) + t23;
t17 = pkin(3) * t29 + t19;
t16 = -rSges(5,2) * t11 - rSges(5,3) * t14;
t13 = cos(qJ(5));
t10 = sin(qJ(5));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t15 - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t15 + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (-rSges(3,2) * t15 + t12 * rSges(3,3) + t22) + g(2) * (-t12 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t15 + t24) + g(3) * (rSges(3,1) + t25)) - m(4) * (g(1) * (rSges(4,1) * t29 + rSges(4,2) * t28 + t19) + g(2) * (t12 * rSges(4,3) + t21) + g(3) * (rSges(4,1) * t14 - rSges(4,2) * t11 + t23) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t11 - rSges(4,2) * t14 - qJ(2))) * t15) - m(5) * (g(1) * t17 + g(2) * t20 + g(3) * (-rSges(5,2) * t14 + rSges(5,3) * t11 + t18) + (g(1) * (t16 - t26) + g(2) * rSges(5,1)) * t12 + (g(1) * rSges(5,1) + g(2) * (-t11 * pkin(3) - qJ(2) - t16)) * t15) - m(6) * (g(1) * (t15 * pkin(4) - t12 * t26 + (-t10 * t28 + t13 * t15) * rSges(6,1) + (-t10 * t15 - t13 * t28) * rSges(6,2) + t17) + g(2) * (t12 * pkin(4) - t15 * qJ(2) + (t10 * t27 + t12 * t13) * rSges(6,1) + (-t12 * t10 + t13 * t27) * rSges(6,2) + t20) + g(3) * (t30 * t14 + t18) + (g(3) * (rSges(6,1) * t10 + rSges(6,2) * t13) + g(1) * t30 * t12 + g(2) * (-pkin(3) - t30) * t15) * t11);
U = t1;
