% Calculate potential energy for
% S5RPPRR12
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:53
% EndTime: 2019-12-31 18:06:53
% DurationCPUTime: 0.34s
% Computational Cost: add. (126->71), mult. (123->69), div. (0->0), fcn. (103->8), ass. (0->28)
t12 = sin(qJ(5));
t14 = cos(qJ(5));
t18 = rSges(6,1) * t14 - rSges(6,2) * t12 + pkin(4);
t8 = pkin(8) + qJ(4);
t2 = sin(t8);
t3 = cos(t8);
t31 = rSges(6,3) + pkin(7);
t33 = t18 * t2 - t31 * t3;
t9 = sin(pkin(8));
t30 = t9 * pkin(3);
t11 = -pkin(6) - qJ(3);
t29 = rSges(5,3) - t11;
t28 = rSges(4,3) + qJ(3);
t27 = pkin(5) + r_base(3);
t13 = sin(qJ(1));
t26 = t13 * pkin(1) + r_base(2);
t25 = pkin(2) + t27;
t15 = cos(qJ(1));
t24 = t15 * pkin(1) + t13 * qJ(2) + r_base(1);
t23 = -qJ(2) - t30;
t22 = g(2) * t26;
t10 = cos(pkin(8));
t21 = t10 * pkin(3) + t25;
t20 = rSges(5,1) * t2 + rSges(5,2) * t3;
t19 = rSges(4,1) * t9 + rSges(4,2) * t10;
t17 = t12 * rSges(6,1) + t14 * rSges(6,2) - t11;
t16 = g(1) * (t13 * t30 + t24) + t22;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t15 - t13 * rSges(2,2) + r_base(1)) + g(2) * (t13 * rSges(2,1) + rSges(2,2) * t15 + r_base(2)) + g(3) * (rSges(2,3) + t27)) - m(3) * (g(1) * (-rSges(3,2) * t15 + t13 * rSges(3,3) + t24) + g(2) * (-t13 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t15 + t26) + g(3) * (rSges(3,1) + t27)) - m(4) * (g(1) * t24 + t22 + g(3) * (rSges(4,1) * t10 - rSges(4,2) * t9 + t25) + (g(1) * t19 + g(2) * t28) * t13 + (g(1) * t28 + g(2) * (-qJ(2) - t19)) * t15) - m(5) * (g(3) * (rSges(5,1) * t3 - rSges(5,2) * t2 + t21) + (g(1) * t20 + g(2) * t29) * t13 + (g(1) * t29 + g(2) * (-t20 + t23)) * t15 + t16) - m(6) * ((g(1) * t33 + g(2) * t17) * t13 + (g(1) * t17 + (t23 - t33) * g(2)) * t15 + t16 + (t18 * t3 + t31 * t2 + t21) * g(3));
U = t1;
