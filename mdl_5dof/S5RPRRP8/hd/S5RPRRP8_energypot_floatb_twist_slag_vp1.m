% Calculate potential energy for
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:04
% EndTime: 2019-12-31 18:47:04
% DurationCPUTime: 0.27s
% Computational Cost: add. (123->65), mult. (156->66), div. (0->0), fcn. (160->6), ass. (0->22)
t31 = -rSges(6,1) - pkin(4);
t30 = rSges(6,3) + qJ(5);
t29 = rSges(6,2) + pkin(7);
t28 = rSges(5,3) + pkin(7);
t27 = cos(qJ(3));
t26 = sin(qJ(1));
t25 = sin(qJ(3));
t24 = pkin(5) + r_base(3);
t23 = t26 * pkin(1) + r_base(2);
t22 = -pkin(6) + t24;
t15 = cos(qJ(1));
t21 = t15 * pkin(1) + t26 * qJ(2) + r_base(1);
t20 = t15 * pkin(2) + t21;
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t19 = -rSges(5,1) * t14 + rSges(5,2) * t13;
t18 = t26 * pkin(2) - t15 * qJ(2) + t23;
t17 = -t30 * t13 + t31 * t14;
t3 = -t15 * t27 - t26 * t25;
t4 = t15 * t25 - t26 * t27;
t16 = g(1) * (-t3 * pkin(3) + t20) + g(2) * (-t4 * pkin(3) + t18);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t15 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t15 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (t15 * rSges(3,1) + t26 * rSges(3,3) + t21) + g(2) * (t26 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t15 + t23) + g(3) * (rSges(3,2) + t24)) - m(4) * (g(1) * (-t3 * rSges(4,1) - t4 * rSges(4,2) + t20) + g(2) * (-t4 * rSges(4,1) + t3 * rSges(4,2) + t18) + g(3) * (-rSges(4,3) + t22)) - m(5) * (g(3) * (-t13 * rSges(5,1) - t14 * rSges(5,2) + t22) + (g(1) * t28 + g(2) * t19) * t4 + (g(1) * t19 - g(2) * t28) * t3 + t16) - m(6) * (g(3) * (t31 * t13 + t30 * t14 + t22) + (g(1) * t29 + g(2) * t17) * t4 + (g(1) * t17 - g(2) * t29) * t3 + t16);
U = t1;
