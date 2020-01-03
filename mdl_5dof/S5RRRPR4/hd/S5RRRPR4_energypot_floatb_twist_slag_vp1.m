% Calculate potential energy for
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:10:55
% EndTime: 2019-12-31 21:10:55
% DurationCPUTime: 0.33s
% Computational Cost: add. (156->73), mult. (131->77), div. (0->0), fcn. (113->8), ass. (0->26)
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t41 = pkin(3) * t18 + qJ(4) * t15;
t40 = rSges(5,1) * t18 + rSges(5,3) * t15;
t39 = rSges(4,1) * t18 - rSges(4,2) * t15;
t37 = -pkin(8) - rSges(6,3);
t31 = pkin(5) + r_base(3);
t16 = sin(qJ(1));
t30 = t16 * pkin(1) + r_base(2);
t19 = cos(qJ(1));
t29 = t19 * pkin(1) + r_base(1);
t28 = pkin(6) + t31;
t13 = qJ(1) + qJ(2);
t8 = sin(t13);
t27 = t8 * pkin(2) + t30;
t26 = t15 * pkin(3) + t28;
t9 = cos(t13);
t25 = t9 * pkin(2) + t8 * pkin(7) + t29;
t24 = t41 * t8 + t27;
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t23 = -t18 * t14 + t15 * t17;
t22 = t15 * t14 + t18 * t17;
t21 = t41 * t9 + t25;
t20 = t22 * rSges(6,1) + t23 * rSges(6,2) + t18 * pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t16 * rSges(2,2) + r_base(1)) + g(2) * (t16 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(1) * (t9 * rSges(3,1) - t8 * rSges(3,2) + t29) + g(2) * (t8 * rSges(3,1) + t9 * rSges(3,2) + t30) + g(3) * (rSges(3,3) + t28)) - m(4) * (g(1) * (t8 * rSges(4,3) + t25) + g(2) * (t39 * t8 + t27) + g(3) * (t15 * rSges(4,1) + t18 * rSges(4,2) + t28) + (g(1) * t39 + g(2) * (-rSges(4,3) - pkin(7))) * t9) - m(5) * (g(1) * (t8 * rSges(5,2) + t21) + g(2) * (t40 * t8 + t24) + g(3) * (t15 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t18 + t26) + (g(1) * t40 + g(2) * (-rSges(5,2) - pkin(7))) * t9) - m(6) * (g(1) * t21 + g(2) * t24 + g(3) * (t23 * rSges(6,1) - t22 * rSges(6,2) + t15 * pkin(4) - t18 * qJ(4) + t26) + (g(1) * t37 + g(2) * t20) * t8 + (g(1) * t20 + g(2) * (-pkin(7) - t37)) * t9);
U = t1;
