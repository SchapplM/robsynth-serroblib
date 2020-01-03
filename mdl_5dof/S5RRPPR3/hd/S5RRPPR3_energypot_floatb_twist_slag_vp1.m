% Calculate potential energy for
% S5RRPPR3
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:25
% EndTime: 2019-12-31 19:26:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (143->63), mult. (76->57), div. (0->0), fcn. (52->8), ass. (0->22)
t28 = rSges(6,3) + pkin(7);
t13 = qJ(1) + qJ(2);
t27 = pkin(5) + r_base(3);
t15 = sin(qJ(1));
t26 = t15 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + r_base(1);
t24 = pkin(6) + t27;
t9 = sin(t13);
t23 = pkin(2) * t9 + t26;
t10 = cos(t13);
t22 = pkin(2) * t10 + t25;
t8 = pkin(8) + t13;
t4 = sin(t8);
t21 = t4 * pkin(3) + t23;
t20 = qJ(3) + t24;
t5 = cos(t8);
t19 = t5 * pkin(3) + t4 * qJ(4) + t22;
t14 = sin(qJ(5));
t16 = cos(qJ(5));
t18 = rSges(6,1) * t14 + rSges(6,2) * t16;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t17 - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + rSges(2,2) * t17 + r_base(2)) + g(3) * (rSges(2,3) + t27)) - m(3) * (g(1) * (rSges(3,1) * t10 - t9 * rSges(3,2) + t25) + g(2) * (rSges(3,1) * t9 + rSges(3,2) * t10 + t26) + g(3) * (rSges(3,3) + t24)) - m(4) * (g(1) * (rSges(4,1) * t5 - rSges(4,2) * t4 + t22) + g(2) * (rSges(4,1) * t4 + rSges(4,2) * t5 + t23) + g(3) * (rSges(4,3) + t20)) - m(5) * (g(1) * (-rSges(5,2) * t5 + rSges(5,3) * t4 + t19) + g(2) * (-rSges(5,2) * t4 + (-rSges(5,3) - qJ(4)) * t5 + t21) + g(3) * (rSges(5,1) + t20)) - m(6) * (g(1) * t19 + g(2) * t21 + g(3) * (rSges(6,1) * t16 - t14 * rSges(6,2) + pkin(4) + t20) + (g(1) * t18 + g(2) * t28) * t4 + (g(1) * t28 + g(2) * (-qJ(4) - t18)) * t5);
U = t1;
