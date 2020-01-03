% Calculate potential energy for
% S5RRRPR3
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:46
% EndTime: 2020-01-03 12:08:47
% DurationCPUTime: 0.22s
% Computational Cost: add. (153->67), mult. (97->64), div. (0->0), fcn. (73->10), ass. (0->28)
t32 = -rSges(4,3) - pkin(7);
t19 = cos(qJ(3));
t4 = t19 * pkin(3) + pkin(2);
t16 = -qJ(4) - pkin(7);
t31 = -rSges(5,3) + t16;
t30 = -rSges(6,3) - pkin(8) + t16;
t29 = pkin(5) + r_base(1);
t14 = qJ(3) + pkin(9);
t18 = sin(qJ(1));
t28 = t18 * pkin(1) + r_base(2);
t27 = pkin(6) + t29;
t20 = cos(qJ(1));
t26 = -pkin(1) * t20 + r_base(3);
t17 = sin(qJ(3));
t25 = t17 * pkin(3) + t27;
t5 = sin(t14);
t6 = cos(t14);
t24 = rSges(5,1) * t6 - rSges(5,2) * t5 + t4;
t7 = qJ(5) + t14;
t2 = sin(t7);
t3 = cos(t7);
t23 = rSges(6,1) * t3 - rSges(6,2) * t2 + pkin(4) * t6 + t4;
t22 = rSges(4,1) * t19 - rSges(4,2) * t17 + pkin(2);
t21 = g(2) * t28 + g(3) * t26;
t15 = qJ(1) + qJ(2);
t9 = cos(t15);
t8 = sin(t15);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t29) + g(2) * (t18 * rSges(2,1) + rSges(2,2) * t20 + r_base(2)) + g(3) * (-rSges(2,1) * t20 + t18 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t27) + g(2) * (rSges(3,1) * t8 + rSges(3,2) * t9 + t28) + g(3) * (-t9 * rSges(3,1) + t8 * rSges(3,2) + t26)) - m(4) * (g(1) * (rSges(4,1) * t17 + rSges(4,2) * t19 + t27) + (g(2) * t32 - g(3) * t22) * t9 + (g(2) * t22 + g(3) * t32) * t8 + t21) - m(5) * (g(1) * (rSges(5,1) * t5 + rSges(5,2) * t6 + t25) + (g(2) * t31 - g(3) * t24) * t9 + (g(2) * t24 + g(3) * t31) * t8 + t21) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t3 + pkin(4) * t5 + t25) + (g(2) * t30 - g(3) * t23) * t9 + (g(2) * t23 + g(3) * t30) * t8 + t21);
U = t1;
