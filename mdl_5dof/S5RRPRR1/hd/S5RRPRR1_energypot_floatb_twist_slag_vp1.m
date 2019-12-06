% Calculate potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:46
% EndTime: 2019-12-05 18:23:46
% DurationCPUTime: 0.32s
% Computational Cost: add. (108->71), mult. (115->81), div. (0->0), fcn. (95->8), ass. (0->26)
t32 = rSges(4,1) + pkin(1);
t31 = rSges(6,3) + pkin(4);
t7 = qJ(2) + qJ(4);
t5 = sin(t7);
t6 = cos(t7);
t30 = rSges(5,1) * t6 - rSges(5,2) * t5;
t11 = sin(qJ(1));
t9 = sin(qJ(5));
t27 = t11 * t9;
t14 = cos(qJ(1));
t26 = t14 * t9;
t12 = cos(qJ(5));
t24 = t11 * t12;
t13 = cos(qJ(2));
t15 = pkin(2) + pkin(1);
t23 = t13 * t15;
t22 = t14 * t12;
t21 = rSges(4,3) + qJ(3);
t20 = t14 * t23 + r_base(1);
t10 = sin(qJ(2));
t19 = t10 * t15 + r_base(3);
t8 = -pkin(3) - qJ(3);
t18 = t11 * t23 + t14 * t8 + r_base(2);
t17 = rSges(3,1) * t13 - rSges(3,2) * t10;
t16 = -rSges(4,2) * t10 + t32 * t13;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t14 * rSges(2,1) - t11 * rSges(2,2) + r_base(1)) + g(2) * (t11 * rSges(2,1) + t14 * rSges(2,2) + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (t11 * rSges(3,3) + t17 * t14 + r_base(1)) + g(2) * (-t14 * rSges(3,3) + t17 * t11 + r_base(2)) + g(3) * (t10 * rSges(3,1) + t13 * rSges(3,2) + r_base(3))) - m(4) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t13 * rSges(4,2) + t32 * t10 + r_base(3)) + (g(1) * t16 - g(2) * t21) * t14 + (g(1) * t21 + g(2) * t16) * t11) - m(5) * (g(1) * (t30 * t14 + t20) + g(2) * (-t14 * rSges(5,3) + t18) + g(3) * (t5 * rSges(5,1) + t6 * rSges(5,2) + t19) + (g(1) * (rSges(5,3) - t8) + g(2) * t30) * t11) - m(6) * (g(1) * (-t11 * t8 + (t6 * t22 + t27) * rSges(6,1) + (-t6 * t26 + t24) * rSges(6,2) + t20) + g(2) * ((t6 * t24 - t26) * rSges(6,1) + (-t6 * t27 - t22) * rSges(6,2) + t18) + g(3) * (-t31 * t6 + t19) + (g(3) * (rSges(6,1) * t12 - rSges(6,2) * t9) + (g(1) * t14 + g(2) * t11) * t31) * t5);
U = t1;
