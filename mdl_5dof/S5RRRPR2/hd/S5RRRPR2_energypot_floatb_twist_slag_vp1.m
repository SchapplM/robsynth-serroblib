% Calculate potential energy for
% S5RRRPR2
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:41
% EndTime: 2019-12-05 18:40:41
% DurationCPUTime: 0.17s
% Computational Cost: add. (152->62), mult. (74->56), div. (0->0), fcn. (50->10), ass. (0->26)
t28 = rSges(6,3) + pkin(8);
t12 = qJ(1) + qJ(2);
t27 = pkin(5) + r_base(1);
t16 = cos(qJ(1));
t26 = t16 * pkin(1) + r_base(3);
t25 = pkin(6) + t27;
t10 = qJ(3) + t12;
t9 = cos(t12);
t24 = pkin(2) * t9 + t26;
t23 = pkin(7) + t25;
t14 = sin(qJ(1));
t22 = -pkin(1) * t14 + r_base(2);
t7 = cos(t10);
t21 = pkin(3) * t7 + t24;
t20 = qJ(4) + t23;
t13 = sin(qJ(5));
t15 = cos(qJ(5));
t19 = rSges(6,1) * t15 - rSges(6,2) * t13 + pkin(4);
t8 = sin(t12);
t18 = -pkin(2) * t8 + t22;
t6 = sin(t10);
t17 = -pkin(3) * t6 + t18;
t5 = pkin(9) + t10;
t2 = cos(t5);
t1 = sin(t5);
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t27) + g(2) * (-t14 * rSges(2,1) - rSges(2,2) * t16 + r_base(2)) + g(3) * (rSges(2,1) * t16 - t14 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t25) + g(2) * (-rSges(3,1) * t8 - rSges(3,2) * t9 + t22) + g(3) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t26)) - m(4) * (g(1) * (rSges(4,3) + t23) + g(2) * (-rSges(4,1) * t6 - rSges(4,2) * t7 + t18) + g(3) * (rSges(4,1) * t7 - rSges(4,2) * t6 + t24)) - m(5) * (g(1) * (rSges(5,3) + t20) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2 + t17) + g(3) * (rSges(5,1) * t2 - rSges(5,2) * t1 + t21)) - m(6) * (g(1) * (rSges(6,1) * t13 + rSges(6,2) * t15 + t20) + g(2) * t17 + g(3) * t21 + (g(2) * t28 + g(3) * t19) * t2 + (-g(2) * t19 + g(3) * t28) * t1);
U = t3;
