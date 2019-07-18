% Calculate potential energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:17
% EndTime: 2019-07-18 13:28:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (87->61), mult. (87->67), div. (0->0), fcn. (67->8), ass. (0->23)
t6 = sin(qJ(5));
t8 = sin(qJ(2));
t24 = t8 * t6;
t9 = cos(qJ(5));
t23 = t8 * t9;
t10 = cos(qJ(3));
t22 = pkin(2) * t10;
t11 = cos(qJ(2));
t21 = t11 * t6;
t20 = t11 * t9;
t5 = qJ(3) + qJ(4);
t3 = sin(t5);
t19 = t3 * rSges(6,3);
t18 = pkin(1) + r_base(1);
t17 = qJ(1) + r_base(3);
t16 = t11 * t22 + t18;
t15 = t8 * t22 + t17;
t7 = sin(qJ(3));
t14 = -pkin(2) * t7 + r_base(2);
t4 = cos(t5);
t13 = rSges(5,1) * t4 - rSges(5,2) * t3;
t12 = rSges(4,1) * t10 - rSges(4,2) * t7;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (r_base(1) + rSges(2,1)) + g(2) * (r_base(2) + rSges(2,2)) + g(3) * (rSges(2,3) + t17)) - m(3) * (g(1) * (rSges(3,1) * t11 - t8 * rSges(3,2) + t18) + g(2) * (r_base(2) - rSges(3,3)) + g(3) * (t8 * rSges(3,1) + rSges(3,2) * t11 + t17)) - m(4) * (g(1) * (t8 * rSges(4,3) + t12 * t11 + t18) + g(2) * (-rSges(4,1) * t7 - rSges(4,2) * t10 + r_base(2)) + g(3) * (-rSges(4,3) * t11 + t12 * t8 + t17)) - m(5) * (g(1) * (t8 * rSges(5,3) + t13 * t11 + t16) + g(2) * (-rSges(5,1) * t3 - rSges(5,2) * t4 + t14) + g(3) * (-rSges(5,3) * t11 + t13 * t8 + t15)) - m(6) * (g(1) * ((t4 * t20 + t24) * rSges(6,1) + (-t4 * t21 + t23) * rSges(6,2) + t11 * t19 + t16) + g(2) * (rSges(6,3) * t4 + (-rSges(6,1) * t9 + rSges(6,2) * t6) * t3 + t14) + g(3) * ((t4 * t23 - t21) * rSges(6,1) + (-t4 * t24 - t20) * rSges(6,2) + t8 * t19 + t15));
U  = t1;
