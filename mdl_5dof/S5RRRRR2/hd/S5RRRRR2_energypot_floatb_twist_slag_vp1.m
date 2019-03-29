% Calculate potential energy for
% S5RRRRR2
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:40
% EndTime: 2019-03-29 15:25:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (118->64), mult. (99->73), div. (0->0), fcn. (79->10), ass. (0->26)
t16 = cos(qJ(3));
t30 = pkin(2) * t16;
t12 = sin(qJ(5));
t11 = qJ(1) + qJ(2);
t4 = sin(t11);
t29 = t12 * t4;
t6 = cos(t11);
t28 = t12 * t6;
t15 = cos(qJ(5));
t27 = t15 * t4;
t26 = t15 * t6;
t10 = qJ(3) + qJ(4);
t3 = sin(t10);
t25 = t3 * rSges(6,3);
t13 = sin(qJ(3));
t24 = t13 * pkin(2) + r_base(3);
t14 = sin(qJ(1));
t23 = t14 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t22 = t17 * pkin(1) + r_base(1);
t21 = t4 * t30 + t23;
t20 = t6 * t30 + t22;
t5 = cos(t10);
t19 = rSges(5,1) * t5 - rSges(5,2) * t3;
t18 = rSges(4,1) * t16 - rSges(4,2) * t13;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t17 - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + rSges(2,2) * t17 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t6 - rSges(3,2) * t4 + t22) + g(2) * (rSges(3,1) * t4 + rSges(3,2) * t6 + t23) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,3) * t4 + t18 * t6 + t22) + g(2) * (-rSges(4,3) * t6 + t18 * t4 + t23) + g(3) * (rSges(4,1) * t13 + rSges(4,2) * t16 + r_base(3))) - m(5) * (g(1) * (rSges(5,3) * t4 + t19 * t6 + t20) + g(2) * (-rSges(5,3) * t6 + t19 * t4 + t21) + g(3) * (rSges(5,1) * t3 + rSges(5,2) * t5 + t24)) - m(6) * (g(1) * ((t5 * t26 + t29) * rSges(6,1) + (-t5 * t28 + t27) * rSges(6,2) + t6 * t25 + t20) + g(2) * ((t5 * t27 - t28) * rSges(6,1) + (-t5 * t29 - t26) * rSges(6,2) + t4 * t25 + t21) + g(3) * (-rSges(6,3) * t5 + (rSges(6,1) * t15 - rSges(6,2) * t12) * t3 + t24));
U  = t1;
