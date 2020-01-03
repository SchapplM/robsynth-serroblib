% Calculate potential energy for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:17
% EndTime: 2019-12-31 17:09:18
% DurationCPUTime: 0.34s
% Computational Cost: add. (102->73), mult. (125->85), div. (0->0), fcn. (113->8), ass. (0->23)
t30 = rSges(5,3) + pkin(6) + qJ(3);
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t29 = g(1) * t14 + g(2) * t12;
t28 = rSges(4,3) + qJ(3);
t8 = sin(pkin(7));
t25 = t12 * t8;
t24 = t14 * t8;
t11 = sin(qJ(2));
t23 = rSges(3,2) * t11;
t13 = cos(qJ(2));
t22 = t12 * t13;
t21 = t14 * t13;
t18 = pkin(4) + r_base(3);
t17 = t12 * pkin(1) + r_base(2);
t16 = t14 * pkin(1) + t12 * pkin(5) + r_base(1);
t15 = -t14 * pkin(5) + t17;
t9 = cos(pkin(7));
t7 = pkin(7) + qJ(4);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(3) * t9 + pkin(2);
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t14 - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t14 + r_base(2)) + g(3) * (rSges(2,3) + t18)) - m(3) * (g(1) * (t12 * rSges(3,3) + t16) + g(2) * (rSges(3,1) * t22 - t12 * t23 + t17) + g(3) * (rSges(3,1) * t11 + rSges(3,2) * t13 + t18) + (g(1) * (rSges(3,1) * t13 - t23) + g(2) * (-rSges(3,3) - pkin(5))) * t14) - m(4) * (g(1) * (pkin(2) * t21 + (t9 * t21 + t25) * rSges(4,1) + (t12 * t9 - t8 * t21) * rSges(4,2) + t16) + g(2) * (pkin(2) * t22 + (t9 * t22 - t24) * rSges(4,1) + (-t14 * t9 - t8 * t22) * rSges(4,2) + t15) + g(3) * (-t28 * t13 + t18) + (g(3) * (rSges(4,1) * t9 - rSges(4,2) * t8 + pkin(2)) + t29 * t28) * t11) - m(5) * (g(1) * (t1 * t21 + pkin(3) * t25 + (t12 * t2 + t3 * t21) * rSges(5,1) + (t12 * t3 - t2 * t21) * rSges(5,2) + t16) + g(2) * (t1 * t22 - pkin(3) * t24 + (-t14 * t2 + t3 * t22) * rSges(5,1) + (-t14 * t3 - t2 * t22) * rSges(5,2) + t15) + g(3) * (-t30 * t13 + t18) + (g(3) * (rSges(5,1) * t3 - rSges(5,2) * t2 + t1) + t29 * t30) * t11);
U = t4;
