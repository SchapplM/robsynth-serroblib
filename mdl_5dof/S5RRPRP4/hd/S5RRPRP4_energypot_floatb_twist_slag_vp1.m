% Calculate potential energy for
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:34
% EndTime: 2019-12-31 19:52:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (132->65), mult. (94->62), div. (0->0), fcn. (70->6), ass. (0->21)
t29 = rSges(6,1) + pkin(4);
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t28 = rSges(5,1) * t11 + rSges(5,2) * t13;
t27 = rSges(6,3) + qJ(5);
t24 = pkin(5) + r_base(3);
t12 = sin(qJ(1));
t23 = t12 * pkin(1) + r_base(2);
t14 = cos(qJ(1));
t22 = t14 * pkin(1) + r_base(1);
t21 = pkin(6) + t24;
t10 = qJ(1) + qJ(2);
t6 = sin(t10);
t20 = t6 * pkin(2) + t23;
t19 = pkin(3) + t21;
t7 = cos(t10);
t18 = t7 * pkin(2) + t6 * qJ(3) + t22;
t17 = t6 * pkin(7) + t20;
t16 = t7 * pkin(7) + t18;
t15 = t29 * t11 - t27 * t13;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t14 - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t14 + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (t7 * rSges(3,1) - t6 * rSges(3,2) + t22) + g(2) * (rSges(3,1) * t6 + rSges(3,2) * t7 + t23) + g(3) * (rSges(3,3) + t21)) - m(4) * (g(1) * (-rSges(4,2) * t7 + rSges(4,3) * t6 + t18) + g(2) * (-rSges(4,2) * t6 + (-rSges(4,3) - qJ(3)) * t7 + t20) + g(3) * (rSges(4,1) + t21)) - m(5) * (g(1) * (t28 * t6 + t16) + g(2) * (rSges(5,3) * t6 + t17) + g(3) * (rSges(5,1) * t13 - rSges(5,2) * t11 + t19) + (g(1) * rSges(5,3) + g(2) * (-qJ(3) - t28)) * t7) - m(6) * (g(1) * t16 + g(2) * t17 + g(3) * (t27 * t11 + t29 * t13 + t19) + (g(2) * rSges(6,2) + g(1) * t15) * t6 + (g(1) * rSges(6,2) + g(2) * (-qJ(3) - t15)) * t7);
U = t1;
