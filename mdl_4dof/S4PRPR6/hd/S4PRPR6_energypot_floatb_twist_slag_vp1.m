% Calculate potential energy for
% S4PRPR6
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:16
% EndTime: 2019-12-31 16:24:17
% DurationCPUTime: 0.34s
% Computational Cost: add. (102->72), mult. (125->84), div. (0->0), fcn. (113->8), ass. (0->23)
t32 = rSges(5,3) + pkin(5) + qJ(3);
t13 = sin(qJ(2));
t14 = cos(qJ(2));
t31 = rSges(3,1) * t14 - rSges(3,2) * t13;
t11 = cos(pkin(6));
t9 = sin(pkin(6));
t30 = g(1) * t11 + g(2) * t9;
t29 = rSges(4,3) + qJ(3);
t8 = sin(pkin(7));
t27 = t9 * t8;
t25 = t11 * t8;
t24 = t9 * t14;
t21 = t11 * t14;
t18 = t9 * pkin(1) + r_base(2);
t17 = qJ(1) + r_base(3);
t16 = t11 * pkin(1) + t9 * pkin(4) + r_base(1);
t15 = -t11 * pkin(4) + t18;
t10 = cos(pkin(7));
t7 = pkin(7) + qJ(4);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(3) * t10 + pkin(2);
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t11 - rSges(2,2) * t9 + r_base(1)) + g(2) * (rSges(2,1) * t9 + rSges(2,2) * t11 + r_base(2)) + g(3) * (rSges(2,3) + t17)) - m(3) * (g(1) * (t9 * rSges(3,3) + t16) + g(2) * (t31 * t9 + t18) + g(3) * (t13 * rSges(3,1) + rSges(3,2) * t14 + t17) + (g(1) * t31 + g(2) * (-rSges(3,3) - pkin(4))) * t11) - m(4) * (g(1) * (pkin(2) * t21 + (t10 * t21 + t27) * rSges(4,1) + (t9 * t10 - t8 * t21) * rSges(4,2) + t16) + g(2) * (pkin(2) * t24 + (t10 * t24 - t25) * rSges(4,1) + (-t11 * t10 - t8 * t24) * rSges(4,2) + t15) + g(3) * (-t29 * t14 + t17) + (g(3) * (rSges(4,1) * t10 - rSges(4,2) * t8 + pkin(2)) + t30 * t29) * t13) - m(5) * (g(1) * (t1 * t21 + pkin(3) * t27 + (t9 * t2 + t3 * t21) * rSges(5,1) + (-t2 * t21 + t9 * t3) * rSges(5,2) + t16) + g(2) * (t1 * t24 - pkin(3) * t25 + (-t11 * t2 + t3 * t24) * rSges(5,1) + (-t11 * t3 - t2 * t24) * rSges(5,2) + t15) + g(3) * (-t32 * t14 + t17) + (g(3) * (rSges(5,1) * t3 - rSges(5,2) * t2 + t1) + t30 * t32) * t13);
U = t4;
