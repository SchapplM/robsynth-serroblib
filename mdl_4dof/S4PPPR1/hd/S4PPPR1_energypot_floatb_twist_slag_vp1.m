% Calculate potential energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:52
% EndTime: 2019-03-08 18:08:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (66->52), mult. (62->47), div. (0->0), fcn. (46->4), ass. (0->13)
t8 = sin(pkin(5));
t17 = t8 * pkin(1) + r_base(2);
t16 = qJ(1) + r_base(3);
t15 = t8 * qJ(3) + t17;
t9 = cos(pkin(5));
t14 = t9 * pkin(1) + t8 * qJ(2) + r_base(1);
t13 = pkin(2) + t16;
t12 = t9 * qJ(3) + t14;
t11 = cos(qJ(4));
t10 = sin(qJ(4));
t2 = t9 * t10 + t11 * t8;
t1 = -t8 * t10 + t11 * t9;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t9 - rSges(2,2) * t8 + r_base(1)) + g(2) * (rSges(2,1) * t8 + rSges(2,2) * t9 + r_base(2)) + g(3) * (rSges(2,3) + t16)) - m(3) * (g(1) * (-rSges(3,2) * t9 + rSges(3,3) * t8 + t14) + g(2) * (-rSges(3,2) * t8 + (-rSges(3,3) - qJ(2)) * t9 + t17) + g(3) * (rSges(3,1) + t16)) - m(4) * (g(1) * (rSges(4,1) * t8 + rSges(4,3) * t9 + t12) + g(2) * (rSges(4,3) * t8 + (-rSges(4,1) - qJ(2)) * t9 + t15) + g(3) * (-rSges(4,2) + t13)) - m(5) * (g(1) * (rSges(5,1) * t2 + rSges(5,2) * t1 + pkin(3) * t8 + t12) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2 + (-pkin(3) - qJ(2)) * t9 + t15) + g(3) * (pkin(4) + rSges(5,3) + t13));
U  = t3;
