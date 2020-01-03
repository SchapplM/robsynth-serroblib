% Calculate potential energy for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:26
% EndTime: 2019-12-31 17:03:26
% DurationCPUTime: 0.13s
% Computational Cost: add. (88->51), mult. (60->47), div. (0->0), fcn. (40->6), ass. (0->16)
t20 = rSges(5,3) + pkin(6);
t19 = pkin(4) + r_base(3);
t10 = sin(qJ(1));
t18 = t10 * pkin(1) + r_base(2);
t12 = cos(qJ(1));
t17 = t12 * pkin(1) + r_base(1);
t16 = pkin(5) + t19;
t8 = qJ(1) + qJ(2);
t4 = sin(t8);
t15 = t4 * pkin(2) + t18;
t5 = cos(t8);
t14 = t5 * pkin(2) + t4 * qJ(3) + t17;
t11 = cos(qJ(4));
t9 = sin(qJ(4));
t13 = rSges(5,1) * t9 + rSges(5,2) * t11;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t12 - t10 * rSges(2,2) + r_base(1)) + g(2) * (t10 * rSges(2,1) + rSges(2,2) * t12 + r_base(2)) + g(3) * (rSges(2,3) + t19)) - m(3) * (g(1) * (rSges(3,1) * t5 - rSges(3,2) * t4 + t17) + g(2) * (rSges(3,1) * t4 + rSges(3,2) * t5 + t18) + g(3) * (rSges(3,3) + t16)) - m(4) * (g(1) * (-rSges(4,2) * t5 + rSges(4,3) * t4 + t14) + g(2) * (-rSges(4,2) * t4 + (-rSges(4,3) - qJ(3)) * t5 + t15) + g(3) * (rSges(4,1) + t16)) - m(5) * (g(1) * t14 + g(2) * t15 + g(3) * (rSges(5,1) * t11 - t9 * rSges(5,2) + pkin(3) + t16) + (g(1) * t13 + g(2) * t20) * t4 + (g(1) * t20 + g(2) * (-qJ(3) - t13)) * t5);
U = t1;
