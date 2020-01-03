% Calculate potential energy for
% S4RRRP5
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:35
% EndTime: 2019-12-31 17:16:35
% DurationCPUTime: 0.24s
% Computational Cost: add. (96->57), mult. (88->58), div. (0->0), fcn. (68->6), ass. (0->20)
t25 = rSges(5,1) + pkin(3);
t8 = qJ(2) + qJ(3);
t5 = sin(t8);
t6 = cos(t8);
t24 = rSges(4,1) * t6 - rSges(4,2) * t5;
t23 = rSges(5,3) + qJ(4);
t20 = rSges(3,3) + pkin(5);
t19 = pkin(4) + r_base(3);
t12 = cos(qJ(1));
t11 = cos(qJ(2));
t3 = pkin(2) * t11 + pkin(1);
t18 = t12 * t3 + r_base(1);
t9 = sin(qJ(2));
t17 = t9 * pkin(2) + t19;
t10 = sin(qJ(1));
t13 = -pkin(6) - pkin(5);
t16 = t10 * t3 + t12 * t13 + r_base(2);
t15 = rSges(3,1) * t11 - rSges(3,2) * t9 + pkin(1);
t14 = t23 * t5 + t25 * t6;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t12 - rSges(2,2) * t10 + r_base(1)) + g(2) * (rSges(2,1) * t10 + rSges(2,2) * t12 + r_base(2)) + g(3) * (rSges(2,3) + t19)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t9 + rSges(3,2) * t11 + t19) + (g(1) * t15 - g(2) * t20) * t12 + (g(1) * t20 + g(2) * t15) * t10) - m(4) * (g(1) * (t12 * t24 + t18) + g(2) * (-rSges(4,3) * t12 + t16) + g(3) * (rSges(4,1) * t5 + rSges(4,2) * t6 + t17) + (g(1) * (rSges(4,3) - t13) + g(2) * t24) * t10) - m(5) * (g(1) * t18 + g(2) * t16 + g(3) * (-t23 * t6 + t25 * t5 + t17) + (-g(2) * rSges(5,2) + g(1) * t14) * t12 + (g(1) * (rSges(5,2) - t13) + g(2) * t14) * t10);
U = t1;
