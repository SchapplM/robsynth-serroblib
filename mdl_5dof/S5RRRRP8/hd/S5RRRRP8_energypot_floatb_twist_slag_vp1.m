% Calculate potential energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:34
% EndTime: 2019-12-31 21:59:34
% DurationCPUTime: 0.44s
% Computational Cost: add. (160->93), mult. (180->106), div. (0->0), fcn. (168->8), ass. (0->31)
t40 = rSges(4,3) + pkin(7);
t22 = -pkin(8) - pkin(7);
t39 = rSges(5,3) - t22;
t38 = rSges(6,3) + qJ(5) - t22;
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t37 = g(1) * t21 + g(2) * t18;
t19 = cos(qJ(3));
t7 = t19 * pkin(3) + pkin(2);
t17 = sin(qJ(2));
t33 = rSges(3,2) * t17;
t16 = sin(qJ(3));
t32 = t18 * t16;
t20 = cos(qJ(2));
t31 = t18 * t20;
t30 = t21 * t16;
t29 = t21 * t20;
t26 = pkin(5) + r_base(3);
t25 = t18 * pkin(1) + r_base(2);
t24 = t21 * pkin(1) + t18 * pkin(6) + r_base(1);
t23 = -t21 * pkin(6) + t25;
t15 = qJ(3) + qJ(4);
t9 = cos(t15);
t8 = sin(t15);
t6 = pkin(3) * t16 + pkin(4) * t8;
t5 = pkin(4) * t9 + t7;
t4 = t18 * t8 + t29 * t9;
t3 = t18 * t9 - t29 * t8;
t2 = -t21 * t8 + t31 * t9;
t1 = -t21 * t9 - t31 * t8;
t10 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t21 - rSges(2,2) * t18 + r_base(1)) + g(2) * (rSges(2,1) * t18 + rSges(2,2) * t21 + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (t18 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t31 - t18 * t33 + t25) + g(3) * (rSges(3,1) * t17 + rSges(3,2) * t20 + t26) + (g(1) * (rSges(3,1) * t20 - t33) + g(2) * (-rSges(3,3) - pkin(6))) * t21) - m(4) * (g(1) * (pkin(2) * t29 + (t19 * t29 + t32) * rSges(4,1) + (-t16 * t29 + t18 * t19) * rSges(4,2) + t24) + g(2) * (pkin(2) * t31 + (t19 * t31 - t30) * rSges(4,1) + (-t16 * t31 - t19 * t21) * rSges(4,2) + t23) + g(3) * (-t40 * t20 + t26) + (g(3) * (rSges(4,1) * t19 - rSges(4,2) * t16 + pkin(2)) + t37 * t40) * t17) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + pkin(3) * t32 + t29 * t7 + t24) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) - pkin(3) * t30 + t31 * t7 + t23) + g(3) * (-t39 * t20 + t26) + (g(3) * (rSges(5,1) * t9 - rSges(5,2) * t8 + t7) + t37 * t39) * t17) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t18 * t6 + t29 * t5 + t24) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t21 * t6 + t31 * t5 + t23) + g(3) * (-t38 * t20 + t26) + (g(3) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t5) + t37 * t38) * t17);
U = t10;
