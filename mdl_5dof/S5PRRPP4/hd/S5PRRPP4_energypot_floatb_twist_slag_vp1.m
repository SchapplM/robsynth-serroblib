% Calculate potential energy for
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:43
% DurationCPUTime: 0.28s
% Computational Cost: add. (146->70), mult. (115->70), div. (0->0), fcn. (91->6), ass. (0->23)
t16 = sin(qJ(3));
t17 = cos(qJ(3));
t38 = pkin(3) * t17 + qJ(4) * t16;
t37 = rSges(6,1) + pkin(4);
t36 = rSges(5,1) * t17 + rSges(5,3) * t16;
t35 = rSges(4,1) * t17 - rSges(4,2) * t16;
t28 = -rSges(6,3) - qJ(5);
t14 = sin(pkin(7));
t27 = t14 * pkin(1) + r_base(2);
t15 = cos(pkin(7));
t26 = t15 * pkin(1) + r_base(1);
t25 = qJ(1) + r_base(3);
t13 = pkin(7) + qJ(2);
t8 = sin(t13);
t24 = t8 * pkin(2) + t27;
t23 = pkin(5) + t25;
t9 = cos(t13);
t22 = t9 * pkin(2) + t8 * pkin(6) + t26;
t21 = t16 * pkin(3) + t23;
t20 = t38 * t8 + t24;
t19 = t38 * t9 + t22;
t18 = rSges(6,2) * t16 + t37 * t17;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t15 - rSges(2,2) * t14 + r_base(1)) + g(2) * (rSges(2,1) * t14 + rSges(2,2) * t15 + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t26) + g(2) * (rSges(3,1) * t8 + rSges(3,2) * t9 + t27) + g(3) * (rSges(3,3) + t23)) - m(4) * (g(1) * (t8 * rSges(4,3) + t22) + g(2) * (t35 * t8 + t24) + g(3) * (t16 * rSges(4,1) + rSges(4,2) * t17 + t23) + (g(1) * t35 + g(2) * (-rSges(4,3) - pkin(6))) * t9) - m(5) * (g(1) * (t8 * rSges(5,2) + t19) + g(2) * (t36 * t8 + t20) + g(3) * (t16 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t17 + t21) + (g(1) * t36 + g(2) * (-rSges(5,2) - pkin(6))) * t9) - m(6) * (g(1) * t19 + g(2) * t20 + g(3) * ((-rSges(6,2) - qJ(4)) * t17 + t37 * t16 + t21) + (g(1) * t28 + g(2) * t18) * t8 + (g(1) * t18 + g(2) * (-pkin(6) - t28)) * t9);
U = t1;
