% Calculate potential energy for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:41
% EndTime: 2019-12-31 17:50:41
% DurationCPUTime: 0.24s
% Computational Cost: add. (129->63), mult. (89->59), div. (0->0), fcn. (65->6), ass. (0->21)
t26 = rSges(6,1) + pkin(4);
t25 = rSges(5,3) + pkin(6);
t24 = rSges(6,3) + qJ(5) + pkin(6);
t23 = pkin(5) + r_base(3);
t11 = sin(qJ(1));
t22 = t11 * pkin(1) + r_base(2);
t13 = cos(qJ(1));
t21 = t13 * pkin(1) + r_base(1);
t8 = qJ(1) + pkin(7);
t4 = sin(t8);
t20 = t4 * pkin(2) + t22;
t19 = qJ(2) + t23;
t5 = cos(t8);
t18 = t5 * pkin(2) + t4 * qJ(3) + t21;
t17 = pkin(3) + t19;
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t16 = rSges(5,1) * t10 + rSges(5,2) * t12;
t15 = rSges(6,2) * t12 + t26 * t10;
t14 = g(1) * t18 + g(2) * t20;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t13 - t11 * rSges(2,2) + r_base(1)) + g(2) * (t11 * rSges(2,1) + rSges(2,2) * t13 + r_base(2)) + g(3) * (rSges(2,3) + t23)) - m(3) * (g(1) * (rSges(3,1) * t5 - rSges(3,2) * t4 + t21) + g(2) * (rSges(3,1) * t4 + rSges(3,2) * t5 + t22) + g(3) * (rSges(3,3) + t19)) - m(4) * (g(1) * (-rSges(4,2) * t5 + rSges(4,3) * t4 + t18) + g(2) * (-rSges(4,2) * t4 + (-rSges(4,3) - qJ(3)) * t5 + t20) + g(3) * (rSges(4,1) + t19)) - m(5) * (g(3) * (rSges(5,1) * t12 - rSges(5,2) * t10 + t17) + (g(1) * t16 + g(2) * t25) * t4 + (g(1) * t25 + g(2) * (-qJ(3) - t16)) * t5 + t14) - m(6) * (g(3) * (-rSges(6,2) * t10 + t26 * t12 + t17) + (g(1) * t15 + g(2) * t24) * t4 + (g(1) * t24 + g(2) * (-qJ(3) - t15)) * t5 + t14);
U = t1;
