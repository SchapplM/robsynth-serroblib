% Calculate potential energy for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:34
% EndTime: 2019-12-31 16:40:35
% DurationCPUTime: 0.24s
% Computational Cost: add. (85->64), mult. (115->72), div. (0->0), fcn. (101->6), ass. (0->20)
t12 = sin(qJ(1));
t9 = sin(pkin(6));
t27 = t12 * t9;
t26 = -pkin(5) - rSges(5,3);
t25 = qJ(3) * t9;
t10 = cos(pkin(6));
t24 = t12 * t10;
t23 = pkin(4) + r_base(3);
t22 = t12 * pkin(1) + r_base(2);
t21 = t9 * pkin(2) + t23;
t14 = cos(qJ(1));
t20 = t14 * pkin(1) + t12 * qJ(2) + r_base(1);
t19 = pkin(2) * t24 + t12 * t25 + t22;
t18 = t20 + (pkin(2) * t10 + t25) * t14;
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t17 = t10 * t13 + t11 * t9;
t16 = -t10 * t11 + t9 * t13;
t15 = t17 * rSges(5,1) + t16 * rSges(5,2) + t10 * pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t14 - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t14 + r_base(2)) + g(3) * (rSges(2,3) + t23)) - m(3) * (g(1) * (t12 * rSges(3,3) + t20) + g(2) * (rSges(3,1) * t24 - rSges(3,2) * t27 + t22) + g(3) * (rSges(3,1) * t9 + rSges(3,2) * t10 + t23) + (g(1) * (rSges(3,1) * t10 - rSges(3,2) * t9) + g(2) * (-rSges(3,3) - qJ(2))) * t14) - m(4) * (g(1) * (t12 * rSges(4,2) + t18) + g(2) * (rSges(4,1) * t24 + rSges(4,3) * t27 + t19) + g(3) * (rSges(4,1) * t9 + (-rSges(4,3) - qJ(3)) * t10 + t21) + (g(1) * (rSges(4,1) * t10 + rSges(4,3) * t9) + g(2) * (-rSges(4,2) - qJ(2))) * t14) - m(5) * (g(1) * t18 + g(2) * t19 + g(3) * (t16 * rSges(5,1) - t17 * rSges(5,2) + t9 * pkin(3) - t10 * qJ(3) + t21) + (g(1) * t26 + g(2) * t15) * t12 + (g(1) * t15 + g(2) * (-qJ(2) - t26)) * t14);
U = t1;
