% Calculate potential energy for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (114->67), mult. (103->63), div. (0->0), fcn. (79->6), ass. (0->23)
t28 = rSges(6,1) + pkin(4);
t14 = -pkin(7) - pkin(6);
t10 = sin(qJ(3));
t27 = pkin(3) * t10;
t26 = rSges(4,3) + pkin(6);
t25 = rSges(6,3) + qJ(5) - t14;
t24 = rSges(5,3) - t14;
t23 = pkin(5) + r_base(3);
t11 = sin(qJ(1));
t22 = t11 * pkin(1) + r_base(2);
t21 = pkin(2) + t23;
t13 = cos(qJ(1));
t20 = t13 * pkin(1) + t11 * qJ(2) + r_base(1);
t12 = cos(qJ(3));
t19 = t12 * pkin(3) + t21;
t18 = rSges(4,1) * t10 + rSges(4,2) * t12;
t9 = qJ(3) + qJ(4);
t2 = sin(t9);
t3 = cos(t9);
t17 = rSges(6,2) * t3 + t28 * t2 + t27;
t16 = rSges(5,1) * t2 + rSges(5,2) * t3 + t27;
t15 = g(1) * t20 + g(2) * t22;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t13 - rSges(2,2) * t11 + r_base(1)) + g(2) * (rSges(2,1) * t11 + rSges(2,2) * t13 + r_base(2)) + g(3) * (rSges(2,3) + t23)) - m(3) * (g(1) * (-rSges(3,2) * t13 + rSges(3,3) * t11 + t20) + g(2) * (-rSges(3,2) * t11 + (-rSges(3,3) - qJ(2)) * t13 + t22) + g(3) * (rSges(3,1) + t23)) - m(4) * (g(3) * (rSges(4,1) * t12 - rSges(4,2) * t10 + t21) + (g(1) * t18 + g(2) * t26) * t11 + (g(1) * t26 + g(2) * (-qJ(2) - t18)) * t13 + t15) - m(5) * (g(3) * (rSges(5,1) * t3 - rSges(5,2) * t2 + t19) + (g(1) * t16 + g(2) * t24) * t11 + (g(1) * t24 + g(2) * (-qJ(2) - t16)) * t13 + t15) - m(6) * (g(3) * (-rSges(6,2) * t2 + t28 * t3 + t19) + (g(1) * t17 + g(2) * t25) * t11 + (g(1) * t25 + g(2) * (-qJ(2) - t17)) * t13 + t15);
U = t1;
