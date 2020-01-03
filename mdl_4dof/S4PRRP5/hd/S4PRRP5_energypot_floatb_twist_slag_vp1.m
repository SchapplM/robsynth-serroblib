% Calculate potential energy for
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:48
% EndTime: 2019-12-31 16:28:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (92->67), mult. (125->78), div. (0->0), fcn. (113->6), ass. (0->26)
t35 = rSges(4,3) + pkin(5);
t34 = rSges(5,3) + qJ(4) + pkin(5);
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t33 = rSges(3,1) * t15 - rSges(3,2) * t13;
t10 = cos(pkin(6));
t9 = sin(pkin(6));
t32 = g(1) * t10 + g(2) * t9;
t29 = t15 * t9;
t12 = sin(qJ(3));
t28 = t9 * t12;
t24 = t10 * t12;
t23 = t10 * t15;
t22 = t12 * t15;
t14 = cos(qJ(3));
t21 = t14 * t15;
t19 = t9 * pkin(1) + r_base(2);
t18 = qJ(1) + r_base(3);
t17 = t10 * pkin(1) + t9 * pkin(4) + r_base(1);
t16 = -t10 * pkin(4) + t19;
t5 = pkin(3) * t14 + pkin(2);
t4 = t10 * t21 + t28;
t3 = -t10 * t22 + t9 * t14;
t2 = t9 * t21 - t24;
t1 = -t10 * t14 - t9 * t22;
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t10 - rSges(2,2) * t9 + r_base(1)) + g(2) * (rSges(2,1) * t9 + rSges(2,2) * t10 + r_base(2)) + g(3) * (rSges(2,3) + t18)) - m(3) * (g(1) * (t9 * rSges(3,3) + t17) + g(2) * (t33 * t9 + t19) + g(3) * (t13 * rSges(3,1) + rSges(3,2) * t15 + t18) + (g(1) * t33 + g(2) * (-rSges(3,3) - pkin(4))) * t10) - m(4) * (g(1) * (t4 * rSges(4,1) + t3 * rSges(4,2) + pkin(2) * t23 + t17) + g(2) * (t2 * rSges(4,1) + t1 * rSges(4,2) + pkin(2) * t29 + t16) + g(3) * (-t35 * t15 + t18) + (g(3) * (rSges(4,1) * t14 - rSges(4,2) * t12 + pkin(2)) + t32 * t35) * t13) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + pkin(3) * t28 + t5 * t23 + t17) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) - pkin(3) * t24 + t5 * t29 + t16) + g(3) * (-t34 * t15 + t18) + (g(3) * (rSges(5,1) * t14 - rSges(5,2) * t12 + t5) + t32 * t34) * t13);
U = t6;
