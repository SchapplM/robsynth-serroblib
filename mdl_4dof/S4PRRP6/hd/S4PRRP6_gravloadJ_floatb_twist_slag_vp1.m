% Calculate Gravitation load on the joints for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:13
% DurationCPUTime: 0.27s
% Computational Cost: add. (78->42), mult. (183->68), div. (0->0), fcn. (172->6), ass. (0->20)
t31 = -m(4) - m(5);
t11 = sin(pkin(6));
t12 = cos(pkin(6));
t29 = g(1) * t12 + g(2) * t11;
t30 = rSges(5,1) + pkin(3);
t21 = rSges(5,3) + qJ(4);
t14 = sin(qJ(2));
t25 = g(3) * t14;
t13 = sin(qJ(3));
t16 = cos(qJ(2));
t23 = t13 * t16;
t15 = cos(qJ(3));
t22 = t15 * t16;
t20 = rSges(4,1) * t15 - rSges(4,2) * t13;
t18 = t21 * t13 + t30 * t15;
t4 = t11 * t13 + t12 * t22;
t3 = -t11 * t15 + t12 * t23;
t2 = t11 * t22 - t12 * t13;
t1 = t11 * t23 + t12 * t15;
t5 = [(-m(2) - m(3) + t31) * g(3), (-m(3) * (g(3) * rSges(3,1) - t29 * rSges(3,2)) - m(4) * (t29 * rSges(4,3) + g(3) * t20) - m(5) * (t29 * rSges(5,2) + g(3) * t18)) * t16 + ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * rSges(5,2)) * g(3) + t29 * (m(3) * rSges(3,1) - m(4) * (-pkin(2) - t20) - m(5) * (-pkin(2) - t18))) * t14 + t31 * (g(3) * (t16 * pkin(2) + t14 * pkin(5)) + t29 * pkin(5) * t16), -m(4) * (g(1) * (-rSges(4,1) * t3 - rSges(4,2) * t4) + g(2) * (-rSges(4,1) * t1 - rSges(4,2) * t2)) - m(5) * (g(1) * (t21 * t4 - t3 * t30) + g(2) * (-t1 * t30 + t21 * t2)) + (-m(4) * (-rSges(4,1) * t13 - rSges(4,2) * t15) - m(5) * (-t30 * t13 + t21 * t15)) * t25, -m(5) * (g(1) * t3 + g(2) * t1 + t13 * t25)];
taug = t5(:);
