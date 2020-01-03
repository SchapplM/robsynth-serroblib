% Calculate Gravitation load on the joints for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:11
% DurationCPUTime: 0.35s
% Computational Cost: add. (123->69), mult. (253->102), div. (0->0), fcn. (236->6), ass. (0->26)
t35 = rSges(5,1) + pkin(3);
t34 = rSges(5,3) + qJ(4) + pkin(6);
t12 = sin(qJ(1));
t15 = cos(qJ(1));
t33 = g(1) * t15 + g(2) * t12;
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t24 = rSges(4,3) + pkin(6);
t32 = pkin(2) * t14 + t24 * t11;
t10 = sin(qJ(3));
t13 = cos(qJ(3));
t5 = pkin(3) * t13 + pkin(2);
t31 = m(3) * rSges(3,1) + m(4) * (rSges(4,1) * t13 - rSges(4,2) * t10 + pkin(2)) + m(5) * (rSges(5,1) * t13 - rSges(5,2) * t10 + t5);
t30 = t15 * pkin(1) + t12 * pkin(5);
t27 = pkin(3) * t10;
t23 = rSges(3,2) * t11;
t22 = t12 * t14;
t21 = t14 * t15;
t3 = -t10 * t21 + t12 * t13;
t1 = t10 * t22 + t13 * t15;
t18 = t11 * t34 + t14 * t5;
t17 = m(3) * rSges(3,2) - m(4) * t24 - m(5) * t34;
t7 = t15 * pkin(5);
t4 = t12 * t10 + t13 * t21;
t2 = t10 * t15 - t13 * t22;
t6 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t15 + t7) + g(2) * (rSges(3,1) * t21 - t15 * t23 + t30) + (g(1) * (-rSges(3,1) * t14 - pkin(1) + t23) + g(2) * rSges(3,3)) * t12) - m(4) * ((t4 * rSges(4,1) + t3 * rSges(4,2) + t15 * t32 + t30) * g(2) + (rSges(4,1) * t2 + rSges(4,2) * t1 + t7 + (-pkin(1) - t32) * t12) * g(1)) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t7) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t30) + (g(1) * t27 + g(2) * t18) * t15 + (g(1) * (-pkin(1) - t18) + g(2) * t27) * t12), (t17 * t11 - t14 * t31) * g(3) + t33 * (t11 * t31 + t17 * t14), -m(4) * (g(1) * (rSges(4,1) * t3 - rSges(4,2) * t4) + g(2) * (-rSges(4,1) * t1 + rSges(4,2) * t2)) - m(5) * (g(1) * (-t4 * rSges(5,2) + t3 * t35) + g(2) * (t2 * rSges(5,2) - t1 * t35)) + (-m(4) * (-rSges(4,1) * t10 - rSges(4,2) * t13) - m(5) * (-rSges(5,1) * t10 - rSges(5,2) * t13 - t27)) * g(3) * t11, -m(5) * (-g(3) * t14 + t11 * t33)];
taug = t6(:);
