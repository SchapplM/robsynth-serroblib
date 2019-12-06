% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:24
% DurationCPUTime: 0.27s
% Computational Cost: add. (225->59), mult. (154->75), div. (0->0), fcn. (117->8), ass. (0->31)
t46 = rSges(6,1) + pkin(4);
t45 = rSges(5,3) + pkin(7);
t44 = rSges(6,3) + qJ(5) + pkin(7);
t21 = cos(qJ(4));
t38 = rSges(5,1) * t21;
t43 = -pkin(3) - t38;
t42 = t46 * t21;
t20 = sin(qJ(1));
t40 = pkin(1) * t20;
t22 = cos(qJ(1));
t39 = pkin(1) * t22;
t17 = qJ(1) + pkin(8);
t15 = qJ(3) + t17;
t10 = sin(t15);
t19 = sin(qJ(4));
t36 = t10 * t19;
t11 = cos(t15);
t35 = t11 * t19;
t34 = rSges(5,2) * t36 + t45 * t11;
t33 = -t11 * rSges(4,1) + t10 * rSges(4,2);
t32 = -pkin(3) - t42;
t13 = sin(t17);
t31 = -pkin(2) * t13 - t40;
t14 = cos(t17);
t30 = -pkin(2) * t14 - t39;
t29 = -rSges(4,1) * t10 - rSges(4,2) * t11;
t27 = rSges(5,2) * t35 + t43 * t11;
t25 = rSges(6,2) * t35 - t44 * t10 + t32 * t11;
t24 = rSges(6,2) * t36 + t32 * t10 + t44 * t11;
t23 = (-g(2) * t45 + g(3) * t43) * t10;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t22 + t20 * rSges(2,2)) + g(3) * (-t20 * rSges(2,1) - rSges(2,2) * t22)) - m(3) * (g(2) * (-t14 * rSges(3,1) + t13 * rSges(3,2) - t39) + g(3) * (-rSges(3,1) * t13 - rSges(3,2) * t14 - t40)) - m(4) * (g(2) * (t30 + t33) + g(3) * (t29 + t31)) - m(5) * (g(2) * (t27 + t30) + g(3) * (t31 + t34) + t23) - m(6) * (g(2) * (t25 + t30) + g(3) * (t24 + t31)), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t33 + g(3) * t29) - m(5) * (g(2) * t27 + g(3) * t34 + t23) - m(6) * (g(2) * t25 + g(3) * t24), (-m(5) * (-rSges(5,2) * t19 + t38) - m(6) * (-rSges(6,2) * t19 + t42)) * g(1) + (-g(2) * t10 + g(3) * t11) * (m(5) * (rSges(5,1) * t19 + rSges(5,2) * t21) + m(6) * (rSges(6,2) * t21 + t46 * t19)), -m(6) * (g(2) * t11 + g(3) * t10)];
taug = t1(:);
