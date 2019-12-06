% Calculate Gravitation load on the joints for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:38
% EndTime: 2019-12-05 18:08:39
% DurationCPUTime: 0.38s
% Computational Cost: add. (123->81), mult. (296->135), div. (0->0), fcn. (293->8), ass. (0->33)
t12 = sin(qJ(5));
t16 = cos(qJ(5));
t25 = -rSges(6,1) * t16 + rSges(6,2) * t12;
t13 = sin(qJ(4));
t14 = sin(qJ(3));
t36 = t13 * t14;
t15 = sin(qJ(1));
t35 = t14 * t15;
t34 = t14 * t16;
t17 = cos(qJ(4));
t33 = t14 * t17;
t19 = cos(qJ(1));
t32 = t14 * t19;
t18 = cos(qJ(3));
t31 = t15 * t18;
t30 = t18 * t12;
t29 = t18 * t16;
t28 = t19 * t13;
t27 = t19 * t17;
t26 = rSges(4,1) * t18 - rSges(4,2) * t14;
t4 = t17 * t31 - t28;
t24 = -t12 * t4 + t15 * t34;
t23 = -t12 * t35 - t16 * t4;
t22 = -t16 * t33 + t30;
t21 = t12 * t33 + t29;
t11 = t19 * qJ(2);
t10 = t15 * qJ(2);
t6 = t13 * t15 + t18 * t27;
t5 = -t15 * t17 + t18 * t28;
t3 = t13 * t31 + t27;
t2 = t12 * t32 + t16 * t6;
t1 = -t12 * t6 + t16 * t32;
t7 = [-m(2) * (g(1) * (-rSges(2,1) * t15 - rSges(2,2) * t19) + g(2) * (rSges(2,1) * t19 - rSges(2,2) * t15)) - m(3) * (g(1) * (-rSges(3,1) * t15 + rSges(3,3) * t19 + t11) + g(2) * (rSges(3,1) * t19 + rSges(3,3) * t15 + t10)) - m(4) * (g(1) * (t19 * rSges(4,3) - t15 * t26 + t11) + g(2) * (t15 * rSges(4,3) + t19 * t26 + t10)) - m(5) * (g(1) * (-rSges(5,1) * t4 + rSges(5,2) * t3 - rSges(5,3) * t35 + t11) + g(2) * (rSges(5,1) * t6 - rSges(5,2) * t5 + rSges(5,3) * t32 + t10)) - m(6) * (g(1) * (rSges(6,1) * t23 - rSges(6,2) * t24 - t3 * rSges(6,3) + t11) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t1 + rSges(6,3) * t5 + t10)), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t15 - g(2) * t19), (-m(4) * t26 - m(5) * (t14 * rSges(5,3) + (rSges(5,1) * t17 - rSges(5,2) * t13) * t18) - m(6) * ((t12 * t14 + t17 * t29) * rSges(6,1) + (-t17 * t30 + t34) * rSges(6,2) + t18 * t13 * rSges(6,3))) * g(3) + (g(1) * t19 + g(2) * t15) * (-m(4) * (-rSges(4,1) * t14 - rSges(4,2) * t18) - m(5) * (-rSges(5,1) * t33 + rSges(5,2) * t36 + rSges(5,3) * t18) - m(6) * (rSges(6,1) * t22 + rSges(6,2) * t21 - rSges(6,3) * t36)), -m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t6) + g(2) * (-rSges(5,1) * t3 - rSges(5,2) * t4)) - m(6) * (g(1) * (t6 * rSges(6,3) + t25 * t5) + g(2) * (t4 * rSges(6,3) + t25 * t3)) + (-m(5) * (-rSges(5,1) * t13 - rSges(5,2) * t17) - m(6) * (rSges(6,3) * t17 + t13 * t25)) * g(3) * t14, -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t24 + rSges(6,2) * t23) + g(3) * (-rSges(6,1) * t21 + rSges(6,2) * t22))];
taug = t7(:);
