% Calculate Gravitation load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:19
% DurationCPUTime: 0.20s
% Computational Cost: add. (266->59), mult. (140->68), div. (0->0), fcn. (102->10), ass. (0->34)
t40 = rSges(6,3) + pkin(8);
t18 = cos(qJ(5));
t34 = rSges(6,1) * t18;
t39 = -pkin(4) - t34;
t15 = qJ(1) + pkin(9);
t14 = qJ(3) + t15;
t9 = sin(t14);
t38 = pkin(3) * t9;
t17 = sin(qJ(1));
t37 = pkin(1) * t17;
t19 = cos(qJ(1));
t36 = pkin(1) * t19;
t10 = cos(t14);
t35 = pkin(3) * t10;
t16 = sin(qJ(5));
t33 = rSges(6,2) * t16;
t11 = qJ(4) + t14;
t7 = sin(t11);
t8 = cos(t11);
t32 = t7 * t33 + t40 * t8;
t31 = -t8 * rSges(5,1) + t7 * rSges(5,2);
t30 = -t10 * rSges(4,1) + t9 * rSges(4,2);
t29 = -rSges(5,1) * t7 - rSges(5,2) * t8;
t12 = sin(t15);
t28 = -pkin(2) * t12 - t37;
t13 = cos(t15);
t27 = -pkin(2) * t13 - t36;
t26 = -rSges(4,1) * t9 - rSges(4,2) * t10;
t24 = t31 - t35;
t23 = (t33 + t39) * t8;
t22 = t28 - t38;
t21 = t23 - t35;
t20 = (-g(2) * t40 + g(3) * t39) * t7;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t19 + t17 * rSges(2,2)) + g(3) * (-t17 * rSges(2,1) - rSges(2,2) * t19)) - m(3) * (g(2) * (-t13 * rSges(3,1) + t12 * rSges(3,2) - t36) + g(3) * (-rSges(3,1) * t12 - rSges(3,2) * t13 - t37)) - m(4) * (g(2) * (t27 + t30) + g(3) * (t26 + t28)) - m(5) * (g(2) * (t24 + t27) + g(3) * (t22 + t29)) - m(6) * (g(2) * (t21 + t27) + g(3) * (t22 + t32) + t20), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t30 + g(3) * t26) - m(5) * (g(2) * t24 + g(3) * (t29 - t38)) - m(6) * (g(2) * t21 + g(3) * (t32 - t38) + t20), -m(5) * (g(2) * t31 + g(3) * t29) - m(6) * (g(2) * t23 + g(3) * t32 + t20), -m(6) * (g(1) * (-t33 + t34) + (g(2) * t7 - g(3) * t8) * (rSges(6,1) * t16 + rSges(6,2) * t18))];
taug = t1(:);
