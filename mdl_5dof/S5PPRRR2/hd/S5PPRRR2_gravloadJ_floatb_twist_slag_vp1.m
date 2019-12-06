% Calculate Gravitation load on the joints for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:10
% EndTime: 2019-12-05 15:14:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (169->46), mult. (188->68), div. (0->0), fcn. (167->8), ass. (0->28)
t12 = qJ(4) + qJ(5);
t10 = cos(t12);
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t9 = sin(t12);
t38 = m(5) * (t16 * rSges(5,1) - t15 * rSges(5,2) + pkin(3)) + m(6) * (rSges(6,1) * t10 - rSges(6,2) * t9 + pkin(4) * t16 + pkin(3)) + m(4) * rSges(4,1);
t11 = pkin(9) + qJ(3);
t7 = sin(t11);
t37 = g(3) * t7;
t14 = cos(pkin(8));
t30 = t10 * t14;
t13 = sin(pkin(8));
t31 = t10 * t13;
t32 = t14 * t9;
t33 = t13 * t9;
t8 = cos(t11);
t36 = (-t8 * t33 - t30) * rSges(6,1) + (-t8 * t31 + t32) * rSges(6,2);
t35 = (-t8 * t32 + t31) * rSges(6,1) + (-t8 * t30 - t33) * rSges(6,2);
t29 = t13 * t15;
t28 = t13 * t16;
t27 = t14 * t15;
t26 = t14 * t16;
t25 = -m(3) - m(4) - m(5) - m(6);
t24 = -rSges(6,1) * t9 - rSges(6,2) * t10;
t22 = -t8 * t27 + t28;
t21 = -t8 * t29 - t26;
t19 = m(4) * rSges(4,2) - m(5) * (rSges(5,3) + pkin(6)) - m(6) * (rSges(6,3) + pkin(7) + pkin(6));
t1 = [(-m(2) + t25) * g(3), t25 * (g(1) * t13 - g(2) * t14), (t19 * t7 - t38 * t8) * g(3) + (g(1) * t14 + g(2) * t13) * (t19 * t8 + t38 * t7), -m(5) * (g(1) * (t22 * rSges(5,1) + (-t8 * t26 - t29) * rSges(5,2)) + g(2) * (t21 * rSges(5,1) + (-t8 * t28 + t27) * rSges(5,2))) - m(6) * (g(1) * (t22 * pkin(4) + t35) + g(2) * (t21 * pkin(4) + t36)) + (-m(5) * (-rSges(5,1) * t15 - rSges(5,2) * t16) - m(6) * (-pkin(4) * t15 + t24)) * t37, -m(6) * (g(1) * t35 + g(2) * t36 + t24 * t37)];
taug = t1(:);
