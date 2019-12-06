% Calculate Gravitation load on the joints for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:49
% EndTime: 2019-12-05 14:57:51
% DurationCPUTime: 0.17s
% Computational Cost: add. (114->37), mult. (159->67), div. (0->0), fcn. (159->8), ass. (0->22)
t24 = rSges(6,3) + pkin(6);
t11 = sin(pkin(8));
t15 = sin(qJ(5));
t23 = t11 * t15;
t16 = cos(qJ(5));
t22 = t11 * t16;
t12 = sin(pkin(7));
t13 = cos(pkin(8));
t21 = t12 * t13;
t14 = cos(pkin(7));
t20 = t13 * t14;
t19 = -m(4) - m(5) - m(6);
t18 = -m(3) + t19;
t17 = rSges(6,1) * t16 - rSges(6,2) * t15 + pkin(4);
t10 = pkin(9) + qJ(4);
t9 = cos(t10);
t8 = sin(t10);
t4 = t12 * t8 + t9 * t20;
t3 = t12 * t9 - t8 * t20;
t2 = -t14 * t8 + t9 * t21;
t1 = -t14 * t9 - t8 * t21;
t5 = [(-m(2) + t18) * g(3), t18 * (g(1) * t12 - g(2) * t14), t19 * (-g(3) * t13 + (g(1) * t14 + g(2) * t12) * t11), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (t17 * t3 + t24 * t4) + g(2) * (t17 * t1 + t24 * t2)) + ((m(5) * rSges(5,2) - m(6) * t24) * t9 + (m(5) * rSges(5,1) + m(6) * t17) * t8) * g(3) * t11, -m(6) * (g(1) * ((t14 * t22 - t4 * t15) * rSges(6,1) + (-t14 * t23 - t16 * t4) * rSges(6,2)) + g(2) * ((t12 * t22 - t2 * t15) * rSges(6,1) + (-t12 * t23 - t16 * t2) * rSges(6,2)) + g(3) * ((-t13 * t16 - t9 * t23) * rSges(6,1) + (t13 * t15 - t9 * t22) * rSges(6,2)))];
taug = t5(:);
