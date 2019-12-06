% Calculate Gravitation load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:04:48
% DurationCPUTime: 0.31s
% Computational Cost: add. (131->53), mult. (201->89), div. (0->0), fcn. (202->10), ass. (0->33)
t14 = cos(pkin(8));
t15 = cos(pkin(7));
t17 = sin(qJ(3));
t27 = t15 * t17;
t13 = sin(pkin(7));
t19 = cos(qJ(3));
t29 = t13 * t19;
t37 = -t14 * t27 + t29;
t34 = rSges(6,3) + pkin(6);
t36 = -m(5) - m(6);
t35 = pkin(3) * t17;
t12 = sin(pkin(8));
t16 = sin(qJ(5));
t33 = t12 * t16;
t18 = cos(qJ(5));
t32 = t12 * t18;
t31 = t13 * t14;
t30 = t13 * t17;
t28 = t14 * t15;
t26 = t15 * t19;
t24 = -m(3) - m(4) + t36;
t23 = t37 * pkin(3);
t22 = rSges(6,1) * t18 - rSges(6,2) * t16 + pkin(4);
t21 = -t14 * t30 - t26;
t20 = t21 * pkin(3);
t11 = qJ(3) + pkin(9);
t10 = cos(t11);
t9 = sin(t11);
t4 = t10 * t28 + t13 * t9;
t3 = t10 * t13 - t9 * t28;
t2 = t10 * t31 - t15 * t9;
t1 = -t10 * t15 - t9 * t31;
t5 = [(-m(2) + t24) * g(3), t24 * (g(1) * t13 - g(2) * t15), -m(4) * (g(1) * (t37 * rSges(4,1) + (-t14 * t26 - t30) * rSges(4,2)) + g(2) * (t21 * rSges(4,1) + (-t14 * t29 + t27) * rSges(4,2))) - m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4 + t23) + g(2) * (t1 * rSges(5,1) - t2 * rSges(5,2) + t20)) - m(6) * (g(1) * (t22 * t3 + t34 * t4 + t23) + g(2) * (t22 * t1 + t34 * t2 + t20)) + (-m(4) * (-rSges(4,1) * t17 - rSges(4,2) * t19) - m(5) * (-rSges(5,2) * t10 - t35) - m(6) * (t34 * t10 - t35) + (m(5) * rSges(5,1) + m(6) * t22) * t9) * g(3) * t12, t36 * (-g(3) * t14 + (g(1) * t15 + g(2) * t13) * t12), -m(6) * (g(1) * ((t15 * t32 - t16 * t4) * rSges(6,1) + (-t15 * t33 - t18 * t4) * rSges(6,2)) + g(2) * ((t13 * t32 - t16 * t2) * rSges(6,1) + (-t13 * t33 - t18 * t2) * rSges(6,2)) + g(3) * ((-t10 * t33 - t14 * t18) * rSges(6,1) + (-t10 * t32 + t14 * t16) * rSges(6,2)))];
taug = t5(:);
