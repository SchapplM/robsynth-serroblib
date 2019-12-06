% Calculate Gravitation load on the joints for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:25
% EndTime: 2019-12-05 14:59:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (100->42), mult. (241->78), div. (0->0), fcn. (272->10), ass. (0->30)
t32 = rSges(6,3) + pkin(6);
t14 = sin(pkin(9));
t15 = sin(pkin(8));
t31 = t14 * t15;
t21 = sin(qJ(4));
t30 = t15 * t21;
t23 = cos(qJ(4));
t29 = t15 * t23;
t16 = sin(pkin(7));
t18 = cos(pkin(8));
t28 = t16 * t18;
t19 = cos(pkin(7));
t27 = t18 * t19;
t26 = -m(4) - m(5) - m(6);
t25 = -m(3) + t26;
t20 = sin(qJ(5));
t22 = cos(qJ(5));
t24 = rSges(6,1) * t22 - rSges(6,2) * t20 + pkin(4);
t17 = cos(pkin(9));
t11 = t17 * t29 - t18 * t21;
t10 = -t17 * t30 - t18 * t23;
t9 = t14 * t16 + t17 * t27;
t8 = t14 * t27 - t16 * t17;
t7 = -t14 * t19 + t17 * t28;
t6 = t14 * t28 + t17 * t19;
t4 = t19 * t30 + t9 * t23;
t3 = t19 * t29 - t9 * t21;
t2 = t16 * t30 + t7 * t23;
t1 = t16 * t29 - t7 * t21;
t5 = [(-m(2) + t25) * g(3), t25 * (g(1) * t16 - g(2) * t19), t26 * (-g(3) * t18 + (g(1) * t19 + g(2) * t16) * t15), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (rSges(5,1) * t1 - rSges(5,2) * t2) + g(3) * (rSges(5,1) * t10 - rSges(5,2) * t11)) - m(6) * (g(1) * (t24 * t3 + t32 * t4) + (t24 * t10 + t32 * t11) * g(3) + (t24 * t1 + t32 * t2) * g(2)), -m(6) * (g(1) * ((-t4 * t20 + t8 * t22) * rSges(6,1) + (-t8 * t20 - t4 * t22) * rSges(6,2)) + g(2) * ((-t2 * t20 + t6 * t22) * rSges(6,1) + (-t2 * t22 - t6 * t20) * rSges(6,2)) + g(3) * ((-t11 * t20 + t22 * t31) * rSges(6,1) + (-t11 * t22 - t20 * t31) * rSges(6,2)))];
taug = t5(:);
