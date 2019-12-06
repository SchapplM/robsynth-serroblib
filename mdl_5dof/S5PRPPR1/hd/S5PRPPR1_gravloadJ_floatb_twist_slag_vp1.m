% Calculate Gravitation load on the joints for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:46
% EndTime: 2019-12-05 15:21:48
% DurationCPUTime: 0.27s
% Computational Cost: add. (167->59), mult. (150->83), div. (0->0), fcn. (136->8), ass. (0->26)
t17 = sin(pkin(9));
t18 = sin(pkin(8));
t19 = cos(pkin(9));
t20 = cos(pkin(8));
t34 = (rSges(5,3) + qJ(4)) * t18 + (rSges(5,1) * t19 - rSges(5,2) * t17 + pkin(3)) * t20;
t31 = -m(5) - m(6);
t16 = pkin(7) + qJ(2);
t12 = sin(t16);
t14 = cos(t16);
t30 = t14 * pkin(2) + t12 * qJ(3);
t29 = pkin(4) * t17;
t28 = rSges(4,2) * t18;
t27 = t12 * t20;
t26 = t14 * t20;
t25 = -m(4) + t31;
t24 = t17 * rSges(5,1) + t19 * rSges(5,2);
t22 = (pkin(4) * t19 + pkin(3)) * t20 + (rSges(6,3) + pkin(6) + qJ(4)) * t18;
t15 = pkin(9) + qJ(5);
t13 = cos(t15);
t11 = sin(t15);
t8 = t14 * qJ(3);
t4 = t11 * t12 + t13 * t26;
t3 = -t11 * t26 + t12 * t13;
t2 = t11 * t14 - t13 * t27;
t1 = t11 * t27 + t13 * t14;
t5 = [(-m(2) - m(3) + t25) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t12 - rSges(3,2) * t14) + g(2) * (rSges(3,1) * t14 - rSges(3,2) * t12)) - m(4) * (g(1) * (rSges(4,3) * t14 + t8) + g(2) * (rSges(4,1) * t26 - t14 * t28 + t30) + (g(1) * (-rSges(4,1) * t20 - pkin(2) + t28) + g(2) * rSges(4,3)) * t12) - m(5) * (g(1) * t8 + g(2) * t30 + (g(1) * t24 + t34 * g(2)) * t14 + (g(2) * t24 + (-pkin(2) - t34) * g(1)) * t12) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t30) + (g(1) * t29 + g(2) * t22) * t14 + (g(1) * (-pkin(2) - t22) + g(2) * t29) * t12), t25 * (g(1) * t12 - g(2) * t14), t31 * (-g(3) * t20 + (g(1) * t14 + g(2) * t12) * t18), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t11 - rSges(6,2) * t13) * t18)];
taug = t5(:);
