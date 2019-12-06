% Calculate Gravitation load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:39
% DurationCPUTime: 0.19s
% Computational Cost: add. (175->43), mult. (125->56), div. (0->0), fcn. (93->8), ass. (0->25)
t15 = pkin(9) + qJ(4);
t10 = cos(t15);
t12 = qJ(5) + t15;
t5 = sin(t12);
t6 = cos(t12);
t26 = t6 * rSges(6,1) - rSges(6,2) * t5;
t35 = pkin(4) * t10 + t26;
t16 = pkin(8) + qJ(2);
t11 = cos(t16);
t9 = sin(t16);
t34 = g(1) * t11 + g(2) * t9;
t18 = cos(pkin(9));
t7 = t18 * pkin(3) + pkin(2);
t19 = -pkin(6) - qJ(3);
t30 = rSges(5,3) - t19;
t29 = rSges(6,3) + pkin(7) - t19;
t28 = rSges(4,3) + qJ(3);
t27 = -m(4) - m(5) - m(6);
t25 = -rSges(6,1) * t5 - rSges(6,2) * t6;
t8 = sin(t15);
t24 = rSges(5,1) * t10 - rSges(5,2) * t8;
t23 = t7 + t35;
t22 = t24 + t7;
t21 = rSges(4,1) * t18 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t1 = [(-m(2) - m(3) + t27) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t9 - rSges(3,2) * t11) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t9)) - m(4) * ((-g(1) * t21 + g(2) * t28) * t9 + (g(1) * t28 + g(2) * t21) * t11) - m(5) * ((-g(1) * t22 + g(2) * t30) * t9 + (g(1) * t30 + g(2) * t22) * t11) - m(6) * ((-g(1) * t23 + g(2) * t29) * t9 + (g(1) * t29 + g(2) * t23) * t11), t27 * (g(1) * t9 - g(2) * t11), (-m(5) * t24 - m(6) * t35) * g(3) + t34 * (-m(5) * (-rSges(5,1) * t8 - rSges(5,2) * t10) - m(6) * (-pkin(4) * t8 + t25)), -m(6) * (g(3) * t26 + t34 * t25)];
taug = t1(:);
