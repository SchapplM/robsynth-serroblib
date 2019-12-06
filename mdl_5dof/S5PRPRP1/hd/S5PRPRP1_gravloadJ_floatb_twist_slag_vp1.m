% Calculate Gravitation load on the joints for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:18
% EndTime: 2019-12-05 15:28:20
% DurationCPUTime: 0.23s
% Computational Cost: add. (165->44), mult. (128->55), div. (0->0), fcn. (99->6), ass. (0->22)
t24 = rSges(6,3) + qJ(5);
t26 = rSges(6,1) + pkin(4);
t8 = pkin(8) + qJ(4);
t4 = sin(t8);
t6 = cos(t8);
t14 = t24 * t4 + t26 * t6;
t9 = pkin(7) + qJ(2);
t7 = cos(t9);
t27 = g(2) * t7;
t5 = sin(t9);
t25 = g(1) * t7 + g(2) * t5;
t11 = cos(pkin(8));
t3 = pkin(3) * t11 + pkin(2);
t22 = t3 * t27;
t12 = -pkin(6) - qJ(3);
t20 = rSges(6,2) - t12;
t19 = rSges(5,3) - t12;
t18 = rSges(4,3) + qJ(3);
t17 = -m(4) - m(5) - m(6);
t16 = rSges(5,1) * t6 - rSges(5,2) * t4;
t15 = rSges(4,1) * t11 - rSges(4,2) * sin(pkin(8)) + pkin(2);
t1 = [(-m(2) - m(3) + t17) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t5 - rSges(3,2) * t7) + g(2) * (rSges(3,1) * t7 - rSges(3,2) * t5)) - m(4) * ((g(1) * t18 + g(2) * t15) * t7 + (-g(1) * t15 + g(2) * t18) * t5) - m(5) * (t22 + (g(1) * t19 + g(2) * t16) * t7 + (g(1) * (-t16 - t3) + g(2) * t19) * t5) - m(6) * (t22 + (g(1) * t20 + g(2) * t14) * t7 + (g(1) * (-t14 - t3) + g(2) * t20) * t5), t17 * (g(1) * t5 - t27), (-m(5) * t16 - m(6) * t14) * g(3) + t25 * (-m(5) * (-rSges(5,1) * t4 - rSges(5,2) * t6) - m(6) * (t24 * t6 - t26 * t4)), -m(6) * (-g(3) * t6 + t25 * t4)];
taug = t1(:);
