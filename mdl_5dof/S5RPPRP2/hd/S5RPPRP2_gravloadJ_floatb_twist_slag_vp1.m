% Calculate Gravitation load on the joints for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:06
% DurationCPUTime: 0.25s
% Computational Cost: add. (176->53), mult. (142->67), div. (0->0), fcn. (111->8), ass. (0->26)
t29 = rSges(6,3) + qJ(5);
t31 = rSges(6,1) + pkin(4);
t9 = pkin(8) + qJ(4);
t4 = sin(t9);
t6 = cos(t9);
t17 = t29 * t4 + t31 * t6;
t10 = qJ(1) + pkin(7);
t5 = sin(t10);
t7 = cos(t10);
t30 = g(1) * t7 + g(2) * t5;
t14 = sin(qJ(1));
t26 = pkin(1) * t14;
t13 = -pkin(6) - qJ(3);
t25 = rSges(6,2) - t13;
t24 = rSges(5,3) - t13;
t23 = rSges(4,3) + qJ(3);
t22 = -m(4) - m(5) - m(6);
t21 = g(1) * t26;
t20 = rSges(5,1) * t6 - rSges(5,2) * t4;
t12 = cos(pkin(8));
t19 = rSges(4,1) * t12 - rSges(4,2) * sin(pkin(8)) + pkin(2);
t3 = pkin(3) * t12 + pkin(2);
t15 = cos(qJ(1));
t8 = t15 * pkin(1);
t18 = -t21 + g(2) * (t7 * t3 + t8);
t1 = [-m(2) * (g(1) * (-t14 * rSges(2,1) - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - t14 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t5 - rSges(3,2) * t7 - t26) + g(2) * (rSges(3,1) * t7 - rSges(3,2) * t5 + t8)) - m(4) * (-t21 + g(2) * t8 + (g(1) * t23 + g(2) * t19) * t7 + (-g(1) * t19 + g(2) * t23) * t5) - m(5) * ((g(1) * t24 + g(2) * t20) * t7 + (g(1) * (-t20 - t3) + g(2) * t24) * t5 + t18) - m(6) * ((g(1) * t25 + g(2) * t17) * t7 + (g(1) * (-t17 - t3) + g(2) * t25) * t5 + t18), (-m(3) + t22) * g(3), t22 * (g(1) * t5 - g(2) * t7), (-m(5) * t20 - m(6) * t17) * g(3) + t30 * (-m(5) * (-rSges(5,1) * t4 - rSges(5,2) * t6) - m(6) * (t29 * t6 - t31 * t4)), -m(6) * (-g(3) * t6 + t30 * t4)];
taug = t1(:);
