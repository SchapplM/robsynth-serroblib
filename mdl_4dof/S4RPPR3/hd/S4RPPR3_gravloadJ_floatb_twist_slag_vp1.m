% Calculate Gravitation load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (90->36), mult. (79->48), div. (0->0), fcn. (58->8), ass. (0->19)
t23 = -m(4) - m(5);
t13 = sin(qJ(1));
t22 = pkin(1) * t13;
t21 = rSges(5,3) + pkin(5) + qJ(3);
t20 = rSges(4,3) + qJ(3);
t8 = pkin(7) + qJ(4);
t3 = sin(t8);
t5 = cos(t8);
t19 = rSges(5,1) * t5 - rSges(5,2) * t3;
t14 = cos(qJ(1));
t7 = t14 * pkin(1);
t17 = -g(1) * t22 + g(2) * t7;
t11 = cos(pkin(7));
t16 = pkin(3) * t11 + pkin(2) + t19;
t15 = rSges(4,1) * t11 - rSges(4,2) * sin(pkin(7)) + pkin(2);
t9 = qJ(1) + pkin(6);
t6 = cos(t9);
t4 = sin(t9);
t1 = [-m(2) * (g(1) * (-t13 * rSges(2,1) - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - t13 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t4 - rSges(3,2) * t6 - t22) + g(2) * (rSges(3,1) * t6 - rSges(3,2) * t4 + t7)) - m(4) * ((g(1) * t20 + g(2) * t15) * t6 + (-g(1) * t15 + g(2) * t20) * t4 + t17) - m(5) * ((g(1) * t21 + g(2) * t16) * t6 + (-g(1) * t16 + g(2) * t21) * t4 + t17), (-m(3) + t23) * g(3), t23 * (g(1) * t4 - g(2) * t6), -m(5) * (g(3) * t19 + (g(1) * t6 + g(2) * t4) * (-rSges(5,1) * t3 - rSges(5,2) * t5))];
taug = t1(:);
