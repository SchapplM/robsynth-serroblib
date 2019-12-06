% Calculate Gravitation load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:31
% EndTime: 2019-12-05 16:43:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (193->51), mult. (160->62), div. (0->0), fcn. (122->6), ass. (0->26)
t39 = rSges(6,1) + pkin(4);
t14 = qJ(3) + qJ(4);
t10 = cos(t14);
t9 = sin(t14);
t27 = t10 * rSges(5,1) - rSges(5,2) * t9;
t38 = -rSges(6,2) * t10 - t39 * t9;
t13 = pkin(8) + qJ(2);
t7 = sin(t13);
t8 = cos(t13);
t37 = g(1) * t8 + g(2) * t7;
t26 = -rSges(6,2) * t9 + t39 * t10;
t17 = -pkin(7) - pkin(6);
t15 = sin(qJ(3));
t31 = pkin(3) * t15;
t30 = rSges(4,3) + pkin(6);
t16 = cos(qJ(3));
t11 = t16 * pkin(3);
t6 = t11 + pkin(2);
t29 = rSges(5,3) - t17;
t28 = rSges(6,3) + qJ(5) - t17;
t25 = -rSges(5,1) * t9 - rSges(5,2) * t10;
t23 = rSges(4,1) * t16 - rSges(4,2) * t15;
t22 = t6 + t27;
t21 = t6 + t26;
t20 = pkin(2) + t23;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t7 - rSges(3,2) * t8) + g(2) * (rSges(3,1) * t8 - rSges(3,2) * t7)) - m(4) * ((g(1) * t30 + g(2) * t20) * t8 + (-g(1) * t20 + g(2) * t30) * t7) - m(5) * ((g(1) * t29 + g(2) * t22) * t8 + (-g(1) * t22 + g(2) * t29) * t7) - m(6) * ((g(1) * t28 + g(2) * t21) * t8 + (-g(1) * t21 + g(2) * t28) * t7), (-m(4) * t23 - m(5) * (t11 + t27) - m(6) * (t11 + t26)) * g(3) + t37 * (-m(4) * (-rSges(4,1) * t15 - rSges(4,2) * t16) - m(5) * (t25 - t31) - m(6) * (-t31 + t38)), (-m(5) * t27 - m(6) * t26) * g(3) + t37 * (-m(5) * t25 - m(6) * t38), -m(6) * (g(1) * t7 - g(2) * t8)];
taug = t1(:);
