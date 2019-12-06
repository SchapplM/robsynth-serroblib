% Calculate Gravitation load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:40:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (105->53), mult. (236->78), div. (0->0), fcn. (217->6), ass. (0->24)
t12 = sin(pkin(7));
t13 = cos(pkin(7));
t37 = g(1) * t13 + g(2) * t12;
t32 = rSges(6,1) + pkin(4);
t27 = rSges(6,3) + qJ(5);
t36 = -pkin(2) - pkin(6);
t17 = cos(qJ(2));
t33 = g(3) * t17;
t15 = sin(qJ(2));
t31 = t17 * pkin(2) + t15 * qJ(3);
t14 = sin(qJ(4));
t30 = t14 * t15;
t16 = cos(qJ(4));
t29 = t15 * t16;
t26 = -m(4) - m(5) - m(6);
t23 = t37 * qJ(3) * t17;
t21 = rSges(5,1) * t14 + rSges(5,2) * t16;
t19 = t32 * t14 - t27 * t16;
t18 = g(3) * (t17 * pkin(6) + t31) + t23;
t5 = -t12 * t30 + t13 * t16;
t4 = t12 * t29 + t13 * t14;
t3 = t12 * t16 + t13 * t30;
t2 = t12 * t14 - t13 * t29;
t1 = [(-m(2) - m(3) + t26) * g(3), -m(3) * (g(3) * (rSges(3,1) * t17 - t15 * rSges(3,2)) + t37 * (-rSges(3,1) * t15 - rSges(3,2) * t17)) - m(4) * (g(3) * (-rSges(4,2) * t17 + t15 * rSges(4,3) + t31) + t23 + t37 * (rSges(4,3) * t17 + (rSges(4,2) - pkin(2)) * t15)) - m(5) * ((g(3) * rSges(5,3) + t37 * t21) * t17 + (g(3) * t21 + t37 * (-rSges(5,3) + t36)) * t15 + t18) - m(6) * ((g(3) * rSges(6,2) + t37 * t19) * t17 + (g(3) * t19 + t37 * (-rSges(6,2) + t36)) * t15 + t18), t26 * (t37 * t15 - t33), -m(5) * (g(1) * (-rSges(5,1) * t2 - rSges(5,2) * t3) + g(2) * (rSges(5,1) * t4 + rSges(5,2) * t5)) - m(6) * (g(1) * (-t32 * t2 + t27 * t3) + g(2) * (-t27 * t5 + t32 * t4)) + (-m(5) * (-rSges(5,1) * t16 + rSges(5,2) * t14) - m(6) * (-t27 * t14 - t32 * t16)) * t33, -m(6) * (g(1) * t2 - g(2) * t4 + t16 * t33)];
taug = t1(:);
