% Calculate Gravitation load on the joints for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (71->51), mult. (143->77), div. (0->0), fcn. (132->6), ass. (0->22)
t27 = -rSges(5,3) - pkin(5);
t26 = -m(4) - m(5);
t15 = sin(qJ(1));
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + t15 * qJ(2);
t12 = sin(pkin(6));
t24 = t12 * t17;
t13 = cos(pkin(6));
t23 = t13 * t17;
t22 = qJ(3) * t12;
t21 = pkin(2) * t23 + t17 * t22 + t25;
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t20 = t12 * t16 - t13 * t14;
t19 = t12 * t14 + t13 * t16;
t18 = -pkin(2) * t13 - pkin(1) - t22;
t10 = t17 * qJ(2);
t5 = t19 * t17;
t4 = t20 * t17;
t3 = t19 * t15;
t2 = t20 * t15;
t1 = [-m(2) * (g(1) * (-t15 * rSges(2,1) - rSges(2,2) * t17) + g(2) * (rSges(2,1) * t17 - t15 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t17 + t10) + g(2) * (rSges(3,1) * t23 - rSges(3,2) * t24 + t25) + (g(1) * (-rSges(3,1) * t13 + rSges(3,2) * t12 - pkin(1)) + g(2) * rSges(3,3)) * t15) - m(4) * (g(1) * (rSges(4,2) * t17 + t10) + g(2) * (rSges(4,1) * t23 + rSges(4,3) * t24 + t21) + (g(1) * (-rSges(4,1) * t13 - rSges(4,3) * t12 + t18) + g(2) * rSges(4,2)) * t15) - m(5) * (g(1) * (-t3 * rSges(5,1) - t2 * rSges(5,2) + t27 * t17 + t10) + g(2) * (t5 * rSges(5,1) + t4 * rSges(5,2) + pkin(3) * t23 + t21) + (g(1) * (-pkin(3) * t13 + t18) + g(2) * t27) * t15), (-m(3) + t26) * (g(1) * t15 - g(2) * t17), t26 * (-g(3) * t13 + (g(1) * t17 + g(2) * t15) * t12), -m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t5) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t3) + g(3) * (-t19 * rSges(5,1) - t20 * rSges(5,2)))];
taug = t1(:);
