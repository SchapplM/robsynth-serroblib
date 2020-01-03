% Calculate Gravitation load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:53:43
% DurationCPUTime: 0.21s
% Computational Cost: add. (125->54), mult. (162->79), div. (0->0), fcn. (143->8), ass. (0->28)
t32 = rSges(5,3) + pkin(6);
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t31 = m(4) * rSges(4,1) + m(5) * (rSges(5,1) * t16 - rSges(5,2) * t14 + pkin(3));
t17 = cos(qJ(1));
t29 = t14 * t17;
t15 = sin(qJ(1));
t28 = t15 * t14;
t27 = t15 * t16;
t26 = t16 * t17;
t13 = -pkin(5) - qJ(2);
t25 = rSges(4,3) - t13;
t24 = rSges(3,3) + qJ(2);
t10 = pkin(7) + qJ(3);
t8 = sin(t10);
t9 = cos(t10);
t23 = rSges(4,1) * t9 - rSges(4,2) * t8;
t12 = cos(pkin(7));
t22 = rSges(3,1) * t12 - rSges(3,2) * sin(pkin(7)) + pkin(1);
t20 = m(4) * rSges(4,2) - m(5) * t32;
t19 = pkin(3) * t9 + t32 * t8;
t7 = pkin(2) * t12 + pkin(1);
t5 = t17 * t7;
t4 = t9 * t26 + t28;
t3 = -t9 * t29 + t27;
t2 = -t9 * t27 + t29;
t1 = t9 * t28 + t26;
t6 = [-m(2) * (g(1) * (-t15 * rSges(2,1) - rSges(2,2) * t17) + g(2) * (rSges(2,1) * t17 - t15 * rSges(2,2))) - m(3) * ((g(1) * t24 + g(2) * t22) * t17 + (-g(1) * t22 + g(2) * t24) * t15) - m(4) * (g(2) * t5 + (g(1) * t25 + g(2) * t23) * t17 + (g(1) * (-t23 - t7) + g(2) * t25) * t15) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t5) + (-g(1) * t13 + g(2) * t19) * t17 + (g(1) * (-t19 - t7) - g(2) * t13) * t15), (-m(3) - m(4) - m(5)) * (g(1) * t15 - g(2) * t17), (t20 * t8 - t31 * t9) * g(3) + (g(1) * t17 + g(2) * t15) * (t20 * t9 + t31 * t8), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2) + g(3) * (-rSges(5,1) * t14 - rSges(5,2) * t16) * t8)];
taug = t6(:);
