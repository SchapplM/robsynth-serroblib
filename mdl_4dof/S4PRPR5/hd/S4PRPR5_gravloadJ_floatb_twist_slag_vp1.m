% Calculate Gravitation load on the joints for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:22:59
% DurationCPUTime: 0.17s
% Computational Cost: add. (75->33), mult. (107->49), div. (0->0), fcn. (89->8), ass. (0->19)
t22 = rSges(5,3) + pkin(5);
t10 = cos(qJ(4));
t8 = sin(qJ(4));
t21 = m(4) * rSges(4,1) + m(5) * (rSges(5,1) * t10 - rSges(5,2) * t8 + pkin(3));
t9 = sin(qJ(2));
t20 = pkin(2) * t9;
t6 = sin(pkin(6));
t19 = t6 * t8;
t7 = cos(pkin(6));
t18 = t7 * t8;
t17 = -m(4) - m(5);
t15 = t10 * t6;
t14 = t10 * t7;
t11 = cos(qJ(2));
t5 = qJ(2) + pkin(7);
t4 = t11 * pkin(2);
t3 = cos(t5);
t2 = sin(t5);
t1 = [(-m(2) - m(3) + t17) * g(3), (-m(3) * (rSges(3,1) * t11 - t9 * rSges(3,2)) - m(4) * (-rSges(4,2) * t2 + t4) - m(5) * (t22 * t2 + t4) - t21 * t3) * g(3) + (g(1) * t7 + g(2) * t6) * (-m(3) * (-rSges(3,1) * t9 - rSges(3,2) * t11) - m(4) * (-rSges(4,2) * t3 - t20) - m(5) * (t22 * t3 - t20) + t21 * t2), t17 * (g(1) * t6 - g(2) * t7), -m(5) * (g(1) * ((-t3 * t18 + t15) * rSges(5,1) + (-t3 * t14 - t19) * rSges(5,2)) + g(2) * ((-t3 * t19 - t14) * rSges(5,1) + (-t3 * t15 + t18) * rSges(5,2)) + g(3) * (-rSges(5,1) * t8 - rSges(5,2) * t10) * t2)];
taug = t1(:);
