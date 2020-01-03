% Calculate Gravitation load on the joints for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:43
% DurationCPUTime: 0.14s
% Computational Cost: add. (67->38), mult. (115->52), div. (0->0), fcn. (112->6), ass. (0->17)
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t8 = t14 * qJ(2);
t24 = t8 + (-pkin(1) - pkin(2)) * t12;
t23 = m(4) + m(5);
t22 = rSges(5,3) + pkin(5);
t21 = t14 * pkin(1) + t12 * qJ(2);
t20 = cos(pkin(6));
t19 = sin(pkin(6));
t18 = t14 * pkin(2) + t21;
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t17 = -rSges(5,1) * t13 + rSges(5,2) * t11;
t15 = pkin(3) - t17;
t2 = -t12 * t20 + t14 * t19;
t1 = -t12 * t19 - t14 * t20;
t3 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t14 + t8 + (-rSges(3,1) - pkin(1)) * t12) + g(2) * (rSges(3,1) * t14 + t12 * rSges(3,3) + t21)) - m(4) * (g(1) * (rSges(4,1) * t2 - rSges(4,2) * t1 + t24) + g(2) * (-rSges(4,1) * t1 - rSges(4,2) * t2 + t18)) - m(5) * (g(1) * t24 + g(2) * t18 + (g(1) * t15 + g(2) * t22) * t2 + (g(1) * t22 - g(2) * t15) * t1), (-m(3) - t23) * (g(1) * t12 - g(2) * t14), t23 * g(3), -m(5) * (g(3) * t17 + (g(1) * t1 + g(2) * t2) * (rSges(5,1) * t11 + rSges(5,2) * t13))];
taug = t3(:);
