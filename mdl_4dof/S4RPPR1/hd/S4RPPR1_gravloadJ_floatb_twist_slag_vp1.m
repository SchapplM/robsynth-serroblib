% Calculate Gravitation load on the joints for
% S4RPPR1
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:25
% EndTime: 2019-03-08 18:27:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (86->36), mult. (76->46), div. (0->0), fcn. (64->6), ass. (0->17)
t23 = -m(4) - m(5);
t14 = sin(qJ(1));
t22 = pkin(1) * t14;
t21 = cos(qJ(4));
t20 = sin(qJ(4));
t13 = qJ(1) + pkin(6);
t10 = sin(t13);
t11 = cos(t13);
t15 = cos(qJ(1));
t12 = t15 * pkin(1);
t19 = t11 * pkin(2) + t10 * qJ(3) + t12;
t18 = t11 * qJ(3) - t22;
t1 = -t10 * t20 - t11 * t21;
t2 = -t10 * t21 + t11 * t20;
t17 = -rSges(5,1) * t2 + rSges(5,2) * t1;
t16 = rSges(5,1) * t1 + rSges(5,2) * t2;
t3 = [-m(2) * (g(1) * (-t14 * rSges(2,1) - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - t14 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t10 - rSges(3,2) * t11 - t22) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t10 + t12)) - m(4) * (g(1) * (rSges(4,3) * t11 + (-rSges(4,1) - pkin(2)) * t10 + t18) + g(2) * (rSges(4,1) * t11 + rSges(4,3) * t10 + t19)) - m(5) * (g(1) * ((-pkin(2) - pkin(3)) * t10 - t17 + t18) + g(2) * (pkin(3) * t11 - t16 + t19)) (-m(3) + t23) * g(3), t23 * (g(1) * t10 - g(2) * t11) -m(5) * (g(1) * t17 + g(2) * t16)];
taug  = t3(:);
