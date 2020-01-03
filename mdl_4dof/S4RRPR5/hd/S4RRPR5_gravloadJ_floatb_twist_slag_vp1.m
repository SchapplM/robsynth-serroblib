% Calculate Gravitation load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (121->42), mult. (107->52), div. (0->0), fcn. (80->6), ass. (0->23)
t39 = rSges(5,3) + pkin(6);
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t38 = rSges(5,1) * t19 + rSges(5,2) * t21;
t18 = qJ(1) + qJ(2);
t16 = cos(t18);
t15 = sin(t18);
t35 = g(1) * t15;
t37 = -g(2) * t16 + t35;
t20 = sin(qJ(1));
t36 = pkin(1) * t20;
t33 = t16 * pkin(2) + t15 * qJ(3);
t7 = t16 * qJ(3);
t30 = t38 * t16 + t7;
t29 = t16 * rSges(3,1) - rSges(3,2) * t15;
t28 = (-pkin(2) - t39) * t35;
t27 = t16 * rSges(4,3) + t7 + (rSges(4,2) - pkin(2)) * t15;
t26 = -rSges(3,1) * t15 - rSges(3,2) * t16;
t24 = -rSges(4,2) * t16 + t15 * rSges(4,3) + t33;
t23 = t38 * t15 + t39 * t16 + t33;
t22 = cos(qJ(1));
t17 = t22 * pkin(1);
t1 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - rSges(2,2) * t22) + g(2) * (rSges(2,1) * t22 - t20 * rSges(2,2))) - m(3) * (g(1) * (t26 - t36) + g(2) * (t17 + t29)) - m(4) * (g(1) * (t27 - t36) + g(2) * (t17 + t24)) - m(5) * (g(1) * (t30 - t36) + g(2) * (t17 + t23) + t28), -m(3) * (g(1) * t26 + g(2) * t29) - m(4) * (g(1) * t27 + g(2) * t24) - m(5) * (g(1) * t30 + g(2) * t23 + t28), (-m(4) - m(5)) * t37, -m(5) * (-g(3) * t38 + t37 * (rSges(5,1) * t21 - rSges(5,2) * t19))];
taug = t1(:);
