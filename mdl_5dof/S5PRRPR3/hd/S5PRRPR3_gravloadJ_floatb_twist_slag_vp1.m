% Calculate Gravitation load on the joints for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:02
% EndTime: 2019-12-05 16:19:03
% DurationCPUTime: 0.22s
% Computational Cost: add. (189->50), mult. (144->61), div. (0->0), fcn. (109->8), ass. (0->28)
t17 = qJ(3) + pkin(9);
t12 = cos(t17);
t20 = cos(qJ(3));
t14 = t20 * pkin(3);
t13 = qJ(5) + t17;
t6 = sin(t13);
t7 = cos(t13);
t28 = t7 * rSges(6,1) - rSges(6,2) * t6;
t40 = pkin(4) * t12 + t14 + t28;
t10 = sin(t17);
t39 = rSges(5,1) * t12 - rSges(5,2) * t10 + t14;
t16 = pkin(8) + qJ(2);
t11 = cos(t16);
t9 = sin(t16);
t38 = g(1) * t11 + g(2) * t9;
t36 = -m(5) - m(6);
t19 = sin(qJ(3));
t34 = pkin(3) * t19;
t32 = rSges(4,3) + pkin(6);
t18 = -qJ(4) - pkin(6);
t30 = rSges(5,3) - t18;
t29 = rSges(6,3) + pkin(7) - t18;
t27 = -rSges(6,1) * t6 - rSges(6,2) * t7;
t26 = rSges(4,1) * t20 - t19 * rSges(4,2);
t24 = pkin(2) + t40;
t23 = pkin(2) + t26;
t22 = pkin(2) + t39;
t1 = [(-m(2) - m(3) - m(4) + t36) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t9 - rSges(3,2) * t11) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t9)) - m(4) * ((-g(1) * t23 + g(2) * t32) * t9 + (g(1) * t32 + g(2) * t23) * t11) - m(5) * ((-g(1) * t22 + g(2) * t30) * t9 + (g(1) * t30 + g(2) * t22) * t11) - m(6) * ((-g(1) * t24 + g(2) * t29) * t9 + (g(1) * t29 + g(2) * t24) * t11), (-m(4) * t26 - m(5) * t39 - m(6) * t40) * g(3) + t38 * (-m(4) * (-rSges(4,1) * t19 - rSges(4,2) * t20) - m(5) * (-rSges(5,1) * t10 - rSges(5,2) * t12 - t34) - m(6) * (-pkin(4) * t10 + t27 - t34)), t36 * (g(1) * t9 - g(2) * t11), -m(6) * (g(3) * t28 + t38 * t27)];
taug = t1(:);
