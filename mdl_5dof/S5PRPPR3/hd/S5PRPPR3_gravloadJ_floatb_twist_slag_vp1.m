% Calculate Gravitation load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:18
% EndTime: 2019-12-05 15:26:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (122->47), mult. (157->62), div. (0->0), fcn. (131->8), ass. (0->24)
t11 = sin(pkin(7));
t12 = cos(pkin(7));
t36 = g(1) * t12 + g(2) * t11;
t10 = qJ(2) + pkin(8);
t8 = cos(t10);
t35 = g(3) * t8;
t34 = -m(5) - m(6);
t14 = sin(qJ(2));
t33 = pkin(2) * t14;
t30 = rSges(6,3) + pkin(6);
t13 = sin(qJ(5));
t28 = t11 * t13;
t15 = cos(qJ(5));
t27 = t11 * t15;
t26 = t12 * t13;
t25 = t12 * t15;
t24 = -m(4) + t34;
t7 = sin(t10);
t16 = cos(qJ(2));
t9 = t16 * pkin(2);
t23 = t8 * pkin(3) + t7 * qJ(4) + t9;
t22 = t36 * qJ(4) * t8;
t20 = t13 * rSges(6,1) + t15 * rSges(6,2);
t1 = [(-m(2) - m(3) + t24) * g(3), -m(3) * (g(3) * (rSges(3,1) * t16 - t14 * rSges(3,2)) + t36 * (-rSges(3,1) * t14 - rSges(3,2) * t16)) - m(4) * (g(3) * (rSges(4,1) * t8 - rSges(4,2) * t7 + t9) + t36 * (-rSges(4,1) * t7 - rSges(4,2) * t8 - t33)) - m(5) * (g(3) * (-rSges(5,2) * t8 + rSges(5,3) * t7 + t23) + t22 + t36 * (rSges(5,3) * t8 - t33 + (rSges(5,2) - pkin(3)) * t7)) - m(6) * (g(3) * (t20 * t7 + t30 * t8 + t23) + t22 + t36 * (-t33 + t20 * t8 + (-pkin(3) - t30) * t7)), t24 * (g(1) * t11 - g(2) * t12), t34 * (t36 * t7 - t35), -m(6) * (g(1) * ((t7 * t25 - t28) * rSges(6,1) + (-t7 * t26 - t27) * rSges(6,2)) + g(2) * ((t7 * t27 + t26) * rSges(6,1) + (-t7 * t28 + t25) * rSges(6,2)) + (-rSges(6,1) * t15 + rSges(6,2) * t13) * t35)];
taug = t1(:);
