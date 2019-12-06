% Calculate Gravitation load on the joints for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:43
% EndTime: 2019-12-05 16:08:45
% DurationCPUTime: 0.47s
% Computational Cost: add. (184->71), mult. (285->102), div. (0->0), fcn. (267->8), ass. (0->37)
t14 = sin(pkin(7));
t15 = cos(pkin(7));
t19 = cos(qJ(3));
t17 = sin(qJ(3));
t20 = cos(qJ(2));
t35 = t17 * t20;
t48 = t14 * t19 - t15 * t35;
t40 = rSges(6,1) + pkin(4);
t47 = g(1) * t15 + g(2) * t14;
t31 = rSges(6,3) + qJ(5);
t13 = qJ(3) + pkin(8);
t11 = sin(t13);
t12 = cos(t13);
t46 = t31 * t11 + t40 * t12;
t45 = -m(5) - m(6);
t44 = pkin(3) * t17;
t18 = sin(qJ(2));
t41 = g(3) * t18;
t39 = rSges(4,3) + pkin(6);
t37 = t14 * t20;
t36 = t15 * t20;
t34 = t19 * t20;
t16 = -qJ(4) - pkin(6);
t33 = rSges(6,2) - t16;
t32 = rSges(5,3) - t16;
t28 = rSges(5,1) * t12 - rSges(5,2) * t11;
t27 = t48 * pkin(3);
t26 = t19 * rSges(4,1) - t17 * rSges(4,2) + pkin(2);
t25 = -t14 * t35 - t15 * t19;
t24 = t25 * pkin(3);
t10 = t19 * pkin(3) + pkin(2);
t6 = t20 * t10;
t4 = t14 * t11 + t12 * t36;
t3 = t11 * t36 - t14 * t12;
t2 = -t15 * t11 + t12 * t37;
t1 = t11 * t37 + t15 * t12;
t5 = [(-m(2) - m(3) - m(4) + t45) * g(3), -m(3) * (g(3) * (t20 * rSges(3,1) - t18 * rSges(3,2)) + t47 * (-rSges(3,1) * t18 - rSges(3,2) * t20)) - m(4) * (g(3) * (t39 * t18 + t26 * t20) + t47 * (-t26 * t18 + t39 * t20)) - m(5) * (g(3) * (t32 * t18 + t28 * t20 + t6) + t47 * (t32 * t20 + (-t10 - t28) * t18)) - m(6) * (g(3) * t6 + (g(3) * t46 + t47 * t33) * t20 + (g(3) * t33 + t47 * (-t10 - t46)) * t18), -m(4) * (g(1) * (t48 * rSges(4,1) + (-t14 * t17 - t15 * t34) * rSges(4,2)) + g(2) * (t25 * rSges(4,1) + (-t14 * t34 + t15 * t17) * rSges(4,2))) - m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2) + t27) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2) + t24)) - m(6) * (g(1) * (-t40 * t3 + t31 * t4 + t27) + g(2) * (-t40 * t1 + t31 * t2 + t24)) + (-m(4) * (-rSges(4,1) * t17 - rSges(4,2) * t19) - m(5) * (-rSges(5,1) * t11 - rSges(5,2) * t12 - t44) - m(6) * (-t40 * t11 + t31 * t12 - t44)) * t41, t45 * (-g(3) * t20 + t47 * t18), -m(6) * (g(1) * t3 + g(2) * t1 + t11 * t41)];
taug = t5(:);
