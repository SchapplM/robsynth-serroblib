% Calculate Gravitation load on the joints for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:46
% EndTime: 2019-12-05 18:44:47
% DurationCPUTime: 0.39s
% Computational Cost: add. (292->70), mult. (258->86), div. (0->0), fcn. (206->8), ass. (0->40)
t61 = rSges(6,1) + pkin(4);
t20 = qJ(2) + qJ(3);
t17 = qJ(4) + t20;
t11 = sin(t17);
t12 = cos(t17);
t60 = -rSges(6,2) * t12 - t61 * t11;
t14 = sin(t20);
t56 = pkin(3) * t14;
t59 = -t56 + t60;
t43 = t12 * rSges(5,1) - t11 * rSges(5,2);
t15 = cos(t20);
t44 = t15 * rSges(4,1) - t14 * rSges(4,2);
t37 = -rSges(5,1) * t11 - rSges(5,2) * t12;
t58 = t37 - t56;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t57 = g(1) * t24 + g(2) * t22;
t42 = -t11 * rSges(6,2) + t61 * t12;
t25 = -pkin(7) - pkin(6);
t21 = sin(qJ(2));
t52 = t21 * pkin(2);
t51 = rSges(3,3) + pkin(6);
t23 = cos(qJ(2));
t18 = t23 * pkin(2);
t13 = t18 + pkin(1);
t47 = rSges(4,3) - t25;
t19 = -pkin(8) + t25;
t46 = rSges(5,3) - t19;
t45 = rSges(6,3) + qJ(5) - t19;
t10 = pkin(3) * t15;
t4 = t10 + t13;
t41 = t10 + t43;
t40 = t23 * rSges(3,1) - t21 * rSges(3,2);
t38 = -rSges(4,1) * t14 - rSges(4,2) * t15;
t35 = t10 + t42;
t34 = pkin(1) + t40;
t33 = t4 + t43;
t31 = t4 + t42;
t29 = t13 + t44;
t1 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t24 * rSges(2,2)) + g(2) * (t24 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * ((g(1) * t51 + g(2) * t34) * t24 + (-g(1) * t34 + g(2) * t51) * t22) - m(4) * ((g(1) * t47 + g(2) * t29) * t24 + (-g(1) * t29 + g(2) * t47) * t22) - m(5) * ((g(1) * t46 + g(2) * t33) * t24 + (-g(1) * t33 + g(2) * t46) * t22) - m(6) * ((g(1) * t45 + g(2) * t31) * t24 + (-g(1) * t31 + g(2) * t45) * t22), -m(3) * (g(3) * t40 + t57 * (-rSges(3,1) * t21 - rSges(3,2) * t23)) - m(4) * (g(3) * (t18 + t44) + t57 * (t38 - t52)) - m(5) * (g(3) * (t18 + t41) + t57 * (-t52 + t58)) - m(6) * (g(3) * (t18 + t35) + t57 * (-t52 + t59)), (-m(4) * t44 - m(5) * t41 - m(6) * t35) * g(3) + t57 * (-m(4) * t38 - m(5) * t58 - m(6) * t59), (-m(5) * t43 - m(6) * t42) * g(3) + t57 * (-m(5) * t37 - m(6) * t60), -m(6) * (g(1) * t22 - g(2) * t24)];
taug = t1(:);
