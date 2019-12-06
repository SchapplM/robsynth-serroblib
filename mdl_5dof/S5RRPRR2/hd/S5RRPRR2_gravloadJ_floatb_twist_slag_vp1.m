% Calculate Gravitation load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:04
% EndTime: 2019-12-05 18:27:06
% DurationCPUTime: 0.32s
% Computational Cost: add. (267->67), mult. (223->83), div. (0->0), fcn. (177->10), ass. (0->38)
t22 = qJ(2) + pkin(9);
t18 = qJ(4) + t22;
t13 = cos(t18);
t14 = qJ(5) + t18;
t10 = cos(t14);
t9 = sin(t14);
t43 = t10 * rSges(6,1) - t9 * rSges(6,2);
t41 = pkin(4) * t13 + t43;
t12 = sin(t18);
t42 = t13 * rSges(5,1) - t12 * rSges(5,2);
t16 = sin(t22);
t17 = cos(t22);
t26 = cos(qJ(2));
t20 = t26 * pkin(2);
t57 = rSges(4,1) * t17 - rSges(4,2) * t16 + t20;
t40 = -rSges(6,1) * t9 - rSges(6,2) * t10;
t56 = -pkin(4) * t12 + t40;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t55 = g(1) * t27 + g(2) * t25;
t24 = sin(qJ(2));
t51 = t24 * pkin(2);
t49 = rSges(3,3) + pkin(6);
t23 = -qJ(3) - pkin(6);
t47 = rSges(4,3) - t23;
t21 = -pkin(7) + t23;
t46 = rSges(5,3) - t21;
t45 = rSges(6,3) + pkin(8) - t21;
t44 = pkin(3) * t17 + t20;
t3 = pkin(1) + t44;
t4 = -pkin(3) * t16 - t51;
t39 = rSges(3,1) * t26 - rSges(3,2) * t24;
t36 = -rSges(5,1) * t12 - rSges(5,2) * t13;
t35 = t3 + t41;
t33 = pkin(1) + t39;
t32 = t3 + t42;
t30 = pkin(1) + t57;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t25 - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - rSges(2,2) * t25)) - m(3) * ((g(1) * t49 + g(2) * t33) * t27 + (-g(1) * t33 + g(2) * t49) * t25) - m(4) * ((g(1) * t47 + g(2) * t30) * t27 + (-g(1) * t30 + g(2) * t47) * t25) - m(5) * ((g(1) * t46 + g(2) * t32) * t27 + (-g(1) * t32 + g(2) * t46) * t25) - m(6) * ((g(1) * t45 + g(2) * t35) * t27 + (-g(1) * t35 + g(2) * t45) * t25), -m(3) * (g(3) * t39 + t55 * (-rSges(3,1) * t24 - rSges(3,2) * t26)) - m(4) * (g(3) * t57 + t55 * (-rSges(4,1) * t16 - rSges(4,2) * t17 - t51)) - m(5) * (g(3) * (t42 + t44) + t55 * (t36 + t4)) - m(6) * (g(3) * (t41 + t44) + t55 * (t4 + t56)), (-m(4) - m(5) - m(6)) * (g(1) * t25 - g(2) * t27), (-m(5) * t42 - m(6) * t41) * g(3) + t55 * (-m(5) * t36 - m(6) * t56), -m(6) * (g(3) * t43 + t40 * t55)];
taug = t1(:);
