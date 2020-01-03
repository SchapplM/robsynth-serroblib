% Calculate Gravitation load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:56:48
% DurationCPUTime: 0.44s
% Computational Cost: add. (243->93), mult. (314->125), div. (0->0), fcn. (287->8), ass. (0->41)
t33 = rSges(6,3) + qJ(5) + pkin(7);
t12 = qJ(2) + pkin(8);
t9 = sin(t12);
t48 = t33 * t9;
t39 = rSges(5,3) + pkin(7);
t47 = t39 * t9;
t46 = rSges(6,1) + pkin(4);
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t45 = g(1) * t20 + g(2) * t17;
t16 = sin(qJ(2));
t44 = pkin(2) * t16;
t15 = sin(qJ(4));
t43 = pkin(4) * t15;
t40 = rSges(3,3) + pkin(6);
t38 = t17 * t15;
t18 = cos(qJ(4));
t37 = t17 * t18;
t36 = t20 * t15;
t35 = t20 * t18;
t14 = -qJ(3) - pkin(6);
t34 = rSges(4,3) - t14;
t32 = -t14 + t43;
t10 = cos(t12);
t31 = t10 * rSges(4,1) - t9 * rSges(4,2);
t19 = cos(qJ(2));
t30 = t19 * rSges(3,1) - t16 * rSges(3,2);
t28 = pkin(1) + t30;
t27 = rSges(5,1) * t18 - rSges(5,2) * t15 + pkin(3);
t7 = t18 * pkin(4) + pkin(3);
t26 = rSges(6,1) * t18 - rSges(6,2) * t15 + t7;
t3 = -t10 * t36 + t37;
t1 = t10 * t38 + t35;
t25 = t10 * pkin(3) + t47;
t24 = t10 * t7 + t48;
t11 = t19 * pkin(2);
t8 = t11 + pkin(1);
t6 = t20 * t8;
t4 = t10 * t35 + t38;
t2 = -t10 * t37 + t36;
t5 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - t20 * rSges(2,2)) + g(2) * (t20 * rSges(2,1) - t17 * rSges(2,2))) - m(3) * ((g(1) * t40 + g(2) * t28) * t20 + (-g(1) * t28 + g(2) * t40) * t17) - m(4) * (g(2) * t6 + (g(1) * t34 + g(2) * t31) * t20 + (g(1) * (-t31 - t8) + g(2) * t34) * t17) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t6) + (-g(1) * t14 + g(2) * t25) * t20 + (g(1) * (-t25 - t8) - g(2) * t14) * t17) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t6) + (g(1) * t32 + g(2) * t24) * t20 + (g(1) * (-t24 - t8) + g(2) * t32) * t17), -m(3) * (g(3) * t30 + t45 * (-rSges(3,1) * t16 - rSges(3,2) * t19)) - m(4) * (g(3) * (t11 + t31) + t45 * (-rSges(4,1) * t9 - rSges(4,2) * t10 - t44)) - m(5) * (g(3) * (t10 * t27 + t11 + t47) + t45 * (t10 * t39 - t27 * t9 - t44)) - m(6) * (g(3) * (t10 * t26 + t11 + t48) + t45 * (t10 * t33 - t26 * t9 - t44)), (-m(4) - m(5) - m(6)) * (g(1) * t17 - g(2) * t20), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) + t2 * rSges(5,2))) - m(6) * (g(1) * (-t4 * rSges(6,2) + t46 * t3) + g(2) * (t2 * rSges(6,2) - t46 * t1)) + (-m(5) * (-rSges(5,1) * t15 - rSges(5,2) * t18) - m(6) * (-rSges(6,1) * t15 - rSges(6,2) * t18 - t43)) * g(3) * t9, -m(6) * (-g(3) * t10 + t45 * t9)];
taug = t5(:);
