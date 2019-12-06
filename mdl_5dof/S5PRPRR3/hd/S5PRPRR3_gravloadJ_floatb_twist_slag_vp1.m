% Calculate Gravitation load on the joints for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:37
% EndTime: 2019-12-05 15:46:38
% DurationCPUTime: 0.31s
% Computational Cost: add. (182->60), mult. (214->86), div. (0->0), fcn. (190->10), ass. (0->35)
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t47 = g(1) * t15 + g(2) * t14;
t12 = qJ(2) + pkin(9);
t7 = sin(t12);
t46 = g(3) * t7;
t13 = qJ(4) + qJ(5);
t10 = cos(t13);
t34 = t15 * t10;
t37 = t10 * t14;
t9 = sin(t13);
t39 = t15 * t9;
t40 = t14 * t9;
t8 = cos(t12);
t45 = (-t8 * t40 - t34) * rSges(6,1) + (-t8 * t37 + t39) * rSges(6,2);
t44 = (-t8 * t39 + t37) * rSges(6,1) + (-t8 * t34 - t40) * rSges(6,2);
t17 = sin(qJ(2));
t43 = pkin(2) * t17;
t38 = rSges(5,3) + pkin(6);
t16 = sin(qJ(4));
t36 = t14 * t16;
t18 = cos(qJ(4));
t35 = t14 * t18;
t33 = t15 * t16;
t32 = t15 * t18;
t31 = rSges(6,3) + pkin(7) + pkin(6);
t30 = -m(4) - m(5) - m(6);
t29 = -rSges(6,1) * t9 - rSges(6,2) * t10;
t27 = rSges(6,1) * t10 - rSges(6,2) * t9 + t18 * pkin(4) + pkin(3);
t26 = -t8 * t33 + t35;
t25 = -t8 * t36 - t32;
t24 = t18 * rSges(5,1) - t16 * rSges(5,2) + pkin(3);
t19 = cos(qJ(2));
t11 = t19 * pkin(2);
t1 = [(-m(2) - m(3) + t30) * g(3), -m(3) * (g(3) * (t19 * rSges(3,1) - t17 * rSges(3,2)) + t47 * (-rSges(3,1) * t17 - rSges(3,2) * t19)) - m(4) * (g(3) * (rSges(4,1) * t8 - rSges(4,2) * t7 + t11) + t47 * (-rSges(4,1) * t7 - rSges(4,2) * t8 - t43)) - m(5) * (g(3) * (t24 * t8 + t38 * t7 + t11) + t47 * (-t24 * t7 + t38 * t8 - t43)) - m(6) * (g(3) * (t27 * t8 + t31 * t7 + t11) + t47 * (-t27 * t7 + t31 * t8 - t43)), t30 * (g(1) * t14 - g(2) * t15), -m(5) * (g(1) * (t26 * rSges(5,1) + (-t8 * t32 - t36) * rSges(5,2)) + g(2) * (t25 * rSges(5,1) + (-t8 * t35 + t33) * rSges(5,2))) - m(6) * (g(1) * (t26 * pkin(4) + t44) + g(2) * (t25 * pkin(4) + t45)) + (-m(5) * (-rSges(5,1) * t16 - rSges(5,2) * t18) - m(6) * (-pkin(4) * t16 + t29)) * t46, -m(6) * (g(1) * t44 + g(2) * t45 + t29 * t46)];
taug = t1(:);
