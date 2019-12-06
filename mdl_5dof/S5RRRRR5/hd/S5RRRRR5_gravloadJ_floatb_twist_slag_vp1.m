% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:10
% DurationCPUTime: 0.38s
% Computational Cost: add. (349->71), mult. (217->89), div. (0->0), fcn. (169->10), ass. (0->46)
t27 = cos(qJ(4));
t23 = qJ(4) + qJ(5);
t17 = sin(t23);
t19 = cos(t23);
t44 = t19 * rSges(6,1) - t17 * rSges(6,2);
t66 = t27 * pkin(4) + t44;
t65 = rSges(5,3) + pkin(8);
t64 = rSges(6,3) + pkin(9) + pkin(8);
t63 = -pkin(3) - t66;
t62 = rSges(6,1) * t17 + rSges(6,2) * t19;
t49 = t27 * rSges(5,1);
t61 = -pkin(3) - t49;
t24 = qJ(1) + qJ(2);
t21 = qJ(3) + t24;
t14 = sin(t21);
t60 = t62 * t14;
t18 = sin(t24);
t59 = pkin(2) * t18;
t20 = cos(t24);
t58 = pkin(2) * t20;
t25 = sin(qJ(4));
t57 = pkin(4) * t25;
t15 = cos(t21);
t56 = g(3) * t15;
t26 = sin(qJ(1));
t55 = t26 * pkin(1);
t28 = cos(qJ(1));
t54 = t28 * pkin(1);
t50 = t25 * rSges(5,2);
t48 = t14 * t50 + t65 * t15;
t47 = -t15 * rSges(4,1) + t14 * rSges(4,2);
t46 = -t20 * rSges(3,1) + t18 * rSges(3,2);
t43 = -t18 * rSges(3,1) - t20 * rSges(3,2);
t42 = -t14 * rSges(4,1) - t15 * rSges(4,2);
t41 = rSges(5,1) * t25 + rSges(5,2) * t27;
t39 = t48 - t59;
t38 = t47 - t58;
t37 = (t50 + t61) * t15;
t36 = t42 - t59;
t35 = t37 - t58;
t34 = -t64 * t14 + t63 * t15;
t33 = t63 * t14 + t64 * t15;
t32 = (-g(2) * t65 + g(3) * t61) * t14;
t31 = t33 - t59;
t30 = t34 - t58;
t1 = [-m(2) * (g(2) * (-t28 * rSges(2,1) + t26 * rSges(2,2)) + g(3) * (-t26 * rSges(2,1) - t28 * rSges(2,2))) - m(3) * (g(2) * (t46 - t54) + g(3) * (t43 - t55)) - m(4) * (g(2) * (t38 - t54) + g(3) * (t36 - t55)) - m(5) * (g(2) * (t35 - t54) + g(3) * (t39 - t55) + t32) - m(6) * (g(2) * (t30 - t54) + g(3) * (t31 - t55)), -m(3) * (g(2) * t46 + g(3) * t43) - m(4) * (g(2) * t38 + g(3) * t36) - m(5) * (g(2) * t35 + g(3) * t39 + t32) - m(6) * (g(2) * t30 + g(3) * t31), -m(4) * (g(2) * t47 + g(3) * t42) - m(5) * (g(2) * t37 + g(3) * t48 + t32) - m(6) * (g(2) * t34 + g(3) * t33), -m(5) * (g(1) * (t49 - t50) + g(2) * t41 * t14) - m(6) * (g(1) * t66 + g(2) * (t14 * t57 + t60)) + (m(5) * t41 - m(6) * (-t62 - t57)) * t56, -m(6) * (g(1) * t44 + g(2) * t60 - t56 * t62)];
taug = t1(:);
