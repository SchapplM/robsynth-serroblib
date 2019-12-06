% Calculate Gravitation load on the joints for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:03
% EndTime: 2019-12-05 17:19:07
% DurationCPUTime: 0.66s
% Computational Cost: add. (372->100), mult. (810->159), div. (0->0), fcn. (956->12), ass. (0->53)
t28 = sin(pkin(10));
t32 = sin(qJ(2));
t35 = cos(qJ(2));
t55 = cos(pkin(10));
t56 = cos(pkin(5));
t47 = t56 * t55;
t14 = t28 * t35 + t32 * t47;
t53 = t28 * t56;
t16 = -t32 * t53 + t55 * t35;
t77 = g(1) * t16 + g(2) * t14;
t13 = t28 * t32 - t35 * t47;
t15 = t55 * t32 + t35 * t53;
t76 = g(1) * t15 + g(2) * t13;
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t29 = sin(pkin(5));
t60 = t29 * t32;
t17 = -t31 * t60 + t56 * t34;
t52 = t29 * t55;
t7 = -t14 * t31 - t34 * t52;
t59 = t29 * t34;
t9 = -t16 * t31 + t28 * t59;
t75 = g(1) * t9 + g(2) * t7 + g(3) * t17;
t10 = t28 * t29 * t31 + t16 * t34;
t18 = t56 * t31 + t32 * t59;
t8 = t14 * t34 - t31 * t52;
t74 = g(1) * t10 + g(2) * t8 + g(3) * t18;
t27 = qJ(4) + qJ(5);
t25 = sin(t27);
t26 = cos(t27);
t33 = cos(qJ(4));
t43 = rSges(6,1) * t26 - rSges(6,2) * t25 + pkin(4) * t33 + pkin(3);
t57 = rSges(6,3) + pkin(9) + pkin(8);
t73 = t57 * t31 + t43 * t34;
t30 = sin(qJ(4));
t46 = rSges(5,1) * t33 - rSges(5,2) * t30 + pkin(3);
t61 = rSges(5,3) + pkin(8);
t72 = t61 * t31 + t46 * t34;
t63 = g(3) * t29;
t62 = rSges(4,3) + pkin(7);
t58 = t29 * t35;
t54 = g(3) * (pkin(2) * t58 + pkin(7) * t60);
t51 = t13 * t33 - t30 * t8;
t50 = rSges(4,1) * t34 - rSges(4,2) * t31;
t49 = rSges(5,1) * t30 + rSges(5,2) * t33;
t48 = -t10 * t30 + t15 * t33;
t44 = -t18 * t30 - t33 * t58;
t42 = rSges(6,1) * t25 + rSges(6,2) * t26 + pkin(4) * t30;
t11 = t13 * pkin(2);
t12 = t15 * pkin(2);
t40 = -g(1) * t12 - g(2) * t11 + t54;
t38 = m(6) * (g(1) * ((-t10 * t25 + t15 * t26) * rSges(6,1) + (-t10 * t26 - t15 * t25) * rSges(6,2)) + g(2) * ((t13 * t26 - t25 * t8) * rSges(6,1) + (-t13 * t25 - t26 * t8) * rSges(6,2)) + g(3) * ((-t18 * t25 - t26 * t58) * rSges(6,1) + (-t18 * t26 + t25 * t58) * rSges(6,2)));
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t15 - rSges(3,2) * t16) + g(2) * (-rSges(3,1) * t13 - rSges(3,2) * t14) + (rSges(3,1) * t35 - rSges(3,2) * t32) * t63) - m(4) * (g(1) * (-t50 * t15 + t62 * t16 - t12) + g(2) * (-t50 * t13 + t62 * t14 - t11) + t54 + (rSges(4,3) * t32 + t50 * t35) * t63) - m(5) * ((t49 * t32 + t72 * t35) * t63 + t40 + t77 * (pkin(7) + t49) - t76 * t72) - m(6) * ((t42 * t32 + t73 * t35) * t63 + t40 + t77 * (pkin(7) + t42) - t76 * t73), -m(4) * (g(1) * (rSges(4,1) * t9 - rSges(4,2) * t10) + g(2) * (rSges(4,1) * t7 - rSges(4,2) * t8) + g(3) * (rSges(4,1) * t17 - rSges(4,2) * t18)) - m(5) * (t75 * t46 + t74 * t61) - m(6) * (t75 * t43 + t74 * t57), -m(5) * (g(1) * (t48 * rSges(5,1) + (-t10 * t33 - t15 * t30) * rSges(5,2)) + g(2) * (t51 * rSges(5,1) + (-t13 * t30 - t33 * t8) * rSges(5,2)) + g(3) * (t44 * rSges(5,1) + (-t18 * t33 + t30 * t58) * rSges(5,2))) - t38 - m(6) * (g(1) * t48 + g(2) * t51 + g(3) * t44) * pkin(4), -t38];
taug = t1(:);
