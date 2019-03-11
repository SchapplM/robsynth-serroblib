% Calculate Gravitation load on the joints for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:42
% EndTime: 2019-03-08 19:46:43
% DurationCPUTime: 0.61s
% Computational Cost: add. (366->100), mult. (801->145), div. (0->0), fcn. (929->12), ass. (0->53)
t77 = -rSges(5,3) - pkin(8);
t28 = sin(pkin(10));
t33 = sin(qJ(2));
t35 = cos(qJ(2));
t54 = cos(pkin(10));
t55 = cos(pkin(6));
t47 = t55 * t54;
t12 = t28 * t35 + t33 * t47;
t52 = t28 * t55;
t14 = -t33 * t52 + t54 * t35;
t76 = g(1) * t14 + g(2) * t12;
t11 = t28 * t33 - t35 * t47;
t13 = t54 * t33 + t35 * t52;
t75 = g(1) * t13 + g(2) * t11;
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t29 = sin(pkin(6));
t60 = t29 * t35;
t15 = t55 * t32 + t34 * t60;
t62 = t28 * t29;
t3 = -t13 * t34 + t32 * t62;
t51 = t29 * t54;
t5 = t11 * t34 + t32 * t51;
t74 = -g(1) * t3 + g(2) * t5 - g(3) * t15;
t16 = -t32 * t60 + t55 * t34;
t4 = t13 * t32 + t34 * t62;
t6 = -t11 * t32 + t34 * t51;
t73 = g(1) * t4 - g(2) * t6 + g(3) * t16;
t70 = -m(6) - m(7);
t63 = g(3) * t29;
t61 = t29 * t33;
t59 = rSges(7,3) + pkin(9) + qJ(5);
t58 = pkin(2) * t60 + qJ(3) * t61;
t57 = rSges(4,3) + qJ(3);
t56 = rSges(6,3) + qJ(5);
t53 = -m(4) - m(5) + t70;
t50 = g(3) * (pkin(8) * t60 + t58);
t49 = rSges(5,1) * t32 + rSges(5,2) * t34;
t27 = sin(pkin(11));
t30 = cos(pkin(11));
t48 = rSges(6,1) * t27 + rSges(6,2) * t30;
t46 = rSges(6,1) * t30 - rSges(6,2) * t27 + pkin(4);
t26 = pkin(11) + qJ(6);
t24 = sin(t26);
t25 = cos(t26);
t44 = rSges(7,1) * t25 - rSges(7,2) * t24 + pkin(5) * t30 + pkin(4);
t42 = rSges(7,1) * t24 + rSges(7,2) * t25 + pkin(5) * t27;
t10 = t13 * pkin(2);
t9 = t11 * pkin(2);
t40 = -g(1) * t10 - g(2) * t9 + t50;
t39 = t46 * t32 - t56 * t34;
t38 = t44 * t32 - t59 * t34;
t1 = [(-m(2) - m(3) + t53) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t13 - rSges(3,2) * t14) + g(2) * (-rSges(3,1) * t11 - rSges(3,2) * t12) + (rSges(3,1) * t35 - rSges(3,2) * t33) * t63) - m(4) * (g(1) * (rSges(4,2) * t13 + t57 * t14 - t10) + g(2) * (rSges(4,2) * t11 + t57 * t12 - t9) + g(3) * ((-rSges(4,2) * t35 + rSges(4,3) * t33) * t29 + t58)) - m(5) * (g(1) * (t77 * t13 - t10) + g(2) * (t77 * t11 - t9) + t50 + (rSges(5,3) * t35 + t49 * t33) * t63 + t76 * (qJ(3) + t49)) - m(6) * ((t39 * t33 + t48 * t35) * t63 + t40 + t75 * (-pkin(8) - t48) + t76 * (qJ(3) + t39)) - m(7) * ((t38 * t33 + t42 * t35) * t63 + t40 + t75 * (-pkin(8) - t42) + t76 * (qJ(3) + t38)) t53 * (-g(3) * t60 + t75) -m(5) * (g(1) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (rSges(5,1) * t5 + rSges(5,2) * t6) + g(3) * (-rSges(5,1) * t15 - rSges(5,2) * t16)) - m(6) * (t74 * t46 + t73 * t56) - m(7) * (t74 * t44 + t73 * t59) -t70 * t74, -m(7) * (g(1) * ((t14 * t25 - t24 * t4) * rSges(7,1) + (-t14 * t24 - t25 * t4) * rSges(7,2)) + g(2) * ((t12 * t25 + t24 * t6) * rSges(7,1) + (-t12 * t24 + t25 * t6) * rSges(7,2)) + g(3) * ((-t16 * t24 + t25 * t61) * rSges(7,1) + (-t16 * t25 - t24 * t61) * rSges(7,2)))];
taug  = t1(:);
