% Calculate Gravitation load on the joints for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:25
% EndTime: 2019-03-09 09:28:27
% DurationCPUTime: 0.72s
% Computational Cost: add. (427->148), mult. (982->206), div. (0->0), fcn. (1134->10), ass. (0->60)
t39 = sin(qJ(2));
t40 = sin(qJ(1));
t43 = cos(qJ(2));
t44 = cos(qJ(1));
t67 = cos(pkin(6));
t60 = t44 * t67;
t21 = t39 * t60 + t40 * t43;
t61 = t40 * t67;
t23 = -t39 * t61 + t43 * t44;
t83 = g(1) * t23;
t87 = -g(2) * t21 - t83;
t20 = t39 * t40 - t43 * t60;
t22 = t44 * t39 + t43 * t61;
t86 = g(1) * t22 + g(2) * t20;
t38 = sin(qJ(5));
t42 = cos(qJ(5));
t37 = sin(qJ(6));
t41 = cos(qJ(6));
t54 = t41 * rSges(7,1) - t37 * rSges(7,2) + pkin(5);
t79 = rSges(7,3) + pkin(10);
t85 = t54 * t38 - t79 * t42;
t36 = sin(pkin(6));
t80 = g(3) * t36;
t78 = t36 * t39;
t77 = t36 * t40;
t76 = t36 * t42;
t75 = t36 * t43;
t74 = t36 * t44;
t73 = pkin(9) - qJ(3);
t72 = pkin(2) * t75 + qJ(3) * t78;
t71 = t44 * pkin(1) + pkin(8) * t77;
t70 = rSges(5,2) + qJ(3);
t69 = rSges(4,3) + qJ(3);
t68 = rSges(5,3) + qJ(4);
t66 = -m(5) - m(6) - m(7);
t65 = rSges(6,3) + t73;
t64 = t23 * pkin(2) + t71;
t63 = qJ(4) * t75 + t72;
t62 = -t40 * pkin(1) + pkin(8) * t74;
t59 = pkin(3) * t77 + t64;
t58 = g(2) * t65;
t57 = -t21 * pkin(2) + t62;
t56 = rSges(6,1) * t38 + rSges(6,2) * t42;
t55 = pkin(3) * t74 + t57;
t53 = -t37 * rSges(7,1) - t41 * rSges(7,2) - pkin(9);
t52 = -t21 * t38 + t42 * t74;
t51 = t21 * t42 + t38 * t74;
t49 = -qJ(3) - t53;
t48 = pkin(4) * t77 + qJ(4) * t23 + t59;
t47 = pkin(4) * t74 - t21 * qJ(4) + t55;
t14 = t20 * pkin(2);
t16 = t22 * pkin(2);
t46 = -g(1) * t16 - g(2) * t14 + g(3) * t63;
t19 = t38 * t78 + t67 * t42;
t18 = -t67 * t38 + t39 * t76;
t6 = t23 * t38 + t40 * t76;
t5 = -t23 * t42 + t38 * t77;
t2 = -t22 * t37 + t41 * t6;
t1 = -t22 * t41 - t37 * t6;
t3 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - t40 * rSges(2,2))) - m(3) * (g(1) * (-t21 * rSges(3,1) + t20 * rSges(3,2) + rSges(3,3) * t74 + t62) + g(2) * (rSges(3,1) * t23 - rSges(3,2) * t22 + rSges(3,3) * t77 + t71)) - m(4) * (g(1) * (rSges(4,1) * t74 + t21 * rSges(4,2) - t69 * t20 + t57) + g(2) * (rSges(4,1) * t77 - t23 * rSges(4,2) + t69 * t22 + t64)) - m(5) * (g(1) * (rSges(5,1) * t74 - t70 * t20 - t68 * t21 + t55) + g(2) * (rSges(5,1) * t77 + t70 * t22 + t68 * t23 + t59)) - m(6) * (g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t48) - t22 * t58 + (rSges(6,1) * t52 - rSges(6,2) * t51 + t65 * t20 + t47) * g(1)) - m(7) * (g(1) * (t49 * t20 + t51 * t79 + t52 * t54 + t47) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t6 - t73 * t22 + t79 * t5 + t48)) -m(3) * (g(1) * (-rSges(3,1) * t22 - rSges(3,2) * t23) + g(2) * (-rSges(3,1) * t20 - rSges(3,2) * t21) + (rSges(3,1) * t43 - rSges(3,2) * t39) * t80) - m(4) * (g(1) * (rSges(4,2) * t22 + t69 * t23 - t16) + g(2) * (rSges(4,2) * t20 + t69 * t21 - t14) + g(3) * ((-rSges(4,2) * t43 + rSges(4,3) * t39) * t36 + t72)) - m(5) * (g(1) * (-t68 * t22 + t70 * t23 - t16) + g(2) * (-t68 * t20 + t70 * t21 - t14) + g(3) * ((rSges(5,2) * t39 + rSges(5,3) * t43) * t36 + t63)) - m(6) * (-t65 * t83 - t21 * t58 + (t56 * t43 + (-rSges(6,3) - pkin(9)) * t39) * t80 + t46 + t86 * (-qJ(4) - t56)) - m(7) * ((t53 * t39 + t85 * t43) * t80 + t46 + t87 * t49 + t86 * (-qJ(4) - t85)) (-m(4) + t66) * (-g(3) * t75 + t86) t66 * (g(3) * t78 - t87) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (rSges(6,1) * t51 + rSges(6,2) * t52) + g(3) * (rSges(6,1) * t18 - rSges(6,2) * t19)) - m(7) * (g(2) * (t51 * t54 - t52 * t79) + (t54 * t18 + t79 * t19) * g(3) + (-t54 * t5 + t79 * t6) * g(1)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * ((-t20 * t41 + t37 * t52) * rSges(7,1) + (t20 * t37 + t41 * t52) * rSges(7,2)) + g(3) * ((-t19 * t37 + t41 * t75) * rSges(7,1) + (-t19 * t41 - t37 * t75) * rSges(7,2)))];
taug  = t3(:);
