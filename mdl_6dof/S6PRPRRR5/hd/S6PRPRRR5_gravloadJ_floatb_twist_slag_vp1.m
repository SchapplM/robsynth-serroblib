% Calculate Gravitation load on the joints for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:56
% EndTime: 2019-03-08 20:40:58
% DurationCPUTime: 0.69s
% Computational Cost: add. (464->122), mult. (837->181), div. (0->0), fcn. (961->12), ass. (0->63)
t105 = rSges(7,3) + pkin(10);
t46 = sin(pkin(11));
t48 = cos(pkin(11));
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t76 = cos(pkin(6));
t69 = t54 * t76;
t32 = t46 * t51 - t48 * t69;
t53 = cos(qJ(4));
t47 = sin(pkin(6));
t50 = sin(qJ(4));
t84 = t47 * t50;
t104 = t32 * t53 + t48 * t84;
t34 = t46 * t69 + t48 * t51;
t103 = t34 * t53 - t46 * t84;
t49 = sin(qJ(6));
t52 = cos(qJ(6));
t102 = rSges(7,1) * t52 - rSges(7,2) * t49 + pkin(5);
t101 = rSges(5,3) + pkin(8);
t100 = -rSges(7,1) * t49 - rSges(7,2) * t52;
t70 = t51 * t76;
t33 = t46 * t54 + t48 * t70;
t35 = -t46 * t70 + t48 * t54;
t99 = g(1) * t35 + g(2) * t33;
t45 = qJ(4) + qJ(5);
t43 = sin(t45);
t44 = cos(t45);
t86 = t46 * t47;
t12 = t34 * t44 - t43 * t86;
t13 = t34 * t43 + t44 * t86;
t98 = t12 * rSges(6,1) - t13 * rSges(6,2);
t85 = t47 * t48;
t14 = t32 * t44 + t43 * t85;
t15 = -t32 * t43 + t44 * t85;
t97 = t14 * rSges(6,1) + t15 * rSges(6,2);
t96 = pkin(4) * t50;
t93 = g(3) * t47;
t83 = t47 * t51;
t82 = t47 * t53;
t81 = t47 * t54;
t26 = -t76 * t43 - t44 * t81;
t27 = -t43 * t81 + t76 * t44;
t80 = t26 * rSges(6,1) - t27 * rSges(6,2);
t79 = t104 * pkin(4);
t78 = pkin(2) * t81 + qJ(3) * t83;
t77 = rSges(4,3) + qJ(3);
t30 = t32 * pkin(2);
t55 = -pkin(9) - pkin(8);
t73 = t32 * t55 + t33 * t96 - t30;
t31 = t34 * pkin(2);
t72 = t34 * t55 + t35 * t96 - t31;
t71 = -m(4) - m(5) - m(6) - m(7);
t68 = g(3) * (t83 * t96 + t78);
t67 = rSges(5,1) * t50 + rSges(5,2) * t53;
t66 = rSges(6,1) * t43 + rSges(6,2) * t44;
t65 = t103 * pkin(4);
t62 = -t76 * t50 - t53 * t81;
t61 = t102 * t12 + t105 * t13;
t60 = t102 * t14 - t105 * t15;
t59 = t102 * t26 + t105 * t27;
t58 = t62 * pkin(4);
t57 = t102 * t43 - t105 * t44;
t1 = [(-m(2) - m(3) + t71) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + g(2) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + (rSges(3,1) * t54 - rSges(3,2) * t51) * t93) - m(4) * (g(1) * (rSges(4,2) * t34 + t77 * t35 - t31) + g(2) * (rSges(4,2) * t32 + t77 * t33 - t30) + g(3) * ((-rSges(4,2) * t54 + rSges(4,3) * t51) * t47 + t78)) - m(5) * (g(1) * (-t101 * t34 - t31) + g(2) * (-t101 * t32 - t30) + g(3) * t78 + (t101 * t54 + t67 * t51) * t93 + t99 * (qJ(3) + t67)) - m(6) * (g(1) * (-rSges(6,3) * t34 + t72) + g(2) * (-rSges(6,3) * t32 + t73) + t68 + ((rSges(6,3) - t55) * t54 + t66 * t51) * t93 + t99 * (qJ(3) + t66)) - m(7) * (g(1) * (t100 * t34 + t72) + g(2) * (t100 * t32 + t73) + t68 + ((-t55 - t100) * t54 + t57 * t51) * t93 + t99 * (qJ(3) + t57)) t71 * (g(1) * t34 + g(2) * t32 - g(3) * t81) -m(5) * (g(1) * (t103 * rSges(5,1) + (-t34 * t50 - t46 * t82) * rSges(5,2)) + g(2) * (t104 * rSges(5,1) + (-t32 * t50 + t48 * t82) * rSges(5,2)) + g(3) * (t62 * rSges(5,1) + (t50 * t81 - t76 * t53) * rSges(5,2))) - m(6) * (g(1) * (t65 + t98) + g(2) * (t79 + t97) + g(3) * (t58 + t80)) - m(7) * (g(1) * (t61 + t65) + g(2) * (t60 + t79) + g(3) * (t58 + t59)) -m(6) * (g(1) * t98 + g(2) * t97 + g(3) * t80) - m(7) * (g(1) * t61 + g(2) * t60 + g(3) * t59) -m(7) * (g(1) * ((-t13 * t49 + t35 * t52) * rSges(7,1) + (-t13 * t52 - t35 * t49) * rSges(7,2)) + g(2) * ((t15 * t49 + t33 * t52) * rSges(7,1) + (t15 * t52 - t33 * t49) * rSges(7,2)) + g(3) * ((-t27 * t49 + t52 * t83) * rSges(7,1) + (-t27 * t52 - t49 * t83) * rSges(7,2)))];
taug  = t1(:);
