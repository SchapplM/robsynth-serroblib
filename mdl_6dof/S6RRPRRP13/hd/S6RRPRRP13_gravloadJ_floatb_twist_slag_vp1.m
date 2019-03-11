% Calculate Gravitation load on the joints for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:36
% EndTime: 2019-03-09 12:55:39
% DurationCPUTime: 1.12s
% Computational Cost: add. (555->182), mult. (1296->258), div. (0->0), fcn. (1529->10), ass. (0->71)
t51 = sin(qJ(2));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t78 = cos(pkin(6));
t70 = t56 * t78;
t32 = t51 * t70 + t52 * t55;
t71 = t52 * t78;
t34 = -t51 * t71 + t55 * t56;
t108 = g(1) * t34 + g(2) * t32;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t31 = t51 * t52 - t55 * t70;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t47 = sin(pkin(6));
t88 = t47 * t56;
t66 = -t31 * t50 + t54 * t88;
t5 = t32 * t53 + t49 * t66;
t107 = -t32 * t49 + t53 * t66;
t33 = t56 * t51 + t55 * t71;
t90 = t47 * t52;
t15 = t33 * t50 + t54 * t90;
t89 = t47 * t55;
t30 = -t50 * t89 + t78 * t54;
t106 = g(1) * t15 - g(2) * t66 + g(3) * t30;
t14 = -t33 * t54 + t50 * t90;
t29 = t78 * t50 + t54 * t89;
t65 = t31 * t54 + t50 * t88;
t105 = g(1) * t14 - g(2) * t65 + g(3) * t29;
t98 = g(3) * t47;
t97 = rSges(5,3) + pkin(9);
t96 = rSges(6,3) + pkin(10);
t95 = t31 * t49;
t92 = t33 * t49;
t91 = t47 * t51;
t87 = t49 * t50;
t86 = t49 * t51;
t85 = t49 * t55;
t84 = t50 * t53;
t83 = t51 * t53;
t82 = rSges(7,3) + qJ(6) + pkin(10);
t81 = pkin(2) * t89 + qJ(3) * t91;
t80 = t56 * pkin(1) + pkin(8) * t90;
t79 = rSges(4,3) + qJ(3);
t77 = t34 * pkin(2) + t80;
t76 = pkin(9) * t89 + t81;
t75 = pkin(5) * t49 + pkin(9);
t74 = -t52 * pkin(1) + pkin(8) * t88;
t25 = t31 * pkin(2);
t73 = -pkin(9) * t31 - t25;
t27 = t33 * pkin(2);
t72 = -pkin(9) * t33 - t27;
t69 = -t32 * pkin(2) + t74;
t68 = rSges(5,1) * t50 + rSges(5,2) * t54;
t1 = -t15 * t49 + t34 * t53;
t12 = -t30 * t49 + t47 * t83;
t62 = pkin(4) * t50 - t96 * t54;
t61 = pkin(3) * t90 + qJ(3) * t33 + t77;
t45 = pkin(5) * t53 + pkin(4);
t60 = t45 * t50 - t82 * t54;
t57 = pkin(3) * t88 - t31 * qJ(3) + t69;
t21 = (t50 * t83 + t85) * t47;
t20 = (-t50 * t86 + t53 * t55) * t47;
t13 = -t30 * t53 - t47 * t86;
t11 = t34 * t84 - t92;
t10 = -t33 * t53 - t34 * t87;
t9 = t32 * t84 - t95;
t8 = -t31 * t53 - t32 * t87;
t2 = t15 * t53 + t34 * t49;
t3 = [-m(2) * (g(1) * (-t52 * rSges(2,1) - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 - t52 * rSges(2,2))) - m(3) * (g(1) * (-t32 * rSges(3,1) + t31 * rSges(3,2) + rSges(3,3) * t88 + t74) + g(2) * (rSges(3,1) * t34 - rSges(3,2) * t33 + rSges(3,3) * t90 + t80)) - m(4) * (g(1) * (rSges(4,1) * t88 + t32 * rSges(4,2) - t79 * t31 + t69) + g(2) * (rSges(4,1) * t90 - rSges(4,2) * t34 + t79 * t33 + t77)) - m(5) * (g(1) * (rSges(5,1) * t66 - rSges(5,2) * t65 - t97 * t32 + t57) + g(2) * (rSges(5,1) * t15 - rSges(5,2) * t14 + t97 * t34 + t61)) - m(6) * (g(1) * (rSges(6,1) * t107 - rSges(6,2) * t5 + pkin(4) * t66 - pkin(9) * t32 + t65 * t96 + t57) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t1 + pkin(4) * t15 + pkin(9) * t34 + t96 * t14 + t61)) - m(7) * (g(1) * (rSges(7,1) * t107 - rSges(7,2) * t5 - t75 * t32 + t45 * t66 + t65 * t82 + t57) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + t82 * t14 + t15 * t45 + t75 * t34 + t61)) -m(3) * (g(1) * (-rSges(3,1) * t33 - rSges(3,2) * t34) + g(2) * (-rSges(3,1) * t31 - rSges(3,2) * t32) + (rSges(3,1) * t55 - rSges(3,2) * t51) * t98) - m(4) * (g(1) * (rSges(4,2) * t33 + t79 * t34 - t27) + g(2) * (rSges(4,2) * t31 + t79 * t32 - t25) + g(3) * ((-rSges(4,2) * t55 + rSges(4,3) * t51) * t47 + t81)) - m(5) * (g(1) * (-rSges(5,3) * t33 + t72) + g(2) * (-rSges(5,3) * t31 + t73) + g(3) * t76 + (rSges(5,3) * t55 + t68 * t51) * t98 + t108 * (qJ(3) + t68)) - m(6) * (g(1) * (rSges(6,1) * t11 + rSges(6,2) * t10 + t72) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + t73) + t108 * (qJ(3) + t62) + (rSges(6,1) * t21 + rSges(6,2) * t20 + t62 * t91 + t76) * g(3)) - m(7) * (g(1) * (rSges(7,1) * t11 + rSges(7,2) * t10 - pkin(5) * t92 + t72) + g(2) * (rSges(7,1) * t9 + rSges(7,2) * t8 - pkin(5) * t95 + t73) + g(3) * (rSges(7,1) * t21 + rSges(7,2) * t20 + t76) + (pkin(5) * t85 + t60 * t51) * t98 + t108 * (qJ(3) + t60)) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t33 + g(2) * t31 - g(3) * t89) -m(5) * (g(1) * (-rSges(5,1) * t14 - rSges(5,2) * t15) + g(2) * (rSges(5,1) * t65 + rSges(5,2) * t66) + g(3) * (-rSges(5,1) * t29 - rSges(5,2) * t30)) - m(6) * (t106 * t96 - t105 * (rSges(6,1) * t53 - rSges(6,2) * t49 + pkin(4))) - m(7) * (t106 * t82 - t105 * (rSges(7,1) * t53 - rSges(7,2) * t49 + t45)) -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t5 + rSges(6,2) * t107) + g(3) * (rSges(6,1) * t12 + rSges(6,2) * t13)) + (-g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) - g(2) * (rSges(7,1) * t5 + rSges(7,2) * t107) - g(3) * (rSges(7,1) * t12 + rSges(7,2) * t13) - (g(1) * t1 + g(2) * t5 + g(3) * t12) * pkin(5)) * m(7), -m(7) * t105];
taug  = t3(:);
