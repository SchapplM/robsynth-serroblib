% Calculate Gravitation load on the joints for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:34
% EndTime: 2019-03-09 15:56:37
% DurationCPUTime: 1.13s
% Computational Cost: add. (451->200), mult. (913->261), div. (0->0), fcn. (986->10), ass. (0->72)
t41 = sin(qJ(1));
t39 = sin(qJ(3));
t44 = cos(qJ(1));
t78 = t44 * t39;
t42 = cos(qJ(3));
t43 = cos(qJ(2));
t80 = t42 * t43;
t15 = t41 * t80 - t78;
t35 = pkin(10) + qJ(6);
t28 = sin(t35);
t81 = t41 * t43;
t14 = t39 * t81 + t42 * t44;
t29 = cos(t35);
t6 = t14 * t29;
t96 = -t15 * t28 + t6;
t91 = g(2) * t41;
t95 = g(1) * t44 + t91;
t94 = -m(6) - m(7);
t93 = g(1) * t41;
t40 = sin(qJ(2));
t18 = t40 * t42 * qJ(4);
t90 = g(3) * t18;
t89 = g(3) * t40;
t32 = t43 * pkin(2);
t88 = -rSges(5,1) - pkin(3);
t36 = sin(pkin(10));
t87 = t14 * t36;
t85 = t36 * t39;
t84 = t36 * t42;
t83 = t39 * t43;
t82 = t40 * t44;
t79 = t43 * t44;
t77 = t40 * pkin(8) + t32;
t76 = t44 * pkin(1) + t41 * pkin(7);
t75 = -pkin(9) - qJ(5) - rSges(7,3);
t74 = rSges(5,3) + qJ(4);
t73 = -qJ(5) - rSges(6,3);
t72 = -pkin(1) - t32;
t37 = cos(pkin(10));
t71 = -t14 * t37 + t15 * t36;
t16 = -t41 * t42 + t43 * t78;
t17 = t41 * t39 + t42 * t79;
t70 = t16 * t36 + t17 * t37;
t69 = t75 * t44;
t25 = pkin(5) * t37 + pkin(4);
t68 = t28 * rSges(7,2) - t25;
t67 = pkin(5) * t36 + qJ(4);
t66 = t73 * t44;
t65 = pkin(3) * t80 + qJ(4) * t83 + t77;
t64 = pkin(2) * t79 + pkin(8) * t82 + t76;
t33 = t44 * pkin(7);
t63 = -t15 * pkin(3) - t14 * qJ(4) + t33;
t62 = t17 * pkin(3) + t64;
t53 = t28 * t39 + t29 * t42;
t54 = t28 * t42 - t29 * t39;
t61 = (-rSges(7,1) * t54 - rSges(7,2) * t53) * t40;
t60 = rSges(3,1) * t43 - rSges(3,2) * t40;
t58 = t36 * rSges(6,1) + t37 * rSges(6,2);
t57 = -t14 * t28 - t15 * t29;
t56 = t15 * t37 + t87;
t55 = t16 * t37 - t17 * t36;
t51 = -t29 * rSges(7,1) + t68;
t50 = -rSges(6,1) * t37 + rSges(6,2) * t36 - pkin(3) - pkin(4);
t49 = t28 * rSges(7,1) + t29 * rSges(7,2) + t67;
t20 = pkin(8) * t81;
t23 = pkin(8) * t79;
t47 = g(1) * t23 + g(2) * t20 + g(3) * t65;
t12 = t16 * pkin(3);
t10 = t14 * pkin(3);
t2 = t16 * t28 + t17 * t29;
t1 = t16 * t29 - t17 * t28;
t3 = [-m(2) * (g(1) * (-t41 * rSges(2,1) - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - t41 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t44 + t33) + g(2) * (rSges(3,1) * t79 - rSges(3,2) * t82 + t76) + (g(1) * (-pkin(1) - t60) + g(2) * rSges(3,3)) * t41) - m(4) * (g(1) * (-rSges(4,1) * t15 + rSges(4,2) * t14 + t33) + g(2) * (t17 * rSges(4,1) - t16 * rSges(4,2) + rSges(4,3) * t82 + t64) + ((-rSges(4,3) - pkin(8)) * t40 + t72) * t93) - m(5) * (g(1) * (-rSges(5,1) * t15 - rSges(5,3) * t14 + t63) + g(2) * (t17 * rSges(5,1) + rSges(5,2) * t82 + t74 * t16 + t62) + ((-rSges(5,2) - pkin(8)) * t40 + t72) * t93) - m(6) * (g(1) * (-t56 * rSges(6,1) + t71 * rSges(6,2) - t15 * pkin(4) + t63) + g(2) * (t70 * rSges(6,1) + t55 * rSges(6,2) + t17 * pkin(4) + t16 * qJ(4) + t40 * t66 + t62) + ((-pkin(8) - t73) * t40 + t72) * t93) - m(7) * (g(1) * (t57 * rSges(7,1) - rSges(7,2) * t96 - pkin(5) * t87 - t15 * t25 + t63) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t67 * t16 + t17 * t25 + t40 * t69 + t62) + ((-pkin(8) - t75) * t40 + t72) * t93) -m(3) * (g(3) * t60 + t95 * (-rSges(3,1) * t40 - rSges(3,2) * t43)) - m(4) * (g(1) * (rSges(4,3) * t79 + t23) + g(2) * (rSges(4,3) * t81 + t20) + g(3) * (rSges(4,1) * t80 - rSges(4,2) * t83 + t77) + (g(3) * rSges(4,3) + t95 * (-rSges(4,1) * t42 + rSges(4,2) * t39 - pkin(2))) * t40) - m(5) * (g(1) * (rSges(5,2) * t79 + t23) + g(2) * (rSges(5,2) * t81 + t20) + g(3) * (rSges(5,1) * t80 + rSges(5,3) * t83 + t65) + (g(3) * rSges(5,2) + t95 * (-t74 * t39 + t88 * t42 - pkin(2))) * t40) - m(6) * ((g(3) * (t42 * pkin(4) + (t37 * t42 + t85) * rSges(6,1) + (t37 * t39 - t84) * rSges(6,2)) + g(1) * t66 + t73 * t91) * t43 + (g(3) * t73 + t95 * (-pkin(2) + (-qJ(4) - t58) * t39 + t50 * t42)) * t40 + t47) - m(7) * ((g(3) * (rSges(7,1) * t53 - rSges(7,2) * t54 + pkin(5) * t85 + t42 * t25) + g(1) * t69 + t75 * t91) * t43 + (g(3) * t75 + t95 * (-pkin(2) + (-pkin(3) + t51) * t42 - t49 * t39)) * t40 + t47) -m(4) * (g(1) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(2) * (-rSges(4,1) * t14 - rSges(4,2) * t15) + (-rSges(4,1) * t39 - rSges(4,2) * t42) * t89) - m(5) * (g(1) * (-rSges(5,1) * t16 + t74 * t17 - t12) + g(2) * (-rSges(5,1) * t14 + t74 * t15 - t10) + t90 + (rSges(5,3) * t42 + t88 * t39) * t89) - m(6) * (g(1) * (-t55 * rSges(6,1) + t70 * rSges(6,2) - t16 * pkin(4) + t17 * qJ(4) - t12) + g(2) * (t71 * rSges(6,1) + t56 * rSges(6,2) - t14 * pkin(4) + t15 * qJ(4) - t10) + t90 + (t39 * t50 + t42 * t58) * t89) - m(7) * (g(1) * (t16 * t51 + t17 * t49 - t12) + g(2) * (-t6 * rSges(7,1) + t14 * t68 + t15 * t49 - t10) + g(3) * (t18 - t61) + (pkin(5) * t84 + (-pkin(3) - t25) * t39) * t89) (-m(5) + t94) * (g(1) * t16 + g(2) * t14 + t39 * t89) t94 * (g(3) * t43 - t40 * t95) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (rSges(7,1) * t96 + t57 * rSges(7,2)) + g(3) * t61)];
taug  = t3(:);
