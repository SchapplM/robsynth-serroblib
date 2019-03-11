% Calculate Gravitation load on the joints for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:17
% EndTime: 2019-03-09 14:29:19
% DurationCPUTime: 0.79s
% Computational Cost: add. (428->163), mult. (578->223), div. (0->0), fcn. (548->10), ass. (0->74)
t84 = rSges(5,3) + pkin(8);
t52 = -pkin(9) - pkin(8);
t74 = rSges(6,3) - t52;
t73 = rSges(7,3) + pkin(10) - t52;
t51 = cos(qJ(1));
t48 = sin(qJ(1));
t87 = g(2) * t48;
t92 = g(1) * t51 + t87;
t45 = qJ(4) + qJ(5);
t38 = qJ(6) + t45;
t32 = sin(t38);
t33 = cos(t38);
t47 = sin(qJ(2));
t78 = t47 * t51;
t5 = -t32 * t48 + t33 * t78;
t6 = t32 * t78 + t33 * t48;
t91 = t5 * rSges(7,1) - t6 * rSges(7,2);
t79 = t47 * t48;
t7 = t32 * t51 + t33 * t79;
t8 = -t32 * t79 + t33 * t51;
t90 = t7 * rSges(7,1) + t8 * rSges(7,2);
t46 = sin(qJ(4));
t89 = pkin(4) * t46;
t50 = cos(qJ(2));
t86 = g(3) * t50;
t41 = t50 * pkin(2);
t83 = rSges(7,1) * t33;
t82 = rSges(4,2) * t50;
t35 = sin(t45);
t23 = pkin(5) * t35 + t89;
t81 = t23 * t47;
t36 = cos(t45);
t80 = t36 * t50;
t49 = cos(qJ(4));
t77 = t48 * t49;
t76 = t49 * t51;
t75 = t50 * t51;
t40 = t49 * pkin(4);
t24 = pkin(5) * t36 + t40;
t37 = t47 * qJ(3);
t72 = t37 + t41;
t71 = t51 * pkin(1) + t48 * pkin(7);
t70 = qJ(3) * t50;
t69 = -pkin(2) - t84;
t68 = t47 * t89;
t67 = -pkin(2) - t74;
t66 = -pkin(2) - t73;
t65 = -pkin(1) - t37;
t64 = pkin(2) * t75 + t51 * t37 + t71;
t63 = g(1) * t69;
t62 = g(1) * t67;
t61 = g(1) * t66;
t60 = rSges(3,1) * t50 - rSges(3,2) * t47;
t58 = rSges(5,1) * t46 + rSges(5,2) * t49;
t13 = -t35 * t48 + t36 * t78;
t15 = t35 * t51 + t36 * t79;
t18 = -t46 * t48 + t47 * t76;
t20 = t46 * t51 + t47 * t77;
t57 = rSges(7,1) * t32 + rSges(7,2) * t33 + t23;
t56 = rSges(6,1) * t35 + rSges(6,2) * t36 + t89;
t27 = t48 * t70;
t29 = t51 * t70;
t55 = g(1) * t29 + g(2) * t27 + g(3) * t72;
t25 = t50 * t32 * rSges(7,2);
t54 = g(1) * t91 + g(2) * t90 + g(3) * (-t50 * t83 + t25);
t14 = t35 * t78 + t36 * t48;
t16 = -t35 * t79 + t36 * t51;
t53 = g(1) * (t13 * rSges(6,1) - t14 * rSges(6,2)) + g(2) * (t15 * rSges(6,1) + t16 * rSges(6,2)) + g(3) * (t50 * t35 * rSges(6,2) - rSges(6,1) * t80);
t42 = t51 * pkin(7);
t34 = t40 + pkin(3);
t22 = pkin(3) + t24;
t21 = -t46 * t79 + t76;
t19 = t46 * t78 + t77;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t48 - rSges(2,2) * t51) + g(2) * (rSges(2,1) * t51 - rSges(2,2) * t48)) - m(3) * (g(1) * (rSges(3,3) * t51 + t42) + g(2) * (rSges(3,1) * t75 - rSges(3,2) * t78 + t71) + (g(1) * (-pkin(1) - t60) + g(2) * rSges(3,3)) * t48) - m(4) * (g(1) * (rSges(4,1) * t51 + t42) + g(2) * (-rSges(4,2) * t75 + rSges(4,3) * t78 + t64) + (g(1) * (-rSges(4,3) * t47 - t41 + t65 + t82) + g(2) * rSges(4,1)) * t48) - m(5) * (g(1) * (rSges(5,1) * t21 - rSges(5,2) * t20 + pkin(3) * t51 + t42) + g(2) * (rSges(5,1) * t19 + rSges(5,2) * t18 + t84 * t75 + t64) + (g(2) * pkin(3) + g(1) * t65 + t50 * t63) * t48) - m(6) * (g(1) * (t16 * rSges(6,1) - t15 * rSges(6,2) + t42) + g(2) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t64) + (g(1) * t34 + g(2) * (t74 * t50 + t68)) * t51 + (g(1) * (t65 - t68) + g(2) * t34 + t50 * t62) * t48) - m(7) * (g(1) * (rSges(7,1) * t8 - rSges(7,2) * t7 + t42) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t64) + (g(1) * t22 + g(2) * (t73 * t50 + t81)) * t51 + (g(1) * (t65 - t81) + g(2) * t22 + t50 * t61) * t48) -m(3) * (g(3) * t60 + t92 * (-rSges(3,1) * t47 - rSges(3,2) * t50)) - m(4) * (g(1) * (rSges(4,3) * t75 + t29) + g(2) * (rSges(4,3) * t48 * t50 + t27) + g(3) * (t72 - t82) + (g(3) * rSges(4,3) + t92 * (rSges(4,2) - pkin(2))) * t47) - m(5) * ((g(3) * t84 + t92 * t58) * t50 + (g(3) * t58 + t51 * t63 + t69 * t87) * t47 + t55) - m(6) * ((g(3) * t74 + t92 * t56) * t50 + (g(3) * t56 + t51 * t62 + t67 * t87) * t47 + t55) - m(7) * ((g(3) * t73 + t92 * t57) * t50 + (g(3) * t57 + t51 * t61 + t66 * t87) * t47 + t55) (-m(4) - m(5) - m(6) - m(7)) * (t92 * t47 - t86) -m(5) * (g(1) * (rSges(5,1) * t18 - rSges(5,2) * t19) + g(2) * (rSges(5,1) * t20 + rSges(5,2) * t21) + (-rSges(5,1) * t49 + rSges(5,2) * t46) * t86) - m(6) * ((g(1) * t18 + g(2) * t20 - t49 * t86) * pkin(4) + t53) - m(7) * (g(1) * (-t23 * t48 + t24 * t78 + t91) + g(2) * (t23 * t51 + t24 * t79 + t90) + g(3) * (t25 + (-t24 - t83) * t50)) -m(6) * t53 - m(7) * ((g(1) * t13 + g(2) * t15 - g(3) * t80) * pkin(5) + t54) -m(7) * t54];
taug  = t1(:);
