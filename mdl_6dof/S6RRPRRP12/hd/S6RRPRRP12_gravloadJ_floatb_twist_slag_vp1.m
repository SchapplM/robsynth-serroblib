% Calculate Gravitation load on the joints for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:25
% EndTime: 2019-03-09 12:50:27
% DurationCPUTime: 0.86s
% Computational Cost: add. (403->160), mult. (613->210), div. (0->0), fcn. (596->8), ass. (0->75)
t94 = rSges(7,1) + pkin(5);
t46 = qJ(4) + qJ(5);
t39 = sin(t46);
t74 = rSges(7,3) + qJ(6);
t105 = t74 * t39;
t51 = cos(qJ(2));
t96 = g(3) * t51;
t91 = rSges(5,3) + pkin(8);
t52 = cos(qJ(1));
t49 = sin(qJ(1));
t98 = g(2) * t49;
t104 = g(1) * t52 + t98;
t40 = cos(t46);
t48 = sin(qJ(2));
t85 = t48 * t52;
t11 = t39 * t49 - t40 * t85;
t12 = t39 * t85 + t40 * t49;
t103 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t86 = t48 * t49;
t13 = t39 * t52 + t40 * t86;
t14 = -t39 * t86 + t40 * t52;
t102 = t13 * rSges(6,1) + t14 * rSges(6,2);
t47 = sin(qJ(4));
t101 = pkin(4) * t47;
t100 = g(1) * t49;
t97 = t39 * rSges(6,2) * t96;
t50 = cos(qJ(4));
t95 = t50 * pkin(4);
t43 = t51 * pkin(2);
t92 = -rSges(7,2) - pkin(2);
t90 = -rSges(6,3) - pkin(2);
t89 = rSges(4,2) * t51;
t88 = t47 * t49;
t87 = t47 * t52;
t84 = t49 * t50;
t53 = -pkin(9) - pkin(8);
t83 = t49 * t53;
t82 = t50 * t52;
t81 = t51 * t52;
t80 = rSges(7,2) - t53;
t79 = rSges(6,3) - t53;
t36 = pkin(4) * t87;
t70 = t48 * t84;
t78 = pkin(4) * t70 + t36;
t41 = t48 * qJ(3);
t77 = t41 + t43;
t76 = t52 * pkin(1) + t49 * pkin(7);
t75 = qJ(3) * t51;
t73 = -pkin(2) - t91;
t72 = pkin(4) * t88;
t71 = t47 * t85;
t69 = t48 * t82;
t38 = pkin(3) + t95;
t44 = t52 * pkin(7);
t68 = t52 * t38 + t51 * t83 + t44;
t67 = -pkin(1) - t41;
t66 = pkin(2) * t81 + t52 * t41 + t76;
t65 = g(1) * t73;
t64 = pkin(4) * t69 - t72;
t63 = rSges(3,1) * t51 - rSges(3,2) * t48;
t61 = rSges(5,1) * t47 + rSges(5,2) * t50;
t60 = rSges(6,1) * t39 + rSges(6,2) * t40;
t59 = -t11 * t94 + t74 * t12;
t58 = t13 * t94 - t74 * t14;
t57 = pkin(4) * t71 + t49 * t38 + t66;
t56 = (-qJ(3) - t101) * t48 - pkin(1);
t55 = t39 * t94 - t40 * t74;
t32 = t49 * t75;
t34 = t52 * t75;
t54 = g(1) * (t51 * t36 + t53 * t85 + t34) + g(2) * (t48 * t83 + t51 * t72 + t32) + g(3) * (t48 * t101 + t77);
t19 = -t47 * t86 + t82;
t18 = t70 + t87;
t17 = t71 + t84;
t16 = t69 - t88;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t49 - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 - t49 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t52 + t44) + g(2) * (rSges(3,1) * t81 - rSges(3,2) * t85 + t76) + (g(1) * (-pkin(1) - t63) + g(2) * rSges(3,3)) * t49) - m(4) * (g(1) * (rSges(4,1) * t52 + t44) + g(2) * (-rSges(4,2) * t81 + rSges(4,3) * t85 + t66) + (g(1) * (-rSges(4,3) * t48 - t43 + t67 + t89) + g(2) * rSges(4,1)) * t49) - m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t18 + pkin(3) * t52 + t44) + g(2) * (rSges(5,1) * t17 + rSges(5,2) * t16 + t91 * t81 + t66) + (g(2) * pkin(3) + g(1) * t67 + t51 * t65) * t49) - m(6) * (g(1) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t68) + g(2) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t79 * t81 + t57) + (t51 * t90 + t56) * t100) - m(7) * (g(1) * (t74 * t13 + t94 * t14 + t68) + g(2) * (t11 * t74 + t12 * t94 + t80 * t81 + t57) + (t51 * t92 + t56) * t100) -m(3) * (g(3) * t63 + t104 * (-rSges(3,1) * t48 - rSges(3,2) * t51)) - m(4) * (g(1) * (rSges(4,3) * t81 + t34) + g(2) * (rSges(4,3) * t49 * t51 + t32) + g(3) * (t77 - t89) + (g(3) * rSges(4,3) + t104 * (rSges(4,2) - pkin(2))) * t48) - m(5) * (g(1) * t34 + g(2) * t32 + g(3) * t77 + (g(3) * t91 + t104 * t61) * t51 + (g(3) * t61 + t52 * t65 + t73 * t98) * t48) - m(6) * ((g(3) * t79 + t104 * t60) * t51 + (g(3) * t60 + t104 * t90) * t48 + t54) - m(7) * ((g(3) * t55 + t104 * t92) * t48 + (g(3) * t80 + t104 * t55) * t51 + t54) (-m(4) - m(5) - m(6) - m(7)) * (t104 * t48 - t96) -m(5) * (g(1) * (rSges(5,1) * t16 - rSges(5,2) * t17) + g(2) * (rSges(5,1) * t18 + rSges(5,2) * t19)) - m(6) * (g(1) * (t64 + t103) + g(2) * (t78 + t102) + t97) - m(7) * (g(1) * (t59 + t64) + g(2) * (t58 + t78)) + (-m(5) * (-rSges(5,1) * t50 + rSges(5,2) * t47) - m(6) * (-rSges(6,1) * t40 - t95) - m(7) * (-t94 * t40 - t105 - t95)) * t96, -m(6) * (g(1) * t103 + g(2) * t102 + t97) - m(7) * (g(1) * t59 + g(2) * t58) + (m(7) * t105 + (m(6) * rSges(6,1) + m(7) * t94) * t40) * t96, -m(7) * (g(1) * t11 - g(2) * t13 + t40 * t96)];
taug  = t1(:);
