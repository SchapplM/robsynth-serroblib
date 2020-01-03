% Calculate Gravitation load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:23:23
% DurationCPUTime: 0.64s
% Computational Cost: add. (263->122), mult. (664->170), div. (0->0), fcn. (733->8), ass. (0->63)
t80 = rSges(6,1) + pkin(4);
t43 = sin(qJ(2));
t40 = sin(pkin(5));
t63 = qJ(3) * t40;
t30 = t43 * t63;
t45 = cos(qJ(2));
t66 = t45 * pkin(2) + t30;
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t87 = g(1) * t46 + g(2) * t44;
t86 = rSges(6,2) + qJ(4);
t85 = rSges(6,3) + qJ(5);
t84 = -m(5) - m(6);
t83 = pkin(2) * t43;
t39 = sin(pkin(8));
t41 = cos(pkin(8));
t42 = cos(pkin(5));
t73 = t42 * t45;
t48 = t39 * t73 + t41 * t43;
t12 = t48 * t44;
t56 = t45 * t63;
t24 = t44 * t56;
t79 = -t12 * pkin(3) + t24;
t78 = t40 * t43;
t77 = t40 * t44;
t76 = t40 * t45;
t75 = t40 * t46;
t74 = t42 * t43;
t72 = t43 * t44;
t71 = t43 * t46;
t70 = t44 * t42;
t69 = t44 * t45;
t68 = t45 * t46;
t14 = t48 * t46;
t26 = t46 * t56;
t67 = -t14 * pkin(3) + t26;
t37 = t46 * pkin(7);
t62 = qJ(3) * t42;
t65 = t46 * t62 + t37;
t64 = t46 * pkin(1) + t44 * pkin(7);
t61 = rSges(5,3) + qJ(4);
t60 = t41 * t69;
t59 = t41 * t68;
t58 = t39 * t74;
t17 = t41 * t45 - t58;
t57 = t17 * pkin(3) + t66;
t55 = pkin(2) * t68 + t46 * t30 + t44 * t62 + t64;
t54 = rSges(3,1) * t45 - rSges(3,2) * t43;
t5 = t39 * t69 + (t43 * t70 + t75) * t41;
t6 = -t39 * t75 - t44 * t58 + t60;
t52 = -t6 * pkin(3) - qJ(4) * t5 + t65;
t8 = t39 * t77 - t46 * t58 + t59;
t51 = t8 * pkin(3) + t55;
t50 = rSges(5,1) * t76 - t83;
t49 = rSges(4,3) * t76 - t83;
t16 = t39 * t45 + t41 * t74;
t47 = g(1) * (-pkin(1) - t66) * t44;
t19 = t40 * t71 + t70;
t18 = t40 * t72 - t46 * t42;
t13 = t39 * t71 - t42 * t59;
t11 = t39 * t72 - t42 * t60;
t7 = t16 * t46 - t41 * t77;
t1 = [-m(2) * (g(1) * (-t44 * rSges(2,1) - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - t44 * rSges(2,2))) - m(3) * (g(1) * (t46 * rSges(3,3) + t37) + g(2) * (rSges(3,1) * t68 - rSges(3,2) * t71 + t64) + (g(1) * (-pkin(1) - t54) + g(2) * rSges(3,3)) * t44) - m(4) * (g(1) * (-rSges(4,1) * t6 + rSges(4,2) * t5 - rSges(4,3) * t18 + t65) + g(2) * (rSges(4,1) * t8 - rSges(4,2) * t7 + rSges(4,3) * t19 + t55) + t47) - m(5) * (g(1) * (-rSges(5,1) * t18 + rSges(5,2) * t6 - rSges(5,3) * t5 + t52) + g(2) * (rSges(5,1) * t19 - rSges(5,2) * t8 + t61 * t7 + t51) + t47) - m(6) * (g(1) * (-rSges(6,2) * t5 - t80 * t18 - t85 * t6 + t52) + g(2) * (t80 * t19 + t86 * t7 + t85 * t8 + t51) + t47), -m(3) * (g(3) * t54 + t87 * (-rSges(3,1) * t43 - rSges(3,2) * t45)) - m(4) * (g(1) * (-t14 * rSges(4,1) + t13 * rSges(4,2) + t49 * t46 + t26) + g(2) * (-t12 * rSges(4,1) + t11 * rSges(4,2) + t49 * t44 + t24) + g(3) * (rSges(4,1) * t17 - rSges(4,2) * t16 + rSges(4,3) * t78 + t66)) - m(5) * (g(1) * (t14 * rSges(5,2) - t61 * t13 + t50 * t46 + t67) + g(2) * (t12 * rSges(5,2) - t61 * t11 + t50 * t44 + t79) + g(3) * (rSges(5,1) * t78 - t17 * rSges(5,2) + t61 * t16 + t57)) - m(6) * (g(1) * (-pkin(2) * t71 - t86 * t13 - t85 * t14 + t67) + g(2) * (-pkin(2) * t72 - t86 * t11 - t85 * t12 + t79) + g(3) * (t86 * t16 + t85 * t17 + t57) + (g(3) * t43 + t87 * t45) * t80 * t40), (-m(4) + t84) * (g(1) * t19 + g(2) * t18 - g(3) * t76), t84 * (g(1) * t7 + g(2) * t5 + g(3) * (t39 * t43 - t41 * t73)), -m(6) * (g(1) * t8 + g(2) * t6 + g(3) * t48)];
taug = t1(:);
