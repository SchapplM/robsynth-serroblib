% Calculate Gravitation load on the joints for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:10
% EndTime: 2019-07-18 17:17:12
% DurationCPUTime: 0.55s
% Computational Cost: add. (335->115), mult. (400->161), div. (0->0), fcn. (365->10), ass. (0->70)
t45 = cos(qJ(2));
t38 = t45 * pkin(1);
t40 = qJ(2) + qJ(3);
t35 = sin(t40);
t37 = cos(t40);
t91 = t37 * rSges(4,1) - t35 * rSges(4,2);
t97 = -t91 - t38;
t43 = sin(qJ(1));
t39 = qJ(4) + qJ(5);
t34 = sin(t39);
t77 = rSges(6,2) * t34;
t62 = t35 * t77;
t74 = t37 * t43;
t96 = rSges(6,3) * t74 + t43 * t62;
t46 = cos(qJ(1));
t73 = t37 * t46;
t95 = rSges(6,3) * t73 + t46 * t62;
t41 = sin(qJ(4));
t78 = rSges(5,2) * t41;
t63 = t35 * t78;
t94 = rSges(5,3) * t74 + t43 * t63;
t93 = rSges(5,3) * t73 + t46 * t63;
t32 = t37 * pkin(2);
t92 = t35 * rSges(5,3) + t32;
t90 = g(1) * t46 + g(2) * t43;
t31 = t35 * pkin(5);
t44 = cos(qJ(4));
t33 = t44 * pkin(3) + pkin(2);
t89 = -t35 * rSges(6,3) - t37 * t33 - t31;
t88 = t35 * t90;
t36 = cos(t39);
t67 = t46 * t36;
t72 = t43 * t34;
t5 = t37 * t72 + t67;
t68 = t46 * t34;
t71 = t43 * t36;
t6 = -t37 * t71 + t68;
t87 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = -t37 * t68 + t71;
t8 = t37 * t67 + t72;
t86 = t7 * rSges(6,1) - t8 * rSges(6,2);
t42 = sin(qJ(2));
t85 = pkin(1) * t42;
t84 = pkin(3) * t41;
t81 = g(3) * t35;
t80 = rSges(5,1) * t44;
t79 = rSges(6,1) * t36;
t75 = t35 * t46;
t70 = t43 * t41;
t69 = t43 * t44;
t66 = t46 * t41;
t65 = t46 * t44;
t27 = t46 * t38;
t64 = pkin(5) * t75 + t27;
t24 = pkin(5) * t74;
t58 = -t43 * t85 + t24;
t26 = pkin(5) * t73;
t57 = -t46 * t85 + t26;
t56 = t45 * rSges(3,1) - t42 * rSges(3,2);
t53 = -rSges(4,1) * t35 - rSges(4,2) * t37;
t52 = -rSges(6,1) * t34 - rSges(6,2) * t36;
t11 = -t37 * t66 + t69;
t9 = t37 * t70 + t65;
t51 = t31 + t92 + (-t78 + t80) * t37;
t50 = -t89 + (-t77 + t79) * t37;
t48 = (-pkin(2) - t80) * t88;
t47 = (-t33 - t79) * t88;
t12 = t37 * t65 + t70;
t10 = -t37 * t69 + t66;
t1 = [-m(2) * (g(1) * (-t43 * rSges(2,1) - t46 * rSges(2,2)) + g(2) * (t46 * rSges(2,1) - t43 * rSges(2,2))) - m(3) * (g(1) * (t46 * rSges(3,3) - t56 * t43) + g(2) * (t43 * rSges(3,3) + t56 * t46)) - m(4) * (g(2) * t27 + (g(1) * rSges(4,3) + g(2) * t91) * t46 + (g(2) * rSges(4,3) + g(1) * t97) * t43) - m(5) * (g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t92 * t46 + t64) + (t10 * rSges(5,1) + t9 * rSges(5,2) + (-t38 - t32 + (-rSges(5,3) - pkin(5)) * t35) * t43) * g(1)) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + pkin(3) * t66) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + rSges(6,3) * t75 + t33 * t73 + t64) + (g(1) * (-t38 + t89) + g(2) * t84) * t43), -m(3) * (g(3) * t56 + t90 * (-rSges(3,1) * t42 - rSges(3,2) * t45)) - m(4) * (-g(3) * t97 + t90 * (t53 - t85)) - m(5) * (g(1) * (t57 + t93) + g(2) * (t58 + t94) + g(3) * (t38 + t51) + t48) - m(6) * (g(1) * (t57 + t95) + g(2) * (t58 + t96) + g(3) * (t38 + t50) + t47), -m(4) * (g(3) * t91 + t90 * t53) - m(5) * (g(1) * (t26 + t93) + g(2) * (t24 + t94) + g(3) * t51 + t48) - m(6) * (g(1) * (t26 + t95) + g(2) * (t24 + t96) + g(3) * t50 + t47), -m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2))) - m(6) * (g(1) * (t11 * pkin(3) + t86) + g(2) * (-t9 * pkin(3) + t87)) + (-m(5) * (-rSges(5,1) * t41 - rSges(5,2) * t44) - m(6) * (t52 - t84)) * t81, -m(6) * (g(1) * t86 + g(2) * t87 + t52 * t81)];
taug  = t1(:);
