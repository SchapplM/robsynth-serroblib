% Calculate Gravitation load on the joints for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:23
% EndTime: 2019-12-31 22:24:26
% DurationCPUTime: 0.68s
% Computational Cost: add. (365->116), mult. (410->166), div. (0->0), fcn. (375->10), ass. (0->66)
t97 = rSges(5,3) + pkin(8);
t40 = qJ(2) + qJ(3);
t35 = sin(t40);
t37 = cos(t40);
t96 = t37 * rSges(4,1) - t35 * rSges(4,2);
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t95 = g(1) * t46 + g(2) * t43;
t44 = cos(qJ(4));
t32 = t44 * pkin(4) + pkin(3);
t47 = -pkin(9) - pkin(8);
t53 = t37 * t32 + (rSges(6,3) - t47) * t35;
t56 = t37 * pkin(3) + t97 * t35;
t39 = qJ(4) + qJ(5);
t36 = cos(t39);
t73 = t46 * t36;
t34 = sin(t39);
t78 = t43 * t34;
t5 = t37 * t78 + t73;
t74 = t46 * t34;
t77 = t43 * t36;
t6 = -t37 * t77 + t74;
t94 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = -t37 * t74 + t77;
t8 = t37 * t73 + t78;
t93 = t7 * rSges(6,1) - t8 * rSges(6,2);
t42 = sin(qJ(2));
t92 = pkin(2) * t42;
t41 = sin(qJ(4));
t91 = pkin(4) * t41;
t88 = g(3) * t35;
t87 = rSges(3,3) + pkin(6);
t86 = rSges(5,1) * t44;
t85 = rSges(6,1) * t36;
t84 = rSges(5,2) * t41;
t83 = rSges(6,2) * t34;
t80 = t37 * t43;
t79 = t37 * t46;
t76 = t43 * t41;
t75 = t43 * t44;
t72 = t46 * t41;
t71 = t46 * t44;
t48 = -pkin(7) - pkin(6);
t70 = rSges(4,3) - t48;
t69 = t35 * t84;
t68 = t35 * t83;
t67 = t43 * t69 + t97 * t80;
t66 = t46 * t69 + t97 * t79;
t65 = -pkin(3) - t86;
t64 = -t48 + t91;
t63 = -t32 - t85;
t45 = cos(qJ(2));
t61 = t45 * rSges(3,1) - t42 * rSges(3,2);
t58 = -rSges(6,1) * t34 - rSges(6,2) * t36;
t57 = pkin(1) + t61;
t11 = -t37 * t72 + t75;
t9 = t37 * t76 + t71;
t55 = t56 + (-t84 + t86) * t37;
t52 = g(1) * (rSges(6,3) * t79 + t46 * t68) + g(2) * (rSges(6,3) * t80 + t43 * t68);
t51 = t53 + (-t83 + t85) * t37;
t38 = t45 * pkin(2);
t33 = t38 + pkin(1);
t19 = t46 * t33;
t12 = t37 * t71 + t76;
t10 = -t37 * t75 + t72;
t1 = [-m(2) * (g(1) * (-t43 * rSges(2,1) - t46 * rSges(2,2)) + g(2) * (t46 * rSges(2,1) - t43 * rSges(2,2))) - m(3) * ((g(1) * t87 + g(2) * t57) * t46 + (-g(1) * t57 + g(2) * t87) * t43) - m(4) * (g(2) * t19 + (g(1) * t70 + g(2) * t96) * t46 + (g(1) * (-t33 - t96) + g(2) * t70) * t43) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2)) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t19) + (-g(1) * t48 + g(2) * t56) * t46 + (g(1) * (-t33 - t56) - g(2) * t48) * t43) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2)) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t19) + (g(1) * t64 + g(2) * t53) * t46 + (g(1) * (-t33 - t53) + g(2) * t64) * t43), -m(3) * (g(3) * t61 + t95 * (-rSges(3,1) * t42 - rSges(3,2) * t45)) - m(4) * (g(3) * (t38 + t96) + t95 * (-rSges(4,1) * t35 - rSges(4,2) * t37 - t92)) - m(5) * (g(1) * (-t46 * t92 + t66) + g(2) * (-t43 * t92 + t67) + g(3) * (t38 + t55) + t95 * t35 * t65) - m(6) * (g(3) * (t38 + t51) + t52 + t95 * (t63 * t35 - t37 * t47 - t92)), -m(4) * g(3) * t96 - m(5) * (g(1) * t66 + g(2) * t67 + g(3) * t55) - m(6) * (g(3) * t51 + t52) + t95 * ((m(4) * rSges(4,2) + m(6) * t47) * t37 + (m(4) * rSges(4,1) - m(5) * t65 - m(6) * t63) * t35), -m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2))) - m(6) * (g(1) * (t11 * pkin(4) + t93) + g(2) * (-t9 * pkin(4) + t94)) + (-m(5) * (-rSges(5,1) * t41 - rSges(5,2) * t44) - m(6) * (t58 - t91)) * t88, -m(6) * (g(1) * t93 + g(2) * t94 + t58 * t88)];
taug = t1(:);
