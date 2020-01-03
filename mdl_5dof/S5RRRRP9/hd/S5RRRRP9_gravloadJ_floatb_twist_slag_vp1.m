% Calculate Gravitation load on the joints for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:48
% EndTime: 2019-12-31 22:03:50
% DurationCPUTime: 0.69s
% Computational Cost: add. (345->119), mult. (493->167), div. (0->0), fcn. (487->8), ass. (0->59)
t75 = rSges(6,1) + pkin(4);
t36 = sin(qJ(3));
t38 = sin(qJ(1));
t39 = cos(qJ(3));
t40 = cos(qJ(2));
t41 = cos(qJ(1));
t66 = t40 * t41;
t17 = -t36 * t66 + t38 * t39;
t60 = rSges(6,3) + qJ(5);
t79 = g(2) * t38;
t87 = g(1) * t41 + t79;
t37 = sin(qJ(2));
t74 = rSges(4,3) + pkin(7);
t86 = t40 * pkin(2) + t74 * t37;
t35 = qJ(3) + qJ(4);
t30 = sin(t35);
t31 = cos(t35);
t85 = t60 * t30 + t75 * t31;
t67 = t38 * t40;
t11 = t30 * t67 + t31 * t41;
t65 = t41 * t30;
t12 = t31 * t67 - t65;
t84 = -t11 * rSges(5,1) - t12 * rSges(5,2);
t13 = -t38 * t31 + t40 * t65;
t14 = t30 * t38 + t31 * t66;
t83 = -t13 * rSges(5,1) - t14 * rSges(5,2);
t82 = pkin(3) * t36;
t81 = g(1) * t38;
t29 = pkin(3) * t39 + pkin(2);
t21 = t40 * t29;
t78 = g(3) * t21;
t77 = g(3) * t37;
t73 = rSges(3,2) * t37;
t72 = t30 * t37;
t70 = t36 * t38;
t69 = t36 * t41;
t42 = -pkin(8) - pkin(7);
t64 = rSges(6,2) - t42;
t63 = rSges(5,3) - t42;
t62 = t60 * t31 * t37;
t61 = t41 * pkin(1) + t38 * pkin(6);
t33 = t41 * pkin(6);
t58 = t38 * t37 * t42 + pkin(3) * t69 + t33;
t57 = -pkin(1) - t21;
t56 = t64 * t41;
t55 = t63 * t41;
t54 = pkin(3) * t70 + t29 * t66 + t61;
t53 = rSges(3,1) * t40 - t73;
t51 = rSges(5,1) * t31 - rSges(5,2) * t30;
t50 = -rSges(5,1) * t30 - rSges(5,2) * t31;
t49 = -t75 * t11 + t60 * t12;
t48 = t17 * pkin(3);
t47 = -t75 * t13 + t60 * t14;
t46 = rSges(4,1) * t39 - rSges(4,2) * t36 + pkin(2);
t15 = t36 * t67 + t39 * t41;
t44 = t15 * pkin(3);
t18 = t39 * t66 + t70;
t16 = -t39 * t67 + t69;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t38 - rSges(2,2) * t41) + g(2) * (rSges(2,1) * t41 - rSges(2,2) * t38)) - m(3) * (g(1) * (rSges(3,3) * t41 + t33) + g(2) * (rSges(3,1) * t66 - t41 * t73 + t61) + (g(1) * (-pkin(1) - t53) + g(2) * rSges(3,3)) * t38) - m(4) * (g(1) * (rSges(4,1) * t16 + rSges(4,2) * t15 + t33) + (-pkin(1) - t86) * t81 + (rSges(4,1) * t18 + rSges(4,2) * t17 + t86 * t41 + t61) * g(2)) - m(5) * (g(1) * (-rSges(5,1) * t12 + rSges(5,2) * t11 + t58) + g(2) * (t14 * rSges(5,1) - t13 * rSges(5,2) + t37 * t55 + t54) + (-t37 * rSges(5,3) + t57) * t81) - m(6) * (g(1) * (-t60 * t11 - t75 * t12 + t58) + g(2) * (t60 * t13 + t75 * t14 + t37 * t56 + t54) + (-t37 * rSges(6,2) + t57) * t81), -m(3) * (g(3) * t53 + t87 * (-rSges(3,1) * t37 - rSges(3,2) * t40)) - m(4) * ((g(3) * t46 + t87 * t74) * t40 + (g(3) * t74 - t87 * t46) * t37) - m(5) * (t78 + (g(1) * t55 + g(3) * t51 + t63 * t79) * t40 + (g(3) * t63 + t87 * (-t29 - t51)) * t37) - m(6) * (t78 + (g(1) * t56 + g(3) * t85 + t64 * t79) * t40 + (g(3) * t64 + t87 * (-t29 - t85)) * t37), -m(4) * (g(1) * (rSges(4,1) * t17 - rSges(4,2) * t18) + g(2) * (-rSges(4,1) * t15 + rSges(4,2) * t16)) - m(5) * (g(1) * (t48 + t83) + g(2) * (-t44 + t84)) - m(6) * (g(1) * (t47 + t48) + g(2) * (-t44 + t49) + g(3) * t62) + (-m(4) * (-rSges(4,1) * t36 - rSges(4,2) * t39) - m(5) * (t50 - t82) - m(6) * (-t75 * t30 - t82)) * t77, -m(5) * (g(1) * t83 + g(2) * t84 + t50 * t77) - m(6) * (g(1) * t47 + g(2) * t49 + g(3) * (-t75 * t72 + t62)), -m(6) * (g(1) * t13 + g(2) * t11 + g(3) * t72)];
taug = t1(:);
