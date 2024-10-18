% Calculate Gravitation load on the joints for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:14
% DurationCPUTime: 0.34s
% Computational Cost: add. (543->108), mult. (365->154), div. (0->0), fcn. (322->20), ass. (0->62)
t83 = pkin(7) + pkin(8);
t52 = qJ(3) + qJ(4);
t68 = pkin(5) - t52;
t67 = -qJ(5) + t68;
t61 = sin(t67);
t46 = pkin(5) + t52;
t39 = qJ(5) + t46;
t72 = sin(t39) / 0.2e1;
t22 = t72 - t61 / 0.2e1;
t49 = qJ(5) + t52;
t41 = cos(t49);
t51 = pkin(10) + qJ(2);
t44 = sin(t51);
t45 = cos(t51);
t33 = cos(t39) / 0.2e1;
t37 = cos(t67);
t23 = t37 / 0.2e1 + t33;
t40 = sin(t49);
t5 = -t23 * t45 + t40 * t44;
t82 = -t5 * rSges(6,1) + (-t22 * t45 - t41 * t44) * rSges(6,2);
t6 = -t23 * t44 - t40 * t45;
t81 = t6 * rSges(6,1) + (t22 * t44 - t41 * t45) * rSges(6,2);
t53 = sin(pkin(5));
t80 = g(3) * t53;
t58 = cos(qJ(3));
t43 = t58 * pkin(3) + pkin(2);
t47 = sin(t52);
t78 = t44 * t47;
t77 = t45 * t47;
t54 = cos(pkin(5));
t56 = sin(qJ(3));
t76 = t54 * t56;
t75 = t54 * t58;
t55 = sin(qJ(4));
t74 = t55 * t56;
t73 = (t72 + t61 / 0.2e1) * rSges(6,1) + (t33 - t37 / 0.2e1) * rSges(6,2);
t71 = sin(t46) / 0.2e1;
t48 = cos(t52);
t70 = rSges(5,1) * t48 + t43;
t69 = rSges(6,1) * t41 + pkin(4) * t48 + t43;
t66 = sin(t68);
t57 = cos(qJ(4));
t42 = pkin(4) * t57 + pkin(3);
t65 = -pkin(4) * t74 + t42 * t58;
t64 = -t44 * t56 + t45 * t75;
t18 = -t44 * t75 - t45 * t56;
t26 = t71 - t66 / 0.2e1;
t63 = -t26 * rSges(5,1) - pkin(3) * t76 + (rSges(5,3) + t83) * t53;
t62 = -t22 * rSges(6,1) - (pkin(4) * t55 * t58 + t42 * t56) * t54 + (rSges(6,3) + pkin(9) + t83) * t53;
t60 = pkin(4) * (t57 * t58 - t74);
t36 = cos(t46) / 0.2e1;
t38 = cos(t68);
t27 = t38 / 0.2e1 + t36;
t11 = -t27 * t45 + t78;
t12 = -t27 * t44 - t77;
t59 = g(1) * (t12 * rSges(5,1) + (t26 * t44 - t45 * t48) * rSges(5,2)) + g(2) * (-t11 * rSges(5,1) + (-t26 * t45 - t44 * t48) * rSges(5,2)) + g(3) * ((t71 + t66 / 0.2e1) * rSges(5,1) + (t36 - t38 / 0.2e1) * rSges(5,2));
t30 = -pkin(3) * t56 - pkin(4) * t47;
t19 = -t44 * t76 + t45 * t58;
t17 = -t44 * t58 - t45 * t76;
t15 = t54 * t60;
t14 = t65 * t54;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t44 - rSges(3,2) * t45) + g(2) * (rSges(3,1) * t45 - rSges(3,2) * t44)) - m(4) * (g(1) * (t17 * rSges(4,1) - rSges(4,2) * t64 - t44 * pkin(2)) + g(2) * (t19 * rSges(4,1) + t18 * rSges(4,2) + t45 * pkin(2)) + (g(1) * t45 + g(2) * t44) * t53 * (rSges(4,3) + pkin(7))) - m(5) * (g(1) * (t11 * rSges(5,2) - t44 * t70 + t45 * t63) + g(2) * (t12 * rSges(5,2) + t44 * t63 + t45 * t70)) - m(6) * (g(1) * (t5 * rSges(6,2) - t44 * t69 + t45 * t62) + g(2) * (t6 * rSges(6,2) + t44 * t62 + t45 * t69)), -m(4) * (g(1) * (rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (rSges(4,1) * t64 + rSges(4,2) * t17) + (rSges(4,1) * t58 - rSges(4,2) * t56) * t80) - m(5) * ((g(1) * t18 + g(2) * t64 + t58 * t80) * pkin(3) + t59) - m(6) * (g(1) * (-t14 * t44 + t30 * t45 + t81) + g(2) * (t14 * t45 + t30 * t44 + t82) + g(3) * (t53 * t65 + t73)), -m(5) * t59 - m(6) * (g(1) * (-pkin(4) * t77 - t15 * t44 + t81) + g(2) * (-pkin(4) * t78 + t15 * t45 + t82) + g(3) * (t53 * t60 + t73)), -m(6) * (g(1) * t81 + g(2) * t82 + g(3) * t73)];
taug = t1(:);
