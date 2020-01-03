% Calculate joint inertia matrix for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:06
% DurationCPUTime: 0.46s
% Computational Cost: add. (1085->121), mult. (898->176), div. (0->0), fcn. (794->8), ass. (0->67)
t58 = qJ(1) + pkin(7);
t53 = sin(t58);
t54 = cos(t58);
t97 = t53 * t54;
t59 = qJ(3) + qJ(4);
t55 = sin(t59);
t56 = cos(t59);
t77 = rSges(5,1) * t56 - rSges(5,2) * t55;
t51 = t53 ^ 2;
t52 = t54 ^ 2;
t96 = t53 / 0.2e1;
t95 = -t54 / 0.2e1;
t94 = t53 * pkin(5);
t61 = sin(qJ(1));
t93 = t61 * pkin(1);
t62 = cos(qJ(3));
t92 = rSges(4,1) * t62;
t60 = sin(qJ(3));
t90 = rSges(4,2) * t60;
t88 = t54 * rSges(4,3);
t65 = t53 * rSges(5,3) + t77 * t54;
t8 = t53 * (-t54 * rSges(5,3) + t77 * t53) + t54 * t65;
t87 = t53 * rSges(4,3) + t54 * t92;
t86 = t51 + t52;
t85 = Icges(4,4) * t60;
t84 = Icges(4,4) * t62;
t83 = Icges(5,4) * t55;
t82 = Icges(5,4) * t56;
t33 = Icges(5,5) * t55 + Icges(5,6) * t56;
t68 = -Icges(5,2) * t55 + t82;
t70 = Icges(5,1) * t56 - t83;
t34 = Icges(5,2) * t56 + t83;
t35 = Icges(5,1) * t55 + t82;
t72 = -t34 * t55 + t35 * t56;
t81 = (t56 * (Icges(5,6) * t53 + t68 * t54) + t55 * (Icges(5,5) * t53 + t70 * t54) + t53 * t33 + t72 * t54) * t96 + (t56 * (-Icges(5,6) * t54 + t68 * t53) + t55 * (-Icges(5,5) * t54 + t70 * t53) - t54 * t33 + t72 * t53) * t95;
t66 = Icges(5,5) * t56 - Icges(5,6) * t55;
t16 = -Icges(5,3) * t54 + t66 * t53;
t17 = Icges(5,3) * t53 + t66 * t54;
t80 = -t54 * (t52 * t16 - t17 * t97) + t53 * (-t16 * t97 + t51 * t17);
t36 = t55 * rSges(5,1) + t56 * rSges(5,2);
t79 = -pkin(3) * t60 - t36;
t78 = -t90 + t92;
t71 = Icges(4,1) * t62 - t85;
t69 = -Icges(4,2) * t60 + t84;
t67 = Icges(4,5) * t62 - Icges(4,6) * t60;
t64 = -pkin(6) - pkin(5);
t63 = cos(qJ(1));
t57 = t63 * pkin(1);
t50 = t62 * pkin(3) + pkin(2);
t49 = t54 * pkin(5);
t45 = t63 * rSges(2,1) - t61 * rSges(2,2);
t44 = -t61 * rSges(2,1) - t63 * rSges(2,2);
t43 = t60 * rSges(4,1) + t62 * rSges(4,2);
t37 = t54 * t50;
t31 = t54 * rSges(3,1) - t53 * rSges(3,2) + t57;
t30 = -t53 * rSges(3,1) - t54 * rSges(3,2) - t93;
t25 = Icges(4,3) * t53 + t67 * t54;
t24 = -Icges(4,3) * t54 + t67 * t53;
t23 = t79 * t54;
t22 = t79 * t53;
t13 = t94 + t57 + (pkin(2) - t90) * t54 + t87;
t12 = t88 - t93 + t49 + (-pkin(2) - t78) * t53;
t11 = -t53 * t64 + t37 + t57 + t65;
t10 = -t93 + (rSges(5,3) - t64) * t54 + (-t50 - t77) * t53;
t9 = t54 * (-t54 * t90 + t87) + (t78 * t53 - t88) * t53;
t3 = (t49 + (-pkin(2) + t50) * t53) * t53 + (-t54 * pkin(2) + t37 - t94) * t54 + t8;
t1 = [t56 * t34 + t55 * t35 + t62 * (Icges(4,2) * t62 + t85) + t60 * (Icges(4,1) * t60 + t84) + Icges(2,3) + Icges(3,3) + m(2) * (t44 ^ 2 + t45 ^ 2) + m(3) * (t30 ^ 2 + t31 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2); 0; m(3) + m(4) + m(5); (t62 * (Icges(4,6) * t53 + t69 * t54) + t60 * (Icges(4,5) * t53 + t71 * t54)) * t96 + (t62 * (-Icges(4,6) * t54 + t69 * t53) + t60 * (-Icges(4,5) * t54 + t71 * t53)) * t95 + m(5) * (t23 * t10 + t22 * t11) + m(4) * (-t12 * t54 - t13 * t53) * t43 + (t51 / 0.2e1 + t52 / 0.2e1) * (Icges(4,5) * t60 + Icges(4,6) * t62) + t81; m(4) * t9 + m(5) * t3; m(4) * (t86 * t43 ^ 2 + t9 ^ 2) + t53 * (-t24 * t97 + t51 * t25) - t54 * (t52 * t24 - t25 * t97) + m(5) * (t22 ^ 2 + t23 ^ 2 + t3 ^ 2) + t80; m(5) * (-t10 * t54 - t11 * t53) * t36 + t81; m(5) * t8; m(5) * (t8 * t3 + (-t22 * t53 - t23 * t54) * t36) + t80; m(5) * (t86 * t36 ^ 2 + t8 ^ 2) + t80;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
