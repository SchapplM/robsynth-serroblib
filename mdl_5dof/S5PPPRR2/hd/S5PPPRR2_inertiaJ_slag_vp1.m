% Calculate joint inertia matrix for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (2704->177), mult. (7027->286), div. (0->0), fcn. (8955->10), ass. (0->86)
t95 = cos(qJ(4));
t75 = sin(pkin(9));
t76 = sin(pkin(8));
t94 = t76 * t75;
t82 = sin(qJ(4));
t93 = t76 * t82;
t77 = sin(pkin(7));
t79 = cos(pkin(8));
t92 = t77 * t79;
t80 = cos(pkin(7));
t91 = t79 * t80;
t78 = cos(pkin(9));
t64 = -t75 * t80 + t78 * t92;
t57 = t64 * t95 + t77 * t93;
t63 = t75 * t92 + t78 * t80;
t81 = sin(qJ(5));
t83 = cos(qJ(5));
t47 = -t57 * t81 + t63 * t83;
t48 = t57 * t83 + t63 * t81;
t85 = t76 * t95;
t56 = t64 * t82 - t77 * t85;
t31 = rSges(6,1) * t48 + rSges(6,2) * t47 + rSges(6,3) * t56;
t90 = pkin(4) * t57 + pkin(6) * t56 + t31;
t66 = t75 * t77 + t78 * t91;
t59 = t66 * t95 + t80 * t93;
t65 = t75 * t91 - t77 * t78;
t49 = -t59 * t81 + t65 * t83;
t50 = t59 * t83 + t65 * t81;
t58 = t66 * t82 - t80 * t85;
t32 = rSges(6,1) * t50 + rSges(6,2) * t49 + rSges(6,3) * t58;
t89 = pkin(4) * t59 + pkin(6) * t58 + t32;
t68 = t78 * t85 - t79 * t82;
t60 = -t68 * t81 + t83 * t94;
t61 = t68 * t83 + t81 * t94;
t67 = t78 * t93 + t79 * t95;
t44 = rSges(6,1) * t61 + rSges(6,2) * t60 + rSges(6,3) * t67;
t88 = pkin(4) * t68 + pkin(6) * t67 + t44;
t87 = t77 ^ 2 + t80 ^ 2;
t86 = -m(4) - m(5) - m(6);
t84 = m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1;
t54 = rSges(5,1) * t68 - rSges(5,2) * t67 + rSges(5,3) * t94;
t53 = Icges(5,1) * t68 - Icges(5,4) * t67 + Icges(5,5) * t94;
t52 = Icges(5,4) * t68 - Icges(5,2) * t67 + Icges(5,6) * t94;
t51 = Icges(5,5) * t68 - Icges(5,6) * t67 + Icges(5,3) * t94;
t43 = Icges(6,1) * t61 + Icges(6,4) * t60 + Icges(6,5) * t67;
t42 = Icges(6,4) * t61 + Icges(6,2) * t60 + Icges(6,6) * t67;
t41 = Icges(6,5) * t61 + Icges(6,6) * t60 + Icges(6,3) * t67;
t40 = rSges(5,1) * t59 - rSges(5,2) * t58 + rSges(5,3) * t65;
t39 = rSges(5,1) * t57 - rSges(5,2) * t56 + rSges(5,3) * t63;
t38 = Icges(5,1) * t59 - Icges(5,4) * t58 + Icges(5,5) * t65;
t37 = Icges(5,1) * t57 - Icges(5,4) * t56 + Icges(5,5) * t63;
t36 = Icges(5,4) * t59 - Icges(5,2) * t58 + Icges(5,6) * t65;
t35 = Icges(5,4) * t57 - Icges(5,2) * t56 + Icges(5,6) * t63;
t34 = Icges(5,5) * t59 - Icges(5,6) * t58 + Icges(5,3) * t65;
t33 = Icges(5,5) * t57 - Icges(5,6) * t56 + Icges(5,3) * t63;
t30 = Icges(6,1) * t50 + Icges(6,4) * t49 + Icges(6,5) * t58;
t29 = Icges(6,1) * t48 + Icges(6,4) * t47 + Icges(6,5) * t56;
t28 = Icges(6,4) * t50 + Icges(6,2) * t49 + Icges(6,6) * t58;
t27 = Icges(6,4) * t48 + Icges(6,2) * t47 + Icges(6,6) * t56;
t26 = Icges(6,5) * t50 + Icges(6,6) * t49 + Icges(6,3) * t58;
t25 = Icges(6,5) * t48 + Icges(6,6) * t47 + Icges(6,3) * t56;
t24 = t40 * t94 - t54 * t65;
t23 = -t39 * t94 + t54 * t63;
t22 = t39 * t65 - t40 * t63;
t21 = t32 * t67 - t44 * t58;
t20 = -t31 * t67 + t44 * t56;
t19 = t41 * t67 + t42 * t60 + t43 * t61;
t18 = t31 * t58 - t32 * t56;
t17 = -t65 * t88 + t89 * t94;
t16 = t63 * t88 - t90 * t94;
t15 = t41 * t58 + t42 * t49 + t43 * t50;
t14 = t41 * t56 + t42 * t47 + t43 * t48;
t13 = -t63 * t89 + t65 * t90;
t12 = t26 * t67 + t28 * t60 + t30 * t61;
t11 = t25 * t67 + t27 * t60 + t29 * t61;
t10 = t26 * t58 + t28 * t49 + t30 * t50;
t9 = t25 * t58 + t27 * t49 + t29 * t50;
t8 = t26 * t56 + t28 * t47 + t30 * t48;
t7 = t25 * t56 + t27 * t47 + t29 * t48;
t6 = t11 * t63 + t12 * t65 + t19 * t94;
t5 = t11 * t56 + t12 * t58 + t19 * t67;
t4 = t10 * t65 + t15 * t94 + t63 * t9;
t3 = t14 * t94 + t63 * t7 + t65 * t8;
t2 = t10 * t58 + t15 * t67 + t56 * t9;
t1 = t14 * t67 + t56 * t7 + t58 * t8;
t45 = [m(2) + m(3) - t86; 0; 0.2e1 * (m(3) / 0.2e1 + t84) * t87; t86 * t79; 0; 0.2e1 * t84 * (t76 ^ 2 * t87 + t79 ^ 2); m(5) * t22 + m(6) * t13; m(5) * (t23 * t77 - t24 * t80) + m(6) * (t16 * t77 - t17 * t80); m(5) * (-t22 * t79 + (t23 * t80 + t24 * t77) * t76) + m(6) * (-t13 * t79 + (t16 * t80 + t17 * t77) * t76); m(5) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t13 ^ 2 + t16 ^ 2 + t17 ^ 2) + (t6 + (t51 * t94 - t52 * t67 + t53 * t68) * t94) * t94 + (t4 + (t34 * t65 - t36 * t58 + t38 * t59) * t65 + (t34 * t94 - t36 * t67 + t38 * t68 + t51 * t65 - t52 * t58 + t53 * t59) * t94) * t65 + (t3 + (t33 * t63 - t35 * t56 + t37 * t57) * t63 + (t33 * t94 - t35 * t67 + t37 * t68 + t51 * t63 - t52 * t56 + t53 * t57) * t94 + (t33 * t65 + t34 * t63 - t35 * t58 - t36 * t56 + t37 * t59 + t38 * t57) * t65) * t63; m(6) * t18; m(6) * (t20 * t77 - t21 * t80); m(6) * (-t18 * t79 + (t20 * t80 + t21 * t77) * t76); m(6) * (t13 * t18 + t16 * t20 + t17 * t21) + t65 * t2 / 0.2e1 + t58 * t4 / 0.2e1 + t63 * t1 / 0.2e1 + t56 * t3 / 0.2e1 + t5 * t94 / 0.2e1 + t67 * t6 / 0.2e1; m(6) * (t18 ^ 2 + t20 ^ 2 + t21 ^ 2) + t58 * t2 + t56 * t1 + t67 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t45(1), t45(2), t45(4), t45(7), t45(11); t45(2), t45(3), t45(5), t45(8), t45(12); t45(4), t45(5), t45(6), t45(9), t45(13); t45(7), t45(8), t45(9), t45(10), t45(14); t45(11), t45(12), t45(13), t45(14), t45(15);];
Mq = res;
