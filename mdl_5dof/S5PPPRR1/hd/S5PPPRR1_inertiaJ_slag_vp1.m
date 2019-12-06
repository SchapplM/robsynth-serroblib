% Calculate joint inertia matrix for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:50
% EndTime: 2019-12-05 14:57:54
% DurationCPUTime: 0.81s
% Computational Cost: add. (2864->172), mult. (4047->292), div. (0->0), fcn. (4835->8), ass. (0->86)
t70 = pkin(9) + qJ(4);
t65 = sin(t70);
t66 = cos(t70);
t71 = sin(pkin(8));
t73 = cos(pkin(8));
t51 = -Icges(5,3) * t73 + (Icges(5,5) * t66 - Icges(5,6) * t65) * t71;
t94 = t73 * t51;
t93 = t73 ^ 2;
t92 = t65 * t71;
t91 = t66 * t71;
t72 = sin(pkin(7));
t90 = t72 * t71;
t89 = t72 * t73;
t74 = cos(pkin(7));
t88 = t74 * t65;
t87 = t74 * t66;
t86 = t74 * t71;
t56 = t66 * t89 - t88;
t75 = sin(qJ(5));
t76 = cos(qJ(5));
t47 = -t56 * t75 + t76 * t90;
t48 = t56 * t76 + t75 * t90;
t55 = t65 * t89 + t87;
t29 = t48 * rSges(6,1) + t47 * rSges(6,2) + t55 * rSges(6,3);
t85 = t56 * pkin(4) + t55 * pkin(6) + t29;
t58 = t72 * t65 + t73 * t87;
t49 = -t58 * t75 + t76 * t86;
t50 = t58 * t76 + t75 * t86;
t57 = -t72 * t66 + t73 * t88;
t30 = t50 * rSges(6,1) + t49 * rSges(6,2) + t57 * rSges(6,3);
t84 = -t58 * pkin(4) - t57 * pkin(6) - t30;
t60 = -t73 * t76 - t75 * t91;
t61 = -t73 * t75 + t76 * t91;
t44 = t61 * rSges(6,1) + t60 * rSges(6,2) + rSges(6,3) * t92;
t83 = t44 + (pkin(4) * t66 + pkin(6) * t65) * t71;
t82 = t72 ^ 2 + t74 ^ 2;
t81 = Icges(5,5) * t71;
t80 = Icges(5,6) * t71;
t79 = Icges(5,3) * t71;
t78 = -m(4) - m(5) - m(6);
t77 = m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1;
t54 = -t73 * rSges(5,3) + (rSges(5,1) * t66 - rSges(5,2) * t65) * t71;
t53 = -Icges(5,5) * t73 + (Icges(5,1) * t66 - Icges(5,4) * t65) * t71;
t52 = -Icges(5,6) * t73 + (Icges(5,4) * t66 - Icges(5,2) * t65) * t71;
t43 = Icges(6,1) * t61 + Icges(6,4) * t60 + Icges(6,5) * t92;
t42 = Icges(6,4) * t61 + Icges(6,2) * t60 + Icges(6,6) * t92;
t41 = Icges(6,5) * t61 + Icges(6,6) * t60 + Icges(6,3) * t92;
t40 = t58 * rSges(5,1) - t57 * rSges(5,2) + rSges(5,3) * t86;
t39 = t56 * rSges(5,1) - t55 * rSges(5,2) + rSges(5,3) * t90;
t38 = Icges(5,1) * t58 - Icges(5,4) * t57 + t74 * t81;
t37 = Icges(5,1) * t56 - Icges(5,4) * t55 + t72 * t81;
t36 = Icges(5,4) * t58 - Icges(5,2) * t57 + t74 * t80;
t35 = Icges(5,4) * t56 - Icges(5,2) * t55 + t72 * t80;
t34 = Icges(5,5) * t58 - Icges(5,6) * t57 + t74 * t79;
t33 = Icges(5,5) * t56 - Icges(5,6) * t55 + t72 * t79;
t32 = -t73 * t40 - t54 * t86;
t31 = t73 * t39 + t54 * t90;
t28 = Icges(6,1) * t50 + Icges(6,4) * t49 + Icges(6,5) * t57;
t27 = Icges(6,1) * t48 + Icges(6,4) * t47 + Icges(6,5) * t55;
t26 = Icges(6,4) * t50 + Icges(6,2) * t49 + Icges(6,6) * t57;
t25 = Icges(6,4) * t48 + Icges(6,2) * t47 + Icges(6,6) * t55;
t24 = Icges(6,5) * t50 + Icges(6,6) * t49 + Icges(6,3) * t57;
t23 = Icges(6,5) * t48 + Icges(6,6) * t47 + Icges(6,3) * t55;
t22 = (t39 * t74 - t40 * t72) * t71;
t21 = t30 * t92 - t57 * t44;
t20 = -t29 * t92 + t55 * t44;
t19 = t41 * t92 + t60 * t42 + t61 * t43;
t18 = t84 * t73 - t83 * t86;
t17 = t85 * t73 + t83 * t90;
t16 = t57 * t29 - t55 * t30;
t15 = t57 * t41 + t49 * t42 + t50 * t43;
t14 = t55 * t41 + t47 * t42 + t48 * t43;
t13 = (t84 * t72 + t85 * t74) * t71;
t12 = t24 * t92 + t60 * t26 + t61 * t28;
t11 = t23 * t92 + t60 * t25 + t61 * t27;
t10 = t57 * t24 + t49 * t26 + t50 * t28;
t9 = t57 * t23 + t49 * t25 + t50 * t27;
t8 = t55 * t24 + t47 * t26 + t48 * t28;
t7 = t55 * t23 + t47 * t25 + t48 * t27;
t6 = -t19 * t73 + (t11 * t72 + t12 * t74) * t71;
t5 = t11 * t55 + t12 * t57 + t19 * t92;
t4 = -t15 * t73 + (t10 * t74 + t72 * t9) * t71;
t3 = -t14 * t73 + (t7 * t72 + t74 * t8) * t71;
t2 = t10 * t57 + t15 * t92 + t9 * t55;
t1 = t14 * t92 + t7 * t55 + t8 * t57;
t45 = [m(2) + m(3) - t78; 0; 0.2e1 * (m(3) / 0.2e1 + t77) * t82; t78 * t73; 0; 0.2e1 * t77 * (t82 * t71 ^ 2 + t93); m(5) * t22 + m(6) * t13; m(5) * (t31 * t72 - t32 * t74) + m(6) * (t17 * t72 - t18 * t74); m(5) * (-t22 * t73 + (t31 * t74 + t32 * t72) * t71) + m(6) * (-t13 * t73 + (t17 * t74 + t18 * t72) * t71); m(5) * (t22 ^ 2 + t31 ^ 2 + t32 ^ 2) - t73 * (t93 * t51 + (((-t36 * t65 + t38 * t66) * t74 + (-t35 * t65 + t37 * t66) * t72) * t71 + (-t33 * t72 - t34 * t74 + t52 * t65 - t53 * t66) * t73) * t71) + m(6) * (t13 ^ 2 + t17 ^ 2 + t18 ^ 2) - t73 * t6 + (-(-t55 * t52 + t56 * t53) * t73 + t3 + (t33 * t90 - t55 * t35 + t56 * t37 - t94) * t90) * t90 + (-(-t57 * t52 + t58 * t53) * t73 + t4 + (t34 * t86 - t57 * t36 + t58 * t38 - t94) * t86 + (t33 * t86 + t34 * t90 - t57 * t35 - t55 * t36 + t58 * t37 + t56 * t38) * t90) * t86; m(6) * t16; m(6) * (t20 * t72 - t21 * t74); m(6) * (-t16 * t73 + (t20 * t74 + t21 * t72) * t71); m(6) * (t16 * t13 + t20 * t17 + t21 * t18) + t57 * t4 / 0.2e1 + t55 * t3 / 0.2e1 - t73 * t5 / 0.2e1 + (t74 * t2 / 0.2e1 + t72 * t1 / 0.2e1 + t65 * t6 / 0.2e1) * t71; m(6) * (t16 ^ 2 + t20 ^ 2 + t21 ^ 2) + t57 * t2 + t55 * t1 + t5 * t92;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t45(1), t45(2), t45(4), t45(7), t45(11); t45(2), t45(3), t45(5), t45(8), t45(12); t45(4), t45(5), t45(6), t45(9), t45(13); t45(7), t45(8), t45(9), t45(10), t45(14); t45(11), t45(12), t45(13), t45(14), t45(15);];
Mq = res;
